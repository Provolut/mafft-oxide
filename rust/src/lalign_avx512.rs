//! The AVX-512 kernel for `L__align11`'s DP row, transliterated from
//! `core/Lalign11.c`'s `LALIGN11_AVX512` block.
//!
//! This is the only module in the port that contains `unsafe`, and it exists
//! because measurement left no alternative. `l_align11` is 86% of `tbfast` and
//! `tbfast` is 80% of the `mafft --auto` a GEMME/ESCOTT realign runs, so this
//! one loop is roughly 69% of the aligner's cost; the port's scalar version is
//! already within 8% of the C's scalar version, and neither GCC nor LLVM will
//! vectorise it. `bench/OPTIMISATION-LOG.md` records both measurements, and the
//! C carries a hand-written kernel for exactly the same reason.
//!
//! **Why it cannot be written safely.** Portable SIMD (`std::simd`) is still
//! nightly-only, so explicit vectorisation on stable Rust means
//! `core::arch::x86_64` intrinsics, which are `unsafe` because calling one on a
//! CPU lacking the feature is undefined. That is the *whole* of the unsafety
//! here, and [`row`] discharges it the way the standard library intends: it is
//! `#[target_feature]`-gated and reached only through [`row_checked`], which
//! asks the CPU at runtime and asserts the lengths. Every memory access inside
//! is a `loadu`/`storeu` within a slice the wrapper has already measured -- the
//! loop stops at `jj + 8 <= lastj`.
//!
//! This is better than the C's arrangement, not merely equal to it: the C
//! compiles the kernel in behind `#if defined(__AVX512F__)`, so a binary built
//! with `-march=native` crashes on an older machine and a portable build never
//! gets the kernel at all. The dispatch here is per-run, so one binary is fast
//! where the hardware allows and correct everywhere else.
//!
//! **Why it is bit-exact.** Every lane performs the scalar computation on its
//! own cell: the same comparisons (`_CMP_GT_OQ` / `_CMP_LT_OQ`) against the
//! same operands in the same order, and the same additions -- nothing is
//! reassociated, and there is no horizontal sum. The one place lanes interact
//! is the running maximum, and the tie rule restores scan order exactly: within
//! a lane the strict `>` keeps the earliest `j`, and across lanes the smallest
//! `j` wins an equal score, which is the cell the scalar loop would have
//! stopped on.
#![allow(unsafe_code)]

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256i, __m512d, _mm256_add_epi32, _mm256_loadu_si256, _mm256_mask_blend_epi32,
    _mm256_mask_permutexvar_epi32, _mm256_set1_epi32, _mm256_setr_epi32, _mm256_setzero_si256,
    _mm256_storeu_si256, _mm256_sub_epi32, _mm512_add_pd, _mm512_cmp_pd_mask, _mm512_cvtsd_f64,
    _mm512_loadu_pd, _mm512_mask_blend_pd, _mm512_mask_permutexvar_pd, _mm512_permutexvar_pd,
    _mm512_set1_pd, _mm512_setr_epi64, _mm512_setr_pd, _mm512_storeu_pd, _CMP_GE_OQ, _CMP_GT_OQ,
    _CMP_LT_OQ,
};

/// Does this CPU have the two feature sets the kernel needs?
///
/// The C tests `__AVX512F__` and `__AVX512VL__` at compile time; this asks the
/// running processor instead. `is_x86_feature_detected!` caches its answer, so
/// the check costs a load per row after the first call.
///
/// Returns:
///     True when [`row`] may be called.
#[cfg(target_arch = "x86_64")]
pub fn available() -> bool {
    // Decided once. `is_x86_feature_detected!` caches on its own, but the
    // environment lookup does not, and this is asked once per DP row.
    static AVAILABLE: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *AVAILABLE.get_or_init(|| {
        // An escape hatch for bisecting a suspected kernel bug, and for proving
        // the two paths agree: `verify_lalign_kernel.sh` runs the same input
        // twice, once with this set, and requires identical output.
        if std::env::var_os("MAFFT_RS_NO_AVX512").is_some() {
            return false;
        }
        is_x86_feature_detected!("avx512f") && is_x86_feature_detected!("avx512vl")
    })
}

/// Always false off x86-64, where the kernel does not exist.
#[cfg(not(target_arch = "x86_64"))]
pub fn available() -> bool {
    false
}

/// The best cell a row's vectorised part found.
pub struct RowBest {
    /// How far along the row the kernel got; the caller's scalar tail resumes
    /// here.
    pub jj: usize,
    /// The highest score seen, and the column it was at. Meaningful only when
    /// `jj` moved past the caller's starting index.
    pub value: f64,
    pub j: i32,
}

/// Eight cells of one DP row at a time.
///
/// Mirrors the C block statement for statement. `mm`/`mpp`, `ijpi` and `cw` are
/// updated in place exactly as the scalar loop updates them.
///
/// Args:
///     jj0: First column to process (the C enters with `jj == 1`).
///     lastj: One past the last column.
///     pw: The previous row (`previousw`); read at `jj - 1`.
///     miv: The prefix-scanned row maxima.
///     mpiv: Their columns.
///     mm: Running column scores, updated.
///     mpp: Their rows, updated.
///     ijpi: This row of the traceback matrix, written.
///     cw: The current row, accumulated into.
///     i: Which row this is.
///     fpenalty: Gap-open penalty.
///     fpenalty_ex: Gap-extension penalty.
///     localthr: Score below which a cell is cut.
///     localthr2: What a cut cell contributes instead.
///     localstop: The traceback marker for a cut cell.
///
/// Returns:
///     Where the tail resumes, and the best cell the kernel saw.
///
/// # Safety
///
/// The caller must have checked [`available`]; every load and store stays
/// inside the slices, since the loop runs only while `jj + 8 <= lastj` and each
/// slice is at least `lastj` long.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512vl")]
#[allow(clippy::too_many_arguments)]
pub unsafe fn row(
    jj0: usize,
    lastj: usize,
    pw: &[f64],
    miv: &[f64],
    mpiv: &[i32],
    mm: &mut [f64],
    mpp: &mut [i32],
    ijpi: &mut [i32],
    cw: &mut [f64],
    i: i32,
    fpenalty: f64,
    fpenalty_ex: f64,
    localthr: f64,
    localthr2: f64,
    localstop: i32,
) -> RowBest {
    let vpen = _mm512_set1_pd(fpenalty);
    let vpenex = _mm512_set1_pd(fpenalty_ex);
    let vthr = _mm512_set1_pd(localthr);
    let vthr2 = _mm512_set1_pd(localthr2);
    let vstop = _mm256_set1_epi32(localstop);
    let vi = _mm256_set1_epi32(i);
    let vim1 = _mm256_set1_epi32(i - 1);
    let vstep = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut vbest: __m512d = _mm512_set1_pd(-1.0e308);
    let mut vbestj: __m256i = _mm256_setzero_si256();

    let mut jj = jj0;
    while jj + 8 <= lastj {
        let vjj = _mm256_add_epi32(_mm256_set1_epi32(jj as i32), vstep);
        let diag = _mm512_loadu_pd(pw.as_ptr().add(jj - 1));
        let gi = _mm512_add_pd(_mm512_loadu_pd(miv.as_ptr().add(jj)), vpen);
        let mmv = _mm512_loadu_pd(mm.as_ptr().add(jj));
        let gj = _mm512_add_pd(mmv, vpen);
        let mppv = _mm256_loadu_si256(mpp.as_ptr().add(jj).cast::<__m256i>());

        let ki = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(gi, diag);
        let w1v = _mm512_mask_blend_pd(ki, diag, gi);
        let ip1 = _mm256_mask_blend_epi32(
            ki,
            _mm256_setzero_si256(),
            _mm256_sub_epi32(
                _mm256_loadu_si256(mpiv.as_ptr().add(jj).cast::<__m256i>()),
                vjj,
            ),
        );

        let kj = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(gj, w1v);
        let w2v = _mm512_mask_blend_pd(kj, w1v, gj);
        let ip2 = _mm256_mask_blend_epi32(kj, ip1, _mm256_sub_epi32(vi, mppv));

        let km = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(diag, mmv);
        _mm256_storeu_si256(
            mpp.as_mut_ptr().add(jj).cast::<__m256i>(),
            _mm256_mask_blend_epi32(km, mppv, vim1),
        );
        _mm512_storeu_pd(
            mm.as_mut_ptr().add(jj),
            _mm512_add_pd(_mm512_mask_blend_pd(km, mmv, diag), vpenex),
        );

        let kbest = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(w2v, vbest);
        vbest = _mm512_mask_blend_pd(kbest, vbest, w2v);
        vbestj = _mm256_mask_blend_epi32(kbest, vbestj, vjj);

        let kstop = _mm512_cmp_pd_mask::<_CMP_LT_OQ>(w2v, vthr);
        _mm256_storeu_si256(
            ijpi.as_mut_ptr().add(jj).cast::<__m256i>(),
            _mm256_mask_blend_epi32(kstop, ip2, vstop),
        );
        _mm512_storeu_pd(
            cw.as_mut_ptr().add(jj),
            _mm512_add_pd(
                _mm512_loadu_pd(cw.as_ptr().add(jj)),
                _mm512_mask_blend_pd(kstop, w2v, vthr2),
            ),
        );

        jj += 8;
    }

    // One reduction per row rather than a second pass over it. Within a lane
    // the strict `>` above already kept the earliest `j`; across lanes the
    // smallest index wins a tie, which is the same cell the scalar loop would
    // have stopped at.
    let mut lane = [0.0f64; 8];
    let mut lanej = [0i32; 8];
    _mm512_storeu_pd(lane.as_mut_ptr(), vbest);
    _mm256_storeu_si256(lanej.as_mut_ptr().cast::<__m256i>(), vbestj);
    let mut bl = 0usize;
    for l in 1..8 {
        if lane[l] > lane[bl] || (lane[l] == lane[bl] && lanej[l] < lanej[bl]) {
            bl = l;
        }
    }

    RowBest {
        jj,
        value: lane[bl],
        j: lanej[bl],
    }
}

/// The kernel when this CPU has it, nothing when it does not.
///
/// The only way in: it performs the feature check itself, so no caller needs an
/// `unsafe` block and the whole of the port's unsafety stays in this file.
///
/// Args:
///     See [`row`]; the arguments are forwarded unchanged.
///
/// Returns:
///     Where the caller's scalar tail should resume and the best cell seen, or
///     `None` when the kernel is unavailable and the caller must run the whole
///     row itself.
///
/// Panics:
///     When `jj0` is zero, or when a buffer is shorter than the row -- the
///     bounds [`row`] would otherwise read past. Being safe to call is what
///     lets the rest of the crate stay `deny(unsafe_code)`, so these are
///     checked, not documented.
#[allow(clippy::too_many_arguments)]
pub fn row_checked(
    jj0: usize,
    lastj: usize,
    pw: &[f64],
    miv: &[f64],
    mpiv: &[i32],
    mm: &mut [f64],
    mpp: &mut [i32],
    ijpi: &mut [i32],
    cw: &mut [f64],
    i: i32,
    fpenalty: f64,
    fpenalty_ex: f64,
    localthr: f64,
    localthr2: f64,
    localstop: i32,
) -> Option<RowBest> {
    #[cfg(target_arch = "x86_64")]
    {
        if !available() {
            return None;
        }
        // `row`'s preconditions, enforced rather than assumed. This is a safe
        // function, so nothing may stop it being called with arguments that
        // would make the kernel read out of bounds -- a debug-only check would
        // leave exactly that hole in a release build. One comparison per row
        // against a loop that does eight cells at a time.
        assert!(jj0 >= 1, "the kernel reads `pw[jj - 1]`");
        assert!(
            pw.len() + 1 >= lastj,
            "previousw is read up to `lastj - 2`"
        );
        assert!(
            miv.len() >= lastj
                && mpiv.len() >= lastj
                && mm.len() >= lastj
                && mpp.len() >= lastj
                && ijpi.len() >= lastj
                && cw.len() >= lastj,
            "every other buffer is read and written up to `lastj - 1`"
        );
        // SAFETY: `available()` just confirmed avx512f and avx512vl on this
        // CPU, which is the whole of `row`'s target-feature requirement; the
        // asserts above establish the rest of it, so every load and store stays
        // inside its slice -- the loop stops at `jj + 8 <= lastj`.
        Some(unsafe {
            row(
                jj0, lastj, pw, miv, mpiv, mm, mpp, ijpi, cw, i, fpenalty, fpenalty_ex,
                localthr, localthr2, localstop,
            )
        })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (
            jj0, lastj, pw, miv, mpiv, mm, mpp, ijpi, cw, i, fpenalty, fpenalty_ex, localthr,
            localthr2, localstop,
        );
        None
    }
}

/// Eight positions of the `mi` prefix scan at a time.
///
/// `l_align11` runs a serial recurrence down each DP row before the row itself:
///
/// ```text
/// miv[j] = m;  mpiv[j] = mp;
/// if pw[j - 1] > m { m = pw[j - 1]; mp = j - 1; }
/// m += e;
/// ```
///
/// Each step waits for the previous one, so it costs a compare plus an add of
/// latency per position -- about eight cycles -- while doing almost no work.
/// With the DP row vectorised it became the largest single item in `tbfast`:
/// breaking the chain (and the answer with it) took a staged payload from
/// 0.294s to 0.167s, so 43% of the binary was this recurrence stalling.
///
/// It is a max-plus scan, and those parallelise. Writing `d_s = p_s + (8 - s)e`
/// puts every candidate on a common footing, so a prefix maximum over `d`
/// answers the whole block at once and only the carry between blocks stays
/// serial.
///
/// # Bit-exactness, and why it is conditional
///
/// Two rewrites happen here. The first is free: `max(a, b) + e == max(a + e,
/// b + e)` for every IEEE `a`, `b` and `e`, because rounding is monotone -- so
/// pushing the `+ e` inside the maximum cannot move a value.
///
/// The second is not. The scan reaches a position with one multiply, `m + t*e`,
/// where the loop got there by adding `e` to itself `t` times, and those agree
/// only when no rounding happened along the way. They do agree when every value
/// involved is an integer and the sums stay under 2^53, which is the ordinary
/// case here -- MAFFT's score matrices and gap penalties are integers -- but it
/// is a property of the input, not of the code, so [`scan_checked`] is handed
/// the answer and the caller falls back to the serial loop when it does not
/// hold.
///
/// That condition also settles the indices, which is the subtler half. Ties go
/// to the earliest position, and the serial loop and the scan compare the same
/// pair of candidates at different points along the decay: if rounding could
/// make two unequal values equal later on, the two would disagree about which
/// came first. Exact arithmetic keeps a difference a difference forever, so
/// they cannot.
///
/// Args:
///     pw: The previous DP row; position `s` reads `pw[s]`.
///     n: How many positions to scan.
///     m0: The running maximum entering the scan.
///     mp0: Its position.
///     e: The gap-extension penalty added after every position.
///     miv: Receives the running maximum before each position, at `1..=n`.
///     mpiv: Receives its position.
///
/// Returns:
///     How many positions were done, and the carry to resume from.
///
/// # Safety
///
/// The caller must have checked [`available`]. `pw` must hold at least `n`
/// values and `miv`/`mpiv` at least `n + 1`; the loop stops at `done + 8 <= n`.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx512f,avx512vl")]
pub unsafe fn prefix_scan(
    pw: &[f64],
    n: usize,
    m0: f64,
    mp0: i32,
    e: f64,
    miv: &mut [f64],
    mpiv: &mut [i32],
) -> (usize, f64, i32) {
    // `e * k` for a small integer `k` is exact whenever `e` is, which is the
    // precondition the caller checked.
    let v8ms_e = _mm512_setr_pd(
        8.0 * e,
        7.0 * e,
        6.0 * e,
        5.0 * e,
        4.0 * e,
        3.0 * e,
        2.0 * e,
        e,
    );
    let vt_e = _mm512_setr_pd(
        // `-0.0`, not `0.0 * e`. Position 0 must reproduce the incumbent
        // untouched, and adding `-0.0` is the only identity that does: with a
        // zero `e` the product is `+0.0`, which turns a `-0.0` incumbent
        // positive, and the serial loop has not added anything to it yet.
        -0.0,
        e,
        2.0 * e,
        3.0 * e,
        4.0 * e,
        5.0 * e,
        6.0 * e,
        7.0 * e,
    );
    let vtm8_e = _mm512_setr_pd(
        -8.0 * e,
        -7.0 * e,
        -6.0 * e,
        -5.0 * e,
        -4.0 * e,
        -3.0 * e,
        -2.0 * e,
        -e,
    );
    let e8 = 8.0 * e;
    let vstep = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let neg_inf = _mm512_set1_pd(f64::NEG_INFINITY);
    // Lane `t` takes lane `t - k`; the low `k` lanes are filled with -inf,
    // which loses every comparison against a real candidate.
    // The same shift, in both lane widths: `_mm512_permutexvar_pd` reads 64-bit
    // indices and `_mm256_permutexvar_epi32` reads 32-bit ones.
    let sh1 = _mm512_setr_epi64(0, 0, 1, 2, 3, 4, 5, 6);
    let sh2 = _mm512_setr_epi64(0, 0, 0, 1, 2, 3, 4, 5);
    let sh4 = _mm512_setr_epi64(0, 0, 0, 0, 0, 1, 2, 3);
    let sh1i = _mm256_setr_epi32(0, 0, 1, 2, 3, 4, 5, 6);
    let sh2i = _mm256_setr_epi32(0, 0, 0, 1, 2, 3, 4, 5);
    let sh4i = _mm256_setr_epi32(0, 0, 0, 0, 0, 1, 2, 3);
    let last = _mm512_setr_epi64(7, 7, 7, 7, 7, 7, 7, 7);

    let mut m = m0;
    let mut mp = mp0;
    let mut done = 0usize;
    while done + 8 <= n {
        let vidx = _mm256_add_epi32(_mm256_set1_epi32(done as i32), vstep);
        let p = _mm512_loadu_pd(pw.as_ptr().add(done));
        // Every candidate in the block, decayed to a common position.
        let d = _mm512_add_pd(p, v8ms_e);

        // Inclusive prefix maximum of `d`, earliest position winning a tie.
        let mut iv = d;
        let mut ii = vidx;
        for (shift, shifti, mask) in [(sh1, sh1i, 0xfeu8), (sh2, sh2i, 0xfcu8), (sh4, sh4i, 0xf0u8)]
        {
            let sv = _mm512_mask_permutexvar_pd(neg_inf, mask, shift, iv);
            let si = _mm256_mask_permutexvar_epi32(ii, mask, shifti, ii);
            // `>=`, not `>`: the shifted candidate is the earlier one, so it
            // takes a tie.
            let k = _mm512_cmp_pd_mask::<_CMP_GE_OQ>(sv, iv);
            iv = _mm512_mask_blend_pd(k, iv, sv);
            ii = _mm256_mask_blend_epi32(k, ii, si);
        }

        // Exclusive: position `t` may only use candidates before it.
        let ev = _mm512_mask_permutexvar_pd(neg_inf, 0xfe, sh1, iv);
        let ei = _mm256_mask_permutexvar_epi32(ii, 0xfe, sh1i, ii);

        let cand = _mm512_add_pd(ev, vtm8_e);
        let inc = _mm512_add_pd(_mm512_set1_pd(m), vt_e);
        // Strict `>`: the incumbent entered the block from an earlier position.
        let k = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(cand, inc);
        _mm512_storeu_pd(
            miv.as_mut_ptr().add(done + 1),
            _mm512_mask_blend_pd(k, inc, cand),
        );
        _mm256_storeu_si256(
            mpiv.as_mut_ptr().add(done + 1).cast::<__m256i>(),
            _mm256_mask_blend_epi32(k, _mm256_set1_epi32(mp), ei),
        );

        // The carry: the block's best against the incumbent decayed across all
        // eight positions. Neither term needs the block's own output, so this
        // is the only thing the next block waits for.
        let best = _mm512_cvtsd_f64(_mm512_permutexvar_pd(last, iv));
        let mut bestj = [0i32; 8];
        _mm256_storeu_si256(bestj.as_mut_ptr().cast::<__m256i>(), ii);
        let incumbent = m + e8;
        if best > incumbent {
            m = best;
            mp = bestj[7];
        } else {
            m = incumbent;
        }
        done += 8;
    }
    (done, m, mp)
}

/// [`prefix_scan`] when this CPU has the kernel and the arithmetic is exact.
///
/// Args:
///     exact: Whether every value the scan will touch is an integer small
///         enough that the sums cannot round. See [`prefix_scan`]; the caller
///         establishes this from the score matrix and the penalties.
///     Others: forwarded to [`prefix_scan`].
///
/// Returns:
///     How far the kernel got and the carry to resume from, or `None` when the
///     caller must run the whole scan itself.
///
/// Panics:
///     When a buffer is shorter than the scan -- the bounds [`prefix_scan`]
///     would otherwise write past.
#[allow(clippy::too_many_arguments)]
pub fn prefix_scan_checked(
    exact: bool,
    pw: &[f64],
    n: usize,
    m0: f64,
    mp0: i32,
    e: f64,
    miv: &mut [f64],
    mpiv: &mut [i32],
) -> Option<(usize, f64, i32)> {
    #[cfg(target_arch = "x86_64")]
    {
        if !exact || !available() {
            return None;
        }
        assert!(pw.len() >= n, "the scan reads `pw[..n]`");
        assert!(
            miv.len() > n && mpiv.len() > n,
            "the scan writes `1..=n`"
        );
        // SAFETY: `available()` just confirmed avx512f and avx512vl, and the
        // asserts cover the rest of `prefix_scan`'s contract -- it stops at
        // `done + 8 <= n`, so every load stays inside `pw[..n]` and every store
        // inside `1..=n`.
        Some(unsafe { prefix_scan(pw, n, m0, mp0, e, miv, mpiv) })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (exact, pw, n, m0, mp0, e, miv, mpiv);
        None
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;

    /// A row's worth of every buffer the kernel touches, long enough that only
    /// the argument under test is wrong.
    struct Row {
        pw: Vec<f64>,
        miv: Vec<f64>,
        mpiv: Vec<i32>,
        m: Vec<f64>,
        mp: Vec<i32>,
        ijp: Vec<i32>,
        cw: Vec<f64>,
    }

    impl Row {
        fn new() -> Self {
            Self {
                pw: vec![0.0; 8],
                miv: vec![0.0; 8],
                mpiv: vec![0; 8],
                m: vec![0.0; 8],
                mp: vec![0; 8],
                ijp: vec![0; 8],
                cw: vec![0.0; 8],
            }
        }

        fn call(&mut self, jj0: usize) {
            row_checked(
                jj0,
                8,
                &self.pw,
                &self.miv,
                &self.mpiv,
                &mut self.m,
                &mut self.mp,
                &mut self.ijp,
                &mut self.cw,
                1,
                0.0,
                0.0,
                0.0,
                0.0,
                0,
            );
        }
    }

    /// `pw.as_ptr().add(jj - 1)` wraps at `jj == 0`, so the wrapper has to stop
    /// it before the kernel is entered.
    #[test]
    #[should_panic(expected = "pw[jj - 1]")]
    fn rejects_a_zero_start() {
        assert!(available(), "pw[jj - 1]: skipped, no AVX-512 on this CPU");
        Row::new().call(0);
    }

    /// Half a row of `cw` -- the kernel would store eight doubles into it.
    #[test]
    #[should_panic(expected = "read and written up to")]
    fn rejects_a_short_buffer() {
        assert!(
            available(),
            "read and written up to: skipped, no AVX-512 on this CPU"
        );
        let mut row = Row::new();
        row.cw.truncate(4);
        row.call(1);
    }

    /// The serial recurrence [`prefix_scan`] replaces.
    fn scan_serially(pw: &[f64], n: usize, m0: f64, mp0: i32, e: f64) -> (Vec<f64>, Vec<i32>) {
        let mut miv = vec![0.0f64; n + 1];
        let mut mpiv = vec![0i32; n + 1];
        let (mut m, mut mp) = (m0, mp0);
        for idx in 0..n {
            miv[idx + 1] = m;
            mpiv[idx + 1] = mp;
            if pw[idx] > m {
                m = pw[idx];
                mp = idx as i32;
            }
            m += e;
        }
        (miv, mpiv)
    }

    /// The kernel plus the serial tail it leaves, as `l_align11` drives it.
    fn scan_with_kernel(pw: &[f64], n: usize, m0: f64, mp0: i32, e: f64) -> (Vec<f64>, Vec<i32>) {
        let mut miv = vec![0.0f64; n + 1];
        let mut mpiv = vec![0i32; n + 1];
        let (mut m, mut mp) = (m0, mp0);
        let mut start = 0usize;
        if let Some((done, cm, cmp)) =
            prefix_scan_checked(true, pw, n, m0, mp0, e, &mut miv, &mut mpiv)
        {
            start = done;
            m = cm;
            mp = cmp;
        }
        for idx in start..n {
            miv[idx + 1] = m;
            mpiv[idx + 1] = mp;
            if pw[idx] > m {
                m = pw[idx];
                mp = idx as i32;
            }
            m += e;
        }
        (miv, mpiv)
    }

    /// Every value and every position, bit for bit, over lengths that exercise
    /// the block loop, the tail, and a scan too short for either.
    ///
    /// The values are integers because that is the kernel's precondition: it
    /// reaches a position with one multiply where the loop adds `e` repeatedly.
    #[test]
    fn the_prefix_scan_matches_the_serial_one() {
        if !available() {
            return;
        }
        let mut x = 0x2545_f491_4f6c_dd1du64;
        let mut next = move || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        for n in [0usize, 1, 5, 7, 8, 9, 15, 16, 17, 31, 64, 101, 400] {
            for &e in &[-1.0f64, -3.0, -12.0, 0.0] {
                // A run of nearly-equal values is where the tie rule bites, so
                // the spread is deliberately narrow.
                let pw: Vec<f64> = (0..n.max(1))
                    .map(|_| ((next() % 41) as f64) - 20.0)
                    .collect();
                for &m0 in &[-100.0f64, 0.0, -0.0, 7.0, 20.0] {
                    let a = scan_serially(&pw, n, m0, -3, e);
                    let b = scan_with_kernel(&pw, n, m0, -3, e);
                    assert_eq!(
                        a.0.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                        b.0.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                        "values differ: n={n} e={e} m0={m0}"
                    );
                    assert_eq!(a.1, b.1, "positions differ: n={n} e={e} m0={m0}");
                }
            }
        }
    }

    /// The lengths the aligner actually passes are accepted.
    #[test]
    fn accepts_a_whole_row() {
        if !available() {
            return;
        }
        let mut row = Row::new();
        row.call(1);
    }
}
