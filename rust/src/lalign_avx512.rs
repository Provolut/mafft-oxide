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
    _mm256_set1_epi32, _mm256_setr_epi32, _mm256_setzero_si256, _mm256_storeu_si256,
    _mm256_sub_epi32, _mm512_add_pd, _mm512_cmp_pd_mask, _mm512_loadu_pd, _mm512_mask_blend_pd,
    _mm512_set1_pd, _mm512_storeu_pd, _CMP_GT_OQ, _CMP_LT_OQ,
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
