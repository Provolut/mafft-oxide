//! The AVX-512 kernel for `partA__align`'s DP row, transliterated from the
//! scalar loop in `core/partSalignmm.c`.
//!
//! The second of the port's two kernels, and the second place it contains
//! `unsafe`; [`crate::lalign_avx512`] carries the argument for why explicit
//! vectorisation on stable Rust means `core::arch` intrinsics at all, and this
//! module discharges the obligation the same way -- `#[target_feature]`-gated,
//! reached only through a safe wrapper that asks the CPU at runtime and asserts
//! the lengths.
//!
//! **Why this loop, and why now.** `perf record` on a staged ESCOTT payload
//! puts `part_a_align` at 20.4% of `dvtditr`, which is itself 41% of a cold
//! GEMME/ESCOTT run's aligner cost. The C has no kernel here either, so unlike
//! `L__align11` this is not catching up with the C -- it is ahead of it.
//!
//! **Why it is bit-exact.** The row is element-wise in `j`: every lane runs the
//! scalar computation on its own cell, with the same comparisons
//! (`_CMP_GT_OQ`) against the same operands in the same order and the same
//! multiplications and additions. Nothing is reassociated. Unlike
//! [`crate::lalign_avx512::row`] there is no running maximum, so no lane ever
//! reads another lane's value and there is no tie rule to preserve -- the only
//! cross-lane operation in this file does not exist.
//!
//! The `mi` recurrence that feeds the row is *not* here: it is the same
//! max-plus scan `l_align11` runs, so [`crate::lalign_avx512::prefix_scan`]
//! does it, over a candidate row the caller builds. See
//! [`crate::partalign::PartAlignState::part_a_align`].
#![allow(unsafe_code)]

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::{
    __m256i, _mm256_add_epi32, _mm256_loadu_si256, _mm256_mask_blend_epi32, _mm256_set1_epi32,
    _mm256_setr_epi32, _mm256_setzero_si256, _mm256_storeu_si256, _mm256_sub_epi32, _mm512_add_pd,
    _mm512_cmp_pd_mask, _mm512_loadu_pd, _mm512_mask_blend_pd, _mm512_mul_pd, _mm512_set1_pd,
    _mm512_storeu_pd, _CMP_GT_OQ,
};

/// Eight cells of one `partA__align` DP row at a time.
///
/// Mirrors the scalar loop statement for statement. `mrow`/`mprow`, `ijpr` and
/// `curw` are updated in place exactly as it updates them.
///
/// Args:
///     jj0: First column to process (the caller enters with `jj == 1`).
///     lastj: One past the last column.
///     prevw: The previous row, read at `jj - 1`.
///     miv: The prefix-scanned `mi` values.
///     mpiv: Their columns.
///     fgcp2: Final-gap counts for group 2, read at `jj - 1`.
///     gapfreq2: Gap frequencies for group 2, read at `jj` and `jj - 1`.
///     mrow: Running column scores, updated.
///     mprow: Their rows, updated.
///     ijpr: This row of the traceback matrix, written.
///     curw: The current row, accumulated into.
///     i: Which row this is.
///     gf1va: `gapfreq1[i]`.
///     fgcp1va: `fgcp1[i - 1]`.
///     ogcp1va: `ogcp1[i]`.
///     fpenalty_ex: Gap-extension penalty, added to `mrow` after each cell.
///
/// Returns:
///     Where the caller's scalar tail resumes, and the last cell's `w` -- the
///     running score the caller returns when the tail is empty.
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
    prevw: &[f64],
    miv: &[f64],
    mpiv: &[i32],
    fgcp2: &[f64],
    gapfreq2: &[f64],
    mrow: &mut [f64],
    mprow: &mut [i32],
    ijpr: &mut [i32],
    curw: &mut [f64],
    i: i32,
    gf1va: f64,
    fgcp1va: f64,
    ogcp1va: f64,
    fpenalty_ex: f64,
) -> (usize, f64) {
    let vgf1va = _mm512_set1_pd(gf1va);
    let vfgcp1va = _mm512_set1_pd(fgcp1va);
    let vogcp1va = _mm512_set1_pd(ogcp1va);
    let vpenex = _mm512_set1_pd(fpenalty_ex);
    let vi = _mm256_set1_epi32(i);
    let vim1 = _mm256_set1_epi32(i - 1);
    let vstep = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    let mut lastw = [0.0f64; 8];

    let mut jj = jj0;
    while jj + 8 <= lastj {
        let vjj = _mm256_add_epi32(_mm256_set1_epi32(jj as i32), vstep);
        let prev = _mm512_loadu_pd(prevw.as_ptr().add(jj - 1));

        // `w = prevw[j - 1]; ijpr[j] = 0;`
        let mut w = prev;
        let mut ijp = _mm256_setzero_si256();

        // `g = miv[j] + fgcp2[j - 1] * gf1va; if g > w { .. }`
        let g = _mm512_add_pd(
            _mm512_loadu_pd(miv.as_ptr().add(jj)),
            _mm512_mul_pd(_mm512_loadu_pd(fgcp2.as_ptr().add(jj - 1)), vgf1va),
        );
        let k = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(g, w);
        w = _mm512_mask_blend_pd(k, w, g);
        // `-(j - mpiv[j])`, which is `mpiv[j] - j`.
        ijp = _mm256_mask_blend_epi32(
            k,
            ijp,
            _mm256_sub_epi32(
                _mm256_loadu_si256(mpiv.as_ptr().add(jj).cast::<__m256i>()),
                vjj,
            ),
        );

        // `g = mrow[j] + fgcp1va * gapfreq2[j]; if g > w { .. }`
        let mrowv = _mm512_loadu_pd(mrow.as_ptr().add(jj));
        let mprowv = _mm256_loadu_si256(mprow.as_ptr().add(jj).cast::<__m256i>());
        let g = _mm512_add_pd(
            mrowv,
            _mm512_mul_pd(vfgcp1va, _mm512_loadu_pd(gapfreq2.as_ptr().add(jj))),
        );
        let k = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(g, w);
        w = _mm512_mask_blend_pd(k, w, g);
        ijp = _mm256_mask_blend_epi32(k, ijp, _mm256_sub_epi32(vi, mprowv));

        // `g = prevw[j - 1] + ogcp1va * gapfreq2[j - 1]; if g > mrow[j] { .. }`
        let g = _mm512_add_pd(
            prev,
            _mm512_mul_pd(vogcp1va, _mm512_loadu_pd(gapfreq2.as_ptr().add(jj - 1))),
        );
        let k = _mm512_cmp_pd_mask::<_CMP_GT_OQ>(g, mrowv);
        _mm256_storeu_si256(
            mprow.as_mut_ptr().add(jj).cast::<__m256i>(),
            _mm256_mask_blend_epi32(k, mprowv, vim1),
        );
        _mm512_storeu_pd(
            mrow.as_mut_ptr().add(jj),
            _mm512_add_pd(_mm512_mask_blend_pd(k, mrowv, g), vpenex),
        );

        _mm256_storeu_si256(ijpr.as_mut_ptr().add(jj).cast::<__m256i>(), ijp);
        _mm512_storeu_pd(
            curw.as_mut_ptr().add(jj),
            _mm512_add_pd(_mm512_loadu_pd(curw.as_ptr().add(jj)), w),
        );

        // The caller returns the last cell's `w` as the alignment score, and
        // the scalar tail is empty whenever the row divides by eight.
        _mm512_storeu_pd(lastw.as_mut_ptr(), w);
        jj += 8;
    }

    (jj, lastw[7])
}

/// The kernel when this CPU has it, nothing when it does not.
///
/// The only way in: it performs the feature check itself, so no caller needs an
/// `unsafe` block and the whole of this module's unsafety stays here.
///
/// Args:
///     See [`row`]; the arguments are forwarded unchanged.
///
/// Returns:
///     Where the caller's scalar tail should resume and the last cell's `w`, or
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
    prevw: &[f64],
    miv: &[f64],
    mpiv: &[i32],
    fgcp2: &[f64],
    gapfreq2: &[f64],
    mrow: &mut [f64],
    mprow: &mut [i32],
    ijpr: &mut [i32],
    curw: &mut [f64],
    i: i32,
    gf1va: f64,
    fgcp1va: f64,
    ogcp1va: f64,
    fpenalty_ex: f64,
) -> Option<(usize, f64)> {
    #[cfg(target_arch = "x86_64")]
    {
        if !crate::lalign_avx512::available() {
            return None;
        }
        assert!(jj0 >= 1, "the kernel reads `prevw[jj - 1]`");
        assert!(
            prevw.len() + 1 >= lastj && fgcp2.len() + 1 >= lastj,
            "both are read up to `lastj - 2`"
        );
        assert!(
            miv.len() >= lastj
                && mpiv.len() >= lastj
                && gapfreq2.len() >= lastj
                && mrow.len() >= lastj
                && mprow.len() >= lastj
                && ijpr.len() >= lastj
                && curw.len() >= lastj,
            "every other buffer is read and written up to `lastj - 1`"
        );
        // SAFETY: `available()` just confirmed avx512f and avx512vl on this
        // CPU, which is the whole of `row`'s target-feature requirement; the
        // asserts above establish the rest of it, so every load and store stays
        // inside its slice -- the loop stops at `jj + 8 <= lastj`.
        Some(unsafe {
            row(
                jj0, lastj, prevw, miv, mpiv, fgcp2, gapfreq2, mrow, mprow, ijpr, curw, i, gf1va,
                fgcp1va, ogcp1va, fpenalty_ex,
            )
        })
    }
    #[cfg(not(target_arch = "x86_64"))]
    {
        let _ = (
            jj0, lastj, prevw, miv, mpiv, fgcp2, gapfreq2, mrow, mprow, ijpr, curw, i, gf1va,
            fgcp1va, ogcp1va, fpenalty_ex,
        );
        None
    }
}

#[cfg(all(test, target_arch = "x86_64"))]
mod tests {
    use super::*;

    /// Reproducible pseudo-random `f64`s, so a failure is the same failure
    /// twice. Values land in roughly `-4..4`, the range MAFFT's scaled scores
    /// and gap terms occupy, and every eighth is made negative zero -- the one
    /// value that distinguishes `x` from `x + 0.0` and so the one a
    /// bit-for-bit comparison must be given the chance to catch.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
            let x = ((self.0 >> 11) as f64) / ((1u64 << 53) as f64);
            if self.0 % 8 == 0 {
                -0.0
            } else {
                x * 8.0 - 4.0
            }
        }

        fn next_i32(&mut self) -> i32 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((self.0 >> 33) % 64) as i32
        }
    }

    /// Every buffer one row of `partA__align` reads or writes.
    #[derive(Clone)]
    struct Row {
        prevw: Vec<f64>,
        miv: Vec<f64>,
        mpiv: Vec<i32>,
        fgcp2: Vec<f64>,
        gapfreq2: Vec<f64>,
        mrow: Vec<f64>,
        mprow: Vec<i32>,
        ijpr: Vec<i32>,
        curw: Vec<f64>,
    }

    impl Row {
        fn random(len: usize, seed: u64) -> Self {
            let mut r = Lcg(seed);
            Self {
                prevw: (0..len).map(|_| r.next_f64()).collect(),
                miv: (0..len).map(|_| r.next_f64()).collect(),
                mpiv: (0..len).map(|_| r.next_i32()).collect(),
                fgcp2: (0..len).map(|_| r.next_f64()).collect(),
                gapfreq2: (0..len).map(|_| r.next_f64()).collect(),
                mrow: (0..len).map(|_| r.next_f64()).collect(),
                mprow: (0..len).map(|_| r.next_i32()).collect(),
                ijpr: (0..len).map(|_| r.next_i32()).collect(),
                curw: (0..len).map(|_| r.next_f64()).collect(),
            }
        }

        /// The scalar loop the kernel replaces, transcribed from
        /// `PartAlignState::part_a_align`. Runs the half-open range the kernel
        /// covered, so the two are comparable.
        fn scalar(
            &mut self,
            jj_end: usize,
            i: i32,
            gf1va: f64,
            fgcp1va: f64,
            ogcp1va: f64,
            fpenalty_ex: f64,
        ) -> f64 {
            let mut last = 0.0;
            for j in 1..jj_end {
                let mut w = self.prevw[j - 1];
                self.ijpr[j] = 0;

                let g = self.miv[j] + self.fgcp2[j - 1] * gf1va;
                if g > w {
                    w = g;
                    self.ijpr[j] = -(j as i32 - self.mpiv[j]);
                }

                let g = self.mrow[j] + fgcp1va * self.gapfreq2[j];
                if g > w {
                    w = g;
                    self.ijpr[j] = i - self.mprow[j];
                }

                let g = self.prevw[j - 1] + ogcp1va * self.gapfreq2[j - 1];
                if g > self.mrow[j] {
                    self.mrow[j] = g;
                    self.mprow[j] = i - 1;
                }
                self.mrow[j] += fpenalty_ex;

                self.curw[j] += w;
                last = w;
            }
            last
        }

        fn kernel(
            &mut self,
            lastj: usize,
            i: i32,
            gf1va: f64,
            fgcp1va: f64,
            ogcp1va: f64,
            fpenalty_ex: f64,
        ) -> Option<(usize, f64)> {
            row_checked(
                1,
                lastj,
                &self.prevw,
                &self.miv,
                &self.mpiv,
                &self.fgcp2,
                &self.gapfreq2,
                &mut self.mrow,
                &mut self.mprow,
                &mut self.ijpr,
                &mut self.curw,
                i,
                gf1va,
                fgcp1va,
                ogcp1va,
                fpenalty_ex,
            )
        }

        /// Bit patterns, not values: `-0.0 == 0.0` compares true, and telling
        /// those two apart is most of the point of this comparison.
        fn assert_same_as(&self, other: &Self, what: &str) {
            let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            assert_eq!(bits(&self.mrow), bits(&other.mrow), "{what}: mrow");
            assert_eq!(bits(&self.curw), bits(&other.curw), "{what}: curw");
            assert_eq!(self.mprow, other.mprow, "{what}: mprow");
            assert_eq!(self.ijpr, other.ijpr, "{what}: ijpr");
        }
    }

    impl Row {
        /// A row drawn from a handful of small integers, so the three
        /// comparisons tie constantly.
        ///
        /// Random `f64`s never collide, which makes them blind to the one
        /// mistake this kernel is most exposed to: `partA__align` uses `>`
        /// where `A__align` uses `>=`, upstream flipped it deliberately in
        /// 2018, and the difference is visible only on a tie at identical
        /// score. A mutation flipping the first comparison to `_CMP_GE_OQ`
        /// survives the random rows above and fails here.
        fn tied(len: usize, seed: u64) -> Self {
            let mut r = Lcg(seed);
            let small = |r: &mut Lcg, n: i32| f64::from(r.next_i32() % n);
            Self {
                prevw: (0..len).map(|_| small(&mut r, 3)).collect(),
                miv: (0..len).map(|_| small(&mut r, 3)).collect(),
                mpiv: (0..len).map(|_| r.next_i32() % 5).collect(),
                fgcp2: (0..len).map(|_| small(&mut r, 2)).collect(),
                gapfreq2: (0..len).map(|_| small(&mut r, 2)).collect(),
                mrow: (0..len).map(|_| small(&mut r, 3)).collect(),
                mprow: (0..len).map(|_| r.next_i32() % 5).collect(),
                ijpr: vec![0; len],
                curw: (0..len).map(|_| small(&mut r, 3)).collect(),
            }
        }
    }

    /// The same comparison, on rows where the values tie.
    ///
    /// `1.0` multipliers keep every product on the same small integer grid, so
    /// `g` lands exactly on `w` often rather than never.
    #[test]
    fn kernel_matches_the_scalar_row_on_ties() {
        if !crate::lalign_avx512::available() {
            return;
        }
        for lastj in 9..40usize {
            for seed in [2u64, 23, 2027] {
                let mut k = Row::tied(lastj, seed);
                let mut s = k.clone();
                let (jj, last_w) = k
                    .kernel(lastj, 2, 1.0, 1.0, 1.0, 0.0)
                    .expect("the kernel is available");
                let want = s.scalar(jj, 2, 1.0, 1.0, 1.0, 0.0);
                k.assert_same_as(&s, &format!("tied lastj={lastj} seed={seed}"));
                assert_eq!(last_w.to_bits(), want.to_bits(), "tied: the last cell's w");
            }
        }
    }

    /// The kernel and the loop it replaces agree bit for bit, over row lengths
    /// that exercise every remainder the block loop can stop on.
    #[test]
    fn kernel_matches_the_scalar_row() {
        if !crate::lalign_avx512::available() {
            return;
        }
        for lastj in 9..40usize {
            for seed in [1u64, 7, 99] {
                let mut k = Row::random(lastj, seed);
                let mut s = k.clone();
                let (jj, last_w) = k
                    .kernel(lastj, 5, 0.25, -1.5, 0.75, 0.0)
                    .expect("the kernel is available");
                assert!(jj > 1 && jj <= lastj, "lastj={lastj}: kernel stopped at {jj}");
                assert_eq!((jj - 1) % 8, 0, "lastj={lastj}: {jj} is not a whole block");
                let want = s.scalar(jj, 5, 0.25, -1.5, 0.75, 0.0);
                k.assert_same_as(&s, &format!("lastj={lastj} seed={seed}"));
                assert_eq!(
                    last_w.to_bits(),
                    want.to_bits(),
                    "lastj={lastj} seed={seed}: the last cell's w"
                );
            }
        }
    }

    /// A non-zero extension penalty is not the caller's ordinary case, but the
    /// kernel takes one and must still agree.
    #[test]
    fn kernel_matches_the_scalar_row_with_an_extension_penalty() {
        if !crate::lalign_avx512::available() {
            return;
        }
        let mut k = Row::random(25, 4242);
        let mut s = k.clone();
        let (jj, _) = k
            .kernel(25, 3, -0.5, 2.0, -0.125, -0.059)
            .expect("the kernel is available");
        s.scalar(jj, 3, -0.5, 2.0, -0.125, -0.059);
        k.assert_same_as(&s, "extension penalty");
    }

    /// A row too short for one block leaves everything to the caller, rather
    /// than reading past the ends of the buffers.
    #[test]
    fn a_row_shorter_than_a_block_does_nothing() {
        if !crate::lalign_avx512::available() {
            return;
        }
        let mut k = Row::random(8, 11);
        let before = k.clone();
        let (jj, _) = k.kernel(8, 1, 1.0, 1.0, 1.0, 0.0).expect("available");
        assert_eq!(jj, 1, "nothing should have been processed");
        k.assert_same_as(&before, "short row");
    }

    /// The bounds the kernel would otherwise read past are checked, not
    /// documented -- a release build must refuse rather than reach outside.
    #[test]
    #[should_panic(expected = "prevw[jj - 1]")]
    fn a_zero_start_column_is_refused() {
        if !crate::lalign_avx512::available() {
            panic!("the kernel reads `prevw[jj - 1]`");
        }
        let mut r = Row::random(16, 3);
        row_checked(
            0,
            16,
            &r.prevw,
            &r.miv,
            &r.mpiv,
            &r.fgcp2,
            &r.gapfreq2,
            &mut r.mrow,
            &mut r.mprow,
            &mut r.ijpr,
            &mut r.curw,
            1,
            0.0,
            0.0,
            0.0,
            0.0,
        );
    }

    /// The same, for a buffer that is simply too short for the row.
    #[test]
    #[should_panic(expected = "read and written up to")]
    fn a_short_buffer_is_refused() {
        if !crate::lalign_avx512::available() {
            panic!("every other buffer is read and written up to `lastj - 1`");
        }
        let mut r = Row::random(16, 3);
        r.curw.truncate(9);
        row_checked(
            1,
            16,
            &r.prevw,
            &r.miv,
            &r.mpiv,
            &r.fgcp2,
            &r.gapfreq2,
            &mut r.mrow,
            &mut r.mprow,
            &mut r.ijpr,
            &mut r.curw,
            1,
            0.0,
            0.0,
            0.0,
            0.0,
        );
    }
}
