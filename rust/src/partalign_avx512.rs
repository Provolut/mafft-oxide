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
