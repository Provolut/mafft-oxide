//! The partial group-to-group aligner (`partSalignmm.c`).
//!
//! `dvtditr`'s iterative refinement pulls a group of rows out of the finished
//! alignment and re-aligns it against the rest — but only across a *window* of
//! the alignment, not end to end. That is what makes this a separate aligner
//! from `align.rs`'s `A__align` rather than a call to it: every position it
//! touches is expressed as an offset into a window (`start1..end1`,
//! `start2..end2`) and mapped back through a gap map.
//!
//! It is also *simpler* than `A__align` in one respect worth knowing before
//! reading either: it carries no `cpmxhist` profile cache. That cache is what
//! makes `A__align`'s summation order subtle, and its absence here means the
//! profile is rebuilt per call, with nothing to keep consistent across calls.
//!
//! This module currently holds the constraint ("localhom") half. The DP itself
//! is not ported yet, so `dvtditr` still runs the C.

use crate::localhom::{fillimp, LocalHom};

/// The constraint matrix `partA__align`'s DP adds to its match scores.
///
/// A struct rather than the C's file-scope `static TLS double **impmtx`: the C
/// keeps it between calls only to avoid reallocating, and passes `seq1 == NULL`
/// to free it. An owned value expresses the same lifetime without the sentinel.
pub struct PartImp {
    mtx: Vec<Vec<f64>>,
    alloclen: usize,
}

impl Default for PartImp {
    fn default() -> Self {
        Self::new()
    }
}

impl PartImp {
    /// An empty matrix; the first `init_strict` sizes it.
    pub fn new() -> PartImp {
        PartImp {
            mtx: Vec::new(),
            alloclen: 0,
        }
    }

    /// `part_imp_match_init_strict()`: size the matrix and fill it.
    ///
    /// The growth rule is the C's -- reallocate to `max(lgth1, lgth2) + 2`
    /// whenever either dimension no longer fits, and keep it otherwise, so a
    /// run of calls on shrinking windows reuses one allocation.
    ///
    /// Args:
    ///     clus1: Members in the first group.
    ///     clus2: Members in the second group.
    ///     lgth1: Window width on the first group.
    ///     lgth2: Window width on the second group.
    ///     seq1: The first group's rows.
    ///     seq2: The second group's rows.
    ///     eff1: The first group's member weights.
    ///     eff2: The second group's member weights.
    ///     shrink: Per-member-pair constraint lists.
    ///     orinum1: Original indices of the first group's members.
    ///     orinum2: Original indices of the second group's members.
    ///     fastathreshold: The `effijx` weight `fillimp` applies.
    #[allow(clippy::too_many_arguments)]
    pub fn init_strict(
        &mut self,
        clus1: usize,
        clus2: usize,
        lgth1: usize,
        lgth2: usize,
        seq1: &[&[u8]],
        seq2: &[&[u8]],
        eff1: &[f64],
        eff2: &[f64],
        shrink: &[Vec<Option<&Vec<LocalHom>>>],
        orinum1: &[i32],
        orinum2: &[i32],
        fastathreshold: f64,
    ) {
        if self.alloclen < lgth1 + 2 || self.alloclen < lgth2 + 2 {
            self.alloclen = lgth1.max(lgth2) + 2;
            self.mtx = vec![vec![0.0; self.alloclen]; self.alloclen];
        }
        fillimp(
            &mut self.mtx,
            clus1,
            clus2,
            lgth1,
            lgth2,
            seq1,
            seq2,
            eff1,
            eff2,
            shrink,
            orinum1,
            orinum2,
            fastathreshold,
        );
    }

    /// `part_imp_match_out_sc()`: one cell, as the traceback reads it.
    pub fn score(&self, i1: usize, j1: usize) -> f64 {
        self.mtx[i1][j1]
    }

    /// `part_imp_match_out_vead_gapmap()`: add row `i1` of the constraint
    /// matrix onto `imp`, walking the second group's gap map.
    ///
    /// The `#if FASTMACHCALC` branch beside this one in the C is unreachable --
    /// that macro is misspelled and defined nowhere (upstream issue #11) -- so
    /// the accumulating `#else` is what the shipped binary runs, and what this
    /// reproduces.
    pub fn add_row(&self, imp: &mut [f64], i1: usize, start2: usize, gapmap2: &[i32]) {
        for (dst, &g) in imp.iter_mut().zip(gapmap2) {
            *dst += self.mtx[i1][start2 + g as usize];
        }
    }

    /// `part_imp_match_out_vead_tate_gapmap()`: the same down a column.
    ///
    /// **This one accumulates**, and its unreachable twin assigns. Fixing the
    /// misspelled macro without also fixing that `=` would change alignments;
    /// see upstream issue #11. The accumulating form is the live one.
    pub fn add_column(&self, imp: &mut [f64], j1: usize, start1: usize, gapmap1: &[i32]) {
        for (dst, &g) in imp.iter_mut().zip(gapmap1) {
            *dst += self.mtx[start1 + g as usize][j1];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A `PartImp` with a known matrix, bypassing `fillimp` so the accessors
    /// can be checked against hand-computed values.
    fn seeded(rows: usize, cols: usize) -> PartImp {
        let mut p = PartImp::new();
        p.alloclen = rows.max(cols);
        p.mtx = (0..rows)
            .map(|i| (0..cols).map(|j| (i * 10 + j) as f64).collect())
            .collect();
        p
    }

    #[test]
    fn score_reads_the_cell_it_is_asked_for() {
        let p = seeded(4, 4);
        assert_eq!(p.score(2, 3), 23.0);
    }

    #[test]
    fn add_row_walks_the_gap_map_and_accumulates() {
        let p = seeded(4, 8);
        let mut imp = vec![1.0, 1.0, 1.0];
        // Columns 0, 2, 3 of the window starting at 1 -> matrix columns 1, 3, 4.
        p.add_row(&mut imp, 2, 1, &[0, 2, 3]);
        assert_eq!(imp, vec![1.0 + 21.0, 1.0 + 23.0, 1.0 + 24.0]);
    }

    #[test]
    fn add_column_accumulates_rather_than_assigning() {
        // The distinction upstream issue #11 is about: the live branch adds to
        // whatever `imp` already holds. A port that assigned would pass a test
        // starting from zeros, so this one deliberately does not.
        let p = seeded(8, 4);
        let mut imp = vec![100.0, 200.0];
        p.add_column(&mut imp, 3, 1, &[0, 2]);
        assert_eq!(imp, vec![100.0 + 13.0, 200.0 + 33.0]);
    }

    #[test]
    fn the_matrix_grows_only_when_a_window_no_longer_fits() {
        let mut p = PartImp::new();
        p.alloclen = 10;
        p.mtx = vec![vec![7.0; 10]; 10];
        // Fits: 8 + 2 == 10, so the existing allocation is kept and its
        // contents with it.
        assert!(!(p.alloclen < 8 + 2));
        // Does not fit: 9 + 2 > 10.
        assert!(p.alloclen < 9 + 2);
    }
}
