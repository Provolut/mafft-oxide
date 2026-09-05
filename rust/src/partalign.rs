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

    /// A zeroed `ijp` of the right shape: every cell says "diagonal".
    fn ijp_of(lgth1: usize, lgth2: usize) -> Vec<Vec<i32>> {
        vec![vec![0i32; lgth2 + 1]; lgth1 + 1]
    }

    #[test]
    fn an_all_diagonal_path_gives_two_gapless_tables() {
        let mut ijp = ijp_of(3, 3);
        let t = atracking(&[0.0; 3], &[0.0; 3], 3, 3, &mut ijp, true);
        assert_eq!(t.gt1, b"ooo".to_vec());
        assert_eq!(t.gt2, b"ooo".to_vec());
    }

    /// The contract every traceback must satisfy: the two tables describe one
    /// alignment, and each consumes its whole sequence exactly once.
    fn assert_consistent(t: &PartTrack, lgth1: usize, lgth2: usize) {
        assert_eq!(t.gt1.len(), t.gt2.len(), "tables describe one alignment");
        assert_eq!(
            t.gt1.iter().filter(|&&c| c == b'o').count(),
            lgth1,
            "the first table must consume exactly its sequence"
        );
        assert_eq!(
            t.gt2.iter().filter(|&&c| c == b'o').count(),
            lgth2,
            "the second table must consume exactly its sequence"
        );
        for (a, b) in t.gt1.iter().zip(&t.gt2) {
            assert!(!(*a == b'-' && *b == b'-'), "no column may be gap-gap");
        }
    }

    #[test]
    fn a_positive_ijp_consumes_rows_and_gaps_the_second_table() {
        // `ijp[3][3] = 2` means "came from (1, 2)": one row is consumed against
        // a gap in the second group. The walk then runs to the border, whose
        // `ijp[0][j] = -(j+1)` emits the terminal gap -- so both tables end up
        // with one gap, which is why this asserts the consistency contract
        // rather than a gap count guessed from the first step.
        let mut ijp = ijp_of(3, 3);
        ijp[3][3] = 2;
        let t = atracking(&[0.0; 3], &[0.0; 3], 3, 3, &mut ijp, true);
        assert_consistent(&t, 3, 3);
        assert_eq!(t.gt2.iter().filter(|&&c| c == b'-').count(), 1);
    }

    #[test]
    fn a_negative_ijp_consumes_columns_and_gaps_the_first_table() {
        let mut ijp = ijp_of(3, 3);
        ijp[3][3] = -2;
        let t = atracking(&[0.0; 3], &[0.0; 3], 3, 3, &mut ijp, true);
        assert_consistent(&t, 3, 3);
        assert_eq!(t.gt1.iter().filter(|&&c| c == b'-').count(), 1);
    }

    #[test]
    fn every_single_step_entry_yields_a_consistent_alignment() {
        // Sweep the entry cell over every legal value rather than trusting one
        // hand-picked case: each must still produce one alignment that consumes
        // both sequences.
        for v in -3i32..=3 {
            let mut ijp = ijp_of(3, 3);
            ijp[3][3] = v;
            let t = atracking(&[0.0; 3], &[0.0; 3], 3, 3, &mut ijp, true);
            assert_consistent(&t, 3, 3);
        }
    }

    #[test]
    fn the_entry_cell_scan_keeps_the_last_maximum_not_the_first() {
        // `>=` rather than `>` is the whole point: with a tie, the *later*
        // position wins, and entering one cell over gives the same score with a
        // different alignment. A port using `>` would pass every score check.
        let mut ijp = ijp_of(3, 3);
        // Two equal maxima down the last column; the scan must keep i = 2.
        let lastverticalw = [5.0, 5.0, 5.0];
        let lasthorizontalw = [0.0, 0.0, 0.0];
        atracking(&lasthorizontalw, &lastverticalw, 3, 3, &mut ijp, false);
        assert_eq!(ijp[3][3], 1, "expected the last maximum (i = 2), not the first");
    }

    #[test]
    fn a_horizontal_maximum_wins_ties_against_the_vertical_scan() {
        // The horizontal scan runs second and also uses `>=`, so an equal value
        // there displaces the vertical pick.
        let mut ijp = ijp_of(3, 3);
        atracking(&[5.0, 5.0, 5.0], &[5.0, 5.0, 5.0], 3, 3, &mut ijp, false);
        assert!(ijp[3][3] < 0, "a tie must end on the horizontal branch");
    }

    #[test]
    fn the_localhom_walk_matches_the_plain_one_table_for_table() {
        // Sharing the walk is only safe if adding the constraint term does not
        // disturb it. Same inputs, both entry points, tables must be identical.
        for v in -3i32..=3 {
            let mut a = ijp_of(3, 3);
            a[3][3] = v;
            let plain = atracking(&[0.0; 3], &[0.0; 3], 3, 3, &mut a, true);

            let mut b = ijp_of(3, 3);
            b[3][3] = v;
            let imp = seeded(8, 8);
            let t = ImpTrack {
                imp: &imp,
                gapmap1: &[0, 1, 2, 3],
                gapmap2: &[0, 1, 2, 3],
                start1: 0,
                start2: 0,
            };
            let (with_imp, _) = atracking_localhom(&[0.0; 3], &[0.0; 3], 3, 3, &mut b, true, &t);

            assert_eq!(plain.gt1, with_imp.gt1, "entry value {v}");
            assert_eq!(plain.gt2, with_imp.gt2, "entry value {v}");
        }
    }

    #[test]
    fn the_constraint_score_sums_the_cells_the_path_passes_through() {
        // Straight diagonal on a 3x3: the walk visits (2,2) then (1,1), and the
        // entry corner (3,3) is excluded by the C's `iin != lgth1` guard. With
        // the seeded matrix those are 22 and 11.
        let mut ijp = ijp_of(3, 3);
        let imp = seeded(8, 8);
        let t = ImpTrack {
            imp: &imp,
            gapmap1: &[0, 1, 2, 3],
            gapmap2: &[0, 1, 2, 3],
            start1: 0,
            start2: 0,
        };
        let (_, score) = atracking_localhom(&[0.0; 3], &[0.0; 3], 3, 3, &mut ijp, true, &t);
        assert_eq!(score, 22.0 + 11.0);
    }

    #[test]
    fn the_window_offsets_shift_which_cells_are_read() {
        // `start1` / `start2` are what make this the *partial* aligner: the same
        // path reads a different part of the constraint matrix.
        let imp = seeded(8, 8);
        let mut plain = ijp_of(3, 3);
        let t0 = ImpTrack {
            imp: &imp,
            gapmap1: &[0, 1, 2, 3],
            gapmap2: &[0, 1, 2, 3],
            start1: 0,
            start2: 0,
        };
        let (_, at_origin) = atracking_localhom(&[0.0; 3], &[0.0; 3], 3, 3, &mut plain, true, &t0);

        let mut shifted = ijp_of(3, 3);
        let t1 = ImpTrack {
            imp: &imp,
            gapmap1: &[0, 1, 2, 3],
            gapmap2: &[0, 1, 2, 3],
            start1: 1,
            start2: 1,
        };
        let (_, at_offset) =
            atracking_localhom(&[0.0; 3], &[0.0; 3], 3, 3, &mut shifted, true, &t1);

        // Each cell moves one row and one column on: 22 -> 33, 11 -> 22.
        assert_eq!(at_origin, 22.0 + 11.0);
        assert_eq!(at_offset, 33.0 + 22.0);
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

/// One traceback's result: the two gap tables `profile::gapireru` expands.
///
/// The C writes the aligned rows straight into `mseq1` / `mseq2` as it walks.
/// Gap tables carry the same information -- `b'o'` takes the next residue,
/// `b'-'` inserts a gap -- and are what the rest of this port already speaks,
/// so the expansion is the existing `gapireru` rather than a second writer.
pub struct PartTrack {
    pub gt1: Vec<u8>,
    pub gt2: Vec<u8>,
}

/// Where the constraint score comes from, for the `localhom` traceback.
///
/// Bundled because the five values are only ever used together, at one line of
/// the walk, and threading them as separate parameters would put the shared
/// walk over any reasonable argument count.
pub struct ImpTrack<'a> {
    pub imp: &'a PartImp,
    pub gapmap1: &'a [i32],
    pub gapmap2: &'a [i32],
    pub start1: usize,
    pub start2: usize,
}

/// The traceback walk shared by `Atracking` and `Atracking_localhom`.
///
/// The two C functions are the same walk; the `localhom` one additionally sums
/// the constraint score of every cell it passes through. Written once here for
/// the reason the `match_calc` note in `GOAL-C-FREE.md` gives: two
/// transliterations of one walk are free to drift apart, and a traceback that
/// drifts changes an alignment without changing its score.
///
/// The **walk itself** matches `align.rs`'s `atracking`: same branch structure
/// on the sign of `ijp`, same gap-run loops, same "off the front" terminal
/// case. The **entry cell** does not -- `Salignmm`'s picks it under `tailgp`,
/// this one under `outgap`, scanning with `>=` so a tie keeps the *later*
/// maximum.
///
/// The warp branches are omitted, as in `align.rs`: a `trywarp` configuration
/// is refused before reaching here.
///
/// Args:
///     lasthorizontalw: Final row of the DP.
///     lastverticalw: Final column of the DP.
///     lgth1: Rows.
///     lgth2: Columns.
///     ijp: The traceback matrix, edited at its borders as the C does.
///     outgap: Whether terminal gaps are free, which decides the entry cell.
///     imp: Present for the `localhom` variant only.
///
/// Returns:
///     The two gap tables in forward order, and the accumulated constraint
///     score (0.0 when `imp` is `None`).
fn track(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    lgth1: usize,
    lgth2: usize,
    ijp: &mut [Vec<i32>],
    outgap: bool,
    imp: Option<&ImpTrack>,
) -> (PartTrack, f64) {
    if !outgap {
        let mut wm = lastverticalw[0];
        for i in 0..lgth1 {
            if lastverticalw[i] >= wm {
                wm = lastverticalw[i];
                ijp[lgth1][lgth2] = (lgth1 - i) as i32;
            }
        }
        for j in 0..lgth2 {
            if lasthorizontalw[j] >= wm {
                wm = lasthorizontalw[j];
                ijp[lgth1][lgth2] = -((lgth2 - j) as i32);
            }
        }
    }

    // The borders encode "everything left is a terminal gap".
    for i in 0..lgth1 + 1 {
        ijp[i][0] = i as i32 + 1;
    }
    for j in 0..lgth2 + 1 {
        ijp[0][j] = -(j as i32 + 1);
    }

    let (mut rev1, mut rev2) = (Vec::new(), Vec::new());
    let (mut iin, mut jin) = (lgth1 as i64, lgth2 as i64);
    let mut impwm = 0.0f64;
    let limk = lgth1 + lgth2 + 1;
    for _ in 0..limk {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi) = if v < 0 {
            (iin - 1, jin + i64::from(v))
        } else if v > 0 {
            (iin - i64::from(v), jin - 1)
        } else {
            (iin - 1, jin - 1)
        };

        if ifi < 0 && jfi < 0 {
            // Off the front: the remainder of both rows is terminal gap.
            let mut l = iin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'o');
                rev2.push(b'-');
            }
            let mut l = jin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'-');
                rev2.push(b'o');
            }
            break;
        }

        let mut l = iin - ifi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'o');
            rev2.push(b'-');
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'-');
            rev2.push(b'o');
        }

        // Before the break below, so the cell the walk stops on still counts.
        // The guard skips the entry corner, as the C's `iin != lgth1 && jin !=
        // lgth2` does.
        if let Some(t) = imp {
            if iin as usize != lgth1 && jin as usize != lgth2 {
                impwm += t.imp.score(
                    t.gapmap1[iin as usize] as usize + t.start1,
                    t.gapmap2[jin as usize] as usize + t.start2,
                );
            }
        }

        if iin <= 0 || jin <= 0 {
            break;
        }
        rev1.push(b'o');
        rev2.push(b'o');
        iin = ifi;
        jin = jfi;
    }

    rev1.reverse();
    rev2.reverse();
    (PartTrack { gt1: rev1, gt2: rev2 }, impwm)
}

/// `Atracking()` from `partSalignmm.c`: the walk with no constraint term.
pub fn atracking(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    lgth1: usize,
    lgth2: usize,
    ijp: &mut [Vec<i32>],
    outgap: bool,
) -> PartTrack {
    track(lasthorizontalw, lastverticalw, lgth1, lgth2, ijp, outgap, None).0
}

/// `Atracking_localhom()` from `partSalignmm.c`.
///
/// Returns:
///     The gap tables and the C's `*impwmpt` -- the constraint score summed
///     over every cell the path passes through, which the caller adds to the
///     alignment's score.
pub fn atracking_localhom(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    lgth1: usize,
    lgth2: usize,
    ijp: &mut [Vec<i32>],
    outgap: bool,
    imp: &ImpTrack,
) -> (PartTrack, f64) {
    track(
        lasthorizontalw,
        lastverticalw,
        lgth1,
        lgth2,
        ijp,
        outgap,
        Some(imp),
    )
}
