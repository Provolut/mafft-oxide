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
//! # STATUS: NOT BIT-VERIFIED. Treat everything here as a draft.
//!
//! This module holds `partA__align`, both of `partSalignmm.c`'s tracebacks, and
//! the constraint ("localhom") half. All of it compiles and the tests below
//! pass — but **nothing has ever executed it against the C**, because
//! `dvtditr-rs` does not exist yet and no other binary calls it. It is
//! unreachable from every shipped binary, so it cannot currently affect any
//! output; it is groundwork for `dvtditr`.
//!
//! What the tests here do and do not establish:
//!
//! * They **do** show the DP terminates, preserves every residue, emits
//!   equal-width rows, and gives the same answer when its state is reused
//!   (which would catch stale scratch buffers).
//! * They **do not** show bit-exactness, and no unit test can. A single
//!   reordered addition or a flipped tie comparison changes an alignment while
//!   leaving its score identical, and nothing in `cargo test` would notice.
//!
//! What will decide it, in order: `dvtditr-rs` calling this so it runs at all;
//! then `verify_tbfast.sh` and `verify_e2e.sh` against a built `core/`; then
//! the ProteinGym corpus through `mafft --auto`, byte-compared; then
//! provolut-ml's ESCOTT scorer, SHA-256 compared. `README.md`'s "Verification"
//! has the
//! commands.
//!
//! **The likeliest place for this to be wrong** is the pair of gap-open tie
//! comparisons in `part_a_align`, which use `>` where `align.rs`'s `A__align`
//! uses `>=` — deliberately, because the two C files disagree there. See the
//! note on that function, and `README.md`'s "Hazards in the C".

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
            // Freed before the replacement is built: an assignment evaluates its
            // right-hand side first, so the old buffer and the new one are both
            // live for the duration -- and a growth happens at the moment the
            // working set is largest, so the overlap lands on the peak. Worth
            // 2.8 MB where `commonip.rs` had the same shape.
            self.mtx = Vec::new();
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
    /// that macro is misspelled and defined nowhere (this repo's issue #11) -- so
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
    /// see this repo's issue #11. The accumulating form is the live one.
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
        // The distinction this repo's issue #11 is about: the live branch adds to
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

        // Each cell moves one row and one column on: 22 -> 33, 11 -> 22, and
        // the border cell 00 -> 11. The walk ends *on* (0, 0) and counts it,
        // as the C does -- at the origin it contributes 0.0 and is invisible,
        // shifted it contributes 11.0. Confirmed against the C's own
        // `Atracking_localhom` (verify_flh.sh, 98 of 98 cases bit-identical).
        assert_eq!(at_origin, 22.0 + 11.0 + 0.0);
        assert_eq!(at_offset, 33.0 + 22.0 + 11.0);
    }

    /// The progressive stage's parameters, matching `tests_unit::tree_constants`.
    fn dp_constants() -> crate::constants::Constants {
        let a = crate::params::Args {
            nblosum: 62,
            ppenalty: -1530,
            ppenalty_ex: 0,
            poffset: 0,
            penalty_shift_factor: 100.0,
            ..Default::default()
        };
        crate::constants::Constants::new(&a)
    }

    /// Run the DP on two groups and return the aligned rows.
    fn run_dp(g1: &[&str], g2: &[&str]) -> (Vec<Vec<u8>>, Vec<Vec<u8>>, f64) {
        let cn = dp_constants();
        let mut group1: Vec<Vec<u8>> = g1.iter().map(|s| s.as_bytes().to_vec()).collect();
        let mut group2: Vec<Vec<u8>> = g2.iter().map(|s| s.as_bytes().to_vec()).collect();
        let eff1 = vec![1.0 / group1.len() as f64; group1.len()];
        let eff2 = vec![1.0 / group2.len() as f64; group2.len()];
        let mut st = PartAlignState::default();
        let (score, imp) = st.part_a_align(
            &cn,
            cn.penalty,
            cn.penalty_ex,
            &mut group1,
            &mut group2,
            &eff1,
            &eff2,
            true,
            false,
            None,
            None,
            b'-',
            usize::MAX,
        );
        assert_eq!(imp, 0.0, "unconstrained runs carry no constraint score");
        (group1, group2, score)
    }

    /// Every row of an alignment is the same width, and dropping its gaps gives
    /// back the input.
    fn assert_is_alignment(out: &[Vec<u8>], input: &[&str]) {
        let w = out[0].len();
        for (row, orig) in out.iter().zip(input) {
            assert_eq!(row.len(), w, "ragged alignment");
            let ungapped: Vec<u8> = row.iter().copied().filter(|&c| c != b'-').collect();
            assert_eq!(
                ungapped,
                orig.as_bytes().to_vec(),
                "the DP altered residues, not just gaps"
            );
        }
    }

    #[test]
    fn the_dp_runs_and_returns_a_well_formed_alignment() {
        // Not a bit-exactness check -- only `verify_tbfast.sh` against a built
        // `core/` can be that. This is the weaker claim that the DP terminates,
        // preserves every residue, and emits equal-width rows, which is what
        // separates "compiles" from "runs at all".
        let g1 = ["MKVLTA", "MKVITA"];
        let g2 = ["MKVLTA"];
        let (a1, a2, _) = run_dp(&g1, &g2);
        assert_is_alignment(&a1, &g1);
        assert_is_alignment(&a2, &g2);
        assert_eq!(a1[0].len(), a2[0].len(), "both groups share one width");
    }

    #[test]
    fn the_dp_aligns_sequences_of_different_lengths() {
        let g1 = ["MKVLTAAG"];
        let g2 = ["MKVTA"];
        let (a1, a2, _) = run_dp(&g1, &g2);
        assert_is_alignment(&a1, &g1);
        assert_is_alignment(&a2, &g2);
        assert_eq!(a1[0].len(), a2[0].len());
        assert!(a1[0].len() >= 8, "the alignment cannot be shorter than its longer input");
    }

    #[test]
    fn an_empty_group_becomes_all_gap() {
        // The C's degenerate branches: one side empty, the other all gap.
        let cn = dp_constants();
        let mut group1: Vec<Vec<u8>> = vec![b"".to_vec()];
        let mut group2: Vec<Vec<u8>> = vec![b"MKVL".to_vec()];
        let mut st = PartAlignState::default();
        let (score, _) = st.part_a_align(
            &cn, cn.penalty, cn.penalty_ex, &mut group1, &mut group2,
            &[1.0], &[1.0], true, false, None, None, b'-', usize::MAX,
        );
        assert_eq!(score, 0.0);
        assert_eq!(group1[0], b"----".to_vec());
    }

    #[test]
    fn the_scratch_survives_a_second_call() {
        // `ensure` keeps its buffers when a later window fits, so a stale value
        // left in them would show up as a different answer on the second run.
        let cn = dp_constants();
        let mut st = PartAlignState::default();
        let mut first = None;
        for _ in 0..2 {
            let mut g1: Vec<Vec<u8>> = vec![b"MKVLTA".to_vec()];
            let mut g2: Vec<Vec<u8>> = vec![b"MKVITA".to_vec()];
            let (score, _) = st.part_a_align(
                &cn, cn.penalty, cn.penalty_ex, &mut g1, &mut g2,
                &[1.0], &[1.0], true, false, None, None, b'-', usize::MAX,
            );
            match first {
                None => first = Some((score, g1, g2)),
                Some((s0, ref r1, ref r2)) => {
                    assert_eq!(score, s0, "a reused state changed the score");
                    assert_eq!(&g1, r1);
                    assert_eq!(&g2, r2);
                }
            }
        }
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
/// the reason the `match_calc` note in `README.md` gives: two
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

        // No "off the front" branch here. The C has one, but it fires on
        // `ifi == -warpbase && jfi == -warpbase` -- a warp entry, and warp is
        // refused before reaching this port -- not on any negative pair. The
        // border cells reach it with `ifi == -1, jfi == -1`, and the ordinary
        // path below already handles those: both gap-run loops run zero or
        // `jin` times exactly as the C's do, and the cell still contributes its
        // constraint score. Taking a shortcut branch there instead drops that
        // last `part_imp_match_out_sc` from `*impwmpt`, which leaves the rows
        // and the DP score identical and only the constraint score wrong --
        // invisible to everything except a differential on `impmatch`, which is
        // how it was found (verify_flh.sh, 3 of 98 cases).

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

/// Scratch buffers for [`part_a_align`], reused across calls.
///
/// The C keeps these as `static TLS` and grows them on demand; the growth rule
/// is the contract, the storage class is not. No `cpmxhist` here -- see the
/// module docs.
#[derive(Default)]
pub struct PartAlignState {
    orlgth1: usize,
    orlgth2: usize,
    w1: Vec<f64>,
    w2: Vec<f64>,
    initverticalw: Vec<f64>,
    lastverticalw: Vec<f64>,
    m: Vec<f64>,
    mp: Vec<i32>,
    /// The `mi`/`mpi` recurrence, computed as a prefix scan ahead of the row
    /// that consumes it. The C runs it inline; see the note in
    /// [`PartAlignState::part_a_align`].
    miv: Vec<f64>,
    mpiv: Vec<i32>,
    /// The scan's candidates, written as a row so `l_align11`'s scan kernel
    /// can consume them unchanged.
    cand: Vec<f64>,
    ogcp1: Vec<f64>,
    ogcp2: Vec<f64>,
    fgcp1: Vec<f64>,
    fgcp2: Vec<f64>,
    cpmx1: Vec<Vec<f64>>,
    cpmx2: Vec<Vec<f64>>,
    gapfreq1: Vec<f64>,
    gapfreq2: Vec<f64>,
    doublework: Vec<Vec<f64>>,
    intwork: Vec<Vec<i32>>,
    ijp: Vec<Vec<i32>>,
}

/// Everything the constraint path needs, absent when running unconstrained.
pub struct PartConstraint<'a> {
    pub imp: &'a PartImp,
    pub gapmap1: &'a [i32],
    pub gapmap2: &'a [i32],
    pub start1: usize,
    pub start2: usize,
}

/// The four terminal-gap patterns `new_*GapCount` needs, when the caller has
/// them. `None` selects the C's `st_*` path and zero head/tail frequencies.
pub type PartGaps<'a> = (&'a [u8], &'a [u8], &'a [u8], &'a [u8]);

impl PartAlignState {
    /// `partA__align()` from `partSalignmm.c`.
    ///
    /// Aligns two groups across the window the caller has already cut, writing
    /// the aligned rows back into `group1` / `group2`.
    ///
    /// **The two gap-open comparisons use `>`, not `>=`.** `Salignmm.c`'s
    /// `A__align` uses `>=` there and this one does not; upstream flipped them
    /// here in April 2018 and left the other file alone. Copying `align.rs`'s
    /// loop would change the alignment on every tie at identical score. See
    /// `README.md`, "Hazards in the C".
    ///
    /// The `trywarp` and `chudanpt` paths are absent: this port refuses the
    /// former before reaching here, and nothing in this repo sets the latter.
    ///
    /// Returns:
    ///     The alignment score, and the accumulated constraint score (the C's
    ///     `*impmatch`, 0.0 when unconstrained).
    ///
    /// `alloclen` is the caller's row capacity; the C aborts when the result
    /// outgrows it.
    #[allow(clippy::too_many_arguments)]
    pub fn part_a_align(
        &mut self,
        cn: &crate::constants::Constants,
        penalty: i32,
        penalty_ex: i32,
        group1: &mut [Vec<u8>],
        group2: &mut [Vec<u8>],
        eff1: &[f64],
        eff2: &[f64],
        outgap: bool,
        legacygapcost: bool,
        gaps: Option<PartGaps>,
        constraint: Option<&PartConstraint>,
        newgapstr: u8,
        alloclen: usize,
    ) -> (f64, f64) {
        use crate::profile;

        let (icyc, jcyc) = (group1.len(), group2.len());
        let lgth1 = group1[0].len();
        let lgth2 = group2[0].len();

        // Degenerate windows: one side is empty, so the other is all gap.
        if lgth1 == 0 && lgth2 == 0 {
            return (0.0, 0.0);
        }
        if lgth1 == 0 {
            for row in group1.iter_mut() {
                *row = vec![newgapstr; lgth2];
            }
            return (0.0, 0.0);
        }
        if lgth2 == 0 {
            for row in group2.iter_mut() {
                *row = vec![newgapstr; lgth1];
            }
            return (0.0, 0.0);
        }

        self.ensure(lgth1, lgth2, cn.charsize);

        let fpenalty = f64::from(penalty);
        let fpenalty_ex = f64::from(penalty_ex);
        let s1: Vec<&[u8]> = group1.iter().map(std::vec::Vec::as_slice).collect();
        let s2: Vec<&[u8]> = group2.iter().map(std::vec::Vec::as_slice).collect();

        profile::cpmx_calc_new(cn, &s1, &mut self.cpmx1, eff1, lgth1, icyc);
        profile::cpmx_calc_new(cn, &s2, &mut self.cpmx2, eff2, lgth2, jcyc);

        let (headgapfreq1, headgapfreq2);
        if let Some((sgap1, sgap2, egap1, egap2)) = gaps {
            profile::new_opening_gap_count(&mut self.ogcp1, icyc, &s1, eff1, lgth1, sgap1);
            profile::new_opening_gap_count(&mut self.ogcp2, jcyc, &s2, eff2, lgth2, sgap2);
            profile::new_final_gap_count(&mut self.fgcp1, icyc, &s1, eff1, lgth1);
            // `new_FinalGapCount` takes an `egappat` the C never reads (see
            // `profile::new_final_gap_count`); `egap1` / `egap2` earn their
            // keep just below, in the tail gap frequencies.
            profile::new_final_gap_count(&mut self.fgcp2, jcyc, &s2, eff2, lgth2);
            headgapfreq1 = profile::outgapcount(icyc, sgap1, eff1);
            headgapfreq2 = profile::outgapcount(jcyc, sgap2, eff2);
            self.gapfreq1[lgth1] = profile::outgapcount(icyc, egap1, eff1);
            self.gapfreq2[lgth2] = profile::outgapcount(jcyc, egap2, eff2);
        } else {
            profile::st_opening_gap_count(&mut self.ogcp1, icyc, &s1, eff1, lgth1);
            profile::st_opening_gap_count(&mut self.ogcp2, jcyc, &s2, eff2, lgth2);
            profile::st_final_gap_count(&mut self.fgcp1, icyc, &s1, eff1, lgth1);
            profile::st_final_gap_count(&mut self.fgcp2, jcyc, &s2, eff2, lgth2);
            headgapfreq1 = 0.0;
            headgapfreq2 = 0.0;
            self.gapfreq1[lgth1] = 0.0;
            self.gapfreq2[lgth2] = 0.0;
        }

        let (headgapfreq1, headgapfreq2) = if legacygapcost {
            for v in self.gapfreq1.iter_mut().take(lgth1 + 1) {
                *v = 1.0;
            }
            for v in self.gapfreq2.iter_mut().take(lgth2 + 1) {
                *v = 1.0;
            }
            (1.0, 1.0)
        } else {
            profile::gapcountf(&mut self.gapfreq1, &s1, icyc, eff1, lgth1);
            profile::gapcountf(&mut self.gapfreq2, &s2, jcyc, eff2, lgth2);
            for v in self.gapfreq1.iter_mut().take(lgth1 + 1) {
                *v = 1.0 - *v;
            }
            for v in self.gapfreq2.iter_mut().take(lgth2 + 1) {
                *v = 1.0 - *v;
            }
            (1.0 - headgapfreq1, 1.0 - headgapfreq2)
        };

        for i in 0..lgth1 {
            self.ogcp1[i] = 0.5 * (1.0 - self.ogcp1[i]) * fpenalty * self.gapfreq1[i];
            self.fgcp1[i] = 0.5 * (1.0 - self.fgcp1[i]) * fpenalty * self.gapfreq1[i];
        }
        for i in 0..lgth2 {
            self.ogcp2[i] = 0.5 * (1.0 - self.ogcp2[i]) * fpenalty * self.gapfreq2[i];
            self.fgcp2[i] = 0.5 * (1.0 - self.fgcp2[i]) * fpenalty * self.gapfreq2[i];
        }

        let mut cur_is_w1 = true;
        {
            let (initv, curr) = (&mut self.initverticalw, &mut self.w1);
            crate::align::match_calc(
                &cn.n_dis_consweight_multi,
                initv,
                &self.cpmx2,
                &self.cpmx1,
                0,
                lgth1,
                &mut self.doublework,
                &mut self.intwork,
                true,
            );
            if let Some(c) = constraint {
                c.imp
                    .add_column(&mut initv[..lgth1], c.gapmap2[0] as usize + c.start2, c.start1, c.gapmap1);
            }
            crate::align::match_calc(
                &cn.n_dis_consweight_multi,
                curr,
                &self.cpmx1,
                &self.cpmx2,
                0,
                lgth2,
                &mut self.doublework,
                &mut self.intwork,
                true,
            );
            if let Some(c) = constraint {
                c.imp
                    .add_row(&mut curr[..lgth2], c.gapmap1[0] as usize + c.start1, c.start2, c.gapmap2);
            }
        }

        if outgap {
            for i in 1..lgth1 + 1 {
                self.initverticalw[i] +=
                    self.ogcp1[0] * headgapfreq2 + self.fgcp1[i - 1] * self.gapfreq2[0];
            }
            for j in 1..lgth2 + 1 {
                self.w1[j] += self.ogcp2[0] * headgapfreq1 + self.fgcp2[j - 1] * self.gapfreq1[0];
            }
        }

        for j in 1..lgth2 + 1 {
            self.m[j] = self.w1[j - 1] + self.ogcp1[1] * self.gapfreq2[j - 1];
            self.mp[j] = 0;
        }
        self.lastverticalw[0] = self.w1[lgth2 - 1];

        let lasti = if outgap { lgth1 + 1 } else { lgth1 };
        let lastj = lgth2 + 1;
        let mut wm = 0.0f64;

        for i in 1..lasti {
            cur_is_w1 = !cur_is_w1;
            let (currentw, previousw) = if cur_is_w1 {
                (&mut self.w1, &mut self.w2)
            } else {
                (&mut self.w2, &mut self.w1)
            };
            previousw[0] = self.initverticalw[i - 1];

            crate::align::match_calc(
                &cn.n_dis_consweight_multi,
                currentw,
                &self.cpmx1,
                &self.cpmx2,
                i,
                lgth2,
                &mut self.doublework,
                &mut self.intwork,
                false,
            );
            if let Some(c) = constraint {
                c.imp.add_row(
                    &mut currentw[..lgth2],
                    c.gapmap1[i] as usize + c.start1,
                    c.start2,
                    c.gapmap2,
                );
            }
            currentw[0] = self.initverticalw[i];

            let mut mi = previousw[0] + self.ogcp2[1] * self.gapfreq1[i - 1];
            let mut mpi = 0i32;
            let fgcp1va = self.fgcp1[i - 1];
            let ogcp1va = self.ogcp1[i];
            let gf1va = self.gapfreq1[i];
            let gf1vapre = self.gapfreq1[i - 1];
            let ijprow = &mut self.ijp[i];

            // Every row below is `&self.<field>[..lastj]` before the loop, not
            // `self.<field>[j]` inside it. The values and their order are
            // untouched; what changes is what the compiler can prove. Indexing
            // a `Vec` field while the same loop writes through `&mut self`
            // leaves the length a memory read the optimiser must reload, so
            // each of the ten accesses per cell keeps its bounds check. Bound
            // to local slices first, the lengths are loop-invariant SSA values,
            // every index is provably in range, and the checks fold away.
            // Measured on the L-INS-i replay: this is the single largest
            // instruction saving in the port.
            let prevw = &previousw[..lastj];
            let fgcp2 = &self.fgcp2[..lastj];
            let ogcp2 = &self.ogcp2[..lastj];
            let gapfreq2 = &self.gapfreq2[..lastj];
            let mrow = &mut self.m[..lastj];
            let mprow = &mut self.mp[..lastj];
            let miv = &mut self.miv[..lastj];
            let mpiv = &mut self.mpiv[..lastj];
            let cand = &mut self.cand[..lastj];
            let ijpr = &mut ijprow[..lastj];
            let curw = &mut currentw[..lastj];

            // The `mi` recurrence, lifted out of the row it used to run
            // inline with. The C keeps them in one loop, and every value and
            // comparison here is the C's -- what moves is only when they are
            // computed. `mi` at step `j` depends on nothing the row writes
            // (it reads `previousw` and `ogcp2` alone), and the row reads it
            // before the step that updates it, so the two are separable.
            //
            // Separating them is what makes the row element-wise in `j`, and
            // an element-wise row is the shape a kernel can take eight cells
            // of at a time.
            // The candidates the scan chooses between, written as a row.
            // `l_align11`'s scan takes its candidates that way, so spelling
            // them out lets that kernel serve this recurrence unchanged
            // instead of the port growing a second copy of it.
            for (c, (&pw, &og)) in cand[..lastj - 1]
                .iter_mut()
                .zip(prevw[..lastj - 1].iter().zip(&ogcp2[1..lastj]))
            {
                *c = pw + og * gf1vapre;
            }

            // The recurrence is `mi = max(mi, cand[j - 1]) + fpenalty_ex`,
            // which is the max-plus scan `l_align11` runs. Its vector form
            // reaches a position with one multiply where the loop adds
            // `fpenalty_ex` to itself repeatedly, and the two agree only when
            // nothing rounds -- so the kernel is asked for only when the
            // penalty is zero, where the multiply is exact for every `f64` and
            // no property of the candidates is needed. `partSalignmm.c`
            // compiles the extension penalty out entirely (`USE_PENALTY_EX 0`)
            // and the `mafft` script never passes `-g` to `dvtditr`, so that
            // is every call this repository makes; anything else runs the
            // serial loop below for the whole row.
            let mut scanned = 0usize;
            if let Some((done, carry_m, carry_mp)) = crate::lalign_avx512::prefix_scan_checked(
                fpenalty_ex == 0.0,
                &cand[..lastj - 1],
                lastj - 1,
                mi,
                mpi,
                fpenalty_ex,
                miv,
                mpiv,
            ) {
                scanned = done;
                mi = carry_m;
                mpi = carry_mp;
            }
            for j in 1 + scanned..lastj {
                miv[j] = mi;
                mpiv[j] = mpi;
                // `>`, not `>=`: see the note on this function.
                if cand[j - 1] > mi {
                    mi = cand[j - 1];
                    mpi = j as i32 - 1;
                }
                mi += fpenalty_ex;
            }

            // The vectorised part of the row. It leaves `jj` wherever it
            // stopped and the scalar loop below finishes from there -- and
            // returns the last cell's `w`, because a row whose length divides
            // by eight leaves that loop nothing to run and `wm` is the score
            // this function returns.
            let miv = &*miv;
            let mpiv = &*mpiv;
            let mut jj_start = 1usize;
            if let Some((jj, last_w)) = crate::partalign_avx512::row_checked(
                1,
                lastj,
                prevw,
                miv,
                mpiv,
                fgcp2,
                gapfreq2,
                mrow,
                mprow,
                ijpr,
                curw,
                i as i32,
                gf1va,
                fgcp1va,
                ogcp1va,
                fpenalty_ex,
            ) {
                if jj > 1 {
                    wm = last_w;
                }
                jj_start = jj;
            }

            for j in jj_start..lastj {
                let mut w = prevw[j - 1];
                ijpr[j] = 0;

                let g = miv[j] + fgcp2[j - 1] * gf1va;
                if g > w {
                    w = g;
                    ijpr[j] = -(j as i32 - mpiv[j]);
                }

                let g = mrow[j] + fgcp1va * gapfreq2[j];
                if g > w {
                    w = g;
                    ijpr[j] = i as i32 - mprow[j];
                }
                // `>`, not `>=`: see the note on this function.
                let g = prevw[j - 1] + ogcp1va * gapfreq2[j - 1];
                if g > mrow[j] {
                    mrow[j] = g;
                    mprow[j] = i as i32 - 1;
                }
                mrow[j] += fpenalty_ex;

                curw[j] += w;
                wm = w;
            }
            self.lastverticalw[i] = currentw[lgth2 - 1];
        }

        let currentw: &[f64] = if cur_is_w1 { &self.w1 } else { &self.w2 };
        let (track, impmatch) = match constraint {
            Some(c) => {
                let t = ImpTrack {
                    imp: c.imp,
                    gapmap1: c.gapmap1,
                    gapmap2: c.gapmap2,
                    start1: c.start1,
                    start2: c.start2,
                };
                atracking_localhom(
                    currentw,
                    &self.lastverticalw,
                    lgth1,
                    lgth2,
                    &mut self.ijp,
                    outgap,
                    &t,
                )
            }
            None => (
                atracking(
                    currentw,
                    &self.lastverticalw,
                    lgth1,
                    lgth2,
                    &mut self.ijp,
                    outgap,
                ),
                0.0,
            ),
        };

        for row in group1.iter_mut() {
            *row = profile::gapireru(row, &track.gt1, newgapstr);
        }
        for row in group2.iter_mut() {
            *row = profile::gapireru(row, &track.gt2, newgapstr);
        }

        // `alloclen < resultlen || resultlen > N`, at the end of the C's
        // `partA__align`. It is an error path, but a real one: the refinement
        // loop can grow an alignment past the buffer `dvtditr` sized from the
        // longest input (`nlenmax * 9`), and where the C stops with this
        // message a port without the check would carry on and produce output
        // the C never would.
        const N: usize = 5_000_000;
        let resultlen = group1[0].len();
        if alloclen < resultlen || resultlen > N {
            eprintln!("alloclen={alloclen}, resultlen={resultlen}, N={N}");
            eprintln!("LENGTH OVER!");
            std::process::exit(1);
        }

        (wm, impmatch)
    }

    /// Grow the scratch to fit a window, the way the C's `orlgth1/2` guard does.
    fn ensure(&mut self, lgth1: usize, lgth2: usize, charsize: usize) {
        if lgth1 <= self.orlgth1 && lgth2 <= self.orlgth2 {
            return;
        }
        let ll1 = (lgth1 as f64 * 1.3) as usize + 100;
        let ll2 = (lgth2 as f64 * 1.3) as usize + 100;
        self.w1 = vec![0.0; ll2 + 2];
        self.w2 = vec![0.0; ll2 + 2];
        self.initverticalw = vec![0.0; ll1 + 2];
        self.lastverticalw = vec![0.0; ll1 + 2];
        self.m = vec![0.0; ll2 + 2];
        self.mp = vec![0; ll2 + 2];
        self.miv = vec![0.0; ll2 + 2];
        self.mpiv = vec![0; ll2 + 2];
        self.cand = vec![0.0; ll2 + 2];
        self.ogcp1 = vec![0.0; ll1 + 2];
        self.ogcp2 = vec![0.0; ll2 + 2];
        self.fgcp1 = vec![0.0; ll1 + 2];
        self.fgcp2 = vec![0.0; ll2 + 2];
        self.gapfreq1 = vec![0.0; ll1 + 2];
        self.gapfreq2 = vec![0.0; ll2 + 2];
        self.cpmx1 = vec![vec![0.0; ll1 + 2]; charsize];
        self.cpmx2 = vec![vec![0.0; ll2 + 2]; charsize];
        self.doublework = vec![vec![0.0; charsize]; ll2 + 2];
        self.intwork = vec![vec![0; charsize + 1]; ll2 + 2];
        self.ijp = vec![vec![0; ll2 + 2]; ll1 + 2];
        self.orlgth1 = ll1 - 100;
        self.orlgth2 = ll2 - 100;
    }
}
