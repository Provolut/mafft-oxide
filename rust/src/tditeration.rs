//! Iterative refinement, transliterated from `core/tditeration.c`
//! (`TreeDependentIteration`, sequential branch).
//!
//! The progressive alignment `tbfast` produces is only as good as the order it
//! was built in: an early mistake is carried by every later step. Refinement
//! attacks that directly. Each branch of the guide tree splits the sequences
//! into two groups; the two are pulled apart, realigned against each other, and
//! the result is kept only if it scores better than what was there. Walk every
//! branch, in alternating directions, until nothing improves.
//!
//! **One configuration is implemented, and it is the one `scripts/mafft` runs.**
//! A traced `dvtditr` reports `utree=1 weight=4 score_check=1 nkozo=0 rnakozo=0
//! alg='A' use_fft=1 nthread=0` on L-INS-i, E-INS-i and FFT-NS-i alike, with
//! `constraint` 2 (a `hat3` table present) or 0. Everything else is refused
//! rather than approximated:
//!
//! * `utree=0` -- rebuilding the tree from the current alignment every round
//!   (`substitution_hosei` then `nj`/`spg`/`supg`/`upg2`). The `mafft` script
//!   never passes it, and it is a second tree builder, not a variation on this
//!   one.
//! * `weight` other than 4, `score_check` other than 1, any `nkozo`, `rnakozo`,
//!   `--target`, `specificityconsideration`, `alg` other than `'A'`, `-C n`,
//!   and the BESTFIRST strategy.
//!
//! The refusals are not defensive noise: `README.md` records that a
//! callgrind profile of the real binary found ten functions in this file's
//! dependency list that no strategy ever enters. Porting them would be writing
//! code against a reading of the source rather than against the program.
use crate::constants::Constants;
use crate::falign::{falign, falign_localhom, FalignState};
use crate::localhom::{calcimportance_half, fastshrinklocalhom_half, LocalHomTable};
use crate::partalign::{PartAlignState, PartImp};
use crate::treeop::{
    calc_branch_weight, counteff_simple, dist_from_a_branch, one_cluster_and_the_other_fast,
    tree_cnv, weight_from_a_branch, Topol,
};

/// What `dvtditr` hands the refinement loop, beyond the alignment itself.
pub struct IterParams {
    /// `-I`: how many passes over every branch.
    pub niter: usize,
    /// `-E`: the acceptance margin, as a percentage of the current score.
    pub cut: f64,
    /// Halve `cut` each round. The script leaves it on with `cut` at zero, so
    /// it changes no number -- but it does gate the oscillation check.
    pub cooling: bool,
    /// `-W`: the floor a member's weight is raised to before normalising.
    pub minimumweight: f64,
    /// 2 when a `hat3` table is present, 0 otherwise. 1 is not reachable.
    pub constraint: i32,
    pub weight: i32,
    pub score_check: i32,
    pub use_fft: bool,
    pub utree: bool,
    pub alg: u8,
    pub nkozo: usize,
    pub outgap: bool,
    pub legacygapcost: bool,
    pub newgapstr: u8,
    /// `-l`: `fillimp`'s weight on a constraint.
    pub fastathreshold: f64,
    /// The row capacity `dvtditr` allocated (`nlenmax * 9`). The DP aborts
    /// with the C's `LENGTH OVER` when a refinement outgrows it.
    pub alloclen: usize,
}

/// Scratch shared across every branch of every round, as the C shares its
/// `static` buffers.
#[derive(Default)]
struct Scratch {
    fs: Option<FalignState>,
    ast: crate::align::AlignState,
    past: PartAlignState,
    imp: PartImp,
}

/// Why the loop stopped, as the C's return value.
///
/// `0` converged, `-1` started oscillating, `2` ran out of rounds. `dvtditr`
/// prints nothing different for any of them, but the value is the C's and
/// callers may come to depend on it.
pub type Outcome = i32;

/// `TreeDependentIteration()`, sequential branch.
///
/// Args:
///     cn: Scoring constants.
///     p: The run's parameters.
///     locnjob: Number of sequences.
///     aseq: Working copy; each branch's realignment lands here first.
///     bseq: The best alignment so far, and the result.
///     topol: The guide tree as a merge sequence.
///     len: Its branch lengths; clamped in place by `treeCnv`.
///     skipthisbranch: Branches `dvtditr` has been told not to touch.
///     localhomtable: The `hat3` constraints, when `constraint == 2`.
///     write_pre: Called with the best alignment whenever it improves, which is
///         where the C writes its `pre` file.
///
/// Returns:
///     The C's outcome code.
#[allow(clippy::too_many_arguments)]
pub fn tree_dependent_iteration(
    cn: &Constants,
    p: &IterParams,
    locnjob: usize,
    aseq: &mut [Vec<u8>],
    bseq: &mut [Vec<u8>],
    topol: &Topol,
    len: &mut Vec<[f64; 2]>,
    skipthisbranch: &[[bool; 2]],
    localhomtable: Option<&mut LocalHomTable>,
    write_pre: &mut dyn FnMut(&[Vec<u8>]),
) -> Outcome {
    assert!(p.utree, "dynamic tree rebuilding (utree=0) is not implemented");
    // `weight` is 4 for every real run, and 0 for the two-sequence case, where
    // `dvtditr`'s `main` sets it along with `niter = 1`. Nothing else is
    // implemented.
    assert!(
        p.weight == 4 || p.weight == 0,
        "only weight=4 (and weight=0, which dvtditr sets for two sequences) is implemented"
    );
    assert_eq!(p.score_check, 1, "only score_check=1 is implemented");
    assert_eq!(p.nkozo, 0, "structural (kozo) weighting is not implemented");
    assert!(p.use_fft, "the C refuses !use_fft here itself");
    assert_eq!(p.alg, b'A', "only alg='A' is implemented");
    assert!(
        p.constraint == 0 || p.constraint == 2,
        "constraint must be 0 or 2"
    );

    let mut st = Scratch::default();

    // The C allocates `len` with `locnjob` rows and reads `len[locnjob-2]`
    // when it unroots; a caller that sized it to the `locnjob-1` real steps
    // would index past the end inside `tree_cnv`.
    if len.len() < locnjob {
        len.resize(locnjob, [0.0, 0.0]);
    }

    let mut table = localhomtable;
    if p.constraint != 0 {
        let eff = counteff_simple(locnjob, topol, len);
        calcimportance_half(
            locnjob,
            &eff,
            aseq,
            table.as_mut().expect("constraint without a hat3 table"),
        );
    }

    // `weight == 4` needs the unrooted graph and its branch weights. Two
    // sequences have no branch to weight, and `treeCnv` would index
    // `stopol[-1]`, so the C guards on `locnjob > 2` and so does this.
    let mut stopol = Vec::new();
    let mut branch_weight = vec![[0.0f64; 2]; locnjob];
    // Two sequences have one branch and no weighting to do, and `tree_cnv`
    // would index `stopol[locnjob-3]`; the C guards this the same way.
    if locnjob > 2 && (p.weight == 4 || p.weight == 0) {
        stopol = tree_cnv(locnjob, topol, len);
        calc_branch_weight(&mut branch_weight, locnjob, &stopol, topol);
    }

    for i in 0..locnjob {
        bseq[i] = aseq[i].clone();
    }
    write_pre(bseq);

    let mut history: Vec<Vec<[f64; 2]>> = vec![vec![[0.0; 2]; locnjob]; p.niter.max(1)];
    let mut converged = 0usize;
    let mut cut = p.cut;
    if p.cooling {
        cut *= 2.0;
    }

    for iterate in 0..p.niter {
        if p.cooling {
            cut *= 0.5;
        }

        // Alternate direction each round: the branch a pass reaches last is the
        // one the next pass reaches first.
        let (lin, ldf): (isize, isize) = if iterate % 2 == 0 {
            (0, 1)
        } else {
            (locnjob as isize - 2, -1)
        };

        let mut l = lin;
        while l < locnjob as isize - 1 && l >= 0 {
            let li = l as usize;
            // The C writes `for( k=0; k<2; k++ ) { if( l == locnjob-2 ) k = 1;`
            // -- on the last step only the second side is a real branch, and
            // assigning to the loop variable is how it skips the first.
            let ks: &[usize] = if li == locnjob - 2 { &[1] } else { &[0, 1] };
            for &k in ks {
                if skipthisbranch[li][k] {
                    continue;
                }

                let mut distarr = vec![0.0f64; locnjob];
                dist_from_a_branch(locnjob, &mut distarr, &stopol, topol, len, li, k);
                let (mem1, mem2) = one_cluster_and_the_other_fast(locnjob, topol, li, k);
                let m1: Vec<i32> = mem1.iter().copied().take_while(|&v| v != -1).collect();
                let m2: Vec<i32> = mem2.iter().copied().take_while(|&v| v != -1).collect();
                let (s1, s2) = (m1[0] as usize, m2[0] as usize);

                let mut effarr = vec![1.0f64; locnjob];
                if p.weight == 4 {
                    weight_from_a_branch(
                        locnjob,
                        &mut effarr,
                        &stopol,
                        &branch_weight,
                        topol,
                        li,
                        k,
                    );
                }

                for i in 0..locnjob {
                    aseq[i] = bseq[i].clone();
                }
                let length = aseq[0].len();

                let (effarr1, _) = fastconjuction_noname(&m1, &effarr, p.minimumweight);
                let (effarr2, _) = fastconjuction_noname(&m2, &effarr, p.minimumweight);
                let clus1 = m1.len();
                let clus2 = m2.len();
                let mut mseq1: Vec<Vec<u8>> = m1.iter().map(|&i| aseq[i as usize].clone()).collect();
                let mut mseq2: Vec<Vec<u8>> = m2.iter().map(|&i| aseq[i as usize].clone()).collect();

                // --- the score to beat, on the alignment as it stands ---
                let tmpdouble = crate::score::intergroup_score(
                    cn, &mseq1, &mseq2, &effarr1, &effarr2, clus1, clus2, length,
                );
                let mut oimpmatchdouble = 0.0f64;
                if p.constraint != 0 {
                    let t = table.as_ref().expect("constraint without a hat3 table");
                    let shrink = fastshrinklocalhom_half(&m1, &m2, t);
                    let r1: Vec<&[u8]> = mseq1.iter().map(std::vec::Vec::as_slice).collect();
                    let r2: Vec<&[u8]> = mseq2.iter().map(std::vec::Vec::as_slice).collect();
                    st.imp.init_strict(
                        clus1,
                        clus2,
                        length,
                        length,
                        &r1,
                        &r2,
                        &effarr1,
                        &effarr2,
                        &shrink,
                        &m1,
                        &m2,
                        p.fastathreshold,
                    );
                    // Summed downwards, as the C does: the order decides the
                    // last bits of a float sum.
                    for i in (0..length).rev() {
                        oimpmatchdouble += st.imp.score(i, i);
                    }
                }
                let mscore = if p.constraint != 0 {
                    oimpmatchdouble + tmpdouble
                } else {
                    tmpdouble
                };

                // --- realign the two groups against each other ---
                let fs = st.fs.get_or_insert_with(FalignState::new);
                let mut impmatchdouble = 0.0f64;
                if p.constraint == 2 {
                    let (_, imp) = falign_localhom(
                        fs,
                        &mut st.past,
                        &st.imp,
                        cn,
                        &mut mseq1,
                        &mut mseq2,
                        &effarr1,
                        &effarr2,
                        p.outgap,
                        p.legacygapcost,
                        p.newgapstr,
                        p.alloclen,
                    );
                    impmatchdouble = imp;
                } else {
                    let mut fftlog = 0i32;
                    falign(
                        fs,
                        &mut st.ast,
                        cn,
                        &cn.n_dis_consweight_multi,
                        &mut mseq1,
                        &mut mseq2,
                        &effarr1,
                        &effarr2,
                        p.outgap,
                        &mut fftlog,
                        p.newgapstr,
                        p.alg,
                        // `dvtditr` sets `kobetsubunkatsu`, which turns off the
                        // FFT anchoring inside `Falign` exactly as it does in
                        // `Falign_localhom`.
                        true,
                    );
                }

                // The C aligns rows of `aseq` in place, because `mseq1[i]` is
                // `aseq[member]`. Here the rows were copied out, so they go
                // back; the two groups are disjoint, so nothing is written
                // twice.
                for (i, &m) in m1.iter().enumerate() {
                    aseq[m as usize] = mseq1[i].clone();
                }
                for (i, &m) in m2.iter().enumerate() {
                    aseq[m as usize] = mseq2[i].clone();
                }

                // --- keep it or throw it away ---
                let identity = aseq[s1] == bseq[s1] && aseq[s2] == bseq[s2];
                let length = mseq1[0].len();
                let mut tscore;
                if identity {
                    tscore = mscore;
                    converged += 1;
                } else {
                    let tmpdouble = crate::score::intergroup_score(
                        cn, &mseq1, &mseq2, &effarr1, &effarr2, clus1, clus2, length,
                    );
                    tscore = if p.constraint == 2 {
                        impmatchdouble + tmpdouble
                    } else {
                        tmpdouble
                    };

                    assert!(!mscore.is_nan(), "mscore became NaN");
                    assert!(!tscore.is_nan(), "tscore became NaN");

                    if tscore > mscore - cut / 100.0 * mscore {
                        write_pre(aseq);
                        for i in 0..locnjob {
                            bseq[i] = aseq[i].clone();
                        }
                        converged = 0;
                    } else {
                        tscore = mscore;
                        converged += 1;
                    }
                }

                history[iterate][li][k] = tscore;

                if converged >= locnjob * 2 {
                    return 0;
                }
                if iterate >= 1 {
                    // Oscillation: this branch has produced this exact score
                    // before, on a round of the same parity. Two rounds back is
                    // the same walk direction, which is why the step is 2.
                    let mut oscillating = false;
                    let mut ii = iterate as isize - 2;
                    while ii >= 0 {
                        if tscore == history[ii as usize][li][k] {
                            oscillating = true;
                            break;
                        }
                        ii -= 2;
                    }
                    if (oscillating && !p.cooling) || (oscillating && cut < 0.001 && p.cooling) {
                        return -1;
                    }
                }
            }
            l += ldf;
        }
    }
    2
}

/// `fastconjuction_noname()`: one group's member weights, floored and
/// normalised.
///
/// `treebase.rs` has this too, but privately and for its own caller. Duplicated
/// rather than shared for now; the two should become one when `dvtditr-rs`
/// lands and both callers are in the tree.
fn fastconjuction_noname(memlist: &[i32], eff: &[f64], mineff: f64) -> (Vec<f64>, f64) {
    let mut peff: Vec<f64> = memlist
        .iter()
        .map(|&m| {
            let e = eff[m as usize];
            if e < mineff {
                mineff
            } else {
                e
            }
        })
        .collect();
    let mut total = 0.0f64;
    for v in peff.iter() {
        total += *v;
    }
    for v in peff.iter_mut() {
        *v /= total;
    }
    (peff, total)
}
