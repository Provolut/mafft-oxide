//! The all-pairs pairwise stage, transliterated from `pairalign()` in
//! `core/pairlocalalign.c` (single-threaded, `alg == 'L'` or `'A'`).
//!
//! For every pair of input sequences it runs one pairwise alignment — local
//! (`L__align11`) for L-INS-i, global (`G__align11`) for G-INS-i — and keeps
//! two things: a distance for the guide tree, and the ungapped blocks of that
//! alignment as constraints for the progressive stage.

use crate::align::AlignState;
use crate::constants::{dist2offset, makedynamicmtx, Constants, NALPHABETS};
use crate::genalign::GenalignState;
use crate::lalign::LalignState;
use crate::localhom::{putlocalhom2, LocalHomTable};
use crate::score::naivepairscore11;

/// `score2dist()`.
fn score2dist(pscore: f64, selfscore1: f64, selfscore2: f64) -> f64 {
    let bunbo = if selfscore1 < selfscore2 {
        selfscore1
    } else {
        selfscore2
    };
    if bunbo == 0.0 {
        2.0
    } else if bunbo < pscore {
        0.0
    } else {
        (1.0 - pscore / bunbo) * 2.0
    }
}

/// `removex()`: strip `X`/`x`, returning the stripped copy and how many went.
fn removex(m: &[u8]) -> (Vec<u8>, i32) {
    let mut d = Vec::with_capacity(m.len());
    let mut val = 0;
    for &c in m {
        if c == b'X' || c == b'x' {
            val += 1;
        } else {
            d.push(c);
        }
    }
    (d, val)
}

/// Which pairwise aligner the strategy selects.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum PairAlg {
    /// `-L`: local pairwise alignment (L-INS-i).
    Local,
    /// `-A`: global pairwise alignment (G-INS-i).
    Global,
    /// `-N`: generalised affine-gap local alignment (E-INS-i).
    Generalized,
}

/// `pairalign()`: fill the distance half-matrix and the constraint table.
///
/// `seq` are the ungapped input sequences. `outgap` is the global flag the C
/// reads inside `G__align11`.
#[allow(clippy::too_many_arguments)]
pub fn pairalign(
    cn: &Constants,
    alg: PairAlg,
    njob: usize,
    seq: &[Vec<u8>],
    distancemtx: &mut [Vec<f64>],
    table: &mut LocalHomTable,
    store_localhom: bool,
    outgap: bool,
    lst: &mut LalignState,
    ast: &mut AlignState,
    gst: &mut GenalignState,
    // `-Z`: score the alignment with `naivepairscore11` rather than keeping
    // the DP's own score. The wrapper turns this on for `--genafpair`.
    usenaivescore: bool,
) {
    // `dseq`: the X-free copy used for the score-only variants. `X` is only
    // stripped for protein against a real matrix -- for nucleotides the C
    // copies the sequence through and sets `thereisx` to zero.
    let mut dseq: Vec<Vec<u8>> = Vec::with_capacity(njob);
    let mut thereisxineachseq = vec![0i32; njob];
    for i in 0..njob {
        if cn.dorp == b'p' {
            let (d, n) = removex(&seq[i]);
            dseq.push(d);
            thereisxineachseq[i] = n;
        } else {
            dseq.push(seq[i].clone());
            thereisxineachseq[i] = 0;
        }
    }

    let mut dynmtx: Vec<Vec<f64>> = vec![vec![0.0; NALPHABETS]; NALPHABETS];
    let spec = cn.specificityconsideration;

    let mut selfscore = vec![0.0f64; njob];
    for i in 0..njob {
        let mut pscore = 0.0f64;
        for &c in &seq[i] {
            pscore += cn.amino_dis[c as usize][c as usize] as f64;
        }
        selfscore[i] = pscore;
    }

    for i in 0..njob - 1 {
        for j in i + 1..njob {
            if seq[i].is_empty() || seq[j].is_empty() {
                distancemtx[i][j - i] = 3.0;
                continue;
            }

            let thereisx = thereisxineachseq[i] + thereisxineachseq[j];
            let mut a1: Vec<u8>;
            let mut a2: Vec<u8>;
            let mut off1 = 0usize;
            let mut off2 = 0usize;
            let mut pscore;

            match alg {
                PairAlg::Local => {
                    let r = lst.l_align11(cn, &cn.n_dis_consweight_multi, 0.0, &seq[i], &seq[j]);
                    pscore = r.score;
                    a1 = r.a1;
                    a2 = r.a2;
                    off1 = r.off1;
                    off2 = r.off2;
                    if thereisx != 0 {
                        pscore = lst.l_align11_noalign(
                            cn,
                            &cn.n_dis_consweight_multi,
                            &dseq[i],
                            &dseq[j],
                        );
                    }
                    // --allowshift: redo the alignment with a matrix flattened
                    // in proportion to how close the pair is. The score keeps
                    // the first alignment's value; only the layout changes.
                    if spec > 0.0 {
                        let dist = score2dist(pscore, selfscore[i], selfscore[j]);
                        let scoreoffset = dist2offset(dist, spec);
                        if scoreoffset < 0.0 {
                            makedynamicmtx(
                                &mut dynmtx,
                                &cn.n_dis_consweight_multi,
                                0.5 * dist,
                                &cn.amino,
                                spec,
                            );
                            let r = lst.l_align11(cn, &dynmtx, scoreoffset, &seq[i], &seq[j]);
                            a1 = r.a1;
                            a2 = r.a2;
                            off1 = r.off1;
                            off2 = r.off2;
                        }
                    }
                }
                PairAlg::Generalized => {
                    let r = gst.gen_l_align11(cn, &cn.n_dis_consweight_multi, &seq[i], &seq[j]);
                    a1 = r.a1;
                    a2 = r.a2;
                    off1 = r.off1;
                    off2 = r.off2;
                    if usenaivescore {
                        // `-Z`: the DP's score is overwritten, and the
                        // `thereisx` / `--allowshift` re-runs below are skipped
                        // entirely -- they live in the other branch of the C.
                        pscore = naivepairscore11(cn, &a1, &a2, 0);
                    } else {
                        pscore = r.score;
                        if thereisx != 0 {
                            let r2 = gst.gen_l_align11(
                                cn,
                                &cn.n_dis_consweight_multi,
                                &dseq[i],
                                &dseq[j],
                            );
                            pscore = r2.score;
                        }
                        if spec > 0.0 {
                            let dist = score2dist(pscore, selfscore[i], selfscore[j]);
                            if dist2offset(dist, spec) < 0.0 {
                                makedynamicmtx(
                                    &mut dynmtx,
                                    &cn.n_dis_consweight_multi,
                                    0.5 * dist,
                                    &cn.amino,
                                    spec,
                                );
                                let r3 = gst.gen_l_align11(cn, &dynmtx, &seq[i], &seq[j]);
                                a1 = r3.a1;
                                a2 = r3.a2;
                                off1 = r3.off1;
                                off2 = r3.off2;
                            }
                        }
                    }
                }
                PairAlg::Global => {
                    a1 = seq[i].clone();
                    a2 = seq[j].clone();
                    pscore = ast.g_align11(
                        cn,
                        &cn.n_dis_consweight_multi,
                        cn.penalty,
                        cn.penalty_ex,
                        &mut a1,
                        &mut a2,
                        outgap,
                        outgap,
                        b'-',
                    );
                    if thereisx != 0 {
                        pscore = lst.g_align11_noalign(
                            cn,
                            &cn.n_dis_consweight_multi,
                            cn.penalty,
                            cn.penalty_ex,
                            &dseq[i],
                            &dseq[j],
                        );
                    }
                    if spec > 0.0 {
                        let dist = score2dist(pscore, selfscore[i], selfscore[j]);
                        if dist2offset(dist, spec) < 0.0 {
                            makedynamicmtx(
                                &mut dynmtx,
                                &cn.n_dis_consweight_multi,
                                0.5 * dist,
                                &cn.amino,
                                spec,
                            );
                            a1 = seq[i].clone();
                            a2 = seq[j].clone();
                            ast.g_align11(
                                cn,
                                &dynmtx,
                                cn.penalty,
                                cn.penalty_ex,
                                &mut a1,
                                &mut a2,
                                outgap,
                                outgap,
                                b'-',
                            );
                        }
                    }
                }
            }

            // C: if( mseq1[0][0] != 0 && mseq2[0][0] != 0 )
            if !a1.is_empty() && !a2.is_empty() {
                if store_localhom {
                    putlocalhom2(cn, &a1, &a2, &mut table[i][j - i], off1, off2, b'h');
                }
                pscore = score2dist(pscore, selfscore[i], selfscore[j]);
            } else {
                pscore = 2.0;
            }
            distancemtx[i][j - i] = pscore;
        }
    }
}
