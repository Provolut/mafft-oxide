//! `Aalign()` (`core/SAalignmm.c`): the plain group-to-group DP, reached by
//! `-&` (`alg == 'a'`).
//!
//! This is the oldest of MAFFT's profile aligners and much the simplest: one
//! flat `penalty * 0.5` for every gap operation, no opening/final gap counts,
//! no gap frequencies, no extension penalty, and no cached profiles. It scores
//! against `n_dis` -- the static matrix -- rather than the dynamic one the
//! other aligners take as an argument.
//!
//! Nothing in the `mafft` wrapper selects it; `--sp` sets `-@` and no option
//! sets `-&`. It is ported because the goal is the whole option surface of the
//! two binaries, not only what the wrapper reaches.

use crate::constants::Constants;

/// `cpmx_calc()`: composition profile, normalised by the total weight.
///
/// Not `cpmx_calc_new`: this one divides by `totaleff` and accumulates
/// column-outer, member-inner. Both differences are visible in the last bits
/// of the result, so the two are not interchangeable.
fn cpmx_calc(cn: &Constants, seq: &[&[u8]], cpmx: &mut [Vec<f64>], eff: &[f64], lgth: usize, clus: usize) {
    let totaleff: f64 = eff.iter().take(clus).sum();
    for row in cpmx.iter_mut().take(cn.nalphabets) {
        for v in row.iter_mut().take(lgth) {
            *v = 0.0;
        }
    }
    for j in 0..lgth {
        for k in 0..clus {
            cpmx[cn.amino_n[seq[k][j] as usize] as usize][j] += eff[k] / totaleff;
        }
    }
}

/// `match_calc()` as `SAalignmm.c` defines it.
///
/// The work arrays are indexed `[k][j]`, the transpose of the layout
/// `Salignmm.c` uses under `FASTMATCHCALC`, and the scores come from `n_dis`.
#[allow(clippy::too_many_arguments)]
fn match_calc(
    cn: &Constants,
    out: &mut [f64],
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    i1: usize,
    lgth2: usize,
    doublework: &mut [Vec<f64>],
    intwork: &mut [Vec<i32>],
    initialize: bool,
) {
    let na = cn.nalphabets;
    if initialize {
        for j in 0..lgth2 {
            let mut count = 0usize;
            for l in 0..na {
                if cpmx2[l][j] != 0.0 {
                    doublework[count][j] = cpmx2[l][j];
                    intwork[count][j] = l as i32;
                    count += 1;
                }
            }
            intwork[count][j] = -1;
        }
    }
    let mut scarr = vec![0.0f64; na];
    for (l, sc) in scarr.iter_mut().enumerate() {
        let mut s = 0.0;
        for k in 0..na {
            s += cn.n_dis[k][l] as f64 * cpmx1[k][i1];
        }
        *sc = s;
    }
    for j in 0..lgth2 {
        let mut acc = 0.0f64;
        let mut k = 0usize;
        while intwork[k][j] > -1 {
            acc += scarr[intwork[k][j] as usize] * doublework[k][j];
            k += 1;
        }
        out[j] = acc;
    }
}

/// `Atracking()` as `SAalignmm.c` defines it: walk `ijp` back from the corner,
/// writing both groups' columns backwards into their result buffers.
#[allow(clippy::too_many_arguments)]
fn atracking(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
    ijp: &mut [Vec<i32>],
    icyc: usize,
    jcyc: usize,
    outgap: bool,
    newgapstr: u8,
) -> (Vec<Vec<u8>>, Vec<Vec<u8>>) {
    let lgth1 = seq1[0].len();
    let lgth2 = seq2[0].len();

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
    for (i, row) in ijp.iter_mut().enumerate().take(lgth1 + 1) {
        row[0] = i as i32 + 1;
    }
    for j in 0..lgth2 + 1 {
        ijp[0][j] = -(j as i32 + 1);
    }

    let cap = lgth1 + lgth2;
    let mut m1: Vec<Vec<u8>> = vec![vec![0u8; cap]; icyc];
    let mut m2: Vec<Vec<u8>> = vec![vec![0u8; cap]; jcyc];
    let mut p = cap; // written backwards from the end

    let mut iin = lgth1 as i32;
    let mut jin = lgth2 as i32;
    let mut k = 0i32;
    while k <= (lgth1 + lgth2) as i32 {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi) = if v < 0 {
            (iin - 1, jin + v)
        } else if v > 0 {
            (iin - v, jin - 1)
        } else {
            (iin - 1, jin - 1)
        };
        let mut l = iin - ifi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            p -= 1;
            for i in 0..icyc {
                m1[i][p] = seq1[i][(ifi + l) as usize];
            }
            for j in 0..jcyc {
                m2[j][p] = newgapstr;
            }
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            p -= 1;
            for i in 0..icyc {
                m1[i][p] = newgapstr;
            }
            for j in 0..jcyc {
                m2[j][p] = seq2[j][(jfi + l) as usize];
            }
            k += 1;
        }
        if iin <= 0 || jin <= 0 {
            break;
        }
        p -= 1;
        for i in 0..icyc {
            m1[i][p] = seq1[i][ifi as usize];
        }
        for j in 0..jcyc {
            m2[j][p] = seq2[j][jfi as usize];
        }
        k += 1;
        iin = ifi;
        jin = jfi;
    }

    for s in m1.iter_mut() {
        *s = s[p..].to_vec();
    }
    for s in m2.iter_mut() {
        *s = s[p..].to_vec();
    }
    (m1, m2)
}

/// `Aalign()`.
pub fn aalign(
    cn: &Constants,
    group1: &mut [Vec<u8>],
    group2: &mut [Vec<u8>],
    eff1: &[f64],
    eff2: &[f64],
    outgap: bool,
    newgapstr: u8,
) -> f64 {
    let icyc = group1.len();
    let jcyc = group2.len();
    let lgth1 = group1[0].len();
    let lgth2 = group2[0].len();
    let na = cn.nalphabets;
    let fpenalty_half = cn.penalty as f64 * 0.5;

    let ll1 = (1.1 * lgth1 as f64) as usize + 100;
    let ll2 = (1.1 * lgth2 as f64) as usize + 100;
    let mut currentw = vec![0.0f64; ll2 + 2];
    let mut previousw = vec![0.0f64; ll2 + 2];
    let mut initverticalw = vec![0.0f64; ll1 + 2];
    let mut lastverticalw = vec![0.0f64; ll1 + 2];
    let mut m = vec![0.0f64; ll2 + 2];
    let mut mp = vec![0i32; ll2 + 2];
    let mut cpmx1 = vec![vec![0.0f64; ll1 + 2]; na];
    let mut cpmx2 = vec![vec![0.0f64; ll2 + 2]; na];
    let wide = ll1.max(ll2) + 2;
    let mut doublework = vec![vec![0.0f64; wide]; na];
    let mut intwork = vec![vec![0i32; wide]; na];
    let mut ijp = vec![vec![0i32; ll2 + 10]; ll1 + 10];

    let g1: Vec<&[u8]> = group1.iter().map(|s| s.as_slice()).collect();
    let g2: Vec<&[u8]> = group2.iter().map(|s| s.as_slice()).collect();
    cpmx_calc(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
    cpmx_calc(cn, &g2, &mut cpmx2, eff2, lgth2, jcyc);

    match_calc(cn, &mut initverticalw, &cpmx2, &cpmx1, 0, lgth1, &mut doublework, &mut intwork, true);
    match_calc(cn, &mut currentw, &cpmx1, &cpmx2, 0, lgth2, &mut doublework, &mut intwork, true);

    if outgap {
        for i in 1..lgth1 + 1 {
            initverticalw[i] += fpenalty_half;
        }
        for j in 1..lgth2 + 1 {
            currentw[j] += fpenalty_half;
        }
    }

    // `for( j=0; ... ) m[j] = currentw[j-1] + penalty * 0.5` reads
    // `currentw[-1]` at j == 0 -- one element before the allocation. `m[0]`
    // and `mp[0]` are never read afterwards (the DP loop starts at j == 1),
    // so the value that lands there cannot reach the output; the read is
    // skipped rather than invented.
    for j in 1..lgth2 + 1 {
        m[j] = currentw[j - 1] + fpenalty_half;
        mp[j] = 0;
    }
    // `lastverticalw[0] = currentw[lgth2-1]` reads `currentw[-1]` when the
    // block is empty -- `Falign` can cut a zero-width block and `Aalign`, like
    // `MSalignmm`, has no `lgth == 0` early-out. The value cannot reach the
    // output: `lastverticalw` is only read by the `!outgap` search in
    // `Atracking`, whose only effect is `ijp[lgth1][lgth2]`, and with
    // `lgth2 == 0` that cell is overwritten by the `ijp[i][0] = i+1` loop
    // immediately after. So the read is skipped rather than invented.
    if lgth2 > 0 {
        lastverticalw[0] = currentw[lgth2 - 1];
    }

    let lasti = if outgap { lgth1 + 1 } else { lgth1 };
    let mut wm = 0.0f64;
    for i in 1..lasti {
        // `doublencpy`: a copy, not the pointer swap the other aligners use.
        previousw[..lgth2 + 1].copy_from_slice(&currentw[..lgth2 + 1]);
        previousw[0] = initverticalw[i - 1];
        match_calc(cn, &mut currentw, &cpmx1, &cpmx2, i, lgth2, &mut doublework, &mut intwork, false);
        currentw[0] = initverticalw[i];

        let mut mi = previousw[0] + fpenalty_half;
        let mut mpi = 0i32;
        for j in 1..lgth2 + 1 {
            wm = previousw[j - 1];
            ijp[i][j] = 0;

            let mut x = mi + fpenalty_half;
            if x > wm {
                wm = x;
                ijp[i][j] = -((j as i32) - mpi);
            }
            x = previousw[j - 1] + fpenalty_half;
            if mi <= x {
                mi = x;
                mpi = j as i32 - 1;
            }

            x = m[j] + fpenalty_half;
            if x > wm {
                wm = x;
                ijp[i][j] = i as i32 - mp[j];
            }
            x = previousw[j - 1] + fpenalty_half;
            if m[j] <= x {
                m[j] = x;
                mp[j] = i as i32 - 1;
            }
            currentw[j] += wm;
        }
        if lgth2 > 0 {
            lastverticalw[i] = currentw[lgth2 - 1];
        }
    }

    let (m1, m2) = atracking(
        &currentw, &lastverticalw, group1, group2, &mut ijp, icyc, jcyc, outgap, newgapstr,
    );
    for (i, s) in group1.iter_mut().enumerate() {
        *s = m1[i].clone();
    }
    for (j, s) in group2.iter_mut().enumerate() {
        *s = m2[j].clone();
    }
    wm
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Args;

    fn cn_prot() -> Constants {
        let a = Args {
            nblosum: 62,
            ppenalty: -1530,
            ppenalty_ex: 0,
            poffset: 0,
            penalty_shift_factor: 100.0,
            ..Default::default()
        };
        Constants::new(&a)
    }

    #[test]
    fn identical_sequences_align_without_a_gap() {
        let cn = cn_prot();
        let s = b"MNGTEGDNFYVPFSNKTGLRSPYEYPQYYL".to_vec();
        let mut g1 = vec![s.clone()];
        let mut g2 = vec![s.clone()];
        aalign(&cn, &mut g1, &mut g2, &[1.0], &[1.0], true, b'-');
        assert_eq!(g1[0], s);
        assert_eq!(g2[0], s);
    }

    #[test]
    fn an_insertion_is_matched_by_a_gap_of_the_same_length() {
        let cn = cn_prot();
        let short = b"MNGTEGDNFYVPFSNKTGL".to_vec();
        let long = b"MNGTEGDNFYVWWWPFSNKTGL".to_vec();
        let mut g1 = vec![short.clone()];
        let mut g2 = vec![long.clone()];
        aalign(&cn, &mut g1, &mut g2, &[1.0], &[1.0], true, b'-');
        assert_eq!(g1[0].len(), g2[0].len());
        assert_eq!(g1[0].iter().filter(|&&c| c == b'-').count(), 3);
        assert_eq!(g2[0].iter().filter(|&&c| c == b'-').count(), 0);
        assert_eq!(
            g1[0].iter().filter(|&&c| c != b'-').copied().collect::<Vec<u8>>(),
            short
        );
    }

    #[test]
    fn an_empty_block_becomes_gaps_opposite_the_other_side() {
        // `Falign` cuts zero-width blocks, and `Aalign` has no early-out for
        // them -- it walks the traceback from `ijp[lgth1][0]` instead.
        let cn = cn_prot();
        let mut g1 = vec![b"MNGTE".to_vec()];
        let mut g2 = vec![Vec::new()];
        aalign(&cn, &mut g1, &mut g2, &[1.0], &[1.0], true, b'-');
        assert_eq!(g1[0], b"MNGTE".to_vec());
        assert_eq!(g2[0], b"-----".to_vec());
    }

    #[test]
    fn a_group_keeps_every_row_the_same_width() {
        let cn = cn_prot();
        let mut g1 = vec![
            b"MNGTEGDNFYVPFSNKTGL".to_vec(),
            b"MNGTEGDNFY--FSNKTGL".to_vec(),
        ];
        let mut g2 = vec![b"MNGTEGPNFYVPFSNITGVVRSP".to_vec()];
        aalign(&cn, &mut g1, &mut g2, &[0.5, 0.5], &[1.0], true, b'-');
        let alen = g1[0].len();
        for s in g1.iter().chain(g2.iter()) {
            assert_eq!(s.len(), alen);
        }
    }
}
