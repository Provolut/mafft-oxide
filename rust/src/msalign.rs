//! `MSalignmm()`: the linear-space group-to-group aligner (`core/MSalignmm.c`).
//!
//! This is the same affine profile DP as `A__align`, arranged as Hirschberg's
//! divide and conquer: instead of keeping an `m x n` traceback matrix, each
//! call finds the row where the optimal path crosses the middle column and
//! recurses into the two halves, so the working space is linear in the shorter
//! side. The C reaches it three ways -- `-M`, `--memsave`, and, without being
//! asked, whenever a profile row passes 30000 columns.
//!
//! Two compile-time switches in `MSalignmm.c` change the recurrence relative
//! to `Salignmm.c`, and both are load-bearing:
//!
//! * `USE_PENALTY_EX 0` -- the running gap scores are **not** charged an
//!   extension penalty here, where `A__align` (`USE_PENALTY_EX 1`) charges one
//!   on every cell. This is why the two aligners disagree on inputs with real
//!   ties, not merely in how much memory they use.
//! * `ATO 1` -- the gap-count tables are scaled by the gap frequencies before
//!   the DP rather than inside it.

use crate::constants::Constants;

/// `match_calc_alphabet_seq()`: the match row for profile column `i1`.
///
/// Same as `align.rs`'s `match_calc` except that both profiles are indexed
/// from an offset, since a recursion step works on a sub-block rather than on
/// the whole profile.
#[allow(clippy::too_many_arguments)]
fn match_calc_alphabet_seq(
    n_dynamicmtx: &[Vec<f64>],
    out: &mut [f64],
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    i1: usize,
    start2: usize,
    lgth2: usize,
    doublework: &mut [Vec<f64>],
    intwork: &mut [Vec<i32>],
    initialize: bool,
) {
    let nalphabets = n_dynamicmtx.len();
    if initialize {
        for j in 0..lgth2 {
            let p = start2 + j;
            let mut count = 0usize;
            for l in 0..nalphabets {
                if cpmx2[l][p] != 0.0 {
                    doublework[j][count] = cpmx2[l][p];
                    intwork[j][count] = l as i32;
                    count += 1;
                }
            }
            intwork[j][count] = -1;
        }
    }

    let mut scarr = vec![0.0f64; nalphabets];
    for (l, sc) in scarr.iter_mut().enumerate() {
        let mut s = 0.0;
        for j in 0..nalphabets {
            s += n_dynamicmtx[j][l] * cpmx1[j][i1];
        }
        *sc = s;
    }

    for j in 0..lgth2 {
        let mut acc = 0.0f64;
        let mut k = 0usize;
        while intwork[j][k] > -1 {
            acc += scarr[intwork[j][k] as usize] * doublework[j][k];
            k += 1;
        }
        out[j] = acc;
    }
}

/// The gap tables the DP is scored against, sliced per recursion step.
///
/// `gapinfo` in the C: six rows, indexed from `ist`/`jst` by the callee, so a
/// sub-block sees its own window of each. `gapfreq1[-1]` is read when the
/// block does not start at column 0, which is why these are whole-profile
/// vectors with an offset rather than owned slices.
pub struct GapInfo {
    pub ogcp1: Vec<f64>,
    pub fgcp1: Vec<f64>,
    pub ogcp2: Vec<f64>,
    pub fgcp2: Vec<f64>,
    pub gapfreq1: Vec<f64>,
    pub gapfreq2: Vec<f64>,
}

/// `Atracking()` as `MSalignmm.c` defines it: walk `ijp` back from the corner
/// and write the two gap tables, then apply them to the sequences.
///
/// Returns the two gap tables. The C builds them backwards into the end of a
/// buffer and hands back a pointer to where it stopped; this returns the
/// filled tail directly.
#[allow(clippy::too_many_arguments)]
fn atracking(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
    mseq1: &mut [Vec<u8>],
    mseq2: &mut [Vec<u8>],
    ijp: &mut [Vec<i32>],
    icyc: usize,
    jcyc: usize,
    ist: usize,
    ien: usize,
    jst: usize,
    jen: usize,
    fulllen1: usize,
    fulllen2: usize,
    tailgp: bool,
) -> (Vec<u8>, Vec<u8>) {
    let lgth1 = ien - ist + 1;
    let lgth2 = jen - jst + 1;
    let cap = lgth1 + lgth2;
    let mut gt1 = vec![0u8; cap + 1];
    let mut gt2 = vec![0u8; cap + 1];
    let mut p1 = cap; // one past the last byte written, counting down
    let mut p2 = cap;

    if tailgp {
        // nothing: the path must end in the corner
    } else if ien == fulllen1 - 1 || jen == fulllen2 - 1 {
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

    let mut iin = lgth1 as i32;
    let mut jin = lgth2 as i32;
    let klim = (lgth1 + lgth2) as i32;
    let mut k = 0i32;
    while k <= klim {
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
            p1 -= 1;
            gt1[p1] = b'o';
            p2 -= 1;
            gt2[p2] = b'-';
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            p1 -= 1;
            gt1[p1] = b'-';
            p2 -= 1;
            gt2[p2] = b'o';
            k += 1;
        }
        if iin <= 0 || jin <= 0 {
            break;
        }
        p1 -= 1;
        gt1[p1] = b'o';
        p2 -= 1;
        gt2[p2] = b'o';
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1;
    }

    let gaptable1 = gt1[p1..cap].to_vec();
    let gaptable2 = gt2[p2..cap].to_vec();
    let klim = cap - p1;

    if gaptable1.contains(&b'-') {
        for i in 0..icyc {
            mseq1[i] = crate::profile::gapireru(&seq1[i][ist..], &gaptable1, b'-');
        }
    } else {
        for i in 0..icyc {
            mseq1[i] = seq1[i][ist..(ist + klim).min(seq1[i].len())].to_vec();
        }
    }
    if gaptable2.contains(&b'-') {
        for j in 0..jcyc {
            mseq2[j] = crate::profile::gapireru(&seq2[j][jst..], &gaptable2, b'-');
        }
    } else {
        for j in 0..jcyc {
            mseq2[j] = seq2[j][jst..(jst + klim).min(seq2[j].len())].to_vec();
        }
    }

    (gaptable1, gaptable2)
}

/// `MSalignmm_tanni()`: the base case -- a full DP with a real traceback
/// matrix, once the block is small enough to afford one.
#[allow(clippy::too_many_arguments)]
fn msalignmm_tanni(
    n_dynamicmtx: &[Vec<f64>],
    icyc: usize,
    jcyc: usize,
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    ist: usize,
    ien: usize,
    jst: usize,
    jen: usize,
    fulllen1: usize,
    fulllen2: usize,
    mseq1: &mut [Vec<u8>],
    mseq2: &mut [Vec<u8>],
    gi: &GapInfo,
    headgp: bool,
    tailgp: bool,
    headgapfreq1_g: f64,
    headgapfreq2_g: f64,
) -> (f64, Vec<u8>, Vec<u8>) {
    let og1 = &gi.ogcp1[ist..];
    let fg1 = &gi.fgcp1[ist..];
    let og2 = &gi.ogcp2[jst..];
    let fg2 = &gi.fgcp2[jst..];
    let gf1 = &gi.gapfreq1[ist..];
    let gf2 = &gi.gapfreq2[jst..];

    // `gapfreq1f[-1]`: one column before the block, when there is one.
    let headgapfreq1 = if ist > 0 {
        gi.gapfreq1[ist - 1]
    } else {
        headgapfreq1_g
    };
    let headgapfreq2 = if jst > 0 {
        gi.gapfreq2[jst - 1]
    } else {
        headgapfreq2_g
    };

    let lgth1 = ien - ist + 1;
    let lgth2 = jen - jst + 1;
    let ll1 = lgth1 + 100;
    let ll2 = lgth2 + 100;

    let mut w1 = vec![0.0f64; ll2 + 2];
    let mut w2 = vec![0.0f64; ll2 + 2];
    let mut initverticalw = vec![0.0f64; ll1 + 2];
    let mut lastverticalw = vec![0.0f64; ll1 + 2];
    let mut m = vec![0.0f64; ll2 + 2];
    let mut mp = vec![0i32; ll2 + 2];

    let nalphabets = n_dynamicmtx.len();
    let wide = ll1.max(ll2) + 2;
    let mut doublework = vec![vec![0.0f64; nalphabets + 1]; wide];
    let mut intwork = vec![vec![0i32; nalphabets + 1]; wide];
    let mut ijp = vec![vec![0i32; ll2 + 1]; ll1 + 1];

    let mut wm = 0.0f64;

    match_calc_alphabet_seq(
        n_dynamicmtx,
        &mut initverticalw,
        cpmx2,
        cpmx1,
        jst,
        ist,
        lgth1,
        &mut doublework,
        &mut intwork,
        true,
    );
    match_calc_alphabet_seq(
        n_dynamicmtx,
        &mut w1,
        cpmx1,
        cpmx2,
        ist,
        jst,
        lgth2,
        &mut doublework,
        &mut intwork,
        true,
    );

    if headgp || ist != 0 {
        for i in 1..lgth1 + 1 {
            initverticalw[i] += og1[0] * headgapfreq2 + fg1[i - 1] * gf2[0];
        }
    }
    if headgp || jst != 0 {
        for j in 1..lgth2 + 1 {
            w1[j] += og2[0] * headgapfreq1 + fg2[j - 1] * gf1[0];
        }
    }

    for j in 1..lgth2 + 1 {
        m[j] = w1[j - 1] + og1[1] * gf2[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = w1[lgth2 - 1];

    // `currentw`/`previousw` swap each row; carried as a flag over w1/w2.
    let mut cur_is_w1 = true;
    let lasti = if tailgp || jen != fulllen2 - 1 {
        lgth1 + 1
    } else {
        lgth1
    };
    for i in 1..lasti {
        cur_is_w1 = !cur_is_w1;
        let (currentw, previousw) = if cur_is_w1 {
            (&mut w1, &mut w2)
        } else {
            (&mut w2, &mut w1)
        };
        previousw[0] = initverticalw[i - 1];
        match_calc_alphabet_seq(
            n_dynamicmtx,
            currentw,
            cpmx1,
            cpmx2,
            ist + i,
            jst,
            lgth2,
            &mut doublework,
            &mut intwork,
            false,
        );
        currentw[0] = initverticalw[i];

        let mut mi = previousw[0] + og2[1] * gf1[i - 1];
        let mut mpi = 0i32;

        for j in 1..lgth2 + 1 {
            let mut wmv = previousw[j - 1];
            ijp[i][j] = 0;

            let mut g = mi + fg2[j - 1] * gf1[i];
            if g > wmv {
                wmv = g;
                ijp[i][j] = -((j as i32) - mpi);
            }
            g = previousw[j - 1] + og2[j] * gf1[i - 1];
            if g >= mi {
                mi = g;
                mpi = j as i32 - 1;
            }
            // USE_PENALTY_EX is 0 here: no `mi += fpenalty_ex`.

            g = m[j] + fg1[i - 1] * gf2[j];
            if g > wmv {
                wmv = g;
                ijp[i][j] = i as i32 - mp[j];
            }
            g = previousw[j - 1] + og1[i] * gf2[j - 1];
            if g >= m[j] {
                m[j] = g;
                mp[j] = i as i32 - 1;
            }
            // ... and no `m[j] += fpenalty_ex` either.

            currentw[j] += wmv;
            wm = wmv;
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    let currentw: &[f64] = if cur_is_w1 { &w1 } else { &w2 };
    let (gt1, gt2) = atracking(
        currentw,
        &lastverticalw,
        seq1,
        seq2,
        mseq1,
        mseq2,
        &mut ijp,
        icyc,
        jcyc,
        ist,
        ien,
        jst,
        jen,
        fulllen1,
        fulllen2,
        tailgp,
    );
    (wm, gt1, gt2)
}

/// `DPTANNI`: below this on either side the block is done with a full DP.
const DPTANNI: usize = 100;

/// `MSalignmm_rec()`: split the block at its middle row and recurse.
///
/// A forward DP runs to row `imid = lgth1/2` keeping the three score vectors
/// there (`midw` match, `midm` vertical gap, `midn` horizontal gap), and a
/// backward DP runs from the far corner back to `imid-1` adding into the same
/// vectors. The column where the sum peaks is where the optimal path crosses,
/// and `jumpbacki`/`jumpbackj`/`jumpforwi`/`jumpforwj` say how far it runs
/// straight through -- so the split is a *segment*, `(jumpi,jumpj)` to
/// `(imid,jmid)`, not a single cell. The segment is emitted as gaps between
/// the two recursive halves.
///
/// `USE_PENALTY_EX` is 0 here, so neither running gap score is charged an
/// extension penalty, unlike `A__align`.
#[allow(clippy::too_many_arguments)]
fn msalignmm_rec(
    n_dynamicmtx: &[Vec<f64>],
    icyc: usize,
    jcyc: usize,
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    ist: usize,
    ien: usize,
    jst: usize,
    jen: usize,
    fulllen1: usize,
    fulllen2: usize,
    gi: &GapInfo,
    headgp: bool,
    tailgp: bool,
    headgapfreq1_g: f64,
    headgapfreq2_g: f64,
    newgapstr: u8,
) -> (f64, Vec<Vec<u8>>, Vec<Vec<u8>>, Vec<u8>, Vec<u8>) {
    let og1 = &gi.ogcp1[ist..];
    let fg1 = &gi.fgcp1[ist..];
    let og2 = &gi.ogcp2[jst..];
    let fg2 = &gi.fgcp2[jst..];
    let gf1 = &gi.gapfreq1[ist..];
    let gf2 = &gi.gapfreq2[jst..];
    let headgapfreq1 = if ist > 0 {
        gi.gapfreq1[ist - 1]
    } else {
        headgapfreq1_g
    };
    let headgapfreq2 = if jst > 0 {
        gi.gapfreq2[jst - 1]
    } else {
        headgapfreq2_g
    };

    let lgth1 = ien + 1 - ist;
    let lgth2 = (jen + 1).saturating_sub(jst);

    // `lgth2 <= 0`: nothing left on side 2, so side 1 runs against pure gaps.
    if jen + 1 <= jst {
        let mut m1: Vec<Vec<u8>> = Vec::with_capacity(icyc);
        for s in seq1.iter().take(icyc) {
            m1.push(s[ist..ist + lgth1].to_vec());
        }
        let m2: Vec<Vec<u8>> = vec![vec![newgapstr; lgth1]; jcyc];
        // `mgt2` is filled with '-' whatever `newgapstr` is; the C says so.
        return (0.0, m1, m2, vec![b'o'; lgth1], vec![b'-'; lgth1]);
    }

    let mut mseq1: Vec<Vec<u8>> = vec![Vec::new(); icyc];
    let mut mseq2: Vec<Vec<u8>> = vec![Vec::new(); jcyc];

    if lgth1 < DPTANNI || lgth2 < DPTANNI {
        let (value, gt1, gt2) = msalignmm_tanni(
            n_dynamicmtx,
            icyc,
            jcyc,
            seq1,
            seq2,
            cpmx1,
            cpmx2,
            ist,
            ien,
            jst,
            jen,
            fulllen1,
            fulllen2,
            &mut mseq1,
            &mut mseq2,
            gi,
            headgp,
            tailgp,
            headgapfreq1_g,
            headgapfreq2_g,
        );
        return (value, mseq1, mseq2, gt1, gt2);
    }

    let ll1 = lgth1 + 100;
    let ll2 = lgth2 + 100;
    let mut w1 = vec![0.0f64; ll2 + 2];
    let mut w2 = vec![0.0f64; ll2 + 2];
    let mut midw = vec![0.0f64; ll2 + 2];
    let mut midn = vec![0.0f64; ll2 + 2];
    let mut midm = vec![0.0f64; ll2 + 2];
    let mut jumpbacki = vec![0i32; ll2 + 2];
    let mut jumpbackj = vec![0i32; ll2 + 2];
    let mut jumpforwi = vec![0i32; ll2 + 2];
    let mut jumpforwj = vec![0i32; ll2 + 2];
    let mut initverticalw = vec![0.0f64; ll1 + 2];
    let mut lastverticalw = vec![0.0f64; ll1 + 2];
    let mut m = vec![0.0f64; ll2 + 2];
    let mut mp = vec![0i32; ll2 + 2];
    let nalphabets = n_dynamicmtx.len();
    let wide = ll1.max(ll2) + 2;
    let mut doublework = vec![vec![0.0f64; nalphabets]; wide];
    let mut intwork = vec![vec![0i32; nalphabets]; wide];

    let mut cur_is_w1 = true;
    match_calc_alphabet_seq(
        n_dynamicmtx,
        &mut initverticalw,
        cpmx2,
        cpmx1,
        jst,
        ist,
        lgth1,
        &mut doublework,
        &mut intwork,
        true,
    );
    match_calc_alphabet_seq(
        n_dynamicmtx,
        &mut w1,
        cpmx1,
        cpmx2,
        ist,
        jst,
        lgth2,
        &mut doublework,
        &mut intwork,
        true,
    );

    // Unlike `_tanni`, the recursive step always charges the head gaps.
    for i in 1..lgth1 + 1 {
        initverticalw[i] += og1[0] * headgapfreq2 + fg1[i - 1] * gf2[0];
    }
    for j in 1..lgth2 + 1 {
        w1[j] += og2[0] * headgapfreq1 + fg2[j - 1] * gf1[0];
    }
    for j in 1..lgth2 + 1 {
        m[j] = w1[j - 1] + og1[1] * gf2[j - 1];
        mp[j] = 0;
    }
    lastverticalw[0] = w1[lgth2 - 1];

    // `imid = lgth1 * 0.5` -- a double multiply truncated to int.
    let mut imid = (lgth1 as f64 * 0.5) as i32;
    let mut jumpi = 0i32;
    let mut jmid = 0i32;
    let mut jumpj = 0i32;

    for i in 1..=imid as usize {
        cur_is_w1 = !cur_is_w1;
        let (currentw, previousw) = if cur_is_w1 {
            (&mut w1, &mut w2)
        } else {
            (&mut w2, &mut w1)
        };
        previousw[0] = initverticalw[i - 1];
        match_calc_alphabet_seq(
            n_dynamicmtx,
            currentw,
            cpmx1,
            cpmx2,
            ist + i,
            jst,
            lgth2,
            &mut doublework,
            &mut intwork,
            false,
        );
        currentw[0] = initverticalw[i];

        m[0] = og1[i];
        if i as i32 == imid {
            midm[0] = m[0];
        }
        let mut mi = previousw[0] + og2[1] * gf1[i - 1];
        let mut mpi = 0i32;

        for j in 1..lgth2 + 1 {
            let mut wm = previousw[j - 1];
            let mut g = mi + fg2[j - 1] * gf1[i];
            if g > wm {
                wm = g;
            }
            g = previousw[j - 1] + og2[j] * gf1[i - 1];
            if g >= mi {
                mi = g;
                mpi = j as i32 - 1;
            }
            g = m[j] + fg1[i - 1] * gf2[j];
            if g > wm {
                wm = g;
            }
            g = previousw[j - 1] + og1[i] * gf2[j - 1];
            if g >= m[j] {
                m[j] = g;
                mp[j] = i as i32 - 1;
            }
            currentw[j] += wm;
            if i as i32 == imid {
                jumpbackj[j] = mp[j];
                jumpbacki[j] = mpi;
                midw[j] = currentw[j];
                midm[j] = m[j];
                midn[j] = mi;
            }
        }
        lastverticalw[i] = currentw[lgth2 - 1];
    }

    // ---- the backward DP ----
    match_calc_alphabet_seq(
        n_dynamicmtx,
        &mut initverticalw,
        cpmx2,
        cpmx1,
        jst + lgth2 - 1,
        ist,
        lgth1,
        &mut doublework,
        &mut intwork,
        true,
    );
    {
        let currentw = if cur_is_w1 { &mut w1 } else { &mut w2 };
        match_calc_alphabet_seq(
            n_dynamicmtx,
            currentw,
            cpmx1,
            cpmx2,
            ist + lgth1 - 1,
            jst,
            lgth2,
            &mut doublework,
            &mut intwork,
            true,
        );
        for j in 0..lgth2 - 1 {
            currentw[j] += fg2[lgth2 - 1] * gf1[lgth1] + og2[j + 1] * gf1[lgth1 - 1];
        }
        for j in (0..lgth2).rev() {
            m[j] = currentw[j + 1] + fg1[lgth1 - 2] * gf2[j + 1];
            mp[j] = lgth1 as i32 - 1;
        }
    }
    for i in 0..lgth1 - 1 {
        initverticalw[i] += fg1[lgth1 - 1] * gf2[lgth2] + og1[i + 1] * gf2[lgth2 - 1];
    }

    let mut firstm = -9999999.9f64;
    let mut firstmp = lgth1 as i32;
    let mut i = lgth1 as i32 - 2;
    while i > -1 {
        cur_is_w1 = !cur_is_w1;
        let (currentw, previousw) = if cur_is_w1 {
            (&mut w1, &mut w2)
        } else {
            (&mut w2, &mut w1)
        };
        let iu = i as usize;
        previousw[lgth2 - 1] = initverticalw[iu + 1];
        match_calc_alphabet_seq(
            n_dynamicmtx,
            currentw,
            cpmx1,
            cpmx2,
            ist + iu,
            jst,
            lgth2,
            &mut doublework,
            &mut intwork,
            false,
        );
        currentw[lgth2 - 1] = initverticalw[iu];

        let mut mi = previousw[lgth2 - 1] + fg2[lgth2 - 2] * gf1[iu + 1];
        let mut mpi = lgth2 as i32 - 1;

        for j in (0..lgth2 - 1).rev() {
            let mut wm = previousw[j + 1];
            let mut ijpi = i + 1;
            let mut ijpj = j as i32 + 1;

            let mut g = mi + og2[j + 1] * gf1[iu];
            if g > wm {
                wm = g;
                ijpj = mpi;
                ijpi = i + 1;
            }
            g = previousw[j + 1] + fg2[j] * gf1[iu + 1];
            if g >= mi {
                mi = g;
                mpi = j as i32 + 1;
            }
            g = m[j] + og1[iu + 1] * gf2[j];
            if g > wm {
                wm = g;
                ijpi = mp[j];
                ijpj = j as i32 + 1;
            }
            g = previousw[j + 1] + fg1[iu] * gf2[j + 1];
            if g >= m[j] {
                m[j] = g;
                mp[j] = i + 1;
            }

            if i == jumpi || i == imid - 1 {
                jumpforwi[j] = ijpi;
                jumpforwj[j] = ijpj;
            }
            if i == imid {
                midw[j] += wm;
                midm[j + 1] += m[j];
            }
            if i == imid - 1 {
                midn[j] += mi;
            }
            currentw[j] += wm;
        }
        // `g = *prept + fgcp1[i]` with `prept` left at index 0 by the loop.
        let g = previousw[0] + fg1[iu];
        if firstm < g {
            firstm = g;
            firstmp = i + 1;
        }
        if i == imid {
            midm[0] += firstm;
        }

        if i == imid - 1 {
            let mut maxwm = midw[1];
            jmid = 0;
            for j in 2..lgth2 - 1 {
                if midw[j] > maxwm {
                    jmid = j as i32;
                    maxwm = midw[j];
                }
            }
            for j in 0..lgth2 + 1 {
                if midm[j] > maxwm {
                    jmid = j as i32;
                    maxwm = midm[j];
                }
            }
            let mut wm = midw[jmid as usize];
            jumpi = imid - 1;
            jumpj = jmid - 1;
            if jmid > 0 && midn[jmid as usize - 1] > wm {
                jumpi = imid - 1;
                jumpj = jumpbacki[jmid as usize];
                wm = midn[jmid as usize - 1];
            }
            if midm[jmid as usize] > wm {
                jumpi = jumpbackj[jmid as usize];
                jumpj = jmid - 1;
            }
        }
        if i == jumpi {
            if jmid == 0 {
                jumpj = 0;
                jmid = 1;
                if imid < firstmp - 1 {
                    jumpi = firstmp;
                    imid = firstmp + 1;
                }
            } else if jmid >= lgth2 as i32 {
                jumpi = imid - 1;
                jmid = lgth2 as i32;
                jumpj = lgth2 as i32 - 1;
            } else {
                imid = jumpforwi[jumpj as usize];
                jmid = jumpforwj[jumpj as usize];
                if imid == jumpi {
                    jumpi = imid - 1;
                }
            }
            break;
        }
        i -= 1;
    }

    // ---- the two halves, and the straight segment between them ----
    let (mut value, mut a1, mut a2, mut g1, mut g2) = msalignmm_rec(
        n_dynamicmtx,
        icyc,
        jcyc,
        seq1,
        seq2,
        cpmx1,
        cpmx2,
        ist,
        ist + jumpi as usize,
        jst,
        jst + jumpj as usize,
        fulllen1,
        fulllen2,
        gi,
        headgp,
        tailgp,
        headgapfreq1_g,
        headgapfreq2_g,
        newgapstr,
    );

    let l = jmid - jumpj - 1;
    if l > 0 {
        let l = l as usize;
        for s in a1.iter_mut().take(icyc) {
            s.extend(std::iter::repeat_n(newgapstr, l));
        }
        for (j, s) in a2.iter_mut().enumerate().take(jcyc) {
            let from = jst + jumpj as usize + 1;
            s.extend_from_slice(&seq2[j][from..from + l]);
        }
        g1.extend(std::iter::repeat_n(b'-', l));
        g2.extend(std::iter::repeat_n(b'o', l));
        value += og2[jumpj as usize + 1] + fg2[jmid as usize - 1];
    }
    let l = imid - jumpi - 1;
    if l > 0 {
        let l = l as usize;
        for (i, s) in a1.iter_mut().enumerate().take(icyc) {
            let from = ist + jumpi as usize + 1;
            s.extend_from_slice(&seq1[i][from..from + l]);
        }
        for s in a2.iter_mut().take(jcyc) {
            s.extend(std::iter::repeat_n(newgapstr, l));
        }
        g1.extend(std::iter::repeat_n(b'o', l));
        g2.extend(std::iter::repeat_n(b'-', l));
        value += og1[jumpi as usize + 1] + fg1[imid as usize - 1];
    }

    let (v2, b1, b2, h1, h2) = msalignmm_rec(
        n_dynamicmtx,
        icyc,
        jcyc,
        seq1,
        seq2,
        cpmx1,
        cpmx2,
        ist + imid as usize,
        ien,
        jst + jmid as usize,
        jen,
        fulllen1,
        fulllen2,
        gi,
        headgp,
        tailgp,
        headgapfreq1_g,
        headgapfreq2_g,
        newgapstr,
    );
    value += v2;
    for (i, s) in a1.iter_mut().enumerate().take(icyc) {
        s.extend_from_slice(&b1[i]);
    }
    for (j, s) in a2.iter_mut().enumerate().take(jcyc) {
        s.extend_from_slice(&b2[j]);
    }
    g1.extend_from_slice(&h1);
    g2.extend_from_slice(&h2);

    (value, a1, a2, g1, g2)
}

/// `MSalignmm()`: the entry point the progressive stage calls when `alg == 'M'`.
///
/// Everything before the recursion is the same profile setup `A__align` does
/// -- composition matrix, opening and final gap counts, gap frequencies, and
/// the cached child profiles when the caller wants a result profile back. The
/// difference is all in what runs after: `MSalignmm_rec` instead of a single
/// quadratic DP.
///
/// Unlike `A__align` this keeps no `static` scratch, so nothing carries
/// between calls and there is no growth rule to reproduce.
#[allow(clippy::too_many_arguments)]
pub fn msalignmm(
    cn: &Constants,
    n_dynamicmtx: &[Vec<f64>],
    group1: &mut [Vec<u8>],
    group2: &mut [Vec<u8>],
    eff1: &[f64],
    eff2: &[f64],
    headgp: bool,
    tailgp: bool,
    gaps: Option<(&[u8], &[u8], &[u8], &[u8])>,
    cpmxchild0: Option<crate::align::CpmxHist>,
    cpmxchild1: Option<crate::align::CpmxHist>,
    orieff1: f64,
    orieff2: f64,
    want_cpmxresult: bool,
    newgapstr: u8,
) -> (f64, Option<crate::align::CpmxHist>) {
    let icyc = group1.len();
    let jcyc = group2.len();
    let lgth1 = group1[0].len();
    let lgth2 = group2[0].len();
    let nalphabets = cn.nalphabets;
    let fpenalty = cn.penalty as f64;

    let ll1 = lgth1 + 100;
    let ll2 = lgth2 + 100;
    let mut ogcp1 = vec![0.0f64; ll1 + 2];
    let mut ogcp2 = vec![0.0f64; ll2 + 2];
    let mut fgcp1 = vec![0.0f64; ll1 + 2];
    let mut fgcp2 = vec![0.0f64; ll2 + 2];
    let mut ogcp1o = vec![0.0f64; ll1 + 2];
    let mut ogcp2o = vec![0.0f64; ll2 + 2];
    let mut fgcp1o = vec![0.0f64; ll1 + 2];
    let mut fgcp2o = vec![0.0f64; ll2 + 2];
    // `ogcp1opt` and friends: what the DP is scaled from. They are the cached
    // child rows when there are any, and the local `*o` buffers otherwise --
    // and the C keeps the two apart, which matters at the result block.
    let mut ogcp1opt = vec![0.0f64; ll1 + 2];
    let mut ogcp2opt = vec![0.0f64; ll2 + 2];
    let mut fgcp1opt = vec![0.0f64; ll1 + 2];
    let mut fgcp2opt = vec![0.0f64; ll2 + 2];
    let mut cpmx1 = vec![vec![0.0f64; ll1 + 2]; nalphabets];
    let mut cpmx2 = vec![vec![0.0f64; ll2 + 2]; nalphabets];
    let mut gapfreq1 = vec![0.0f64; ll1 + 2];
    let mut gapfreq2 = vec![0.0f64; ll2 + 2];

    let g1: Vec<&[u8]> = group1.iter().map(|s| s.as_slice()).collect();
    let g2: Vec<&[u8]> = group2.iter().map(|s| s.as_slice()).collect();

    let headgapfreq1: f64;
    let headgapfreq2: f64;
    // The profiles the DP scores against: a child's cached one when there is
    // one, otherwise freshly computed here.
    let hist0 = cpmxchild0;
    let hist1 = cpmxchild1;

    if want_cpmxresult && cn.specificityconsideration == 0.0 {
        assert!(
            gaps.is_none(),
            "The combination of sgap1 and cpmxhist is not supported. See Salignmm.c"
        );
        if let Some(h) = hist0.as_ref() {
            for i in 0..nalphabets {
                cpmx1[i][..h[i].len()].copy_from_slice(&h[i]);
            }
            gapfreq1[..h[nalphabets].len()].copy_from_slice(&h[nalphabets]);
            // `ogcp1opt`/`fgcp1opt` point at the cached rows -- but `ogcp1o`
            // and `fgcp1o`, the local buffers, are deliberately left as the
            // zeros they were allocated with. See the note at the result
            // block below: the C hands *those* to `createogresult`.
            ogcp1opt[..h[nalphabets + 1].len()].copy_from_slice(&h[nalphabets + 1]);
            fgcp1opt[..h[nalphabets + 2].len()].copy_from_slice(&h[nalphabets + 2]);
        } else {
            crate::profile::cpmx_calc_new(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
            crate::profile::gapcountf(&mut gapfreq1, &g1, icyc, eff1, lgth1);
            for v in gapfreq1.iter_mut().take(lgth1) {
                *v = 1.0 - *v;
            }
            crate::profile::st_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1);
            crate::profile::st_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
            ogcp1opt.copy_from_slice(&ogcp1o);
            fgcp1opt.copy_from_slice(&fgcp1o);
        }
        if let Some(h) = hist1.as_ref() {
            for i in 0..nalphabets {
                cpmx2[i][..h[i].len()].copy_from_slice(&h[i]);
            }
            gapfreq2[..h[nalphabets].len()].copy_from_slice(&h[nalphabets]);
            ogcp2opt[..h[nalphabets + 1].len()].copy_from_slice(&h[nalphabets + 1]);
            fgcp2opt[..h[nalphabets + 2].len()].copy_from_slice(&h[nalphabets + 2]);
        } else {
            crate::profile::cpmx_calc_new(cn, &g2, &mut cpmx2, eff2, lgth2, jcyc);
            crate::profile::gapcountf(&mut gapfreq2, &g2, jcyc, eff2, lgth2);
            for v in gapfreq2.iter_mut().take(lgth2) {
                *v = 1.0 - *v;
            }
            crate::profile::st_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2);
            crate::profile::st_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
            ogcp2opt.copy_from_slice(&ogcp2o);
            fgcp2opt.copy_from_slice(&fgcp2o);
        }
        // With a cached profile the C assumes no head or tail gap.
        headgapfreq1 = 1.0;
        headgapfreq2 = 1.0;
        gapfreq1[lgth1] = 1.0;
        gapfreq2[lgth2] = 1.0;
    } else {
        // `ogcp1opt = ogcp1o` and friends (MSalignmm.c:2265-2267): on this
        // branch the two names alias one buffer, so everything filled below
        // reaches both the DP scaling and the result block.
        crate::profile::cpmx_calc_new(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
        crate::profile::cpmx_calc_new(cn, &g2, &mut cpmx2, eff2, lgth2, jcyc);
        let (mut hgf1, mut hgf2);
        if let Some((sgap1, sgap2, egap1, egap2)) = gaps {
            crate::profile::new_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1, sgap1);
            crate::profile::new_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
            crate::profile::new_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2, sgap2);
            crate::profile::new_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
            hgf1 = crate::profile::outgapcount(icyc, sgap1, eff1);
            hgf2 = crate::profile::outgapcount(jcyc, sgap2, eff2);
            gapfreq1[lgth1] = crate::profile::outgapcount(icyc, egap1, eff1);
            gapfreq2[lgth2] = crate::profile::outgapcount(jcyc, egap2, eff2);
        } else {
            crate::profile::st_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1);
            crate::profile::st_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
            crate::profile::st_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2);
            crate::profile::st_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
            hgf1 = 0.0;
            hgf2 = 0.0;
            gapfreq1[lgth1] = 0.0;
            gapfreq2[lgth2] = 0.0;
        }
        ogcp1opt.copy_from_slice(&ogcp1o);
        ogcp2opt.copy_from_slice(&ogcp2o);
        fgcp1opt.copy_from_slice(&fgcp1o);
        fgcp2opt.copy_from_slice(&fgcp2o);
        if !cn.legacygapcost {
            crate::profile::gapcountf(&mut gapfreq1, &g1, icyc, eff1, lgth1);
            crate::profile::gapcountf(&mut gapfreq2, &g2, jcyc, eff2, lgth2);
            for v in gapfreq1.iter_mut().take(lgth1 + 1) {
                *v = 1.0 - *v;
            }
            for v in gapfreq2.iter_mut().take(lgth2 + 1) {
                *v = 1.0 - *v;
            }
            hgf1 = 1.0 - hgf1;
            hgf2 = 1.0 - hgf2;
        } else {
            for v in gapfreq1.iter_mut().take(lgth1 + 1) {
                *v = 1.0;
            }
            for v in gapfreq2.iter_mut().take(lgth2 + 1) {
                *v = 1.0;
            }
            hgf1 = 1.0;
            hgf2 = 1.0;
        }
        headgapfreq1 = hgf1;
        headgapfreq2 = hgf2;
    }

    for i in 0..lgth1 {
        ogcp1[i] = 0.5 * (1.0 - ogcp1opt[i]) * fpenalty * gapfreq1[i];
        fgcp1[i] = 0.5 * (1.0 - fgcp1opt[i]) * fpenalty * gapfreq1[i];
    }
    for i in 0..lgth2 {
        ogcp2[i] = 0.5 * (1.0 - ogcp2opt[i]) * fpenalty * gapfreq2[i];
        fgcp2[i] = 0.5 * (1.0 - fgcp2opt[i]) * fpenalty * gapfreq2[i];
    }

    let gi = GapInfo {
        ogcp1,
        fgcp1,
        ogcp2,
        fgcp2,
        gapfreq1,
        gapfreq2,
    };

    let (wm, m1, m2, mgt1, mgt2) = msalignmm_rec(
        n_dynamicmtx,
        icyc,
        jcyc,
        group1,
        group2,
        &cpmx1,
        &cpmx2,
        0,
        lgth1 - 1,
        0,
        lgth2 - 1,
        lgth1,
        lgth2,
        &gi,
        headgp,
        tailgp,
        headgapfreq1,
        headgapfreq2,
        newgapstr,
    );

    // `if( icyc + jcyc > 20 )`: below that the parent recomputes its own
    // profile more cheaply than carrying this one.
    let mut result: Option<crate::align::CpmxHist> = None;
    if want_cpmxresult && icyc + jcyc > 20 {
        let totaleff1: f64 = eff1.iter().take(icyc).sum();
        let totaleff2: f64 = eff2.iter().take(jcyc).sum();
        if (totaleff1 - 1.0).abs() > 0.001 || (totaleff2 - 1.0).abs() > 0.001 {
            eprintln!(
                "Warning: rounding error may be large.  totaleff1 = {:50.40}",
                totaleff1
            );
            eprintln!(
                "Warning: rounding error may be large.  totaleff2 = {:50.40}",
                totaleff2
            );
            std::process::exit(1);
        }
        let te1 = totaleff1 * orieff1 / (orieff1 + orieff2);
        let te2 = totaleff2 * orieff2 / (orieff1 + orieff2);
        let limk = mgt1.len();
        let mut out =
            crate::align::create_cpmx_result(limk, te1, te2, &cpmx1, &cpmx2, &mgt1, &mgt2);
        out.push(crate::align::create_gapfreq_result(
            limk,
            te1,
            te2,
            &gi.gapfreq1,
            &gi.gapfreq2,
            &mgt1,
            &mgt2,
        ));
        // The C passes `ogcp1o`/`fgcp1o` here, not `ogcp1opt`/`fgcp1opt` --
        // Salignmm.c:2119 passes the `opt` pointers in the same place, and
        // MSalignmm.c:2459 does not. So whenever this node reused a child's
        // cached profile, the opening- and final-gap rows it hands *up* are
        // built from the zeros those local buffers still hold for that side.
        // Reproduced, not fixed: it changes the alignment (measured on 5 of
        // the 86 corpus inputs, wherever a node has more than 20 members).
        out.push(crate::align::create_og_result(
            limk,
            te1,
            te2,
            &ogcp1o,
            &ogcp2o,
            &gi.gapfreq1,
            &gi.gapfreq2,
            &mgt1,
            &mgt2,
        ));
        out.push(crate::align::create_fg_result(
            limk,
            te1,
            te2,
            &fgcp1o,
            &fgcp2o,
            &gi.gapfreq1,
            &gi.gapfreq2,
            &mgt1,
            &mgt2,
        ));
        result = Some(out);
    }

    for (i, s) in group1.iter_mut().enumerate().take(icyc) {
        *s = m1[i].clone();
    }
    for (j, s) in group2.iter_mut().enumerate().take(jcyc) {
        *s = m2[j].clone();
    }

    (wm, result)
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

    /// The memsave aligner must agree with `A__align` wherever the two
    /// recurrences coincide -- which is where no gap is ever extended, since
    /// `MSalignmm.c` compiles with `USE_PENALTY_EX 0`. With `ppenalty_ex` 0
    /// in the constants above, that is every input.
    #[test]
    fn agrees_with_a_align_when_no_extension_penalty_applies() {
        let cn = cn_prot();
        let dyn_mtx = cn.n_dis_consweight_multi.clone();
        let cases: [(&[u8], &[u8]); 4] = [
            (b"MNGTEGDNFYVPFSNKTGL", b"MNGTEGPNFYVPFSNITGV"),
            (b"ACDEFGHIKLMNPQRSTVWY", b"ACDEFGHIKLMNPQRSTVWY"),
            (b"AAAAAAAAAACCCCCCCCCC", b"CCCCCCCCCCAAAAAAAAAA"),
            (b"MKVLAAGIVGLNLGGKW", b"MKVLAAGIVGLLLGGKWMKV"),
        ];
        for (a, b) in cases {
            let mut g1 = vec![a.to_vec()];
            let mut g2 = vec![b.to_vec()];
            let (sm, _) = msalignmm(
                &cn,
                &dyn_mtx,
                &mut g1,
                &mut g2,
                &[1.0],
                &[1.0],
                true,
                true,
                None,
                None,
                None,
                1.0,
                1.0,
                false,
                b'-',
            );
            let mut h1 = vec![a.to_vec()];
            let mut h2 = vec![b.to_vec()];
            let mut st = crate::align::AlignState::new();
            let (sa, _) = st.a_align(
                &cn,
                &dyn_mtx,
                cn.penalty,
                cn.penalty_ex,
                &mut h1,
                &mut h2,
                &[1.0],
                &[1.0],
                true,
                true,
                None,
                None,
                -1,
                0,
                None,
                None,
                1.0,
                1.0,
                false,
                b'-',
            );
            assert_eq!(
                (g1[0].clone(), g2[0].clone()),
                (h1[0].clone(), h2[0].clone()),
                "columns differ for {:?}",
                String::from_utf8_lossy(a)
            );
            assert!((sm - sa).abs() < 1e-9, "score {} vs {}", sm, sa);
        }
    }

    #[test]
    fn output_rows_stay_the_same_length_and_keep_every_residue() {
        let cn = cn_prot();
        let dyn_mtx = cn.n_dis_consweight_multi.clone();
        let mut g1 = vec![
            b"MNGTEGDNFYVPFSNKTGLRSPYEYPQ".to_vec(),
            b"MNGTEGDNFYVPF--KTGLRSPYEYPQ".to_vec(),
        ];
        let mut g2 = vec![b"MNGTEGPNFYVPFSNITGVVRSPFEQPQYY".to_vec()];
        let before1: Vec<usize> = g1
            .iter()
            .map(|s| s.iter().filter(|&&c| c != b'-').count())
            .collect();
        let before2: Vec<usize> = g2
            .iter()
            .map(|s| s.iter().filter(|&&c| c != b'-').count())
            .collect();
        msalignmm(
            &cn,
            &dyn_mtx,
            &mut g1,
            &mut g2,
            &[0.5, 0.5],
            &[1.0],
            true,
            true,
            None,
            None,
            None,
            1.0,
            1.0,
            false,
            b'-',
        );
        let alen = g1[0].len();
        for s in g1.iter().chain(g2.iter()) {
            assert_eq!(s.len(), alen, "ragged output");
        }
        for (i, s) in g1.iter().enumerate() {
            assert_eq!(s.iter().filter(|&&c| c != b'-').count(), before1[i]);
        }
        for (j, s) in g2.iter().enumerate() {
            assert_eq!(s.iter().filter(|&&c| c != b'-').count(), before2[j]);
        }
    }

    /// Above `DPTANNI` on both sides the recursion actually splits; below it
    /// everything goes through the single-block DP. The answer must not
    /// depend on which side of that threshold the input falls, for an input
    /// where the alignment is forced.
    #[test]
    fn splitting_and_not_splitting_agree_on_an_unambiguous_alignment() {
        let cn = cn_prot();
        let dyn_mtx = cn.n_dis_consweight_multi.clone();
        let unit = b"MNGTEGDNFYVPFSNKTGLRSPYEYPQYYLEPWKYSLYMFFLILVGFPVNFLTLFVTVQH";
        for reps in [1usize, 3] {
            let seq: Vec<u8> = unit
                .iter()
                .cycle()
                .take(unit.len() * reps)
                .copied()
                .collect();
            let mut g1 = vec![seq.clone()];
            let mut g2 = vec![seq.clone()];
            msalignmm(
                &cn,
                &dyn_mtx,
                &mut g1,
                &mut g2,
                &[1.0],
                &[1.0],
                true,
                true,
                None,
                None,
                None,
                1.0,
                1.0,
                false,
                b'-',
            );
            // Two identical sequences align without a single gap, whether or
            // not the block was split.
            assert!(!g1[0].contains(&b'-'), "reps={} gapped", reps);
            assert_eq!(g1[0], seq);
            assert_eq!(g2[0], seq);
        }
        // And the longer case really did cross the threshold.
        assert!(unit.len() * 3 > DPTANNI);
        assert!(unit.len() < DPTANNI);
    }
}
