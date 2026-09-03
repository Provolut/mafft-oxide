//! Profile statistics shared by the group-to-group DP.
//!
//! Transliterated from `cpmx_calc_new` (`core/tddis.c`) and
//! `st_OpeningGapCount` / `st_FinalGapCount` / `gapcountf` / `gapireru`
//! (`core/mltaln9.c`).
//!
//! The C also has incremental variants (`cpmx_calc_add`, `st_OpeningGapAdd`,
//! `st_FinalGapAdd`, `gapcountadd`).  They are deliberately absent here: they
//! live on the branch of `A__align` that runs only when no profile-history slot
//! is supplied, and `treebase()` always supplies one, so nothing in the ported
//! path can reach them.

use crate::constants::Constants;

/// `cpmx_calc_new()`: weighted residue composition, `cpmx[letter][column]`.
///
/// Only the first `lgth` columns are touched, exactly as in the C, so any
/// trailing values from a previous call survive in the reused buffer.
pub fn cpmx_calc_new(
    cn: &Constants,
    seq: &[&[u8]],
    cpmx: &mut [Vec<f64>],
    eff: &[f64],
    lgth: usize,
    clus: usize,
) {
    for row in cpmx.iter_mut().take(cn.nalphabets) {
        for v in row.iter_mut().take(lgth) {
            *v = 0.0;
        }
    }
    for k in 0..clus {
        let feff = eff[k];
        let s = seq[k];
        for j in 0..lgth {
            cpmx[cn.amino_n[s[j] as usize] as usize][j] += feff;
        }
    }
}

/// `st_OpeningGapCount()`: weighted frequency of a gap *opening* at each column.
pub fn st_opening_gap_count(ogcp: &mut [f64], clus: usize, seq: &[&[u8]], eff: &[f64], len: usize) {
    for v in ogcp.iter_mut().take(len) {
        *v = 0.0;
    }
    for j in 0..clus {
        let feff = eff[j];
        let s = seq[j];
        let mut gc = false;
        for i in 0..len {
            let gb = gc;
            gc = s[i] == b'-';
            if !gb && gc {
                ogcp[i] += feff;
            }
        }
    }
    ogcp[len] = 0.0;
}

/// `st_FinalGapCount()`: weighted frequency of a gap *closing* at each column.
pub fn st_final_gap_count(fgcp: &mut [f64], clus: usize, seq: &[&[u8]], eff: &[f64], len: usize) {
    for v in fgcp.iter_mut().take(len) {
        *v = 0.0;
    }
    for j in 0..clus {
        let feff = eff[j];
        let s = seq[j];
        // C walks one past each position: gc looks at s[i+1], with s[len] == 0.
        let mut gc = s[0] == b'-';
        for i in 0..len {
            let gb = gc;
            gc = i + 1 < len && s[i + 1] == b'-';
            if gb && !gc {
                fgcp[i] += feff;
            }
        }
    }
}

/// `gapcountf()`: weighted gap frequency per column.
pub fn gapcountf(freq: &mut [f64], seq: &[&[u8]], nseq: usize, eff: &[f64], lgth: usize) {
    for i in 0..lgth {
        let mut fr = 0.0;
        for j in 0..nseq {
            if seq[j][i] == b'-' {
                fr += eff[j];
            }
        }
        freq[i] = fr;
    }
}

/// `new_OpeningGapCount()`: like [`st_opening_gap_count`], but the run state
/// starts from the boundary pattern just before the block, and `ogcp[len]` is
/// deliberately left alone (the C only zeroes `0..len`).
pub fn new_opening_gap_count(
    ogcp: &mut [f64],
    clus: usize,
    seq: &[&[u8]],
    eff: &[f64],
    len: usize,
    sgappat: &[u8],
) {
    for v in ogcp.iter_mut().take(len) {
        *v = 0.0;
    }
    for j in 0..clus {
        let feff = eff[j];
        let s = seq[j];
        let mut gc = sgappat[j] == b'-';
        for i in 0..len {
            let gb = gc;
            gc = s[i] == b'-';
            if !gb && gc {
                ogcp[i] += feff;
            }
        }
    }
}

/// `new_FinalGapCount()`.
///
/// The C takes an `egappat` argument but cannot use it: the extra step after
/// its loop has already read the NUL terminator into `gb`, so the test is
/// always false. The compiled function is therefore identical to
/// [`st_final_gap_count`], and the parameter is dropped here rather than
/// pretending it has an effect.
pub fn new_final_gap_count(fgcp: &mut [f64], clus: usize, seq: &[&[u8]], eff: &[f64], len: usize) {
    st_final_gap_count(fgcp, clus, seq, eff, len)
}

/// `outgapcount()`: weighted frequency of `-` in a boundary pattern.
pub fn outgapcount(nseq: usize, gappat: &[u8], eff: &[f64]) -> f64 {
    let mut fr = 0.0f64;
    for j in 0..nseq {
        if gappat[j] == b'-' {
            fr += eff[j];
        }
    }
    fr
}

/// `cpmx_calc_add()`: incremental update when one member was appended.
///
/// Reachable only under `--allowshift`, where the profile-history cache is
/// off but consecutive tree steps can still grow the same group by one.
pub fn cpmx_calc_add(
    cn: &Constants,
    seq: &[&[u8]],
    cpmx: &mut [Vec<f64>],
    eff: &[f64],
    lgth: usize,
    clus: usize,
) {
    let newmem = clus - 1;
    let neweff = eff[clus - 1];
    let orieff = 1.0 - neweff;
    for j in 0..lgth {
        for row in cpmx.iter_mut().take(cn.nalphabets) {
            row[j] *= orieff;
        }
        cpmx[cn.amino_n[seq[newmem][j] as usize] as usize][j] += neweff;
    }
}

/// `st_OpeningGapAdd()`: incremental form of [`st_opening_gap_count`].
pub fn st_opening_gap_add(ogcp: &mut [f64], clus: usize, seq: &[&[u8]], eff: &[f64], len: usize) {
    let j = clus - 1;
    let feff = eff[j];
    let orieff = 1.0 - feff;
    let s = seq[j];
    let mut gc = false;
    for i in 0..len {
        let gb = gc;
        gc = s[i] == b'-';
        ogcp[i] *= orieff;
        if !gb && gc {
            ogcp[i] += feff;
        }
    }
    ogcp[len] = 0.0;
}

/// `st_FinalGapAdd()`: incremental form of [`st_final_gap_count`].
pub fn st_final_gap_add(fgcp: &mut [f64], clus: usize, seq: &[&[u8]], eff: &[f64], len: usize) {
    let j = clus - 1;
    let feff = eff[j];
    let orieff = 1.0 - feff;
    let s = seq[j];
    let mut gc = s[0] == b'-';
    for i in 0..len {
        let gb = gc;
        gc = i + 1 < len && s[i + 1] == b'-';
        fgcp[i] *= orieff;
        if gb && !gc {
            fgcp[i] += feff;
        }
    }
}

/// `gapcountadd()`: incremental form of [`gapcountf`], operating on the
/// already-complemented values (the C flips them back, updates, and leaves
/// them un-flipped for the caller to flip again).
pub fn gapcountadd(freq: &mut [f64], seq: &[&[u8]], nseq: usize, eff: &[f64], lgth: usize) {
    let j = nseq - 1;
    let newfr = eff[j];
    let orifr = 1.0 - newfr;
    for i in 0..lgth {
        freq[i] = 1.0 - freq[i];
        freq[i] *= orifr;
        if seq[j][i] == b'-' {
            freq[i] += newfr;
        }
    }
}

/// `gapireru()`: expand `ori` according to a gap table of `o`/`-`.
pub fn gapireru(ori: &[u8], gt: &[u8], gapchar: u8) -> Vec<u8> {
    let mut res = Vec::with_capacity(gt.len());
    let mut p = 0usize;
    for &g in gt {
        if g == b'-' {
            res.push(gapchar);
        } else {
            res.push(ori[p]);
            p += 1;
        }
    }
    res
}
