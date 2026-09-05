//! Position-specific constraints ("localhom"), the mechanism that turns
//! L-INS-i's all-pairs local alignments into a bias on the progressive stage.
//!
//! Transliterated from `putlocalhom2` (`core/io.c`), `calcimportance_half`,
//! `movereg` and `fillimp` (`core/mltaln9.c`), and `fastshrinklocalhom_half`
//! (`core/tddis.c`).
//!
//! Each pairwise local alignment is cut at its gaps into ungapped blocks; each
//! block records the residue ranges it spans in both sequences and a score
//! (`opt`). `calcimportance_half` then reweights every block by how many other
//! sequences agree with it, and during the progressive alignment `fillimp`
//! paints those weights onto a matrix that `A__align` adds to its match
//! scores.

use crate::constants::Constants;

/// One ungapped block of a pairwise local alignment (`LocalHom` in mltaln.h).
///
/// The C keeps these as a linked list whose head always exists, with
/// `opt == -1.0` meaning "no constraint here"; the same convention is used
/// here, with the head at index 0 of the vector.
#[derive(Clone, Debug)]
pub struct LocalHom {
    pub start1: i32,
    pub end1: i32,
    pub start2: i32,
    pub end2: i32,
    pub opt: f64,
    pub importance: f64,
    pub rimportance: f64,
    pub overlapaa: i32,
    pub korh: u8,
}

impl Default for LocalHom {
    fn default() -> Self {
        // Matches the initialisation loop in tbfast.c's main().
        LocalHom {
            start1: -1,
            end1: -1,
            start2: -1,
            end2: -1,
            opt: -1.0,
            importance: -1.0,
            rimportance: -1.0,
            overlapaa: 0,
            korh: b'h',
        }
    }
}

/// `localhomtable[i][j-i]` for `j >= i`: the constraint chain for one pair.
pub type LocalHomTable = Vec<Vec<Vec<LocalHom>>>;

/// Allocate the half table the way tbfast.c's main() does.
pub fn new_table(njob: usize) -> LocalHomTable {
    (0..njob)
        .map(|i| (0..njob - i).map(|_| vec![LocalHom::default()]).collect())
        .collect()
}

/// `putlocalhom2()`: split an aligned pair into ungapped blocks.
///
/// `divpairscore` is 0 in every configuration this port accepts, so the score
/// of every block is the alignment-wide average rather than per block.
pub fn putlocalhom2(
    cn: &Constants,
    al1: &[u8],
    al2: &[u8],
    chain: &mut Vec<LocalHom>,
    off1: usize,
    off2: usize,
    korh: u8,
) {
    let mut pos1 = off1 as i32;
    let mut pos2 = off2 as i32;
    let mut isumscore: i64 = 0;
    let mut sumoverlap: i32 = 0;
    let mut start1 = 0i32;
    let mut start2 = 0i32;
    let mut st = 0;
    let mut iscore: i64 = 0;
    let mut nlocalhom = 0usize;

    // The head element already exists; further blocks are appended.
    let push = |chain: &mut Vec<LocalHom>, nlocalhom: &mut usize| -> usize {
        let idx = if *nlocalhom > 0 {
            chain.push(LocalHom::default());
            chain.len() - 1
        } else {
            0
        };
        *nlocalhom += 1;
        idx
    };

    let mut p = 0usize;
    while p < al1.len() {
        if st == 1 && (al1[p] == b'-' || al2[p] == b'-') {
            let end1 = pos1 - 1;
            let end2 = pos2 - 1;
            let idx = push(chain, &mut nlocalhom);
            let e = &mut chain[idx];
            e.start1 = start1;
            e.start2 = start2;
            e.end1 = end1;
            e.end2 = end2;
            e.korh = korh;
            isumscore += iscore;
            sumoverlap += end2 - start2 + 1;
            iscore = 0;
            st = 0;
        } else if al1[p] != b'-' && al2[p] != b'-' {
            if st == 0 {
                start1 = pos1;
                start2 = pos2;
                st = 1;
            }
            iscore += cn.n_dis[cn.amino_n[al1[p] as usize] as usize]
                [cn.amino_n[al2[p] as usize] as usize] as i64;
        }
        if al1[p] != b'-' {
            pos1 += 1;
        }
        if al2[p] != b'-' {
            pos2 += 1;
        }
        p += 1;
    }

    if !al1.is_empty() && al1[al1.len() - 1] != b'-' && al2[al2.len() - 1] != b'-' {
        let idx = push(chain, &mut nlocalhom);
        let end1 = pos1 - 1;
        let end2 = pos2 - 1;
        let e = &mut chain[idx];
        e.start1 = start1;
        e.start2 = start2;
        e.end1 = end1;
        e.end2 = end2;
        e.korh = korh;
        isumscore += iscore;
        sumoverlap += end2 - start2 + 1;
    }

    for e in chain.iter_mut() {
        e.overlapaa = sumoverlap;
        e.opt = if sumoverlap > 0 {
            isumscore as f64 * 5.8 / (600.0 * sumoverlap as f64)
        } else {
            -1.0
        };
    }
}

/// `seqlen()`: length ignoring gaps.
fn seqlen(s: &[u8]) -> usize {
    s.iter().filter(|&&c| c != b'-').count()
}

/// `calcimportance_half()`: weight each block by how much of the rest of the
/// data set supports the same region.
pub fn calcimportance_half(nseq: usize, eff: &[f64], seq: &[Vec<u8>], table: &mut LocalHomTable) {
    let alloclen = seq.iter().map(|s| s.len()).max().unwrap_or(0) + 2;
    let mut ieff = vec![0.0f64; nseq];
    let mut totaleff = 0.0f64;
    for i in 0..nseq {
        let nogaplen = seqlen(&seq[i]);
        ieff[i] = if nogaplen == 0 { 0.0 } else { eff[i] };
        totaleff += ieff[i];
    }
    for v in ieff.iter_mut() {
        *v /= totaleff;
    }

    let mut importance = vec![0.0f64; alloclen];
    for i in 0..nseq {
        for v in importance.iter_mut() {
            *v = 0.0;
        }
        for j in 0..nseq {
            if i == j {
                continue;
            }
            if i < j {
                for e in table[i][j - i].iter() {
                    if e.opt == -1.0 {
                        continue;
                    }
                    for pos in e.start1..=e.end1 {
                        importance[pos as usize] += ieff[j];
                    }
                }
            } else {
                for e in table[j][i - j].iter() {
                    if e.opt == -1.0 {
                        continue;
                    }
                    for pos in e.start2..=e.end2 {
                        importance[pos as usize] += ieff[j];
                    }
                }
            }
        }
        for j in 0..nseq {
            if i == j {
                continue;
            }
            if i < j {
                if table[i][j - i][0].opt == -1.0 {
                    continue;
                }
                for e in table[i][j - i].iter_mut() {
                    if e.opt == -1.0 {
                        continue;
                    }
                    let mut tmp = 0.0f64;
                    let mut len = 0i32;
                    for pos in e.start1..=e.end1 {
                        tmp += importance[pos as usize];
                        len += 1;
                    }
                    tmp /= len as f64;
                    e.importance = tmp * e.opt;
                }
            } else {
                if table[j][i - j][0].opt == -1.0 {
                    continue;
                }
                for e in table[j][i - j].iter_mut() {
                    if e.opt == -1.0 {
                        continue;
                    }
                    let mut tmp = 0.0f64;
                    let mut len = 0i32;
                    for pos in e.start2..=e.end2 {
                        tmp += importance[pos as usize];
                        len += 1;
                    }
                    tmp /= len as f64;
                    e.rimportance = tmp * e.opt;
                }
            }
        }
    }

    // Average the two directions.
    for i in 0..nseq - 1 {
        for j in i + 1..nseq {
            for e in table[i][j - i].iter_mut() {
                if e.opt == -1.0 {
                    continue;
                }
                let imp = 0.5 * (e.importance + e.rimportance);
                e.importance = imp;
                e.rimportance = imp;
            }
        }
    }
}

/// Map one ungapped residue range onto the current (gapped) alignment columns.
fn map_range(s: &[u8], from: i32, to: i32) -> (usize, usize) {
    let mut tmpint: i32 = -1;
    let mut p = 0usize;
    while p < s.len() {
        if s[p] != b'-' {
            tmpint += 1;
        }
        p += 1;
        if tmpint == from {
            break;
        }
    }
    let start = p - 1;
    if from == to {
        return (start, start);
    }
    // The C keeps scanning from where the first loop stopped.
    while p < s.len() {
        if s[p] != b'-' {
            tmpint += 1;
        }
        p += 1;
        if tmpint == to {
            break;
        }
    }
    (start, p - 1)
}

/// `movereg()` / `movereg_swap()`.
///
/// A constraint is stored once, indexed by the lower-numbered sequence. When
/// the member of group 1 has the *higher* original index the two halves of the
/// block have to be read the other way round -- that is `movereg_swap`, and
/// `fillimp` chooses between them per pair. Getting this wrong shifts the
/// constraint bonus by a few columns and silently changes the alignment.
fn movereg(seq1: &[u8], seq2: &[u8], e: &LocalHom, swap: bool) -> (usize, usize, usize, usize) {
    if swap {
        let (start1, end1) = map_range(seq1, e.start2, e.end2);
        let (start2, end2) = map_range(seq2, e.start1, e.end1);
        (start1, start2, end1, end2)
    } else {
        let (start1, end1) = map_range(seq1, e.start1, e.end1);
        let (start2, end2) = map_range(seq2, e.start2, e.end2);
        (start1, start2, end1, end2)
    }
}

/// `fillimp()`: paint the constraint weights onto the match-bonus matrix.
///
/// `shrink[k1][k2]` is the chain for the pair (member k1 of group 1, member k2
/// of group 2), or `None` when there is no constraint for that pair.
#[allow(clippy::too_many_arguments)]
pub fn fillimp(
    impmtx: &mut [Vec<f64>],
    clus1: usize,
    clus2: usize,
    lgth1: usize,
    lgth2: usize,
    seq1: &[&[u8]],
    seq2: &[&[u8]],
    eff1: &[f64],
    eff2: &[f64],
    shrink: &[Vec<Option<&Vec<LocalHom>>>],
    // Original sequence indices of each group's members; they decide which
    // way round a constraint's two halves are read.
    orinum1: &[i32],
    orinum2: &[i32],
    fastathreshold: f64,
) {
    // `fill` lowers to a memset; the element-by-element loop it replaces did
    // not, and this matrix is re-zeroed once per branch of every refinement
    // round. Same values written, so nothing downstream can tell.
    for row in impmtx.iter_mut().take(lgth1) {
        row[..lgth2].fill(0.0);
    }
    let effijx = fastathreshold;
    for i in 0..clus1 {
        for j in 0..clus2 {
            let effij = eff1[i] * eff2[j] * effijx;
            let chain = match shrink[i][j] {
                Some(c) => c,
                None => continue,
            };
            let swap = orinum1[i] > orinum2[j];
            for e in chain.iter() {
                let (start1, start2, end1, end2) = movereg(seq1[i], seq2[j], e, swap);
                let mut k1 = start1;
                let mut k2 = start2;
                // `korh` is always 'h' here, so the kozo weight (zero without
                // structural input) is never selected.
                while k1 < seq1[i].len() && k2 < seq2[j].len() {
                    let c1 = seq1[i][k1];
                    let c2 = seq2[j][k2];
                    if c1 != b'-' && c2 != b'-' {
                        impmtx[k1][k2] += e.importance * effij;
                        k1 += 1;
                        k2 += 1;
                    } else if c1 != b'-' && c2 == b'-' {
                        k2 += 1;
                    } else if c1 == b'-' && c2 != b'-' {
                        k1 += 1;
                    } else {
                        k1 += 1;
                        k2 += 1;
                    }
                    if k1 > end1 || k2 > end2 {
                        break;
                    }
                }
            }
        }
    }
}

/// `fastshrinklocalhom_half()`: the chains for one pair of groups.
pub fn fastshrinklocalhom_half<'a>(
    mem1: &[i32],
    mem2: &[i32],
    table: &'a LocalHomTable,
) -> Vec<Vec<Option<&'a Vec<LocalHom>>>> {
    let mut out = Vec::with_capacity(mem1.len());
    for &a in mem1 {
        let mut row = Vec::with_capacity(mem2.len());
        for &b in mem2 {
            let (lo, hi) = if a < b { (a, b) } else { (b, a) };
            let chain = &table[lo as usize][(hi - lo) as usize];
            row.push(if chain[0].opt == -1.0 {
                None
            } else {
                Some(chain)
            });
        }
        out.push(row);
    }
    out
}
