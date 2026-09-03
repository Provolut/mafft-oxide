//! The 6-mer ("ktuples") distance matrix, transliterated from `disttbfast.c`
//! (`seq_grp`, `makepointtable`, `makecompositiontable_p`) and `mltaln9.c`
//! (`commonsextet_p`).

use crate::constants::Constants;

pub const END_OF_VEC: i32 = -1;

// PLENFAC* from disttbfast.c (the `#else` branch that is actually compiled).
const PLENFACA: f64 = 0.01;
const PLENFACB: f64 = 10000.0;
const PLENFACC: f64 = 10000.0;
const PLENFACD: f64 = 0.1;

/// `seq_grp()`: map residues to their 6 groups, dropping unscored characters.
///
/// Returns the group vector terminated by [`END_OF_VEC`], or a vector whose
/// first element is `-1` when fewer than 6 residues survived (the C sentinel
/// for "too short to have a 6-mer").
pub fn seq_grp(cn: &Constants, seq: &[u8], nunknown: &mut i64) -> Vec<i32> {
    let mut grp: Vec<i32> = Vec::with_capacity(seq.len() + 1);
    for &c in seq {
        let tmp = cn.amino_grp[c as usize];
        if tmp < 6 {
            grp.push(tmp);
        } else {
            *nunknown += 1;
        }
    }
    let n = grp.len();
    grp.push(END_OF_VEC);
    if n < 6 {
        grp[0] = -1;
    }
    grp
}

/// `makepointtable()`: rolling base-6 code of every 6-mer.
pub fn makepointtable(grp: &[i32]) -> Vec<i32> {
    if grp[0] == -1 {
        return vec![-1];
    }
    let mut pointt = Vec::with_capacity(grp.len());
    let mut p = 0usize; // C's `p`, trailing pointer
    let mut n = 0usize; // C's `n`, leading pointer
    let mut point = grp[n] * 7776;
    n += 1;
    point += grp[n] * 1296;
    n += 1;
    point += grp[n] * 216;
    n += 1;
    point += grp[n] * 36;
    n += 1;
    point += grp[n] * 6;
    n += 1;
    point += grp[n];
    n += 1;
    pointt.push(point);

    while grp[n] != END_OF_VEC {
        point -= grp[p] * 7776;
        p += 1;
        point *= 6;
        point += grp[n];
        n += 1;
        pointt.push(point);
    }
    pointt.push(END_OF_VEC);
    pointt
}

/// `makecompositiontable_p()`: 6-mer multiset of one sequence.
pub fn makecompositiontable_p(table: &mut [i32], pointt: &[i32]) {
    for &point in pointt {
        if point == END_OF_VEC {
            break;
        }
        table[point as usize] += 1;
    }
}

/// State behind the `static` scratch arrays of `commonsextet_p()`.
pub struct Sextet {
    memo: Vec<i32>,
    ct: Vec<i32>,
}

impl Sextet {
    pub fn new(tsize: usize, maxl: usize) -> Sextet {
        Sextet {
            memo: vec![0; tsize],
            ct: vec![0; maxl.min(tsize) + 1],
        }
    }

    /// `commonsextet_p()`: size of the 6-mer multiset intersection.
    pub fn common(&mut self, table: &[i32], pointt: &[i32]) -> i32 {
        if pointt[0] == -1 {
            return 0;
        }
        let mut value = 0i32;
        let mut cp = 0usize;
        for &point in pointt {
            if point == END_OF_VEC {
                break;
            }
            let idx = point as usize;
            let tmp = self.memo[idx];
            self.memo[idx] += 1;
            if tmp < table[idx] {
                value += 1;
            }
            if tmp == 0 {
                self.ct[cp] = point;
                cp += 1;
            }
        }
        for k in 0..cp {
            self.memo[self.ct[k] as usize] = 0;
        }
        value
    }
}

/// The whole "Making a distance matrix .." block of `disttbfast()`.
///
/// Returns the half matrix `mtx[i][j-i]` for `j >= i`, plus `nogaplen`.
pub fn distance_matrix(
    cn: &Constants,
    seq: &[Vec<u8>],
    njob: usize,
) -> (Vec<Vec<f64>>, Vec<usize>) {
    let tsize = 6usize.pow(6); // (int)pow( 6, 6 ) for amino acids
    let mut nogaplen = vec![0usize; njob];
    let mut pointt: Vec<Vec<i32>> = Vec::with_capacity(njob);
    let mut maxl = 0usize;
    let mut nunknown = 0i64;

    for i in 0..njob {
        let tmpseq = crate::io::gappick0(&seq[i]);
        nogaplen[i] = tmpseq.len();
        if nogaplen[i] > maxl {
            maxl = nogaplen[i];
        }
        let grpseq = seq_grp(cn, &tmpseq, &mut nunknown);
        pointt.push(makepointtable(&grpseq));
    }

    let mut sextet = Sextet::new(tsize, maxl);
    let mut mtx: Vec<Vec<f64>> = (0..njob).map(|i| vec![0.0; njob - i]).collect();
    let mut table1 = vec![0i32; tsize];
    for i in 0..njob {
        table1.iter_mut().for_each(|v| *v = 0);
        makecompositiontable_p(&mut table1, &pointt[i]);
        for j in i..njob {
            mtx[i][j - i] = sextet.common(&table1, &pointt[j]) as f64;
        }
    }

    let ien = njob - 1;
    for i in 0..ien {
        for j in i + 1..njob {
            let (longer, shorter) = if nogaplen[i] > nogaplen[j] {
                (nogaplen[i] as f64, nogaplen[j] as f64)
            } else {
                (nogaplen[j] as f64, nogaplen[i] as f64)
            };
            let lenfac =
                1.0 / (shorter / longer * PLENFACD + PLENFACB / (longer + PLENFACC) + PLENFACA);
            let bunbo = mtx[i][0].min(mtx[j][0]);
            mtx[i][j - i] = if bunbo == 0.0 {
                2.0
            } else {
                (1.0 - mtx[i][j - i] / bunbo) * lenfac * 2.0
            };
        }
    }

    (mtx, nogaplen)
}
