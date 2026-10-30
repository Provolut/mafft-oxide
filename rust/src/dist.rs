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
// D6LENFAC*, the nucleotide equivalents.
const D6LENFACA: f64 = 0.01;
const D6LENFACB: f64 = 2500.0;
const D6LENFACC: f64 = 2500.0;
const D6LENFACD: f64 = 0.1;
// D10LENFAC*, for `--tuplesize 10`.
const D10LENFACA: f64 = 0.01;
const D10LENFACB: f64 = 1000000.0;
const D10LENFACC: f64 = 1000000.0;
const D10LENFACD: f64 = 0.0;

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

/// `seq_grp_nuc()`: map bases to their 4 groups, dropping everything else.
pub fn seq_grp_nuc(cn: &Constants, seq: &[u8], nunknown: &mut i64, tuplesize: usize) -> Vec<i32> {
    let mut grp: Vec<i32> = Vec::with_capacity(seq.len() + 1);
    for &c in seq {
        let tmp = cn.amino_grp[c as usize];
        if tmp < 4 {
            grp.push(tmp);
        } else {
            *nunknown += 1;
        }
    }
    let n = grp.len();
    grp.push(END_OF_VEC);
    // Unlike `seq_grp`, the nucleotide version compares against `tuplesize`.
    if n < tuplesize {
        grp[0] = -1;
    }
    grp
}

/// `makepointtable_nuc()`: rolling base-4 code of every 6-mer.
pub fn makepointtable_nuc(grp: &[i32]) -> Vec<i32> {
    if grp[0] == -1 {
        return vec![-1];
    }
    let mut pointt = Vec::with_capacity(grp.len());
    let mut p = 0usize;
    let mut n = 0usize;
    let mut point = grp[n] * 1024;
    n += 1;
    point += grp[n] * 256;
    n += 1;
    point += grp[n] * 64;
    n += 1;
    point += grp[n] * 16;
    n += 1;
    point += grp[n] * 4;
    n += 1;
    point += grp[n];
    n += 1;
    pointt.push(point);

    while grp[n] != END_OF_VEC {
        point -= grp[p] * 1024;
        p += 1;
        point *= 4;
        point += grp[n];
        n += 1;
        pointt.push(point);
    }
    pointt.push(END_OF_VEC);
    pointt
}

/// `makepointtable_nuc_dectet()`: rolling base-4 code of every 10-mer.
pub fn makepointtable_nuc_dectet(grp: &[i32]) -> Vec<i32> {
    if grp[0] == -1 {
        return vec![-1];
    }
    let mut pointt = Vec::with_capacity(grp.len());
    let mut p = 0usize;
    let mut n = 0usize;
    let mut point = 0i32;
    for w in [262144, 65536, 16384, 4096, 1024, 256, 64, 16, 4, 1] {
        point += grp[n] * w;
        n += 1;
    }
    pointt.push(point);

    while grp[n] != END_OF_VEC {
        point -= grp[p] * 262144;
        p += 1;
        point *= 4;
        point += grp[n];
        n += 1;
        pointt.push(point);
    }
    pointt.push(END_OF_VEC);
    pointt
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
/// The k-mer tables and length factors that every distance pass starts from:
/// the "Making a distance matrix .." prologue of `disttbfast()`, shared by the
/// full-matrix builder below and by the compact ones in `compacttree.rs`.
pub struct Tuples {
    pub pointt: Vec<Vec<i32>>,
    pub nogaplen: Vec<usize>,
    pub tsize: usize,
    pub maxl: usize,
    pub lenfaca: f64,
    pub lenfacb: f64,
    pub lenfacc: f64,
    pub lenfacd: f64,
}

pub fn build_tuples(cn: &Constants, seq: &[Vec<u8>], njob: usize, tuplesize: usize) -> Tuples {
    let nuc = cn.dorp == b'd';
    // (int)pow( 4, tuplesize ) for nucleotides, (int)pow( 6, 6 ) for amino acids
    let tsize = if nuc {
        4usize.pow(tuplesize as u32)
    } else {
        6usize.pow(6)
    };
    let (lenfaca, lenfacb, lenfacc, lenfacd) = if nuc && tuplesize == 10 {
        (D10LENFACA, D10LENFACB, D10LENFACC, D10LENFACD)
    } else if nuc {
        (D6LENFACA, D6LENFACB, D6LENFACC, D6LENFACD)
    } else {
        (PLENFACA, PLENFACB, PLENFACC, PLENFACD)
    };
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
        if nuc {
            let grpseq = seq_grp_nuc(cn, &tmpseq, &mut nunknown, tuplesize);
            pointt.push(if tuplesize == 10 {
                makepointtable_nuc_dectet(&grpseq)
            } else {
                makepointtable_nuc(&grpseq)
            });
        } else {
            let grpseq = seq_grp(cn, &tmpseq, &mut nunknown);
            pointt.push(makepointtable(&grpseq));
        }
    }
    Tuples {
        pointt,
        nogaplen,
        tsize,
        maxl,
        lenfaca,
        lenfacb,
        lenfacc,
        lenfacd,
    }
}

/// `nthreads` is what `-C` asked for; 1 or less runs inline. Every cell is a
/// pure function of its own `(i, j)`, so the matrix does not depend on it --
/// see `crate::parallel`.
pub fn distance_matrix(
    cn: &Constants,
    seq: &[Vec<u8>],
    njob: usize,
    tuplesize: usize,
    nthreads: usize,
) -> (Vec<Vec<f64>>, Vec<usize>) {
    let Tuples {
        pointt,
        nogaplen,
        tsize,
        maxl,
        lenfaca,
        lenfacb,
        lenfacc,
        lenfacd,
    } = build_tuples(cn, seq, njob, tuplesize);

    let mut mtx: Vec<Vec<f64>> = (0..njob).map(|i| vec![0.0; njob - i]).collect();

    // The `Sextet` and the composition table are per-worker scratch, which is
    // what the C allocates per thread as well. `table1` is rebuilt for every
    // row regardless, so sharing one per worker costs nothing.
    crate::parallel::fill_rows(
        &mut mtx,
        nthreads,
        || (Sextet::new(tsize, maxl), vec![0i32; tsize]),
        |(sextet, table1), i, row| {
            table1.iter_mut().for_each(|v| *v = 0);
            makecompositiontable_p(table1, &pointt[i]);
            for j in i..njob {
                row[j - i] = sextet.common(table1, &pointt[j]) as f64;
            }
        },
    );

    // The second pass reads the diagonal of row `j` while writing row `i`, so
    // it cannot borrow the matrix row-wise with those reads still in it. The
    // diagonals are finished and never written again here -- this pass only
    // touches `j > i`, and `mtx[i][0]` is the `j == i` cell -- so lifting them
    // out first is not a reordering, it is the same values read from a copy.
    let diag: Vec<f64> = mtx.iter().map(|row| row[0]).collect();
    crate::parallel::fill_rows(
        &mut mtx,
        nthreads,
        || (),
        |_, i, row| {
            if i + 1 >= njob {
                return;
            }
            for j in i + 1..njob {
                let (longer, shorter) = if nogaplen[i] > nogaplen[j] {
                    (nogaplen[i] as f64, nogaplen[j] as f64)
                } else {
                    (nogaplen[j] as f64, nogaplen[i] as f64)
                };
                let lenfac =
                    1.0 / (shorter / longer * lenfacd + lenfacb / (longer + lenfacc) + lenfaca);
                let bunbo = diag[i].min(diag[j]);
                row[j - i] = if bunbo == 0.0 {
                    2.0
                } else {
                    (1.0 - row[j - i] / bunbo) * lenfac * 2.0
                };
            }
        },
    );

    (mtx, nogaplen)
}
