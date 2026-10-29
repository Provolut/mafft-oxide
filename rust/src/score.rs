//! Sum-of-pairs scores used to rebuild the distance matrix from an alignment.
//!
//! Transliterated from `naivepairscore11`, `naivepairscorefast`,
//! `commongappickpair`, `commongappickpairfast` and `makeskiptable` in
//! `core/mltaln9.c`.

use crate::constants::Constants;
use std::cell::RefCell;

/// `commongappickpair()`: drop columns that are gaps in both sequences.
fn commongappickpair(i1: &[u8], i2: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(i1.len());
    let mut r2 = Vec::with_capacity(i1.len());
    for k in 0..i1.len() {
        if i1[k] == b'-' && i2[k] == b'-' {
            continue;
        }
        r1.push(i1[k]);
        r2.push(i2[k]);
    }
    (r1, r2)
}

/// `makeskiptable()`: `skip[p]` is the run of gaps just before residue `p`.
pub fn makeskiptable(seq: &[u8]) -> Vec<i32> {
    let nogaplen = seq.iter().filter(|&&c| c != b'-').count();
    let mut skip = vec![0i32; nogaplen + 1];
    let mut posinseq = 0usize;
    for &c in seq {
        if c == b'-' {
            skip[posinseq] += 1;
        } else {
            posinseq += 1;
        }
    }
    skip
}

/// `commongappickpairfast()`: the skip-table accelerated form of the above.
///
/// Writes into `r1`/`r2` rather than returning two fresh `Vec`s, and returns
/// how many columns it wrote. `naivepairscorefast` calls this once per
/// sequence pair -- 311,655 times on a 790-sequence input -- and the two
/// allocations cost 1.53G instructions in `Vec::push` plus 0.42G in allocator
/// plumbing, measured with callgrind against the C's 5.89G for this function
/// and the scoring loop together. The C `calloc`s and `free`s two buffers per
/// call; this reuses one pair per thread.
///
/// The bytes written and the order they are written in are unchanged, so the
/// f64 sum `pairscore` then computes over them is bit-identical.
fn commongappickpairfast(
    i1: &[u8],
    i2: &[u8],
    skip1: &[i32],
    skip2: &[i32],
    r1: &mut Vec<u8>,
    r2: &mut Vec<u8>,
) -> usize {
    // The two aligned rows are the same length by construction; the C walks
    // them with one pair of pointers that it advances together. Slicing both
    // to one hoisted length lets the indexing below drop its bounds checks.
    debug_assert_eq!(i1.len(), i2.len());
    let len = i1.len().min(i2.len());
    let i1 = &i1[..len];
    let i2 = &i2[..len];
    if r1.len() < len {
        r1.resize(len, 0);
    }
    if r2.len() < len {
        r2.resize(len, 0);
    }
    let o1 = &mut r1[..len];
    let o2 = &mut r2[..len];

    let mut n = 0usize;
    // `p1` and `p2` in the C advance by the same amounts and are compared
    // against the same bound, so they are one cursor.
    let mut p = 0usize;
    let mut k1 = 0usize;
    let mut k2 = 0usize; // cursors into skip1 / skip2
    let mut skipped1 = 0i32;
    let mut skipped2 = 0i32;
    loop {
        let cand1 = skip1[k1] - skipped1;
        let cand2 = skip2[k2] - skipped2;
        let skip = if cand1 < cand2 { cand1 } else { cand2 };
        p += skip as usize;
        skipped1 += skip;
        skipped2 += skip;
        if p >= len {
            break;
        }
        // Read before the two updates below, which cannot touch `i1`/`i2`.
        let b1 = i1[p];
        let b2 = i2[p];
        if b1 != b'-' {
            skipped1 = 0;
            k1 += 1;
        } else {
            skipped1 += 1;
        }
        if b2 != b'-' {
            skipped2 = 0;
            k2 += 1;
        } else {
            skipped2 += 1;
        }
        o1[n] = b1;
        o2[n] = b2;
        n += 1;
        p += 1;
    }
    n
}

thread_local! {
    /// The scratch pair `commongappickpairfast` writes into, grown
    /// monotonically and never shrunk -- the same shape as `commonip`'s
    /// matrix pool, and for the same reason.
    static PICKED: RefCell<(Vec<u8>, Vec<u8>)> = const { RefCell::new((Vec::new(), Vec::new())) };
}

/// The shared scoring loop of `naivepairscore11` / `naivepairscorefast`.
fn pairscore(cn: &Constants, s1: &[u8], s2: &[u8], penal: i32) -> f64 {
    // Both callers build the two rows in lockstep, so they are the same
    // length and the C's `while( *p1 )` and its inner `while( *p2 == '-' )`
    // run to the same bound. Hoisting it lets the indexing below drop its
    // bounds checks; the additions, and their order, are untouched.
    debug_assert_eq!(s1.len(), s2.len());
    let n = s1.len().min(s2.len());
    let s1 = &s1[..n];
    let s2 = &s2[..n];
    let mut vali = 0.0f64;
    let mut p = 0usize;
    while p < n {
        if s1[p] == b'-' {
            vali += penal as f64;
            while p < n && s1[p] == b'-' {
                p += 1;
            }
            continue;
        }
        if s2[p] == b'-' {
            vali += penal as f64;
            while p < n && s2[p] == b'-' {
                p += 1;
            }
            continue;
        }
        vali += cn.amino_dis[s1[p] as usize][s2[p] as usize] as f64;
        p += 1;
    }
    vali
}

/// `naivepairscore11()`.
pub fn naivepairscore11(cn: &Constants, seq1: &[u8], seq2: &[u8], penal: i32) -> f64 {
    let (s1, s2) = commongappickpair(seq1, seq2);
    pairscore(cn, &s1, &s2, penal)
}

/// `naivepairscorefast()`.
pub fn naivepairscorefast(
    cn: &Constants,
    seq1: &[u8],
    seq2: &[u8],
    skip1: &[i32],
    skip2: &[i32],
    penal: i32,
) -> f64 {
    PICKED.with(|cell| {
        let mut cell = cell.borrow_mut();
        let (r1, r2) = &mut *cell;
        let n = commongappickpairfast(seq1, seq2, skip1, skip2, r1, r2);
        pairscore(cn, &r1[..n], &r2[..n], penal)
    })
}

/// `intergroup_score()`: the weighted sum-of-pairs score between two aligned
/// groups, with a flat gap-open charge and no extension term.
///
/// Only `dooneiteration` uses it, to decide whether a re-alignment improved
/// on what was there. The gap handling is the load-bearing part: a run of
/// gaps costs `penalty` once, and the run is charged residue by residue
/// against the *gap* character, not the residue that opened it.
pub fn intergroup_score(
    cn: &Constants,
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
    eff1: &[f64],
    eff2: &[f64],
    clus1: usize,
    clus2: usize,
    len: usize,
) -> f64 {
    // The C indexes past the end into the NUL terminator on purpose.
    fn at(s: &[u8], k: usize) -> u8 {
        if k < s.len() {
            s[k]
        } else {
            0
        }
    }
    let m = &cn.amino_dis_consweight_multi;
    let penalty = cn.penalty as f64;
    // `len2 = len - 2` as a plain int: for len < 2 the C's comparison is
    // against a negative number, so use i64 here rather than usize.
    let len2 = len as i64 - 2;

    let mut value = 0.0f64;
    for i in 0..clus1 {
        for j in 0..clus2 {
            let efficient = eff1[i] * eff2[j];
            let mseq1 = &seq1[i];
            let mseq2 = &seq2[j];
            let mut tmpscore = 0.0f64;
            let mut k = 0usize;
            while k < len {
                let ms1 = at(mseq1, k);
                let ms2 = at(mseq2, k);
                if ms1 == b'-' && ms2 == b'-' {
                    k += 1;
                    continue;
                }
                tmpscore += m[ms1 as usize][ms2 as usize];

                if ms1 == b'-' {
                    tmpscore += penalty;
                    tmpscore += m[ms1 as usize][ms2 as usize];
                    loop {
                        k += 1;
                        let c = at(mseq1, k);
                        if c != b'-' {
                            break;
                        }
                        tmpscore += m[c as usize][ms2 as usize];
                    }
                    k -= 1;
                    if k as i64 > len2 {
                        break;
                    }
                    k += 1;
                    continue;
                }
                if ms2 == b'-' {
                    tmpscore += penalty;
                    tmpscore += m[ms1 as usize][ms2 as usize];
                    loop {
                        k += 1;
                        let c = at(mseq2, k);
                        if c != b'-' {
                            break;
                        }
                        tmpscore += m[ms1 as usize][c as usize];
                    }
                    k -= 1;
                    if k as i64 > len2 {
                        break;
                    }
                    k += 1;
                    continue;
                }
                k += 1;
            }
            value += tmpscore * efficient;
        }
    }
    value
}
