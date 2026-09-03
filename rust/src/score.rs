//! Sum-of-pairs scores used to rebuild the distance matrix from an alignment.
//!
//! Transliterated from `naivepairscore11`, `naivepairscorefast`,
//! `commongappickpair`, `commongappickpairfast` and `makeskiptable` in
//! `core/mltaln9.c`.

use crate::constants::Constants;

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
fn commongappickpairfast(i1: &[u8], i2: &[u8], skip1: &[i32], skip2: &[i32]) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(i1.len());
    let mut r2 = Vec::with_capacity(i1.len());
    let mut p1 = 0usize;
    let mut p2 = 0usize; // cursors into i1 / i2 (they advance together)
    let mut k1 = 0usize;
    let mut k2 = 0usize; // cursors into skip1 / skip2
    let mut skipped1 = 0i32;
    let mut skipped2 = 0i32;
    loop {
        let cand1 = skip1[k1] - skipped1;
        let cand2 = skip2[k2] - skipped2;
        let skip = if cand1 < cand2 { cand1 } else { cand2 };
        p1 += skip as usize;
        p2 += skip as usize;
        skipped1 += skip;
        skipped2 += skip;
        if p1 >= i1.len() {
            break;
        }
        if i1[p1] != b'-' {
            skipped1 = 0;
            k1 += 1;
        } else {
            skipped1 += 1;
        }
        if i2[p2] != b'-' {
            skipped2 = 0;
            k2 += 1;
        } else {
            skipped2 += 1;
        }
        r1.push(i1[p1]);
        r2.push(i2[p2]);
        p1 += 1;
        p2 += 1;
    }
    (r1, r2)
}

/// The shared scoring loop of `naivepairscore11` / `naivepairscorefast`.
fn pairscore(cn: &Constants, s1: &[u8], s2: &[u8], penal: i32) -> f64 {
    let mut vali = 0.0f64;
    let mut p = 0usize;
    while p < s1.len() {
        if s1[p] == b'-' {
            vali += penal as f64;
            while p < s1.len() && s1[p] == b'-' {
                p += 1;
            }
            continue;
        }
        if s2[p] == b'-' {
            vali += penal as f64;
            while p < s2.len() && s2[p] == b'-' {
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
    let (s1, s2) = commongappickpairfast(seq1, seq2, skip1, skip2);
    pairscore(cn, &s1, &s2, penal)
}
