//! Sum-of-pairs scoring over two aligned rows (`mltaln9.c`).
//!
//! `dndpre` turns an existing alignment into a distance matrix, and this is
//! how it scores a pair: drop the columns where *both* rows are gapped, then
//! walk what is left, charging `penal` once per gap run and the substitution
//! matrix for every residue pair.
//!
//! Two things here are easy to get wrong and are the reason the C's shape is
//! kept. The gap charge is **once per run, not per column** -- the loop
//! advances to the end of the run without adding again, which is why a long
//! insertion costs the same as a short one. And the run is defined by *one*
//! row: `while( *p1 == '-' )` looks only at `p1`, so a run in row 1 that is
//! interrupted by a column where row 2 is also gapped does not end there.

use crate::io::gappick0;

/// `makeskiptable()`: for each ungapped position of a row, how many gap
/// columns immediately precede it.
///
/// The last entry counts the trailing gaps, so the table has one more entry
/// than the row has residues.
///
/// Args:
///     seq: The aligned rows.
///
/// Returns:
///     One table per row.
pub fn make_skiptable(seq: &[Vec<u8>]) -> Vec<Vec<i32>> {
    let mut out = Vec::with_capacity(seq.len());
    for row in seq {
        let nogaplen = gappick0(row).len();
        let mut skip = vec![0i32; nogaplen + 1];
        let mut posinseq = 0usize;
        for &c in row {
            if c == b'-' {
                skip[posinseq] += 1;
            } else {
                posinseq += 1;
            }
        }
        out.push(skip);
    }
    out
}

/// `commongappickpair()`: drop every column gapped in both rows.
///
/// Returns:
///     The two rows with those columns removed.
fn commongappickpair(i1: &[u8], i2: &[u8]) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(i1.len());
    let mut r2 = Vec::with_capacity(i2.len());
    for (&a, &b) in i1.iter().zip(i2) {
        if a == b'-' && b == b'-' {
            continue;
        }
        r1.push(a);
        r2.push(b);
    }
    (r1, r2)
}

/// `commongappickpairfast()`: the same, using the skip tables to jump whole
/// runs of shared gaps instead of testing each column.
///
/// The jump length is the smaller of the two rows' remaining gap counts, which
/// is exactly the number of columns that are gapped in both.
fn commongappickpairfast(
    i1: &[u8],
    i2: &[u8],
    skip1: &[i32],
    skip2: &[i32],
) -> (Vec<u8>, Vec<u8>) {
    let mut r1 = Vec::with_capacity(i1.len());
    let mut r2 = Vec::with_capacity(i2.len());
    let (mut p1, mut p2) = (0usize, 0usize);
    let (mut s1, mut s2) = (0usize, 0usize);
    let (mut skipped1, mut skipped2) = (0i32, 0i32);
    loop {
        let a = skip1.get(s1).copied().unwrap_or(0) - skipped1;
        let b = skip2.get(s2).copied().unwrap_or(0) - skipped2;
        let skip = a.min(b).max(0) as usize;
        p1 += skip;
        p2 += skip;
        skipped1 += skip as i32;
        skipped2 += skip as i32;
        // The C tests `!*i1`, the NUL terminator; the end of the slice is the
        // same place.
        if p1 >= i1.len() {
            break;
        }
        if i1[p1] != b'-' {
            skipped1 = 0;
            s1 += 1;
        } else {
            skipped1 += 1;
        }
        if p2 < i2.len() && i2[p2] != b'-' {
            skipped2 = 0;
            s2 += 1;
        } else {
            skipped2 += 1;
        }
        r1.push(i1[p1]);
        r2.push(if p2 < i2.len() { i2[p2] } else { b'-' });
        p1 += 1;
        p2 += 1;
    }
    (r1, r2)
}

/// Walk two gap-picked rows, charging one `penal` per gap run.
fn score_picked(s1: &[u8], s2: &[u8], penal: i32, amino_dis: &[Vec<i32>]) -> f64 {
    let mut vali = 0.0f64;
    let mut k = 0usize;
    while k < s1.len() {
        if s1[k] == b'-' {
            vali += f64::from(penal);
            while k < s1.len() && s1[k] == b'-' {
                k += 1;
            }
            continue;
        }
        if k < s2.len() && s2[k] == b'-' {
            vali += f64::from(penal);
            while k < s2.len() && s2[k] == b'-' {
                k += 1;
            }
            continue;
        }
        let (a, b) = (s1[k] as usize, s2[k] as usize);
        vali += f64::from(amino_dis[a][b]);
        k += 1;
    }
    vali
}

/// `naivepairscore11()`: sum-of-pairs score of two aligned rows.
pub fn naive_pairscore11(
    seq1: &[u8],
    seq2: &[u8],
    penal: i32,
    amino_dis: &[Vec<i32>],
) -> f64 {
    let (s1, s2) = commongappickpair(seq1, seq2);
    score_picked(&s1, &s2, penal, amino_dis)
}

/// `naivepairscorefast()`: the same score, reached through the skip tables.
pub fn naive_pairscore_fast(
    seq1: &[u8],
    seq2: &[u8],
    skip1: &[i32],
    skip2: &[i32],
    penal: i32,
    amino_dis: &[Vec<i32>],
) -> f64 {
    let (s1, s2) = commongappickpairfast(seq1, seq2, skip1, skip2);
    score_picked(&s1, &s2, penal, amino_dis)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Identity matrix over the letters used in these tests: 1 for a match,
    /// 0 otherwise, so the arithmetic is checkable by hand. Shaped like
    /// `Constants::amino_dis`, which is what the real callers pass.
    fn ident() -> Vec<Vec<i32>> {
        let mut m = vec![vec![0i32; 0x80]; 0x80];
        for (i, row) in m.iter_mut().enumerate() {
            row[i] = 1;
        }
        m
    }

    #[test]
    fn skiptable_counts_the_gaps_before_each_residue() {
        let seq = vec![b"--AB-C".to_vec()];
        // two gaps before A, none before B, one before C, none trailing
        assert_eq!(make_skiptable(&seq)[0], vec![2, 0, 1, 0]);
    }

    #[test]
    fn skiptable_counts_trailing_gaps_in_its_last_entry() {
        let seq = vec![b"AB--".to_vec()];
        assert_eq!(make_skiptable(&seq)[0], vec![0, 0, 2]);
    }

    #[test]
    fn a_gap_run_is_charged_once_not_per_column() {
        let m = ident();
        // One run of three gaps in row 2 -> one penalty, plus two matches.
        let one_run = naive_pairscore11(b"ABCDE", b"AB---", -5, &m);
        // Two separate runs -> two penalties.
        let two_runs = naive_pairscore11(b"ABCDE", b"A-C-E", -5, &m);
        assert!((one_run - (2.0 - 5.0)).abs() < 1e-12, "got {one_run}");
        assert!((two_runs - (3.0 - 10.0)).abs() < 1e-12, "got {two_runs}");
    }

    #[test]
    fn columns_gapped_in_both_rows_are_dropped_before_scoring() {
        let m = ident();
        // The shared gap column must not become a charge on either row.
        assert_eq!(
            naive_pairscore11(b"A-B", b"A-B", -5, &m),
            naive_pairscore11(b"AB", b"AB", -5, &m)
        );
    }

    #[test]
    fn the_fast_path_agrees_with_the_plain_one() {
        let m = ident();
        let rows: Vec<Vec<u8>> = vec![
            b"--ABC--DE".to_vec(),
            b"AB--C-DE-".to_vec(),
            b"ABCDE----".to_vec(),
            b"---------".to_vec(),
        ];
        let skip = make_skiptable(&rows);
        for i in 0..rows.len() {
            for j in 0..rows.len() {
                let slow = naive_pairscore11(&rows[i], &rows[j], -5, &m);
                let fast =
                    naive_pairscore_fast(&rows[i], &rows[j], &skip[i], &skip[j], -5, &m);
                assert!(
                    (slow - fast).abs() < 1e-12,
                    "rows {i},{j}: slow {slow} fast {fast}"
                );
            }
        }
    }
}
