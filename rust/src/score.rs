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

/// The shared scoring loop of `naivepairscore11`.
///
/// `naivepairscorefast` had the same loop until it was fused into its own
/// column walk; this is still the independent second implementation that
/// `fast_and_plain_pair_scores_agree` checks the fused one against, and it
/// reaches the same answer without consulting a skip table at all.
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

/// `naivepairscorefast()`: the skip-table form, with the column pick and the
/// scoring fused into one pass.
///
/// The C runs this in two passes. `commongappickpairfast` walks the two
/// aligned rows, using the skip tables to jump whole gap runs, and copies the
/// surviving columns into two freshly `calloc`ed buffers; `naivepairscorefast`
/// then walks those buffers and adds up the score. This does both in the one
/// walk, so the buffers do not exist and the columns are read once.
///
/// **The additions, and their order, are untouched** -- which is the only
/// thing that may not move here. The second pass consumed the picked columns
/// strictly left to right, and so does this, with the same value added at the
/// same point of the sequence: `amino_dis[b1][b2]` for a residue pair, and
/// `penal` once per gap run. The f64 accumulator therefore sees an identical
/// chain of addends and the result is bit-identical, not merely close. This is
/// not a reassociation, and nothing here would be sound as one.
///
/// Carrying `run` across columns is what replaces the second pass's inner
/// `while( *p1 == '-' )` loops. That loop existed to charge a gap run once and
/// then walk to its end; charging only when the run *opens* is the same thing
/// expressed forward, and `run` is the only state it needed.
///
/// One ordering detail is load-bearing and is preserved: the C tests row 1 for
/// a gap **before** row 2, so a column gapped in both rows is charged as a
/// row-1 run. `commongappickpairfast` is not supposed to emit such a column
/// and the C has a disabled check saying so, but the two passes agreed on what
/// to do if one appeared, and so do this one and `pairscore`. The fusion
/// therefore does not rest on that invariant holding.
///
/// Worth 1.7G instructions of the 6.81G this function cost at 790 sequences:
/// the two scratch writes per column, and the whole of the second walk.
pub fn naivepairscorefast(
    cn: &Constants,
    seq1: &[u8],
    seq2: &[u8],
    skip1: &[i32],
    skip2: &[i32],
    penal: i32,
) -> f64 {
    // The two aligned rows are the same length by construction; the C walks
    // them with one pair of pointers that it advances together. Slicing both
    // to one hoisted length lets the indexing below drop its bounds checks.
    debug_assert_eq!(seq1.len(), seq2.len());
    let len = seq1.len().min(seq2.len());
    let i1 = &seq1[..len];
    let i2 = &seq2[..len];

    let penal = penal as f64;
    let mut vali = 0.0f64;
    // Which row the previous *kept* column was gapped in: 1, 2, or 0 for
    // neither. Starts at 0 so the first gap column opens a run. This is what
    // replaces the second pass's inner `while( *p1 == '-' )` loops: they
    // existed to charge a run once and then walk to its end, and charging only
    // where a run opens is the same thing expressed forward.
    let mut run = 0u8;

    // `p1` and `p2` in the C advance by the same amounts and are compared
    // against the same bound, so they are one cursor.
    let mut p = 0usize;
    let mut k1 = 0usize;
    let mut k2 = 0usize; // cursors into skip1 / skip2
    let mut skipped1 = 0i32;
    let mut skipped2 = 0i32;
    while p < len {
        let b1 = i1[p];
        let b2 = i2[p];
        // The C reads both skip tables at the top of every column, to ask how
        // far it may jump. That jump is non-zero *exactly* when both rows are
        // inside a gap run: a row sitting on a residue has already consumed
        // every gap recorded before that residue, so its candidate is 0 and
        // the minimum with it is 0. Both rows gapped is also exactly the
        // condition for dropping the column. So testing the two bytes first
        // and reading the tables only then makes the same jumps while taking
        // two loads, two bounds checks and a serial load-to-branch dependency
        // off the common path -- which is nearly every column, since a column
        // gapped in both rows is the minority case.
        //
        // That is where this function's time was. Dropping the tables
        // altogether and walking every column ran the 790-sequence build in
        // 708 ms against the C's 822 and this loop's 1123, so the per-column
        // table traffic cost far more than the columns it let us skip. But it
        // is O(alignment width) per pair where the table is O(columns kept),
        // and those diverge badly for two short rows in a wide alignment --
        // which is the case the table was added for. Hence neither: keep the
        // jump, pay for it only where it can fire.
        if b1 == b'-' && b2 == b'-' {
            let cand1 = skip1[k1] - skipped1;
            let cand2 = skip2[k2] - skipped2;
            let skip = if cand1 < cand2 { cand1 } else { cand2 };
            if skip > 0 {
                p += skip as usize;
                skipped1 += skip;
                skipped2 += skip;
                continue;
            }
            // `skip == 0` on a both-gap column means the tables disagree with
            // the rows they describe. The C cannot spin here because it
            // recomputes nothing: it falls through, emits the column, and its
            // disabled check at `mltaln9.c:13351` says it does not expect to
            // reach this at all. Falling through does the same and keeps the
            // cursor moving, rather than re-reading a table that will give the
            // same answer forever.
        }
        // Per-row bookkeeping, each row independently, exactly as the C.
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
        // This column is one the C would have copied into its scratch pair and
        // then read back out to score. Score it here instead, which is what
        // makes the second pass and the two buffers unnecessary.
        //
        // Row 1 is tested before row 2 because the C's scoring loop does, so a
        // column gapped in both rows is charged as a row-1 run. Such a column
        // should not reach here, but the C's two passes agreed on what to do
        // with one and this agrees with them, so the fusion does not rest on
        // that invariant.
        if b1 == b'-' {
            if run != 1 {
                vali += penal;
                run = 1;
            }
        } else if b2 == b'-' {
            if run != 2 {
                vali += penal;
                run = 2;
            }
        } else {
            vali += cn.amino_dis[b1 as usize][b2 as usize] as f64;
            run = 0;
        }
        p += 1;
    }
    vali
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
