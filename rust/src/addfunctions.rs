//! Gap bookkeeping for `--add`, transliterated from `core/addfunctions.c`.
//!
//! Adding sequences to an existing alignment means the alignment's own columns
//! must survive. The C does that by collapsing each group's all-gap columns
//! before the DP, remembering where they were (`gapmap`), and putting them back
//! afterwards -- while any column the DP newly inserted is marked `'='` rather
//! than `'-'` so it can be told apart and propagated to the sequences that were
//! not part of this step.

use crate::align::AlignState;
use crate::constants::Constants;

/// `commongappick()`: drop every column that is a gap in all of `seq`.
///
/// The C's loop runs to `<= len`, so it also examines the NUL terminator; that
/// position is never all-gap, so it is copied and becomes the new terminator.
pub fn commongappick(seq: &mut [Vec<u8>]) {
    if seq.is_empty() {
        return;
    }
    let len = seq[0].len();
    let mut keep: Vec<usize> = Vec::with_capacity(len + 1);
    for i in 0..=len {
        let all_gap = seq.iter().all(|s| s.get(i).copied().unwrap_or(0) == b'-');
        if !all_gap {
            keep.push(i);
        }
    }
    for s in seq.iter_mut() {
        let old = s.clone();
        for (i, &p) in keep.iter().enumerate() {
            let c = old.get(p).copied().unwrap_or(0);
            if i < s.len() {
                s[i] = c;
            }
        }
        // The C leaves the buffer long and relies on the copied NUL; we shorten.
        let n = keep.iter().take_while(|&&p| p < len).count();
        s.truncate(n);
    }
}

/// `findcommongaps()`: how many all-gap columns sit before each surviving
/// column. Length `len+1`; the last entry counts the trailing run.
pub fn findcommongaps(seq: &[Vec<u8>]) -> Vec<i32> {
    let len = seq[0].len();
    let mut gapmap = vec![0i32; len + 1];
    let mut pos = 0usize;
    for i in 0..len {
        if seq.iter().all(|s| s[i] == b'-') {
            gapmap[pos] += 1;
        } else {
            pos += 1;
        }
    }
    gapmap
}

/// `findnewgaps()`: how many `'='` columns -- the ones the DP just inserted --
/// sit before each real column of `seq[rep]`.
pub fn findnewgaps(seq: &[Vec<u8>], rep: usize) -> Vec<i32> {
    let len = seq[0].len();
    let mut gaplen = vec![0i32; len + 1];
    let mut pos = 0usize;
    for i in 0..len {
        if seq[rep][i] == b'=' {
            gaplen[pos] += 1;
        } else {
            pos += 1;
        }
    }
    gaplen
}

/// `adjustgapmap()`: re-index `gapmap` onto the widened sequence, giving the
/// newly inserted `'='` columns a count of zero.
pub fn adjustgapmap(newlen: usize, gapmap: &mut Vec<i32>, seq: &[u8]) {
    let mut tmpmap = vec![0i32; newlen + 2];
    let mut j = 0usize;
    let mut pos = 0usize;
    for &c in seq.iter() {
        if c == b'=' {
            tmpmap[j] = 0;
        } else {
            tmpmap[j] = gapmap[pos];
            pos += 1;
        }
        j += 1;
    }
    tmpmap[j] = gapmap[pos];
    if gapmap.len() < newlen + 1 {
        gapmap.resize(newlen + 1, 0);
    }
    gapmap[..newlen + 1].copy_from_slice(&tmpmap[..newlen + 1]);
}

/// `isallgap()`.
pub fn isallgap(s: &[u8]) -> bool {
    s.iter().all(|&c| c == b'-')
}

/// `restorecommongaps()`: put the all-gap columns `findcommongaps` removed back
/// into the members of the two groups just aligned.
///
/// The C's inner loop runs to `len+1`, so it copies the NUL terminator too;
/// any gaps `gapmap` records after the last column land before it.
pub fn restorecommongaps(seq: &mut [Vec<u8>], mem: &[i32], gapmap: &mut Vec<i32>, gapchar: u8) {
    if mem.is_empty() {
        return;
    }
    let rep = mem[0] as usize;
    let len = seq[rep].len();
    for &i in mem.iter() {
        let iu = i as usize;
        let mut out: Vec<u8> = Vec::with_capacity(len + 1);
        for j in 0..len + 1 {
            for _ in 0..gapmap[j] {
                out.push(gapchar);
            }
            if j < len {
                out.push(seq[iu][j]);
            }
        }
        seq[iu] = out;
    }

    // The map is widened in place to match: each entry keeps its count and is
    // followed by that many zeros, one per column just reinserted.
    let mut tmp: Vec<i32> = Vec::with_capacity(gapmap.len() * 2);
    for j in 0..len + 1 {
        tmp.push(gapmap[j]);
        tmp.resize(tmp.len() + gapmap[j] as usize, 0);
    }
    if gapmap.len() < tmp.len() {
        gapmap.resize(tmp.len(), 0);
    }
    gapmap[..tmp.len()].copy_from_slice(&tmp);
}

/// `profilealignment()`: realign the three groups around one block of newly
/// inserted columns. A thin wrapper over `A__align` with flat weights.
///
/// `aln1` is the all-gap group; it is rebuilt afterwards, with `'='` written
/// wherever `aln0` is all-gap too, which is how the marker propagates.
#[allow(clippy::too_many_arguments)]
pub fn profilealignment(
    cn: &Constants,
    ast: &mut AlignState,
    aln0: &mut [Vec<u8>],
    aln1: &mut [Vec<u8>],
    aln2: &mut [Vec<u8>],
) {
    let n0 = aln0.len();
    let n1 = aln1.len();
    let n2 = aln2.len();
    // "--allowshift no tokiha..." -- both single columns, nothing to do.
    if aln0[0].len() <= 1 && aln2[0].len() <= 1 {
        return;
    }
    commongappick(aln0);
    commongappick(aln2);
    if aln2[0].is_empty() {
        let newlen = aln0[0].len();
        for s in aln2.iter_mut() {
            *s = vec![b'-'; newlen];
        }
        return;
    }

    // Weights: flat over the members that are not entirely gap.
    let mk = |aln: &[Vec<u8>]| -> Vec<f64> {
        let alcount = aln.iter().filter(|s| !isallgap(s)).count();
        aln.iter()
            .map(|s| {
                if isallgap(s) {
                    0.0
                } else {
                    1.0 / alcount as f64
                }
            })
            .collect()
    };
    let effarr0 = mk(aln0);
    let effarr2 = mk(aln2);

    ast.a_align(
        cn,
        &cn.n_dis_consweight_multi,
        cn.penalty,
        cn.penalty_ex,
        aln0,
        aln2,
        &effarr0,
        &effarr2,
        true,
        true,
        None,
        None,
        -1,
        -1,
        None,
        None,
        0.0,
        0.0,
        false,
        b'-',
    );

    let newlen = aln0[0].len();
    for s in aln1.iter_mut() {
        *s = vec![b'-'; newlen];
    }
    for j in 0..newlen {
        let allgap0 = (0..n0).all(|i| aln0[i][j] == b'-');
        if allgap0 {
            let allgap1 = (0..n1).all(|i| aln1[i][j] == b'-');
            if allgap1 {
                for i in 0..n1 {
                    aln1[i][j] = b'=';
                }
            }
        }
    }
    let _ = n2;
}

/// `countnogaplen()` / `countgapmap()`: how many leading zeros a slice has.
fn count_leading_zeros(v: &[i32]) -> usize {
    v.iter().take_while(|&&x| x == 0).count()
}

/// `insertnewgaps()`: propagate the columns the last step inserted into every
/// sequence that was not part of it.
///
/// Sequences are split three ways: group 1 and group 2 are the two sides just
/// aligned, group 0 is everything already aligned but not involved, and
/// anything not yet aligned at all is excluded.
#[allow(clippy::too_many_arguments)]
pub fn insertnewgaps(
    cn: &Constants,
    ast: &mut AlignState,
    njob: usize,
    alreadyaligned: &[bool],
    seq: &mut [Vec<u8>],
    ex1: &[i32],
    ex2: &[i32],
    gaplen: &[i32],
    gapmap: &[i32],
    gapchar: u8,
) {
    let mut mar = vec![0u8; njob];
    for i in 0..njob {
        if !alreadyaligned[i] {
            mar[i] = 3;
        }
    }
    for &k in ex1.iter() {
        mar[k as usize] = 1;
    }
    for &k in ex2.iter() {
        mar[k as usize] = 2;
    }
    let list =
        |want: u8, mar: &[u8]| -> Vec<usize> { (0..njob).filter(|&i| mar[i] == want).collect() };
    let list2 = list(2, &mar);
    let list1 = list(1, &mar);
    let list0 = list(0, &mar);
    if list0.is_empty() {
        return;
    }
    let rep = list0[0];
    let len = seq[rep].len();
    let len0 = len + 1;

    let mut aseq: Vec<Vec<u8>> = vec![Vec::new(); njob];
    let mut posin12 = 0usize;
    let mut j = 0usize;
    while j < len0 {
        if gaplen[j] != 0 {
            let mut mseq0: Vec<Vec<u8>> = vec![Vec::new(); list0.len()];
            let mut mseq1: Vec<Vec<u8>> = vec![Vec::new(); list1.len()];
            let mut mseq2: Vec<Vec<u8>> = vec![Vec::new(); list2.len()];

            let gapshift = gaplen[j] as usize;
            for s in mseq0.iter_mut() {
                s.extend(std::iter::repeat_n(gapchar, gapshift));
            }
            for (i, s) in mseq1.iter_mut().enumerate() {
                s.extend(take(&seq[list1[i]], posin12, gapshift));
            }
            for (i, s) in mseq2.iter_mut().enumerate() {
                s.extend(take(&seq[list2[i]], posin12, gapshift));
            }
            posin12 += gapshift;

            let gapshift = gapmap[posin12] as usize;
            for (i, s) in mseq0.iter_mut().enumerate() {
                s.extend(take(&seq[list0[i]], j, gapshift));
            }
            for (i, s) in mseq1.iter_mut().enumerate() {
                s.extend(take(&seq[list1[i]], posin12, gapshift));
            }
            for (i, s) in mseq2.iter_mut().enumerate() {
                s.extend(take(&seq[list2[i]], posin12, gapshift));
            }

            if gapshift != 0 {
                profilealignment(cn, ast, &mut mseq0, &mut mseq1, &mut mseq2);
            }
            j += gapshift;
            posin12 += gapshift;

            for (i, &lp) in list0.iter().enumerate() {
                let m = std::mem::take(&mut mseq0[i]);
                aseq[lp].extend(m);
            }
            for (i, &lp) in list1.iter().enumerate() {
                let m = std::mem::take(&mut mseq1[i]);
                aseq[lp].extend(m);
            }
            for (i, &lp) in list2.iter().enumerate() {
                let m = std::mem::take(&mut mseq2[i]);
                aseq[lp].extend(m);
            }
        }
        let blocklen = 1 + count_leading_zeros(&gaplen[(j + 1).min(len0)..len0]);
        for &lp in list0.iter() {
            let t = take(&seq[lp], j, blocklen);
            aseq[lp].extend(t);
        }
        for &lp in list1.iter() {
            let t = take(&seq[lp], posin12, blocklen);
            aseq[lp].extend(t);
        }
        for &lp in list2.iter() {
            let t = take(&seq[lp], posin12, blocklen);
            aseq[lp].extend(t);
        }
        j += blocklen;
        posin12 += blocklen;
    }

    for &lp in list0.iter().chain(list1.iter()).chain(list2.iter()) {
        seq[lp] = std::mem::take(&mut aseq[lp]);
    }
}

/// `strncpy0()`: copy at most `n` bytes, stopping at the end of the source.
fn take(s: &[u8], from: usize, n: usize) -> Vec<u8> {
    if from >= s.len() {
        return Vec::new();
    }
    s[from..(from + n).min(s.len())].to_vec()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(rows: &[&str]) -> Vec<Vec<u8>> {
        rows.iter().map(|r| r.as_bytes().to_vec()).collect()
    }
    fn s(rows: &[Vec<u8>]) -> Vec<String> {
        rows.iter()
            .map(|r| String::from_utf8(r.clone()).unwrap())
            .collect()
    }

    #[test]
    fn commongappick_drops_only_all_gap_columns() {
        let mut a = v(&["A--C-D", "A-GC-D"]);
        commongappick(&mut a);
        // Columns 1 and 4 are gaps in both rows; column 2 is not.
        assert_eq!(s(&a), vec!["A-CD".to_string(), "AGCD".to_string()]);
    }

    #[test]
    fn commongappick_is_a_no_op_without_common_gaps() {
        let mut a = v(&["ACD", "A-D"]);
        commongappick(&mut a);
        assert_eq!(s(&a), vec!["ACD".to_string(), "A-D".to_string()]);
    }

    #[test]
    fn commongappick_collapses_an_all_gap_alignment_to_nothing() {
        let mut a = v(&["---", "---"]);
        commongappick(&mut a);
        assert_eq!(s(&a), vec![String::new(), String::new()]);
    }

    /// `findcommongaps` describes exactly what `commongappick` removes: the
    /// counts sum to the number of dropped columns, and the map has one entry
    /// per surviving column plus one for the trailing run.
    #[test]
    fn findcommongaps_and_commongappick_agree() {
        for rows in [
            vec!["A--C-D", "A-GC-D"],
            vec!["--AC", "--GC"],
            vec!["AC--", "GC--"],
            vec!["ACD", "A-D"],
            vec!["---", "---"],
        ] {
            let orig = v(&rows);
            let map = findcommongaps(&orig);
            let mut picked = orig.clone();
            commongappick(&mut picked);
            let dropped = orig[0].len() - picked[0].len();
            assert_eq!(map.iter().sum::<i32>(), dropped as i32, "{rows:?}");
            assert_eq!(map.len(), orig[0].len() + 1, "{rows:?}");
            // Everything past the surviving columns must be zero.
            for (i, &g) in map.iter().enumerate() {
                if i > picked[0].len() {
                    assert_eq!(g, 0, "{rows:?} at {i}");
                }
            }
        }
    }

    /// `restorecommongaps` is the inverse of `commongappick` when the DP has
    /// not widened anything in between.
    #[test]
    fn restorecommongaps_inverts_commongappick() {
        for rows in [
            vec!["A--C-D", "A-GC-D"],
            vec!["--AC", "--GC"],
            vec!["AC--", "GC--"],
            vec!["A-C--D-", "AG-C-D-"],
            vec!["ACD", "A-D"],
        ] {
            let orig = v(&rows);
            let mut map = findcommongaps(&orig);
            let mut work = orig.clone();
            commongappick(&mut work);
            let mem: Vec<i32> = (0..work.len() as i32).collect();
            restorecommongaps(&mut work, &mem, &mut map, b'-');
            assert_eq!(s(&work), s(&orig), "{rows:?}");
        }
    }

    #[test]
    fn findnewgaps_counts_marker_columns_before_each_real_one() {
        let a = v(&["=A==CD", "=A==CD"]);
        let g = findnewgaps(&a, 0);
        // One marker before column A, two more before C, none after.
        // The tail past the last real column stays zero; the map is sized
        // over the *widened* sequence, as in the C.
        assert_eq!(g, vec![1, 2, 0, 0, 0, 0, 0]);
        assert_eq!(g.len(), a[0].len() + 1);
    }

    #[test]
    fn findnewgaps_is_all_zero_without_markers() {
        let a = v(&["ACD", "A-D"]);
        assert_eq!(findnewgaps(&a, 0), vec![0, 0, 0, 0]);
    }

    #[test]
    fn adjustgapmap_gives_inserted_columns_a_zero_count() {
        // Before: two real columns with 3 then 1 hidden gaps, 2 trailing.
        let mut map = vec![3, 1, 2];
        // The DP inserted a marker column between them.
        let widened = b"A=C";
        adjustgapmap(widened.len(), &mut map, widened);
        assert_eq!(&map[..4], &[3, 0, 1, 2]);
    }

    #[test]
    fn isallgap_only_for_gaps() {
        assert!(isallgap(b"---"));
        assert!(isallgap(b""));
        assert!(!isallgap(b"--A"));
        assert!(!isallgap(b"-=-"));
    }

    #[test]
    fn take_stops_at_the_end_of_the_source() {
        assert_eq!(take(b"ABCDE", 1, 3), b"BCD".to_vec());
        assert_eq!(take(b"ABCDE", 3, 9), b"DE".to_vec());
        assert_eq!(take(b"ABCDE", 5, 2), Vec::<u8>::new());
        assert_eq!(take(b"ABCDE", 9, 2), Vec::<u8>::new());
    }

    #[test]
    fn count_leading_zeros_measures_the_run() {
        assert_eq!(count_leading_zeros(&[0, 0, 1, 0]), 2);
        assert_eq!(count_leading_zeros(&[1, 0]), 0);
        assert_eq!(count_leading_zeros(&[]), 0);
    }
}
