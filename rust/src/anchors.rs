//! `--anchors` (`-l <file>`): user-supplied anchors that pin regions of two
//! sequences to each other, read from `_externalanchors`.
//!
//! An anchor is a pair of intervals, one in each of two named input sequences,
//! plus a score. During the progressive alignment each node's two groups are
//! asked which anchors apply to any of their member pairs; those become fixed
//! cut points, and the DP only runs between them.
//!
//! The C sorts anchors three times with `qsort` and comparators that break
//! ties on the *addresses* of the elements being compared. On glibc that
//! reduces to stability -- `anchcomp` keeps the earlier element first,
//! `anchscorecomp` reverses its operands and so keeps the later one -- which
//! is what the stable sorts below do, with the position as the secondary key.

/// The C's `ExtAnch`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ExtAnch {
    pub i: i32,
    pub j: i32,
    pub starti: i32,
    pub endi: i32,
    pub startj: i32,
    pub endj: i32,
    pub score: i32,
}

impl ExtAnch {
    fn blank() -> ExtAnch {
        ExtAnch {
            i: -1,
            j: -1,
            starti: -1,
            endi: -1,
            startj: -1,
            endj: -1,
            score: 0,
        }
    }
}

/// `readexternalanchors()`: parse `_externalanchors`.
///
/// Seven whitespace-separated fields per line, all 1-origin and converted to
/// 0-origin here: the two sequence ids, then start/end in each, then a score.
/// `#` starts a comment line. The first id must be less than the second, and
/// both must name an existing sequence. The two *start* positions are checked
/// against the ungapped lengths -- with `>`, so one past the last residue
/// passes -- and the two end positions are not checked at all. That gap is
/// load-bearing: see `map_lookup` below and the cut-list check in
/// `falign_givenanchors`.
pub fn read_external_anchors(nseq: usize, nogaplen: &[usize]) -> Vec<ExtAnch> {
    let text = match std::fs::read("_externalanchors") {
        Ok(t) => t,
        Err(_) => {
            eprintln!("Cannot open _externalanchors");
            std::process::exit(1);
        }
    };
    let mut out: Vec<ExtAnch> = Vec::new();
    let mut lineno = 0usize;
    for raw in text.split_inclusive(|&c| c == b'\n') {
        lineno += 1;
        // `fgets` then `feof`: an unterminated last line is not read.
        if raw.last() != Some(&b'\n') {
            break;
        }
        if raw.first() == Some(&b'#') {
            continue;
        }
        let mut f = raw
            .split(|c: &u8| c.is_ascii_whitespace())
            .filter(|t| !t.is_empty())
            .map(|t| {
                std::str::from_utf8(t)
                    .ok()
                    .and_then(|s| s.parse::<i32>().ok())
            });
        // `sscanf` leaves unmatched fields untouched; the C zeroes i and j
        // first and leaves the rest at whatever the previous line put there.
        // Every real file gives all seven, so take them straight.
        let mut a = ExtAnch::blank();
        a.i = f.next().flatten().unwrap_or(0);
        a.j = f.next().flatten().unwrap_or(0);
        a.starti = f.next().flatten().unwrap_or(0);
        a.endi = f.next().flatten().unwrap_or(0);
        a.startj = f.next().flatten().unwrap_or(0);
        a.endj = f.next().flatten().unwrap_or(0);
        a.score = f.next().flatten().unwrap_or(0);
        a.i -= 1;
        a.j -= 1;
        a.starti -= 1;
        a.startj -= 1;
        a.endi -= 1;
        a.endj -= 1;

        if a.i as usize >= nseq || a.j as usize >= nseq {
            eprintln!("\nOut of range?  The input file has {} sequences but pair {}-{} was specified in line {}.\nNote that sequence IDs are counted from 1.", nseq, a.i + 1, a.j + 1, lineno);
            std::process::exit(1);
        }
        if a.i >= a.j {
            eprintln!("\nFormat problem?  \"{} {}\" in line {}.\nThe sequence id of the first column must be less than the second.", a.i + 1, a.j + 1, lineno);
            std::process::exit(1);
        }
        if a.starti > nogaplen[a.i as usize] as i32 {
            eprintln!("\nOut of range?  len(seq{})={}, but anchor={} in line {}.\nNote that position is counted from 1.", a.i + 1, nogaplen[a.i as usize], a.starti + 1, lineno);
            std::process::exit(1);
        }
        if a.startj > nogaplen[a.j as usize] as i32 {
            eprintln!("\nOut of range?  len(seq{})={}, but anchor={} in line {}.\nNote that position is counted from 1.", a.j, nogaplen[a.j as usize] + 1, a.startj + 1, lineno);
            std::process::exit(1);
        }
        out.push(a);
    }
    out
}

/// `indexanchors()`: sort by (i, j) and record where each pair's run starts.
///
/// `anchidcomp` has no tie-break of its own, so equal (i, j) keep their file
/// order -- glibc's mergesort is stable.
pub fn index_anchors(a: &mut [ExtAnch], nseq: usize) -> Vec<Vec<i32>> {
    a.sort_by(|p, q| (p.i, p.j).cmp(&(q.i, q.j)));
    let mut idx = vec![vec![-1i32; nseq]; nseq];
    for (n, e) in a.iter().enumerate() {
        if idx[e.i as usize][e.j as usize] == -1 {
            idx[e.i as usize][e.j as usize] = n as i32;
        }
    }
    idx
}

/// `checkanchors_strongestfirst()`: keep the highest-scoring anchors and drop
/// any that cannot coexist with one already kept.
///
/// Two anchors are consistent when one lies entirely before the other in both
/// sequences. Sharing exactly one endpoint is repaired by moving that endpoint
/// one residue; anything else disqualifies the weaker anchor, which is marked
/// by setting its coordinates to -1 rather than removed.
pub fn check_anchors_strongest_first(a: &mut Vec<ExtAnch>) {
    // `anchscorecomp`: descending by score, and its tie-break compares
    // `q - p`, so among equal scores the *later* element comes first.
    let n = a.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&x, &y| a[y].score.cmp(&a[x].score).then(y.cmp(&x)));
    let sorted: Vec<ExtAnch> = order.iter().map(|&k| a[k]).collect();
    *a = sorted;

    for p in 0..n {
        if a[p].starti == -1 {
            continue;
        }
        for q in p + 1..n {
            if a[q].starti == -1 {
                continue;
            }
            if a[p].endi < a[q].starti && a[p].endj < a[q].startj {
                // consistent
            } else if a[p].endi == a[q].starti && a[p].endj < a[q].startj && a[q].starti < a[q].endi
            {
                a[q].starti += 1;
            } else if a[p].endi < a[q].starti && a[p].endj == a[q].startj && a[q].startj < a[q].endj
            {
                a[q].startj += 1;
            } else if a[q].endi < a[p].starti && a[q].endj < a[p].startj {
                // consistent
            } else if a[q].endi == a[p].starti && a[q].endj < a[p].startj && a[q].starti < a[q].endi
            {
                a[q].endi -= 1;
            } else if a[q].endi < a[p].starti && a[q].endj == a[p].startj && a[q].startj < a[q].endj
            {
                a[q].endj -= 1;
            } else {
                // The C writes `a[q].startj` twice and never clears `starti`'s
                // partner `endi`; kept exactly, because which fields end up at
                // -1 decides what the caller skips.
                a[q].starti = -1;
                a[q].startj = -1;
                a[q].endj = -1;
            }
        }
    }

    // `anchcomp`: ascending by starti, ties keeping the earlier element.
    let n2 = a.len();
    let mut order: Vec<usize> = (0..n2).collect();
    order.sort_by(|&x, &y| a[x].starti.cmp(&a[y].starti).then(x.cmp(&y)));
    let sorted: Vec<ExtAnch> = order.iter().map(|&k| a[k]).collect();
    *a = sorted;
}

/// `recountpositions()`: anchor coordinates arrive as residue numbers and are
/// rewritten as column numbers in the group they belong to.
///
/// The map is refilled for every member `k` and applied to the anchors whose
/// `i` (or `j`) is that member -- so each anchor is translated through its own
/// sequence's gaps, not the group's first row.
///
/// The C allocates `map` once, before the `k` loop, and never clears it
/// between members. Entries above a member's residue count therefore keep
/// whatever an earlier, longer member wrote there, and this reuse is visible
/// whenever an anchor position is out of range for its own row -- which
/// `readexternalanchors` permits, see `map_lookup` below. So the buffer is
/// hoisted here too rather than being allocated per member.
fn recount_positions(pairanch: &mut [ExtAnch], seq1: &[Vec<u8>], seq2: &[Vec<u8>]) {
    // `len = strlen(seq1[0]) + 1`: the map has room for a position one past
    // the last column, which an anchor ending at the sequence end needs.
    let len = seq1[0].len() + 1;
    let mut map = vec![0i32; len];
    for (k, s) in seq1.iter().enumerate() {
        let mut pos = 0usize;
        for i in 0..len {
            if s.get(i).copied().unwrap_or(0) != b'-' {
                map[pos] = i as i32;
                pos += 1;
            }
        }
        for a in pairanch.iter_mut() {
            if a.i == k as i32 {
                a.starti = map_lookup(&map, a.starti, k, 1);
                a.endi = map_lookup(&map, a.endi, k, 1);
            }
        }
    }
    let len = seq2[0].len() + 1;
    let mut map = vec![0i32; len];
    for (k, s) in seq2.iter().enumerate() {
        let mut pos = 0usize;
        for i in 0..len {
            if s.get(i).copied().unwrap_or(0) != b'-' {
                map[pos] = i as i32;
                pos += 1;
            }
        }
        for a in pairanch.iter_mut() {
            if a.j == k as i32 {
                a.startj = map_lookup(&map, a.startj, k, 2);
                a.endj = map_lookup(&map, a.endj, k, 2);
            }
        }
    }
}

/// `map[p]` as the C indexes it, with the one case it cannot survive refused.
///
/// `readexternalanchors` (io.c:6378) checks only `starti` and `startj` against
/// `nogaplen`, and even those with `>` rather than `>=`; `endi` and `endj` are
/// never checked at all. A position past the end of the group's column count
/// is therefore an out-of-bounds read of `map`, which the C `calloc`s over the
/// *column* count, not the residue count. Reproducing whatever the heap held
/// there is not possible, so refuse.
///
/// A position that is merely past this member's *residue* count but still
/// inside `map` is not refused: the C reads the zero or the stale entry left
/// by an earlier member, and the reuse above reproduces exactly that. It is
/// still a defect -- it is what produces a descending cut list in
/// `falign_givenanchors`, see the check there -- but it is a defined read.
fn map_lookup(map: &[i32], p: i32, k: usize, side: u32) -> i32 {
    if p < 0 || p as usize >= map.len() {
        eprintln!(
            "disttbfast-rs: an anchor puts position {} on member {} of group {}, whose\n  \
            alignment is {} columns wide. readexternalanchors never range-checks the\n  \
            end positions (io.c:6378), so the C reads recountpositions' map out of\n  \
            bounds here (disttbfast.c:1373). Refusing rather than guessing what that\n  \
            read returned. Check the end positions in _externalanchors: they count\n  \
            residues that MAFFT recognises, which is fewer than the input has if the\n  \
            sequence contains characters MAFFT drops.",
            p + 1,
            k + 1,
            side,
            map.len() - 1
        );
        std::process::exit(1);
    }
    map[p as usize]
}

/// `pickpairanch()`: the anchors that apply to any member pair of the two
/// groups, renumbered onto the groups' own coordinates and made consistent.
pub fn pick_pair_anch(
    extanch: &[ExtAnch],
    anchindex: &[Vec<i32>],
    m1: &[i32],
    m2: &[i32],
    seq1: &[Vec<u8>],
    seq2: &[Vec<u8>],
) -> Vec<ExtAnch> {
    let mut out: Vec<ExtAnch> = Vec::new();
    for (i, &a) in m1.iter().enumerate() {
        for (j, &b) in m2.iter().enumerate() {
            let (lo, hi, swap) = if a < b { (a, b, false) } else { (b, a, true) };
            let mut k = anchindex[lo as usize][hi as usize];
            while k != -1 {
                let ku = k as usize;
                if ku >= extanch.len() || extanch[ku].i != lo || extanch[ku].j != hi {
                    break;
                }
                let e = extanch[ku];
                out.push(ExtAnch {
                    i: i as i32,
                    j: j as i32,
                    starti: if swap { e.startj } else { e.starti },
                    endi: if swap { e.endj } else { e.endi },
                    startj: if swap { e.starti } else { e.startj },
                    endj: if swap { e.endi } else { e.endj },
                    score: e.score,
                });
                k += 1;
            }
        }
    }
    recount_positions(&mut out, seq1, seq2);
    check_anchors_strongest_first(&mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn anch(i: i32, j: i32, si: i32, ei: i32, sj: i32, ej: i32, score: i32) -> ExtAnch {
        ExtAnch {
            i,
            j,
            starti: si,
            endi: ei,
            startj: sj,
            endj: ej,
            score,
        }
    }

    #[test]
    fn index_anchors_sorts_by_pair_and_indexes_the_first_of_each_run() {
        let mut a = vec![
            anch(1, 2, 0, 4, 0, 4, 10),
            anch(0, 1, 0, 4, 0, 4, 20),
            anch(1, 2, 8, 9, 8, 9, 30),
            anch(0, 2, 0, 4, 0, 4, 40),
        ];
        let idx = index_anchors(&mut a, 3);
        assert_eq!(
            a.iter().map(|e| (e.i, e.j)).collect::<Vec<_>>(),
            vec![(0, 1), (0, 2), (1, 2), (1, 2)]
        );
        // The two 1x2 anchors keep their file order: glibc's qsort is a
        // mergesort, and `anchidcomp` has no tie-break of its own.
        assert_eq!(a[2].score, 10);
        assert_eq!(a[3].score, 30);
        assert_eq!(idx[0][1], 0);
        assert_eq!(idx[0][2], 1);
        assert_eq!(idx[1][2], 2);
        assert_eq!(idx[2][0], -1);
    }

    #[test]
    fn strongest_first_keeps_consistent_anchors_and_reorders_by_start() {
        // Disjoint and in the same order in both sequences, so both survive.
        let mut a = vec![anch(0, 1, 40, 49, 40, 49, 5), anch(0, 1, 0, 9, 0, 9, 90)];
        check_anchors_strongest_first(&mut a);
        assert_eq!(a[0].starti, 0);
        assert_eq!(a[1].starti, 40);
        assert_eq!(a[0].score, 90);
    }

    #[test]
    fn strongest_first_drops_the_weaker_of_a_crossing_pair() {
        // The weak anchor runs backwards relative to the strong one: ahead in
        // sequence i, behind in sequence j. Nothing can repair that.
        let mut a = vec![anch(0, 1, 40, 49, 0, 9, 5), anch(0, 1, 0, 9, 40, 49, 90)];
        check_anchors_strongest_first(&mut a);
        let dead = a.iter().find(|e| e.starti == -1).expect("one is dropped");
        assert_eq!(dead.score, 5);
        // The C's invalidation assigns `startj` twice and never touches
        // `endi`, so the dead anchor keeps its end position in sequence i.
        assert_eq!(dead.endi, 49);
        assert_eq!(dead.startj, -1);
        assert_eq!(dead.endj, -1);
        assert_eq!(a.iter().filter(|e| e.starti != -1).count(), 1);
    }

    #[test]
    fn strongest_first_repairs_a_single_shared_endpoint() {
        // The weak anchor starts exactly where the strong one ends in i, and
        // strictly after it in j: one residue of overlap, so its start moves.
        let mut a = vec![anch(0, 1, 9, 19, 30, 39, 5), anch(0, 1, 0, 9, 0, 9, 90)];
        check_anchors_strongest_first(&mut a);
        let weak = a.iter().find(|e| e.score == 5).unwrap();
        assert_eq!(weak.starti, 10);
        assert_eq!(weak.endi, 19);
    }

    #[test]
    fn strongest_first_breaks_score_ties_towards_the_later_anchor() {
        // `anchscorecomp` compares `q - p` on the addresses, so among equal
        // scores glibc's stable mergesort leaves the later element first.
        // Both anchors here are mutually consistent, so ordering is all that
        // is observable -- and the final sort is by starti, so give them the
        // same start and read the scores back through `endi`.
        let mut a = vec![anch(0, 1, 0, 100, 0, 9, 7), anch(0, 1, 0, 200, 0, 9, 7)];
        // Only the score sort is under test; make them conflict so the second
        // in score order is the one dropped.
        check_anchors_strongest_first(&mut a);
        let dead = a.iter().find(|e| e.starti == -1).expect("one is dropped");
        assert_eq!(dead.endi, 100);
    }

    #[test]
    fn recount_positions_translates_residue_numbers_into_columns() {
        let seq1 = vec![b"AB--CDE".to_vec()];
        let seq2 = vec![b"--ABCDE".to_vec()];
        let mut a = vec![anch(0, 0, 1, 3, 0, 2, 1)];
        recount_positions(&mut a, &seq1, &seq2);
        // residues of seq1 sit at columns 0,1,4,5,6
        assert_eq!(a[0].starti, 1);
        assert_eq!(a[0].endi, 5);
        // residues of seq2 sit at columns 2,3,4,5,6
        assert_eq!(a[0].startj, 2);
        assert_eq!(a[0].endj, 4);
    }

    #[test]
    fn recount_positions_reuses_the_map_across_members() {
        // The C allocates `map` once for the whole group. Member 1 has only
        // two residues, so it never writes map[3] or map[4]; an anchor
        // position that high reads what member 0 left there. Allocating the
        // map per member would give 0 instead.
        let seq1 = vec![b"ABCDE".to_vec(), b"AB---".to_vec()];
        let seq2 = vec![b"ABCDE".to_vec()];
        let mut a = vec![anch(1, 9, 1, 4, 0, 1, 1)];
        recount_positions(&mut a, &seq1, &seq2);
        assert_eq!(a[0].starti, 1);
        assert_eq!(a[0].endi, 4, "stale entry from the longer member 0");
    }

    #[test]
    fn recount_positions_maps_one_past_the_last_residue_to_the_terminator() {
        // `len` is `strlen(seq[0]) + 1`, and the fill loop runs to `len`, so
        // the NUL at the end counts as a residue and gets a map entry. An
        // anchor ending one past the sequence lands on that column.
        let seq1 = vec![b"ABCD".to_vec()];
        let seq2 = vec![b"ABCD".to_vec()];
        let mut a = vec![anch(0, 0, 0, 4, 0, 4, 1)];
        recount_positions(&mut a, &seq1, &seq2);
        assert_eq!(a[0].endi, 4);
        assert_eq!(a[0].endj, 4);
    }
}
