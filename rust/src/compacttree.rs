//! `--treein` with a `compact`, `very compact`, `stepadd` or `youngestlinkage`
//! keyword: guide trees built without ever holding a full distance matrix.
//!
//! The ordinary builder needs an `njob x njob` half matrix. These build the
//! same shape of tree from, per sequence, only its nearest neighbour and the
//! distance to it -- `mindistfrom[i]` and `mindist[i]` -- optionally backed by
//! as many full matrix rows as a stated RAM budget allows (`partmtx`).
//!
//! The C reaches them through `compacttree`, set by `check_guidetreefile`:
//!
//! | keyword           | `compacttree` | rows in `partmtx`  | builder                        |
//! |-------------------|---------------|--------------------|--------------------------------|
//! | `compact <size>`  | 1             | up to the budget   | `compacttreegivendist`         |
//! | `very compact`    | 2             | none               | `compacttreegivendist`         |
//! | `stepadd`         | 2             | none               | `compacttreegivendist`         |
//! | `youngestlinkage` | 4             | none               | `compacttree_memsaveselectable` |
//!
//! `compact` alone falls back to the ordinary builder when the budget is large
//! enough to hold the whole matrix, which is the common case and is already
//! covered elsewhere.

use crate::tree::{Tree, Treedep};

/// `preferenceval()`: a tie-break, tiny enough not to reorder real distances.
///
/// It biases the nearest-neighbour search towards the sequence that follows
/// `ori` in input order, wrapping at `max`. The C calls it "for debug"; it is
/// not, it is what makes the search deterministic.
pub fn preferenceval(ori: i32, pos: i32, max: i32) -> f64 {
    let mut pos = pos - ori;
    if pos < 0 {
        pos += max;
    }
    0.00000000000001 * pos as f64
}

/// `preparepartmtx()`: the rows of the full distance matrix that fit the budget.
///
/// Only `compacttree == 1` caches anything. The C fills rows in order until
/// the running total passes `maxdistmtxsize` and leaves the rest NULL, so the
/// cache covers a prefix of the sequences.
pub fn prepare_partmtx(nseq: usize, compacttree: i32, maxdistmtxsize: f64) -> Vec<Option<Vec<f64>>> {
    let mut val: Vec<Option<Vec<f64>>> = (0..nseq).map(|_| None).collect();
    if compacttree != 1 {
        return val;
    }
    let mut size = 0.0f64;
    let mut i = 0usize;
    while i < nseq {
        size += 8.0 * nseq as f64; // sizeof( double ) * nseq
        if size > maxdistmtxsize {
            eprintln!(
                "\n\nThe size of full distance matrix is estimated to exceed {:.2}GB.",
                maxdistmtxsize / 1000.0 / 1000.0 / 1000.0
            );
            eprintln!("Will try the calculation using a {} x {} matrix.", nseq, i);
            break;
        }
        val[i] = Some(vec![0.0; nseq]);
        i += 1;
    }
    if i == nseq {
        eprintln!("The full matrix will be used.");
    }
    val
}

/// `distcompact()`: one pairwise distance from a 6-mer intersection count.
pub fn distcompact(
    len1: usize,
    len2: usize,
    common: i32,
    ss1: i32,
    ss2: i32,
    lenfaca: f64,
    lenfacb: f64,
    lenfacc: f64,
    lenfacd: f64,
) -> f64 {
    let (longer, shorter) = if len1 > len2 {
        (len1 as f64, len2 as f64)
    } else {
        (len2 as f64, len1 as f64)
    };
    let lenfac = 1.0 / (shorter / longer * lenfacd + lenfacb / (longer + lenfacc) + lenfaca);
    if ss1 == 0 || ss2 == 0 {
        return 2.0;
    }
    (1.0 - common as f64 / ss1.min(ss2) as f64) * lenfac * 2.0
}

/// `distcompact_msa()`: the same, once the sequences are aligned.
pub fn distcompact_msa(naivescore: f64, ss1: i32, ss2: i32) -> f64 {
    let bunbo = ss1.min(ss2);
    if bunbo == 0 {
        return 2.0;
    }
    let value = (1.0 - naivescore / bunbo as f64) * 2.0;
    if value > 10.0 {
        10.0
    } else {
        value
    }
}

/// One node of the C's `Treept`, with array indices where the C has pointers.
///
/// `-1` stands for NULL. Nodes `0..njob` are the leaves, in input order; the
/// internal nodes are appended from `njob` upwards in the order they are made.
#[derive(Clone, Copy)]
struct Treept {
    child0: i32,
    child1: i32,
    parent: i32,
    rep0: i32,
    rep1: i32,
    len0: f64,
    len1: f64,
    height: f64,
}

impl Treept {
    fn leaf(i: i32) -> Treept {
        Treept {
            child0: -1,
            child1: -1,
            parent: -1,
            rep0: i,
            rep1: -1,
            len0: 0.0,
            len1: 0.0,
            height: 0.0,
        }
    }
}

/// `compacttreegivendist()`: build the tree from nearest neighbours alone.
///
/// Sequence `i` is attached to the branch above its recorded neighbour, at
/// half its recorded distance -- walking up from that neighbour until a node
/// is found whose height already exceeds it. Nothing is recomputed, so the
/// whole tree costs one pass.
///
/// `TREE7325` is 0, so the seed is the *first* pair and the loop runs forwards.
pub fn compacttreegivendist(
    njob: usize,
    mindists: &[f64],
    neighbors: &[i32],
    names: Option<&[Vec<u8>]>,
    outnumber: bool,
    treeout: bool,
) -> Tree {
    let mut treept: Vec<Treept> = (0..njob * 2).map(|i| Treept::leaf(i as i32)).collect();

    let mut n = njob;
    let tmpdist = mindists[1];
    treept[0].parent = n as i32;
    treept[1].parent = n as i32;
    treept[n].child0 = 0;
    treept[n].child1 = 1;
    treept[n].height = tmpdist * 0.5;
    treept[n].len0 = tmpdist * 0.5;
    treept[n].len1 = tmpdist * 0.5;
    treept[n].parent = -1;
    treept[n].rep0 = 0;
    treept[n].rep1 = 1;
    let mut root = n;

    for i in 2..njob {
        n += 1;
        let neighbor = neighbors[i];
        let mindist = mindists[i];

        // Climb from the neighbour until the branch is old enough to carry
        // this join, i.e. until a node higher than `mindist * 0.5`.
        let mut b = neighbor;
        let mut p = treept[b as usize].parent;
        while p != -1 {
            if treept[p as usize].height > mindist * 0.5 {
                break;
            }
            b = p;
            p = treept[p as usize].parent;
        }

        if p == -1 {
            treept[n].parent = -1;
            root = n;
        } else if treept[p as usize].child0 == b {
            treept[p as usize].child0 = n as i32;
            treept[p as usize].len0 = treept[p as usize].height - mindist * 0.5;
            treept[n].parent = p;
        } else if treept[p as usize].child1 == b {
            treept[p as usize].child1 = n as i32;
            treept[p as usize].len1 = treept[p as usize].height - mindist * 0.5;
            treept[n].parent = p;
        } else {
            eprintln!("okashii");
            std::process::exit(1);
        }

        treept[i].parent = n as i32;
        treept[b as usize].parent = n as i32;

        let brep0 = treept[b as usize].rep0;
        let bheight = treept[b as usize].height;
        treept[n].child0 = b;
        treept[n].child1 = i as i32;
        treept[n].height = mindist * 0.5;
        treept[n].rep0 = brep0;
        treept[n].rep1 = treept[i].rep0;
        treept[n].len0 = mindist * 0.5 - bheight;
        treept[n].len1 = mindist * 0.5;
    }

    reformattree(root, &treept, njob, names, outnumber, treeout)
}

/// `reformattree()`: turn the linked `Treept` tree into `topol`/`len`/`dep`.
///
/// `YOUNGER0TREE` is 1, so the live traversal is `reformat_younger0_rec`,
/// which orders each node's two sides by representative rather than keeping
/// the tree's own child order. The traversal is post-order, and `lastappear`
/// carries which step each representative was last seen at, which is what
/// `dep[].child0`/`child1` mean.
///
/// The C recurses; this walks with an explicit stack so a comb-shaped tree
/// over a large input cannot overflow. The visit order is unchanged, which is
/// all that reaches the output.
fn reformattree(
    root: usize,
    treept: &[Treept],
    n: usize,
    names: Option<&[Vec<u8>]>,
    outnumber: bool,
    treeout: bool,
) -> Tree {
    if treeout {
        if let Some(names) = names {
            let leaves: Vec<String> = (0..n)
                .map(|i| format!("\n{}_{}\n", i + 1, crate::tree::newick_name(&names[i], outnumber)))
                .collect();
            let newick = reformat_rec_newick(&leaves, treept, root as i32);
            std::fs::write("infile.tree", format!("{};\n", newick))
                .expect("Cannot open infile.tree");
        }
    }

    // `n` rows, not `n-1`: the C allocates one spare, zeroed, and
    // `counteff_simple` walks over it. Only the first `n-1` are ever filled.
    let mut topol: Vec<[i32; 2]> = vec![[-1, -1]; n];
    let mut len: Vec<[f64; 2]> = vec![[0.0, 0.0]; n];
    let mut dep: Vec<Treedep> = vec![Treedep::default(); n];
    let mut lastappear = vec![-1i32; n];
    let mut pos = 0usize;

    // Post-order over the internal nodes; a leaf has `rep1 == -1` and the C
    // returns from it before touching anything.
    let mut stack: Vec<(i32, bool)> = vec![(root as i32, false)];
    while let Some((node, expanded)) = stack.pop() {
        if node == -1 || treept[node as usize].rep1 == -1 {
            continue;
        }
        if !expanded {
            stack.push((node, true));
            // Pushed in reverse so child0 is visited first.
            stack.push((treept[node as usize].child1, false));
            stack.push((treept[node as usize].child0, false));
            continue;
        }
        let pt = &treept[node as usize];
        if pt.rep0 < pt.rep1 {
            topol[pos] = [pt.rep0, pt.rep1];
            len[pos] = [pt.len0, pt.len1];
            dep[pos].child0 = lastappear[pt.rep0 as usize];
            dep[pos].child1 = lastappear[pt.rep1 as usize];
        } else {
            topol[pos] = [pt.rep1, pt.rep0];
            len[pos] = [pt.len1, pt.len0];
            dep[pos].child0 = lastappear[pt.rep1 as usize];
            dep[pos].child1 = lastappear[pt.rep0 as usize];
        }
        lastappear[pt.rep0 as usize] = pos as i32;
        lastappear[pt.rep1 as usize] = pos as i32;
        dep[pos].distfromtip = pt.height;
        pos += 1;
    }

    Tree { topol, len, dep }
}

/// `reformat_rec_newick()`: the Newick string for `infile.tree`.
///
/// Each node prints its lower-numbered representative's side first, which is
/// the same reordering `reformat_younger0_rec` applies to `topol`.
fn reformat_rec_newick(leaves: &[String], treept: &[Treept], node: i32) -> String {
    let pt = &treept[node as usize];
    if pt.rep1 == -1 {
        return leaves[pt.rep0 as usize].clone();
    }
    let n0 = reformat_rec_newick(leaves, treept, pt.child0);
    let n1 = reformat_rec_newick(leaves, treept, pt.child1);
    if pt.rep0 < pt.rep1 {
        format!("({}:{:7.5},{}:{:7.5})", n0, pt.len0, n1, pt.len1)
    } else {
        format!("({}:{:7.5},{}:{:7.5})", n1, pt.len1, n0, pt.len0)
    }
}

/// The `if( compacttree )` arm of "Making a distance matrix .." in
/// `disttbfast()` (disttbfast.c:3638).
///
/// Instead of an `njob x njob` half matrix this keeps, per sequence, only its
/// nearest neighbour and the distance to it -- plus whatever full rows
/// `partmtx` was allowed to hold.
///
/// The traversal direction is what separates `youngestlinkage` from the rest.
/// `compacttree == 4` walks `i` upwards and, for each pair, offers the
/// distance to *both* endpoints; the others walk `i` downwards from `njob-1`
/// and offer it only to `i`. So under `youngestlinkage` a sequence can be
/// claimed by a later one, and otherwise it cannot.
///
/// `mindist` carries `preferenceval`'s tie-break while the search runs and has
/// it subtracted again at the end, exactly as the C does -- so the bias
/// decides ties without reaching the distances the tree is built from.
pub fn compact_distances(
    cn: &crate::constants::Constants,
    seq: &[Vec<u8>],
    njob: usize,
    tuplesize: usize,
    compacttree: i32,
    maxdistmtxsize: f64,
) -> (Vec<f64>, Vec<i32>, Vec<usize>, Vec<Option<Vec<f64>>>, Vec<i32>) {
    let t = crate::dist::build_tuples(cn, seq, njob, tuplesize);

    eprintln!("Compact tree, step 1");
    let mut mindistfrom = vec![-1i32; njob];
    let mut mindist = vec![999.9f64; njob];
    let mut selfscore = vec![0i32; njob];
    let mut partmtx = prepare_partmtx(njob, compacttree, maxdistmtxsize);

    let mut sextet = crate::dist::Sextet::new(t.tsize, t.maxl);
    let mut table1 = vec![0i32; t.tsize];
    for i in 0..njob {
        table1.iter_mut().for_each(|v| *v = 0);
        crate::dist::makecompositiontable_p(&mut table1, &t.pointt[i]);
        selfscore[i] = sextet.common(&table1, &t.pointt[i]);
    }

    let dist = |sx: &mut crate::dist::Sextet, tab: &[i32], i: usize, j: usize| -> f64 {
        distcompact(
            t.nogaplen[i],
            t.nogaplen[j],
            sx.common(tab, &t.pointt[j]),
            selfscore[i],
            selfscore[j],
            t.lenfaca,
            t.lenfacb,
            t.lenfacc,
            t.lenfacd,
        )
    };

    if compacttree == 4 {
        // ylcompactdisthalfmtxthread: i ascending, j > i, both sides offered.
        for i in 0..njob - 1 {
            table1.iter_mut().for_each(|v| *v = 0);
            crate::dist::makecompositiontable_p(&mut table1, &t.pointt[i]);
            for j in i + 1..njob {
                let tmpdist = dist(&mut sextet, &table1, i, j);
                let tmpdistx = tmpdist + preferenceval(i as i32, j as i32, njob as i32);
                if tmpdistx < mindist[i] {
                    mindist[i] = tmpdistx;
                    mindistfrom[i] = j as i32;
                }
                let tmpdisty = tmpdist + preferenceval(j as i32, i as i32, njob as i32);
                if tmpdisty < mindist[j] {
                    mindist[j] = tmpdisty;
                    mindistfrom[j] = i as i32;
                }
                if let Some(row) = partmtx[i].as_mut() {
                    row[j] = tmpdist;
                }
                if let Some(row) = partmtx[j].as_mut() {
                    row[i] = tmpdist;
                }
            }
        }
    } else {
        // compactdisthalfmtxthread: i descending from njob-1, j < i, only i
        // is offered the distance. i == 0 does no work at all, so sequence 0
        // keeps mindist 999.9 and no neighbour.
        for i in (0..njob).rev() {
            table1.iter_mut().for_each(|v| *v = 0);
            crate::dist::makecompositiontable_p(&mut table1, &t.pointt[i]);
            for j in (0..i).rev() {
                let tmpdist = dist(&mut sextet, &table1, i, j);
                let tmpdistx = tmpdist + preferenceval(i as i32, j as i32, njob as i32);
                if tmpdistx < mindist[i] {
                    mindist[i] = tmpdistx;
                    mindistfrom[i] = j as i32;
                }
                if let Some(row) = partmtx[i].as_mut() {
                    row[j] = tmpdist;
                }
                if let Some(row) = partmtx[j].as_mut() {
                    row[i] = tmpdist;
                }
            }
        }
    }

    // `mindist[i] -= preferenceval( i, mindistfrom[i], njob )`, including for
    // the sequence that never found a neighbour: `preferenceval(i, -1, njob)`
    // is a real number, so its 999.9 is nudged too.
    for i in 0..njob {
        mindist[i] -= preferenceval(i as i32, mindistfrom[i], njob as i32);
    }
    eprintln!("\ndone.");

    (mindist, mindistfrom, t.nogaplen, partmtx, selfscore)
}

/// The `else if( calcpairdists && compacttree )` arm of `disttbfast()`
/// (disttbfast.c:4699): the same nearest-neighbour reduction, but taken from
/// the alignment the previous guide-tree pass produced rather than from k-mers.
///
/// The direction rule is the one from `compact_distances`: `youngestlinkage`
/// ascends and offers both sides, the others descend and offer only `i`.
#[allow(clippy::too_many_arguments)]
pub fn msa_compact_distances(
    cn: &crate::constants::Constants,
    bseq: &[Vec<u8>],
    selfscore: &[i32],
    njob: usize,
    compacttree: i32,
    maxdistmtxsize: f64,
) -> (Vec<f64>, Vec<i32>) {
    eprintln!("Making a compact tree from msa, step 1.. ");
    let skiptable: Vec<Vec<i32>> = bseq.iter().map(|s| crate::score::makeskiptable(s)).collect();
    let mut mindistfrom = vec![-1i32; njob];
    let mut mindist = vec![999.9f64; njob];
    let mut partmtx = prepare_partmtx(njob, compacttree, maxdistmtxsize);

    let pairdist = |i: usize, j: usize| -> f64 {
        distcompact_msa(
            crate::score::naivepairscorefast(
                cn,
                &bseq[i],
                &bseq[j],
                &skiptable[i],
                &skiptable[j],
                cn.penalty_dist,
            ),
            selfscore[i],
            selfscore[j],
        )
    };

    if compacttree == 4 {
        for i in 0..njob - 1 {
            for j in i + 1..njob {
                let tmpdist = pairdist(i, j);
                let tmpdistx = tmpdist + preferenceval(i as i32, j as i32, njob as i32);
                if tmpdistx < mindist[i] {
                    mindist[i] = tmpdistx;
                    mindistfrom[i] = j as i32;
                }
                let tmpdisty = tmpdist + preferenceval(j as i32, i as i32, njob as i32);
                if tmpdisty < mindist[j] {
                    mindist[j] = tmpdisty;
                    mindistfrom[j] = i as i32;
                }
                if let Some(row) = partmtx[i].as_mut() {
                    row[j] = tmpdist;
                }
                if let Some(row) = partmtx[j].as_mut() {
                    row[i] = tmpdist;
                }
            }
        }
    } else {
        for i in (0..njob).rev() {
            for j in (0..i).rev() {
                let tmpdist = pairdist(i, j);
                let tmpdistx = tmpdist + preferenceval(i as i32, j as i32, njob as i32);
                if tmpdistx < mindist[i] {
                    mindist[i] = tmpdistx;
                    mindistfrom[i] = j as i32;
                }
                if let Some(row) = partmtx[i].as_mut() {
                    row[j] = tmpdist;
                }
                if let Some(row) = partmtx[j].as_mut() {
                    row[i] = tmpdist;
                }
            }
        }
    }

    for i in 0..njob {
        mindist[i] -= preferenceval(i as i32, mindistfrom[i], njob as i32);
    }
    eprintln!("\ndone.");
    (mindist, mindistfrom)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preferenceval_wraps_and_stays_below_any_real_distance() {
        // The bias is a fraction of a step forward in input order, wrapping.
        assert_eq!(preferenceval(0, 0, 10), 0.0);
        assert_eq!(preferenceval(0, 3, 10), 3.0e-14);
        // 2 -> 1 wraps the long way round: 9 steps, not -1.
        assert_eq!(preferenceval(2, 1, 10), 9.0e-14);
        // Even the largest bias over a big input stays far below the
        // resolution of a k-mer distance.
        assert!(preferenceval(0, 999_999, 1_000_000) < 1e-7);
    }

    #[test]
    fn prepare_partmtx_caches_a_prefix_up_to_the_budget() {
        // Ten sequences, room for three rows of 10 doubles (240 bytes) but
        // not four: the C stops the moment the running total passes.
        let m = prepare_partmtx(10, 1, 250.0);
        assert_eq!(m.iter().filter(|r| r.is_some()).count(), 3);
        assert!(m[0].is_some() && m[2].is_some() && m[3].is_none());
    }

    #[test]
    fn prepare_partmtx_caches_nothing_unless_compacttree_is_one() {
        for ct in [0, 2, 4] {
            let m = prepare_partmtx(10, ct, 1e12);
            assert!(m.iter().all(|r| r.is_none()), "compacttree={}", ct);
        }
    }

    #[test]
    fn distcompact_is_two_when_either_selfscore_is_zero() {
        assert_eq!(distcompact(10, 10, 5, 0, 7, 1.0, 1.0, 1.0, 1.0), 2.0);
        assert_eq!(distcompact(10, 10, 5, 7, 0, 1.0, 1.0, 1.0, 1.0), 2.0);
    }

    #[test]
    fn distcompact_msa_saturates_at_ten() {
        // A strongly negative pair score would otherwise run away.
        assert_eq!(distcompact_msa(-100.0, 10, 10), 10.0);
        assert_eq!(distcompact_msa(0.0, 0, 10), 2.0);
        assert_eq!(distcompact_msa(10.0, 10, 20), 0.0);
    }

    #[test]
    fn givendist_builds_a_ladder_from_increasing_distances() {
        // Each sequence's nearest neighbour is the one before it, at a
        // strictly greater distance, so every join sits above the last: a
        // comb, deepening left to right.
        let njob = 5;
        let mindists = vec![0.0, 0.10, 0.20, 0.30, 0.40];
        let neighbors = vec![-1, 0, 1, 2, 3];
        let t = compacttreegivendist(njob, &mindists, &neighbors, None, false, false);
        assert_eq!(t.topol.len(), njob);
        // n-1 real joins, then the spare zeroed row.
        assert_eq!(t.topol[njob - 1], [-1, -1]);
        // Post-order over a comb visits the deepest join first, and each node
        // is ordered with the smaller representative on side 0.
        assert_eq!(t.topol[0], [0, 1]);
        for k in 1..njob - 1 {
            assert_eq!(t.topol[k][0], 0, "step {} keeps representative 0", k);
            assert_eq!(t.topol[k][1], (k + 1) as i32);
        }
        // Heights are half the recorded distance and increase up the comb.
        for k in 0..njob - 1 {
            assert!((t.dep[k].distfromtip - mindists[k + 1] * 0.5).abs() < 1e-12);
        }
        // Each join's left child is the join before it; the right child is a leaf.
        assert_eq!(t.dep[0].child0, -1);
        assert_eq!(t.dep[0].child1, -1);
        for k in 1..njob - 1 {
            assert_eq!(t.dep[k].child0, (k - 1) as i32);
            assert_eq!(t.dep[k].child1, -1);
        }
    }

    #[test]
    fn givendist_puts_a_closer_late_pair_below_the_seed() {
        // Sequence 2 is much closer to 0 than 0 and 1 are to each other, so
        // its join has to be inserted *under* the seed join rather than above.
        let mindists = vec![0.0, 0.80, 0.10];
        let neighbors = vec![-1, 0, 0];
        let t = compacttreegivendist(3, &mindists, &neighbors, None, false, false);
        // The 0-2 join is made first (it is deeper), then 0-1 above it.
        assert_eq!(t.topol[0], [0, 2]);
        assert!((t.dep[0].distfromtip - 0.05).abs() < 1e-12);
        assert_eq!(t.topol[1], [0, 1]);
        assert!((t.dep[1].distfromtip - 0.40).abs() < 1e-12);
        // The upper join's side 0 is the lower join, side 1 a fresh leaf.
        assert_eq!(t.dep[1].child0, 0);
        assert_eq!(t.dep[1].child1, -1);
    }

    fn test_constants() -> crate::constants::Constants {
        let a = crate::params::Args {
            nblosum: 62,
            ppenalty: -1530,
            ppenalty_ex: 0,
            poffset: 0,
            penalty_shift_factor: 100.0,
            ..Default::default()
        };
        crate::constants::Constants::new(&a)
    }

    #[test]
    fn memsave_picks_the_smallest_mindist_but_never_the_last_active_element() {
        // Sequence 3 is the last live element of the chain, and it holds the
        // globally smallest `mindist`. The C's scan stops *before* the last
        // element, so the first join is 0-1 -- the smallest among those
        // actually examined -- not 2-3.
        let mut md = vec![0.20, 0.30, 0.40, 0.05];
        let mut nr = vec![1i32, 0, 3, 2];
        let seq: Vec<Vec<u8>> = vec![b"AA".to_vec(); 4];
        let skip: Vec<Vec<i32>> = seq.iter().map(|s| crate::score::makeskiptable(s)).collect();
        let cn = test_constants();
        let src = CompactDist::Msa { cn: &cn, seq: &seq, skiptable: &skip };
        let t = compacttree_memsaveselectable(
            4, &mut nr, &mut md, &[100, 100, 100, 100], &src, 0.1, None, false, false,
        );
        assert_eq!(t.topol[0], [0, 1]);
        assert!((t.dep[0].distfromtip - 0.10).abs() < 1e-12);
    }

    #[test]
    fn memsave_names_a_merged_node_by_its_smaller_representative() {
        let mut md = vec![0.20, 0.20, 0.20, 999.9];
        let mut nr = vec![1i32, 0, 1, 2];
        let seq: Vec<Vec<u8>> = vec![b"AA".to_vec(); 4];
        let skip: Vec<Vec<i32>> = seq.iter().map(|s| crate::score::makeskiptable(s)).collect();
        let cn = test_constants();
        let src = CompactDist::Msa { cn: &cn, seq: &seq, skiptable: &skip };
        let t = compacttree_memsaveselectable(
            4, &mut nr, &mut md, &[100, 100, 100, 100], &src, 0.1, None, false, false,
        );
        assert_eq!(t.topol[0], [0, 1]);
        // Step 1's first side stands for the 0-1 node, and is named by that
        // node's smaller member.
        assert_eq!(t.topol[1][0], 0);
        assert_eq!(t.dep[1].child0, 0);
        // The spare row the C allocates is left untouched.
        assert_eq!(t.topol[3], [-1, -1]);
    }

    #[test]
    fn memsave_leaves_the_merged_node_active_and_unlinks_the_larger() {
        // After every step exactly one sequence leaves the chain, so an
        // n-sequence run makes n-1 joins and every one of them is a real pair.
        let mut md = vec![0.10, 0.10, 0.10, 0.10, 0.10];
        let mut nr = vec![1i32, 0, 3, 2, 0];
        let seq: Vec<Vec<u8>> = vec![b"AC".to_vec(); 5];
        let skip: Vec<Vec<i32>> = seq.iter().map(|s| crate::score::makeskiptable(s)).collect();
        let cn = test_constants();
        let src = CompactDist::Msa { cn: &cn, seq: &seq, skiptable: &skip };
        let t = compacttree_memsaveselectable(
            5, &mut nr, &mut md, &[100; 5], &src, 0.1, None, false, false,
        );
        for k in 0..4 {
            assert!(t.topol[k][0] >= 0 && t.topol[k][1] >= 0, "step {} unfilled", k);
            assert_ne!(t.topol[k][0], t.topol[k][1], "step {} joins a node to itself", k);
        }
    }

    #[test]
    fn givendist_orders_each_node_by_representative_not_by_child_order() {
        // reformat_younger0_rec swaps the two sides when rep0 > rep1, and
        // swaps `len` with them. Sequence 1's neighbour is 2, so the join is
        // built as (rep0=2, rep1=1) and must come back as [1, 2].
        let mindists = vec![0.0, 0.50, 0.20];
        let neighbors = vec![-1, 2, 0];
        let t = compacttreegivendist(3, &mindists, &neighbors, None, false, false);
        for step in 0..2 {
            assert!(
                t.topol[step][0] < t.topol[step][1],
                "step {} is {:?}",
                step,
                t.topol[step]
            );
        }
    }
}

/// Where `compacttree_memsaveselectable` gets a pairwise distance when it
/// needs one. It caches nothing under `howcompact == 2`, so every distance in
/// the loop below is recomputed from one of these.
pub enum CompactDist<'a> {
    /// The first guide-tree pass: k-mer tables, before any alignment exists.
    Kmer {
        pointt: &'a [Vec<i32>],
        nlen: &'a [usize],
        tsize: usize,
        lenfaca: f64,
        lenfacb: f64,
        lenfacc: f64,
        lenfacd: f64,
    },
    /// The second pass: the alignment the first one produced.
    Msa {
        cn: &'a crate::constants::Constants,
        seq: &'a [Vec<u8>],
        skiptable: &'a [Vec<i32>],
    },
}

/// `compacttree_memsaveselectable()` as `disttbfast` calls it: `howcompact 2`,
/// `memsave 1`, which is `--treein youngestlinkage`.
///
/// This is UPGMA over a matrix that is never stored. Each step takes the
/// smallest `mindist`, joins it with its `nearest`, and then recomputes the
/// merged node's distance to every sequence still active -- two real distance
/// evaluations per surviving sequence per step, instead of a matrix lookup.
///
/// Under `howcompact == 2` the C frees the recomputed row instead of keeping
/// it (`partmtx` stays empty from `preparepartmtx` on, since `compacttree` is
/// 4), and `continue`s past the whole `resetnearest` section. So `partmtx` is
/// inert here and is not modelled; its only other use is a progress counter.
///
/// Two things the search depends on, both shared with the ordinary UPGMA:
/// the scan for the smallest `mindist` stops at the *last* active element
/// rather than passing it, and `im`/`jm` persist across steps, so a step where
/// nothing beats the running minimum reuses the previous step's pair.
#[allow(clippy::too_many_arguments)]
pub fn compacttree_memsaveselectable(
    nseq: usize,
    nearest: &mut [i32],
    mindist: &mut [f64],
    selfscore: &[i32],
    src: &CompactDist,
    sueff_global: f64,
    names: Option<&[Vec<u8>]>,
    outnumber: bool,
    treeout: bool,
) -> Tree {
    let sueff1 = 1.0 - sueff_global;
    let sueff05 = sueff_global * 0.5;
    // `cluster_mix_double`, the only `treemethod` the C accepts here.
    let cluster_mix = |d1: f64, d2: f64| -> f64 {
        let m = if d1 < d2 { d1 } else { d2 };
        m * sueff1 + (d1 + d2) * sueff05
    };

    let mut sextet = match src {
        CompactDist::Kmer { tsize, nlen, .. } => Some(crate::dist::Sextet::new(
            *tsize,
            nlen.iter().copied().max().unwrap_or(0),
        )),
        CompactDist::Msa { .. } => None,
    };

    // The two composition tables of the pair being merged, rebuilt each step.
    let mut table1: Vec<i32> = Vec::new();
    let mut table2: Vec<i32> = Vec::new();
    if let CompactDist::Kmer { tsize, .. } = src {
        table1 = vec![0i32; *tsize];
        table2 = vec![0i32; *tsize];
    }

    let mut topol: Vec<[i32; 2]> = vec![[-1, -1]; nseq];
    let mut len: Vec<[f64; 2]> = vec![[0.0, 0.0]; nseq];
    let mut dep: Vec<Treedep> = vec![Treedep::default(); nseq];

    // `ac`: the doubly linked list of still-active sequences. Position 0 is
    // never unlinked -- `jm` is always the larger of the pair, so at least 1 --
    // which is why the C can write through `ac[jm].prev` unguarded.
    let mut next: Vec<i32> = (0..nseq as i32).map(|i| i + 1).collect();
    next[nseq - 1] = -1;
    let mut prev: Vec<i32> = (0..nseq as i32).map(|i| i - 1).collect();

    let mut tmptmplen = vec![0.0f64; nseq];
    let mut hist = vec![-1i32; nseq];
    let mut nmemar = vec![1i32; nseq];

    let mut tree: Vec<Option<String>> = match (treeout, names) {
        (true, Some(names)) => (0..nseq)
            .map(|i| Some(format!("\n{}_{}\n", i + 1, crate::tree::newick_name(&names[i], outnumber))))
            .collect(),
        _ => vec![None; nseq],
    };
    let mut treetmp = String::new();

    // `im` is declared outside the loop and initialised to -1: a step that
    // finds nothing below 999.9 reuses the previous step's `im` -- including
    // the value the swap below left there. `jm` is reassigned from
    // `nearest[im]` before every use, so only `im` has to survive the step.
    let mut im: i32 = -1;

    for k in 0..nseq - 1 {
        let mut minscore = 999.9f64;
        // `for( acpti=ac; acpti->next!=NULL; acpti=acpti->next )`: the last
        // live element is never examined.
        let mut p = 0i32;
        while next[p as usize] != -1 {
            if mindist[p as usize] < minscore {
                im = p;
                minscore = mindist[p as usize];
            }
            p = next[p as usize];
        }
        let mut jm = nearest[im as usize];
        if jm < im {
            std::mem::swap(&mut im, &mut jm);
        }
        let (imu, jmu) = (im as usize, jm as usize);

        // `topol[k][side]`, memsave form: a node is named by the smaller of
        // its previous node's two representatives.
        let prev0 = hist[imu];
        dep[k].child0 = prev0;
        let nmemim = nmemar[imu];
        topol[k][0] = if prev0 == -1 {
            im
        } else {
            topol[prev0 as usize][0].min(topol[prev0 as usize][1])
        };

        let prev1 = hist[jmu];
        dep[k].child1 = prev1;
        let nmemjm = nmemar[jmu];
        topol[k][1] = if prev1 == -1 {
            jm
        } else {
            topol[prev1 as usize][0].min(topol[prev1 as usize][1])
        };

        let minscore = minscore * 0.5;
        len[k][0] = minscore - tmptmplen[imu];
        len[k][1] = minscore - tmptmplen[jmu];
        dep[k].distfromtip = minscore;
        tmptmplen[imu] = minscore;
        hist[imu] = k as i32;
        nmemar[imu] = nmemim + nmemjm;
        mindist[imu] = 999.9;

        if let CompactDist::Kmer { pointt, .. } = src {
            table1.iter_mut().for_each(|v| *v = 0);
            table2.iter_mut().for_each(|v| *v = 0);
            crate::dist::makecompositiontable_p(&mut table1, &pointt[imu]);
            crate::dist::makecompositiontable_p(&mut table2, &pointt[jmu]);
        }

        let mut newarr = vec![0.0f64; nseq];
        let joblist: Vec<i32> = {
            let mut v = Vec::new();
            let mut p = 0i32;
            while p != -1 {
                v.push(p);
                p = next[p as usize];
            }
            v
        };

        // `verycompact{kmer,msa}distarrthreadjoblist`, serial.
        for &ji in &joblist {
            let i = ji as usize;
            if ji == im || ji == jm {
                continue;
            }
            let (d1, d2) = match src {
                CompactDist::Kmer {
                    pointt,
                    nlen,
                    lenfaca,
                    lenfacb,
                    lenfacc,
                    lenfacd,
                    ..
                } => {
                    let sx = sextet.as_mut().unwrap();
                    let a = distcompact(
                        nlen[imu], nlen[i], sx.common(&table1, &pointt[i]),
                        selfscore[imu], selfscore[i],
                        *lenfaca, *lenfacb, *lenfacc, *lenfacd,
                    );
                    let b = distcompact(
                        nlen[jmu], nlen[i], sx.common(&table2, &pointt[i]),
                        selfscore[jmu], selfscore[i],
                        *lenfaca, *lenfacb, *lenfacc, *lenfacd,
                    );
                    (a, b)
                }
                CompactDist::Msa { cn, seq, skiptable } => {
                    let a = distcompact_msa(
                        crate::score::naivepairscorefast(
                            cn, &seq[imu], &seq[i], &skiptable[imu], &skiptable[i], cn.penalty_dist,
                        ),
                        selfscore[imu],
                        selfscore[i],
                    );
                    let b = distcompact_msa(
                        crate::score::naivepairscorefast(
                            cn, &seq[jmu], &seq[i], &skiptable[jmu], &skiptable[i], cn.penalty_dist,
                        ),
                        selfscore[jmu],
                        selfscore[i],
                    );
                    (a, b)
                }
            };
            let tmpdouble = cluster_mix(d1, d2);
            newarr[i] = tmpdouble;
            if tmpdouble < mindist[i] {
                mindist[i] = tmpdouble;
                nearest[i] = im;
            }
            if nearest[i] == jm {
                nearest[i] = im;
            }
        }

        // The merged node's own nearest neighbour, in a second pass -- the C
        // notes it cannot be done inside the worker, which is parallel there.
        for &ji in &joblist {
            let i = ji as usize;
            if ji != im && ji != jm && newarr[i] < mindist[imu] {
                mindist[imu] = newarr[i];
                nearest[imu] = ji;
            }
        }

        if treeout && names.is_some() {
            treetmp = format!(
                "({}:{:7.5},{}:{:7.5})",
                tree[imu].as_deref().unwrap_or(""),
                len[k][0],
                tree[jmu].as_deref().unwrap_or(""),
                len[k][1]
            );
            tree[imu] = Some(treetmp.clone());
            tree[jmu] = None;
        }

        // Unlink jm. jm >= 1 always, so `prev[jm]` is a real element.
        let pj = prev[jmu];
        let nj = next[jmu];
        next[pj as usize] = nj;
        if nj != -1 {
            prev[nj as usize] = pj;
        }
    }

    if treeout && names.is_some() {
        std::fs::write("infile.tree", format!("{};\n", treetmp)).expect("Cannot open infile.tree");
    }

    Tree { topol, len, dep }
}
