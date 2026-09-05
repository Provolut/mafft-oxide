//! Branch weighting on the guide tree, transliterated from
//! `core/treeOperation.c`.
//!
//! `tditeration` refines one branch at a time, and to do that it needs the tree
//! as an *unrooted* graph rather than as the merge sequence `topol` records: for
//! a given branch it asks which sequences fall on each side, how far each one is
//! from the branch, and what weight the path to it carries. `treeCnv` builds
//! that graph (`stopol`), `calcBranchWeight` fills in the per-branch weights,
//! and `distFromABranch` / `weightFromABranch` walk it from one branch outwards.
//!
//! The C holds the graph as `Node`s pointing at each other, with `members[]`
//! aliasing rows of `topol` and `weightptr[]` pointing into cells of `bw`. Here
//! the pointers are indices and the member lists are owned copies: nothing
//! mutates a member list while the graph is alive, and `weightptr` is read only
//! after `calcBranchWeight` has written the cell, so the observable behaviour is
//! the same. `tmpChildren` is gone -- it is the C's own index copy of
//! `children`, kept in lockstep with it at every assignment.
//!
//! The C's member lists are `-1`-terminated and that terminator is kept here.
//! It is not decoration: `negativeMember2` writes it, and every walk stops on
//! it.

/// `Node` from `mltaln.h`, as an index graph.
#[derive(Clone)]
pub struct Node {
    /// Index into the node array, `-1` for "no child" (the C's NULL).
    pub children: [i32; 3],
    pub length: [f64; 3],
    /// The C's `weightptr[]`: which cell of `bw` this arm's weight lives in,
    /// as `(step, lr)`.
    pub weightptr: [Option<(usize, usize)>; 3],
    /// The C's `members[]`, `-1`-terminated.
    pub members: [Option<Vec<i32>>; 3],
}

impl Default for Node {
    /// `stopolInit()`: every arm empty.
    fn default() -> Node {
        Node {
            children: [-1, -1, -1],
            length: [0.0, 0.0, 0.0],
            weightptr: [None, None, None],
            members: [None, None, None],
        }
    }
}

/// `NodeInCub` from `mltaln.h`: a branch, named by merge step and side.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NodeInCub {
    pub step: usize,
    pub lorr: usize,
}

/// The full topology: per merge step, the two member lists, each
/// `-1`-terminated. This is the C's `topol`, not the port's memsave `Tree`,
/// which keeps only each side's first member.
pub type Topol = Vec<[Vec<i32>; 2]>;

/// `MINLEN` from `treeOperation.c`.
const MINLEN: f64 = 0.001;
/// `EF_THREEWAY`, `MAXBW` and `MINBW`.
const EF_THREEWAY: f64 = 1.0;
const MAXBW: f64 = 1.0;
const MINBW: f64 = 0.01;

/// `checkMinusLength()`: clamp every branch length up to `MINLEN`.
///
/// The C comments this "uwagaki" -- overwrite -- and means it: `treeCnv` edits
/// the caller's `len` before reading it, so a later call sees the clamped
/// values.
pub fn check_minus_length(len: &mut [[f64; 2]]) {
    for row in len.iter_mut() {
        for v in row.iter_mut() {
            if *v < MINLEN {
                *v = MINLEN;
            }
        }
    }
}

/// `IntExistsInVec()`: is `query` in this `-1`-terminated vector?
fn int_exists_in_vec(query: i32, vector: &[i32]) -> bool {
    for &v in vector {
        if v == -1 {
            return false;
        }
        if v == query {
            return true;
        }
    }
    false
}

/// `negativeMember2()`: the complement of `query` within `0..locnseq`,
/// `-1`-terminated.
fn negative_member2(query: &[i32], locnseq: usize) -> Vec<i32> {
    let mut tmp = vec![false; locnseq];
    for &n in query {
        if n == -1 {
            break;
        }
        tmp[n as usize] = true;
    }
    let mut out: Vec<i32> = (0..locnseq as i32).filter(|&i| !tmp[i as usize]).collect();
    out.push(-1);
    out
}

/// `searchParent()`: the first branch in `Start..End` whose member list
/// contains `top`.
///
/// The C exits the process when nothing matches; a caller that reaches that has
/// handed in a topology that is not a tree.
fn search_parent(top: i32, topol: &Topol, start: usize, end: usize) -> NodeInCub {
    for i in start..end {
        for j in 0..2 {
            if int_exists_in_vec(top, &topol[i][j]) {
                return NodeInCub { step: i, lorr: j };
            }
        }
    }
    eprintln!("ERROR!!!");
    eprintln!("Error in searchParent");
    std::process::exit(1);
}

/// `isLeaf()`: a node with no second child.
fn is_leaf(node: &Node) -> bool {
    node.children[1] < 0
}

/// `treeCnv()`: build the unrooted graph from the merge sequence.
///
/// `len` is clamped in place first (`checkMinusLength`), and `bw` is not read
/// here -- only the cell each arm's weight will live in is recorded, for
/// `calcBranchWeight` to fill and `weightFromABranch` to read.
///
/// Args:
///     locnseq: Number of sequences.
///     topol: The merge sequence, `locnseq - 1` steps.
///     len: Branch lengths, clamped in place.
///
/// Returns:
///     `2 * locnseq` nodes: `0..locnseq-1` internal, `locnseq..` the leaves.
pub fn tree_cnv(locnseq: usize, topol: &Topol, len: &mut [[f64; 2]]) -> Vec<Node> {
    assert!(
        locnseq >= 3,
        "treeCnv needs at least three sequences; the C indexes stopol[locnseq-3] \
         and is undefined below that"
    );
    check_minus_length(len);

    let mut stopol = vec![Node::default(); locnseq * 2];
    let mut count = vec![0usize; 2 * locnseq];

    for i in locnseq..locnseq * 2 {
        let rep = (i - locnseq) as i32;
        let parent = search_parent(rep, topol, 0, locnseq - 1);
        let arm_len = len[parent.step][parent.lorr];
        let members = topol[parent.step][parent.lorr].clone();

        let c = count[parent.step];
        stopol[parent.step].length[c] = arm_len;
        stopol[parent.step].weightptr[c] = Some((parent.step, parent.lorr));
        stopol[parent.step].children[c] = i as i32;
        stopol[parent.step].members[c] = Some(members.clone());
        count[parent.step] += 1;

        let c = count[i];
        stopol[i].length[c] = arm_len;
        stopol[i].weightptr[c] = Some((parent.step, parent.lorr));
        stopol[i].children[c] = parent.step as i32;
        stopol[i].members[c] = Some(members);
        count[i] += 1;
    }

    for i in 0..locnseq.saturating_sub(2) {
        let rep = topol[i][0][0].min(topol[i][1][0]);
        let parent = search_parent(rep, topol, i + 1, locnseq - 1);
        let arm_len = len[parent.step][parent.lorr];

        let c = count[parent.step];
        stopol[parent.step].length[c] = arm_len;
        stopol[parent.step].weightptr[c] = Some((parent.step, parent.lorr));
        stopol[parent.step].children[c] = i as i32;
        stopol[parent.step].members[c] = Some(topol[parent.step][parent.lorr].clone());
        count[parent.step] += 1;

        let c = count[i];
        stopol[i].length[c] = arm_len;
        stopol[i].weightptr[c] = Some((parent.step, parent.lorr));
        stopol[i].children[c] = parent.step as i32;
        // The other side of the same branch: everything the parent's list does
        // not hold. The C allocates one of these per step and frees them all
        // when `treeCnv(NULL, ...)` is called; owning it here retires that
        // bookkeeping.
        stopol[i].members[c] = Some(negative_member2(
            &topol[parent.step][parent.lorr],
            locnseq,
        ));
        count[i] += 1;
    }

    // Unroot: the last merge step has no branch of its own, so its two arms
    // become one arm of `locnseq-3`, carrying the summed length.
    let root = locnseq - 2;
    let far = locnseq - 3;
    let i = if stopol[root].children[0] == far as i32 {
        1
    } else if stopol[root].children[1] == far as i32 {
        0
    } else {
        eprintln!("ERROR in stopol ?");
        std::process::exit(1);
    };
    let joined = len[root][0] + len[root][1];
    let other = stopol[root].children[i];

    stopol[far].length[2] = joined;
    stopol[far].weightptr[2] = Some((root, 0));
    stopol[far].children[2] = other;

    let other = other as usize;
    stopol[other].children[2] = far as i32;
    stopol[other].length[2] = joined;
    stopol[other].weightptr[2] = Some((root, 0));

    stopol
}

/// `syntheticLength()`: the resistance of the subtree hanging off `ob`, seen
/// from `op`.
///
/// Two arms in parallel, as in an electrical network, plus the arm back towards
/// `op`. A zero-length arm anywhere below collapses the whole subtree to zero
/// rather than to infinity -- the C's `nanflag`, which is the guard against
/// dividing by a zero branch and is why the result is not simply the harmonic
/// sum.
fn synthetic_length(stopol: &[Node], ob: usize, op: i32) -> f64 {
    if is_leaf(&stopol[ob]) {
        return stopol[ob].length[0];
    }
    let mut dir_ch = [0usize; 3];
    let mut count = 0usize;
    let mut dir_pa = usize::MAX;
    for i in 0..3 {
        if stopol[ob].children[i] != op {
            dir_ch[count] = i;
            count += 1;
        } else {
            dir_pa = i;
        }
    }
    if count != 2 {
        eprintln!("Invalid call");
        std::process::exit(1);
    }

    let mut t0 = synthetic_length(stopol, stopol[ob].children[dir_ch[0]] as usize, ob as i32);
    let mut t1 = synthetic_length(stopol, stopol[ob].children[dir_ch[1]] as usize, ob as i32);

    let mut nanflag = false;
    if t0 != 0.0 {
        t0 = 1.0 / t0;
    } else {
        nanflag = true;
    }
    if t1 != 0.0 {
        t1 = 1.0 / t1;
    } else {
        nanflag = true;
    }

    let mut value = if nanflag { 0.0 } else { 1.0 / (t0 + t1) };
    value += stopol[ob].length[dir_pa];
    value
}

/// `calcW()`: the three-way weight of the branch between `ob` and `op`.
fn calc_w(stopol: &[Node], ob: usize, op: i32) -> f64 {
    if is_leaf(&stopol[ob]) {
        return 1.0;
    }
    let mut dir_ch = [0usize; 3];
    let mut count = 0usize;
    let mut dir_pa = usize::MAX;
    for i in 0..3 {
        if stopol[ob].children[i] != op {
            dir_ch[count] = i;
            count += 1;
        } else {
            dir_pa = i;
        }
    }
    if count != 2 {
        eprintln!("Invalid call of calcW");
        std::process::exit(1);
    }

    let a = synthetic_length(stopol, stopol[ob].children[dir_ch[0]] as usize, ob as i32);
    let b = synthetic_length(stopol, stopol[ob].children[dir_ch[1]] as usize, ob as i32);
    let c = synthetic_length(stopol, stopol[ob].children[dir_pa] as usize, ob as i32);

    if c == 0.0 {
        return MAXBW;
    }
    if a == 0.0 || b == 0.0 {
        return MINBW;
    }

    let f = EF_THREEWAY;
    let s = b * c + c * a + a * b;
    let value = a * b * (c + a) * (c + b) / (c * (a + b) * f * s);
    value.sqrt()
}

/// `calcBranchWeight()`: fill `bw` with each branch's weight.
///
/// Args:
///     bw: `locnseq - 1` rows of two, written here.
///     locnseq: Number of sequences.
///     stopol: The graph from [`tree_cnv`].
///     topol: The merge sequence.
pub fn calc_branch_weight(bw: &mut [[f64; 2]], locnseq: usize, stopol: &[Node], topol: &Topol) {
    for i in locnseq..locnseq * 2 {
        let rep = (i - locnseq) as i32;
        let parent = search_parent(rep, topol, 0, locnseq - 1);
        if parent.step == locnseq - 2 {
            continue;
        }
        let top = parent.step;
        let btm = i;
        bw[parent.step][parent.lorr] =
            calc_w(stopol, top, btm as i32) * calc_w(stopol, btm, top as i32);
    }
    for i in 0..locnseq.saturating_sub(3) {
        let rep = topol[i][0][0].min(topol[i][1][0]);
        let parent = search_parent(rep, topol, i + 1, locnseq - 1);
        if parent.step == locnseq - 2 {
            continue;
        }
        let top = parent.step;
        bw[parent.step][parent.lorr] =
            calc_w(stopol, top, i as i32) * calc_w(stopol, i, top as i32);
    }

    // `btmNode = stopol + i` reads the loop counter *after* the loop, so this
    // last branch hangs off node `locnseq-3` -- the bound the loop exited on,
    // one past the last index it visited. Rust's `for` does not leak its
    // counter, so the value is written out; reading it as "the last index
    // visited" would pick the wrong node.
    let i = locnseq.saturating_sub(3);
    let top = stopol[locnseq - 3].children[2];
    bw[locnseq - 2][0] = calc_w(stopol, top as usize, i as i32) * calc_w(stopol, i, top);
    bw[locnseq - 2][1] = 1.0;
}

/// `distFromABranch_rec()`: add each arm's length to every member below it.
fn dist_from_a_branch_rec(result: &mut [f64], stopol: &[Node], ob: usize, op: i32) {
    if is_leaf(&stopol[ob]) {
        return;
    }
    let mut dir_ch = [0usize; 3];
    let mut count = 0usize;
    for i in 0..3 {
        if stopol[ob].children[i] != op {
            dir_ch[count] = i;
            count += 1;
        }
    }
    if count != 2 {
        eprintln!("Incorrect call of distFromABranch_rec");
        std::process::exit(1);
    }
    for k in 0..2 {
        let d = dir_ch[k];
        for &n in stopol[ob].members[d].as_ref().expect("member list") {
            if n == -1 {
                break;
            }
            result[n as usize] += stopol[ob].length[d];
        }
        dist_from_a_branch_rec(result, stopol, stopol[ob].children[d] as usize, ob as i32);
    }
}

/// Pick the two nodes a branch joins, given the branch as `(step, lorr)`.
fn branch_ends(nseq: usize, stopol: &[Node], topol: &Topol, step: usize, lorr: usize) -> (usize, usize) {
    if step == nseq - 2 {
        // The unrooted graph has no node for the last merge step; the branch it
        // named is the one `tree_cnv` joined onto `nseq-3`.
        (stopol[nseq - 2].children[0] as usize, nseq - 3)
    } else {
        let mut i = 3;
        for k in 0..3 {
            if let Some(m) = stopol[step].members[k].as_ref() {
                if m[0] == topol[step][lorr][0] {
                    i = k;
                    break;
                }
            }
        }
        if i == 3 {
            eprintln!("Incorrect call of distFromABranch.");
            std::process::exit(1);
        }
        (step, stopol[step].children[i] as usize)
    }
}

/// `distFromABranch()`: distance from one branch to every sequence.
///
/// Args:
///     nseq: Number of sequences.
///     result: `nseq` distances, overwritten.
///     stopol: The graph from [`tree_cnv`].
///     topol: The merge sequence.
///     len: Branch lengths, read for the two-sequence case only.
///     step: Which merge step names the branch.
///     lorr: Which side of it.
pub fn dist_from_a_branch(
    nseq: usize,
    result: &mut [f64],
    stopol: &[Node],
    topol: &Topol,
    len: &[[f64; 2]],
    step: usize,
    lorr: usize,
) {
    if nseq == 2 {
        result[0] = len[0][0];
        result[1] = len[0][1];
        return;
    }
    let (top, btm) = branch_ends(nseq, stopol, topol, step, lorr);
    for v in result.iter_mut().take(nseq) {
        *v = 0.0;
    }
    dist_from_a_branch_rec(result, stopol, btm, top as i32);
    dist_from_a_branch_rec(result, stopol, top, btm as i32);
}

/// `weightFromABranch_rec()`: multiply each member below an arm by its weight.
fn weight_from_a_branch_rec(
    result: &mut [f64],
    stopol: &[Node],
    bw: &[[f64; 2]],
    ob: usize,
    op: i32,
) {
    if is_leaf(&stopol[ob]) {
        return;
    }
    let mut dir_ch = [0usize; 3];
    let mut count = 0usize;
    for i in 0..3 {
        if stopol[ob].children[i] != op {
            dir_ch[count] = i;
            count += 1;
        }
    }
    if count != 2 {
        eprintln!("Incorrect call of weightFromABranch_rec");
        std::process::exit(1);
    }
    for k in 0..2 {
        let d = dir_ch[k];
        let (s, l) = stopol[ob].weightptr[d].expect("weightptr");
        let w = bw[s][l];
        for &n in stopol[ob].members[d].as_ref().expect("member list") {
            if n == -1 {
                break;
            }
            result[n as usize] *= w;
        }
        weight_from_a_branch_rec(result, stopol, bw, stopol[ob].children[d] as usize, ob as i32);
    }
}

/// `weightFromABranch()`: the weight of the path from one branch to every
/// sequence.
///
/// `bw` stands in for the C's `weightptr` dereference; `calc_branch_weight` must
/// have filled it first, exactly as the C requires.
///
/// Args:
///     nseq: Number of sequences.
///     result: `nseq` weights, overwritten.
///     stopol: The graph from [`tree_cnv`].
///     bw: Branch weights from [`calc_branch_weight`].
///     topol: The merge sequence.
///     step: Which merge step names the branch.
///     lorr: Which side of it.
pub fn weight_from_a_branch(
    nseq: usize,
    result: &mut [f64],
    stopol: &[Node],
    bw: &[[f64; 2]],
    topol: &Topol,
    step: usize,
    lorr: usize,
) {
    let (top, btm) = branch_ends(nseq, stopol, topol, step, lorr);
    for v in result.iter_mut().take(nseq) {
        *v = 1.0;
    }
    weight_from_a_branch_rec(result, stopol, bw, btm, top as i32);
    weight_from_a_branch_rec(result, stopol, bw, top, btm as i32);
}

/// `OneClusterAndTheOther_fast()`: split the sequences across one branch.
///
/// One side is the branch's member list from `topol`; the other is everything
/// else, in ascending order. Both come back `-1`-terminated, and the first
/// entry of each is the representative `tditeration` aligns around.
///
/// The C's `smalldistmtx` output is not produced here: `tditeration` allocates
/// it only when `specificityconsideration != 0`, and the `mafft` script runs
/// `dvtditr -s 0.0`. See `GOAL-C-FREE.md` for the other refusals that follow
/// from that flag.
///
/// Args:
///     locnjob: Number of sequences.
///     topol: The merge sequence.
///     step: Which merge step names the branch.
///     branch: Which side of it.
///
/// Returns:
///     The two member lists, each `-1`-terminated.
pub fn one_cluster_and_the_other_fast(
    locnjob: usize,
    topol: &Topol,
    step: usize,
    branch: usize,
) -> (Vec<i32>, Vec<i32>) {
    let mut pair = vec![false; locnjob];
    let mut memlist1 = Vec::with_capacity(locnjob + 1);
    for &r1 in &topol[step][branch] {
        if r1 <= -1 {
            break;
        }
        pair[r1 as usize] = true;
        memlist1.push(r1);
    }
    memlist1.push(-1);

    let mut memlist2: Vec<i32> = (0..locnjob as i32)
        .filter(|&i| !pair[i as usize])
        .collect();
    memlist2.push(-1);

    (memlist1, memlist2)
}

/// `GETA3` from `mltaln.h`.
const GETA3: f64 = 0.001;

/// `counteff_simple()` (`mltaln9.c`), the variant `tditeration` calls.
///
/// Not the one `tree.rs` ports: that is the memsave form, which rebuilds each
/// side's members from `Treedep` and clamps negative branch lengths with a
/// warning first. This one reads the full member lists straight out of `topol`
/// and does no clamping, because by the time `tditeration` calls it `treeCnv`
/// has already raised every length to `MINLEN` in place.
///
/// Args:
///     nseq: Number of sequences.
///     topol: The merge sequence.
///     len: Branch lengths.
///
/// Returns:
///     One weight per sequence, summing to 1.
pub fn counteff_simple(nseq: usize, topol: &Topol, len: &[[f64; 2]]) -> Vec<f64> {
    let mut rootnode = vec![0.0f64; nseq];
    let mut eff = vec![1.0f64; nseq];

    for i in 0..nseq - 1 {
        for side in 0..2 {
            for &s in &topol[i][side] {
                if s <= -1 {
                    break;
                }
                let s = s as usize;
                rootnode[s] += len[i][side] * eff[s];
                eff[s] *= 0.5;
            }
        }
    }
    for v in rootnode.iter_mut() {
        *v += GETA3;
    }
    let total: f64 = rootnode.iter().sum();
    rootnode.iter().map(|v| v / total).collect()
}

/// Expand a memsave `Tree` into the full member lists `tditeration` needs.
///
/// `dvtditr` builds its tree with `fixed_musclesupg_double_treeout`, which
/// records every member of both sides at every step;
/// `tree::upgma_treeout` ports the *memsave* variant, which records each side's
/// representative and a `Treedep` saying which earlier step it came from. The
/// two C functions are separate copies of one algorithm, and they agree: on 40
/// staged distance matrices — 3 to 36 sequences, with ties, zeros, near-equal
/// and sub-micro distances — every branch length matched bit for bit, and the
/// full member lists were exactly this reconstruction. So the port keeps one
/// tree builder and expands its output here, rather than carrying a second
/// 380-line copy.
///
/// Order matters and is the C's: a step's list is its first child's members
/// followed by its second's, each in the order that child recorded them.
///
/// Args:
///     nseq: Number of sequences.
///     tree: The memsave tree.
///
/// Returns:
///     `nseq` rows, the last a placeholder, each side `-1`-terminated.
pub fn expand_topol(nseq: usize, tree: &crate::tree::Tree) -> Topol {
    let mut memhist: Vec<Vec<i32>> = vec![Vec::new(); nseq.saturating_sub(1)];
    let mut out: Topol = Vec::with_capacity(nseq);
    for i in 0..nseq - 1 {
        let m0: Vec<i32> = if tree.dep[i].child0 == -1 {
            vec![tree.topol[i][0]]
        } else {
            memhist[tree.dep[i].child0 as usize].clone()
        };
        let m1: Vec<i32> = if tree.dep[i].child1 == -1 {
            vec![tree.topol[i][1]]
        } else {
            memhist[tree.dep[i].child1 as usize].clone()
        };
        let mut merged = m0.clone();
        merged.extend_from_slice(&m1);
        memhist[i] = merged;

        let mut s0 = m0;
        s0.push(-1);
        let mut s1 = m1;
        s1.push(-1);
        out.push([s0, s1]);
    }
    // The C allocates `nseq` rows and leaves the last unused.
    out.push([vec![-1], vec![-1]]);
    out
}

/// `SEGMENTSIZE` from `fftFunctions.c`, as `searchAnchors` uses it.
const SEGMENTSIZE: i32 = 150;

/// `searchAnchors()` (`mltaln9.c`): where an alignment may be cut in two.
///
/// `dvtditr` refines a long alignment in segments rather than whole, and this
/// picks the cut points: a sliding window over the column-wise mean pairwise
/// score, and every run of columns scoring above the threshold becomes a
/// segment whose centre is a cut. Conserved blocks are what it finds, which is
/// the point — cutting through one would be cutting through the part the
/// alignment is most sure about.
///
/// Only reached unconstrained and with segmentation on: `-l` (which every
/// constrained call passes) turns it off.
///
/// Args:
///     cn: Scoring constants, for `n_dis` and the alphabet map.
///     seq: The alignment.
///     div_win_size: `divWinSize`, 20.
///     div_threshold: `divThreshold`, 65.
///
/// Returns:
///     Each segment's centre column, in order.
pub fn search_anchors(
    cn: &crate::constants::Constants,
    seq: &[Vec<u8>],
    div_win_size: i32,
    div_threshold: i32,
) -> Vec<i32> {
    let nseq = seq.len();
    let len = seq[0].len();
    // The C computes the threshold in ints and then divides, so the truncation
    // is part of the answer.
    let threshold = f64::from(div_threshold) / 100.0 * 600.0 * f64::from(div_win_size);

    let mut stra = vec![0.0f64; len];
    for (i, v) in stra.iter_mut().enumerate() {
        let mut acc = 0.0f64;
        for k in 0..nseq - 1 {
            for j in k + 1..nseq {
                let a = cn.amino_n[seq[k][i] as usize] as usize;
                let b = cn.amino_n[seq[j][i] as usize] as usize;
                acc += f64::from(cn.n_dis[a][b]);
            }
        }
        *v = acc / (nseq as f64 * (nseq as f64 - 1.0) / 2.0);
    }

    let mut centers = Vec::new();
    let mut status = false;
    let mut cumscore = 0.0f64;
    let mut length = 0i32;
    let mut start = 0i32;
    let mut score: f64 = stra.iter().take(div_win_size as usize).sum();

    let mut i = 1i32;
    while i < len as i32 - div_win_size {
        score = score - stra[(i - 1) as usize] + stra[(i + div_win_size - 1) as usize];
        if score > threshold {
            if !status {
                status = true;
                start = i;
                length = 0;
                cumscore = 0.0;
            }
            length += 1;
            cumscore += score;
        }
        // Not an `else`: a run longer than SEGMENTSIZE is closed on the same
        // column that extended it.
        if score <= threshold || length > SEGMENTSIZE {
            if status {
                centers.push((start + i + div_win_size) / 2);
                length = 0;
                cumscore = 0.0;
                status = false;
            }
        }
        i += 1;
    }
    let _ = cumscore;
    if status {
        // `i` here is the loop counter's exit value, as in the C -- the first
        // column the loop did not process, not the last one it did.
        centers.push((start + i + div_win_size) / 2);
    }
    centers
}
