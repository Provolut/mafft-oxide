//! Guide-tree construction and sequence weighting.
//!
//! `upgma()` is `fixed_musclesupg_double_realloc_nobk_halfmtx_memsave()` from
//! `core/mltaln9.c`; `counteff_simple()` is
//! `counteff_simple_double_nostatic_memsave()` from the same file.  The
//! linked-list bookkeeping of `Bchain` is reproduced with index arrays because
//! the *order* in which candidates are visited decides ties.

/// One internal node of the guide tree (`Treedep` in mltaln.h).
#[derive(Clone, Copy, Debug)]
pub struct Treedep {
    pub child0: i32,
    pub child1: i32,
    pub distfromtip: f64,
}

impl Default for Treedep {
    fn default() -> Self {
        Treedep {
            child0: -1,
            child1: -1,
            distfromtip: 0.0,
        }
    }
}

/// The guide tree: for each of the `nseq-1` merge steps, the representative of
/// each child group, the two branch lengths, and the dependency record.
pub struct Tree {
    /// `topol[k][side]` — the C keeps only the first member ("memsave").
    pub topol: Vec<[i32; 2]>,
    pub len: Vec<[f64; 2]>,
    pub dep: Vec<Treedep>,
}

/// A doubly linked list over `0..nseq`, mirroring `Bchain`.
struct Chain {
    next: Vec<i32>,
    prev: Vec<i32>,
}

impl Chain {
    fn new(nseq: usize) -> Chain {
        let mut next = vec![0i32; nseq];
        let mut prev = vec![0i32; nseq];
        for i in 0..nseq {
            next[i] = if i + 1 < nseq { (i + 1) as i32 } else { -1 };
            prev[i] = i as i32 - 1;
        }
        Chain { next, prev }
    }

    /// Unlink `jm`, as the C does with `acjmprev`/`acjmnext`.
    fn unlink(&mut self, jm: usize) {
        let prev = self.prev[jm];
        let next = self.next[jm];
        self.next[prev as usize] = next;
        if next != -1 {
            self.prev[next as usize] = prev;
        }
    }
}

/// `setnearest()`: nearest live neighbour of `pos` in the half matrix.
fn setnearest(ac: &Chain, eff: &[Vec<f64>], pos: usize) -> (f64, i32) {
    let mut mindisfrom = 999.9f64;
    let mut nearest = -1i32;

    let mut j = ac.next[pos];
    while j != -1 {
        let ju = j as usize;
        let tmp = eff[pos][ju - pos];
        if tmp < mindisfrom {
            mindisfrom = tmp;
            nearest = j;
        }
        j = ac.next[ju];
    }

    let mut j = 0i32; // the C always restarts from the list head, ac[0]
    while j != -1 && j as usize != pos {
        let ju = j as usize;
        let tmp = eff[ju][pos - ju];
        if tmp < mindisfrom {
            mindisfrom = tmp;
            nearest = j;
        }
        j = ac.next[ju];
    }

    (mindisfrom, nearest)
}

/// `fixed_musclesupg_double_realloc_nobk_halfmtx_memsave()`.
///
/// `eff` is the half distance matrix and is consumed (the C frees rows as it
/// goes, and overwrites surviving ones in place).
pub fn upgma(nseq: usize, eff: &mut Vec<Vec<f64>>, sueff_global: f64) -> Tree {
    upgma_treeout(nseq, eff, sueff_global, None, false, false)
}

/// `..._treeout_memsave()`: the same UPGMA, plus a Newick string per live
/// cluster, written to `infile.tree` when `names` is given.
///
/// The C keeps this as a separate copy of the function rather than a flag, and
/// the copies differ in one thing besides the output: `im` starts at 1 there
/// and -1 here. It is only read before being assigned if no live element beats
/// the 999.9 seed on the very first step, which needs a degenerate distance
/// matrix -- but the seed is reproduced anyway.
pub fn upgma_treeout(
    nseq: usize,
    eff: &mut Vec<Vec<f64>>,
    sueff_global: f64,
    names: Option<&[Vec<u8>]>,
    outnumber: bool,
    nodeout: bool,
) -> Tree {
    let sueff1 = 1.0 - sueff_global;
    let sueff05 = sueff_global * 0.5;
    let clusterfunc = |d1: f64, d2: f64| -> f64 {
        let m = if d1 < d2 { d1 } else { d2 };
        m * sueff1 + (d1 + d2) * sueff05
    };

    let mut topol = vec![[-1i32, -1i32]; nseq];
    let mut len = vec![[0.0f64, 0.0f64]; nseq];
    let mut dep = vec![Treedep::default(); nseq];

    let mut ac = Chain::new(nseq);
    let mut mindisfrom = vec![0.0f64; nseq];
    let mut nearest = vec![0i32; nseq];
    for i in 0..nseq {
        let (m, n) = setnearest(&ac, eff, i);
        mindisfrom[i] = m;
        nearest[i] = n;
    }

    let mut tmptmplen = vec![0.0f64; nseq];
    let mut hist = vec![-1i32; nseq];
    let mut nmemar = vec![1i32; nseq];

    // The C keeps `im` across iterations, so a step in which nothing beats the
    // 999.9 seed reuses the previous pick.  Kept identical.
    let mut im_g: i32 = if names.is_some() { 1 } else { -1 };

    // One Newick fragment per live cluster; `tree[im]` grows as clusters merge.
    let mut tree: Vec<String> = Vec::new();
    if let Some(name) = names {
        tree = (0..nseq)
            .map(|i| format!("\n{}_{}\n", i + 1, newick_name(&name[i], outnumber)))
            .collect();
    }

    // `-^`: how much of the matrix each sequence is close to, measured once on
    // the full live chain before any merging.
    let density: Vec<f64> = if nodeout {
        (0..nseq).map(|i| setdensity(&ac, eff, i)).collect()
    } else {
        Vec::new()
    };

    for k in 0..nseq - 1 {
        let mut minscore = 999.9f64;
        // C: for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) -- the last
        // live element is skipped, which matters when it holds the minimum.
        let mut i = 0i32;
        while ac.next[i as usize] != -1 {
            if mindisfrom[i as usize] < minscore {
                im_g = i;
                minscore = mindisfrom[i as usize];
            }
            i = ac.next[i as usize];
        }
        let mut jm_g = nearest[im_g as usize];
        if jm_g < im_g {
            std::mem::swap(&mut im_g, &mut jm_g);
        }
        let im = im_g as usize;
        let jm = jm_g as usize;
        let prevnode = hist[im];
        dep[k].child0 = prevnode;
        let nmemim = nmemar[im];
        topol[k][0] = if prevnode == -1 {
            im as i32
        } else {
            let p = prevnode as usize;
            topol[p][0].min(topol[p][1])
        };

        let prevnode = hist[jm];
        dep[k].child1 = prevnode;
        let nmemjm = nmemar[jm];
        topol[k][1] = if prevnode == -1 {
            jm as i32
        } else {
            let p = prevnode as usize;
            topol[p][0].min(topol[p][1])
        };

        let minscore = minscore * 0.5;

        len[k][0] = minscore - tmptmplen[im];
        len[k][1] = minscore - tmptmplen[jm];
        dep[k].distfromtip = minscore;
        tmptmplen[im] = minscore;

        if names.is_some() {
            let merged = format!(
                "({}:{:7.5},{}:{:7.5})",
                tree[im], len[k][0], tree[jm], len[k][1]
            );
            tree[im] = merged;
            tree[jm] = String::new();
        }

        hist[im] = k as i32;
        nmemar[im] = nmemim + nmemjm;

        mindisfrom[im] = 999.9;
        let mut i = 0i32;
        while i != -1 {
            let iu = i as usize;
            if iu != im && iu != jm {
                let (miniim, maxiim, minijm, maxijm) = if iu < im {
                    (iu, im, iu, jm)
                } else if iu < jm {
                    (im, iu, iu, jm)
                } else {
                    (im, iu, jm, iu)
                };
                let eff0 = eff[miniim][maxiim - miniim];
                let eff1 = eff[minijm][maxijm - minijm];
                let tmpdouble = clusterfunc(eff0, eff1);
                eff[miniim][maxiim - miniim] = tmpdouble;
                if tmpdouble < mindisfrom[iu] {
                    mindisfrom[iu] = tmpdouble;
                    nearest[iu] = im as i32;
                }
                if tmpdouble < mindisfrom[im] {
                    mindisfrom[im] = tmpdouble;
                    nearest[im] = iu as i32;
                }
                if nearest[iu] == jm as i32 {
                    nearest[iu] = im as i32;
                }
            }
            i = ac.next[iu];
        }

        // Unlink jm.
        let acjmprev = ac.prev[jm];
        let acjmnext = ac.next[jm];
        ac.next[acjmprev as usize] = acjmnext;
        if acjmnext != -1 {
            ac.prev[acjmnext as usize] = acjmprev;
        }

        // Refresh anyone whose cached nearest was im but is now stale.
        let mut i = 0i32;
        while i != -1 {
            let iu = i as usize;
            if nearest[iu] == im as i32 {
                let (miniim, maxiim) = if iu < im { (iu, im) } else { (im, iu) };
                if eff[miniim][maxiim - miniim] > mindisfrom[iu] {
                    let (m, n) = setnearest(&ac, eff, iu);
                    mindisfrom[iu] = m;
                    nearest[iu] = n;
                }
            }
            i = ac.next[iu];
        }
    }

    if names.is_some() {
        // `treetmp` -- the last merge -- is what gets written, not `tree[0]`.
        let last = tree
            .iter()
            .find(|t| !t.is_empty())
            .cloned()
            .unwrap_or_default();
        let mut body = format!("{};\n", last);
        if nodeout {
            use std::fmt::Write as _;
            let _ = write!(body, "\nDensity:");
            for k in 0..nseq {
                let _ = write!(body, "\nSequence {}, {:7.4}", k + 1, density[k]);
            }
            let _ = write!(body, "\n\nNode info:");
            for k in 0..nseq - 1 {
                let _ = write!(
                    body,
                    "\nNode {}, Height={}\n",
                    k + 1,
                    fmt_f(dep[k].distfromtip)
                );
                for child in 0..2 {
                    let mem = topolorderz(&topol, &dep, k, child);
                    let _ = write!(body, "{}:", getdensest(&mem, &density) + 1);
                    for &m in mem.iter() {
                        let _ = write!(body, " {}", m + 1);
                    }
                    let _ = writeln!(body);
                }
            }
        }
        std::fs::write("infile.tree", body).expect("write infile.tree");
    }

    Tree { topol, len, dep }
}

/// The name sanitisation `..._treeout_memsave()` does before building the
/// Newick string: everything that is not alphanumeric or one of `/=-{}`
/// becomes `_`, the leading `=` the reader added is dropped, and a
/// `_numo_e`/`_oe_` marker moves the start of the name past itself.
pub(crate) fn newick_name(name: &[u8], outnumber: bool) -> String {
    let mut nametmp: Vec<u8> = Vec::with_capacity(name.len());
    for &c in name.iter().take(999) {
        if c == 0 {
            break;
        }
        if c.is_ascii_alphanumeric() || matches!(c, b'/' | b'=' | b'-' | b'{' | b'}') {
            nametmp.push(c);
        } else {
            nametmp.push(b'_');
        }
    }
    let find = |hay: &[u8], needle: &[u8]| hay.windows(needle.len()).position(|w| w == needle);
    let mut start = if outnumber {
        match find(&nametmp, b"_numo_e") {
            Some(p) => p + 8,
            // `strstr` returning NULL and then + 8 is undefined in the C; it
            // cannot happen, because -n always writes the marker.
            None => 1,
        }
    } else {
        1
    };
    if let Some(p) = find(&nametmp[start.min(nametmp.len())..], b"_oe_") {
        start += p + 4;
    }
    let tail = &nametmp[start.min(nametmp.len())..];
    String::from_utf8_lossy(&tail[..tail.len().min(900)]).into_owned()
}

/// `loadtree()`: build the guide tree from `_guidetree` instead of from a
/// distance matrix.
///
/// The file is one line per internal node, `i j len0 len1`, with 1-based
/// sequence numbers and `i < j`; `newick2mafft.rb` produces it from a Newick
/// tree. The bookkeeping is the same as `upgma`'s, only the merge order comes
/// from the file rather than from the minimum. The C's inner loop over the
/// live chain computes four indices and uses none of them -- there is no
/// distance matrix to update here -- so it is not reproduced.
pub fn load_tree(nseq: usize, names: Option<&[Vec<u8>]>, outnumber: bool) -> Result<Tree, String> {
    let text =
        std::fs::read_to_string("_guidetree").map_err(|_| "cannot open _guidetree".to_string())?;
    let mut lines = text.lines();

    let mut topol = vec![[-1i32, -1i32]; nseq];
    let mut len = vec![[0.0f64, 0.0f64]; nseq];
    let mut dep = vec![Treedep::default(); nseq];
    let mut ac = Chain::new(nseq);
    let mut hist = vec![-1i32; nseq];
    let mut nmemar = vec![1i32; nseq];
    let mut height = vec![0.0f64; nseq];

    let mut tree: Vec<String> = Vec::new();
    if let Some(name) = names {
        tree = (0..nseq)
            .map(|i| format!("\n{}_{}\n", i + 1, newick_name(&name[i], outnumber)))
            .collect();
    }
    let mut treetmp = String::new();

    let fmt_err = "\n\nFormat error (1) in the tree?  It has to be a bifurcated and rooted tree.\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n";
    for k in 0..nseq - 1 {
        let line = lines.next().ok_or_else(|| fmt_err.to_string())?;
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 4 {
            return Err("\n\nFormat error (2) in the tree?  It has to be a bifurcated and rooted tree.\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n".into());
        }
        let a0: i64 = f[0].parse().map_err(|_| fmt_err.to_string())?;
        let a1: i64 = f[1].parse().map_err(|_| fmt_err.to_string())?;
        let l0: f64 = f[2].parse().map_err(|_| fmt_err.to_string())?;
        let l1: f64 = f[3].parse().map_err(|_| fmt_err.to_string())?;
        let im = (a0 - 1) as i64;
        let jm = (a1 - 1) as i64;
        if im >= jm {
            return Err("\n\nIncorrect guide tree\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n".into());
        }
        if im > nseq as i64 - 1 || jm > nseq as i64 - 1 {
            return Err(format!(
                "\n\nCheck the guide tree.\nim={}, jm={}\nPlease use newick2mafft.rb to generate a tree file from a newick tree.\n",
                im + 1,
                jm + 1
            ));
        }
        let (im, jm) = (im as usize, jm as usize);
        len[k][0] = if l0 < 0.0 { 0.0 } else { l0 };
        len[k][1] = if l1 < 0.0 { 0.0 } else { l1 };

        let prevnode = hist[im];
        dep[k].child0 = prevnode;
        let nmemim = nmemar[im];
        topol[k][0] = if prevnode == -1 {
            im as i32
        } else {
            let p = prevnode as usize;
            topol[p][0].min(topol[p][1])
        };

        let nmemjm = nmemar[jm];
        let prevnode = hist[jm];
        dep[k].child1 = prevnode;
        topol[k][1] = if prevnode == -1 {
            jm as i32
        } else {
            let p = prevnode as usize;
            topol[p][0].min(topol[p][1])
        };

        hist[im] = k as i32;
        nmemar[im] = nmemim + nmemjm;

        if names.is_some() {
            treetmp = format!(
                "({}:{:7.5},{}:{:7.5})",
                tree[im], len[k][0], tree[jm], len[k][1]
            );
            tree[im] = treetmp.clone();
            tree[jm] = String::new();
        }

        ac.unlink(jm);

        height[im] += len[k][0];
        dep[k].distfromtip = height[im];
    }

    if names.is_some() {
        std::fs::write("infile.tree", format!("{};\n#by loadtree\n", treetmp))
            .expect("write infile.tree");
    }

    Ok(Tree { topol, len, dep })
}

/// `setdensity()`: how many of the still-live pairs `pos` is closer than 1.0
/// to, weighted by how much closer.
fn setdensity(ac: &Chain, eff: &[Vec<f64>], pos: usize) -> f64 {
    let mut tmp = 0.0f64;
    let mut j = ac.next[pos];
    while j != -1 {
        let ju = j as usize;
        if eff[pos][ju - pos] < 1.0 {
            tmp += 2.0 - eff[pos][ju - pos];
        }
        j = ac.next[ju];
    }
    // The second pass walks from the head until it reaches `pos`.
    let mut j = 0i32;
    while j != -1 && j as usize != pos {
        let ju = j as usize;
        if eff[ju][pos - ju] < 1.0 {
            tmp += 2.0 - eff[ju][pos - ju];
        }
        j = ac.next[ju];
    }
    tmp
}

/// `getdensest()`: the member of `m` with the highest density.
fn getdensest(m: &[i32], d: &[f64]) -> i32 {
    let mut dmax = -100.0f64;
    let mut pmax = -1i32;
    for &i in m.iter() {
        if i < 0 {
            break;
        }
        if d[i as usize] > dmax {
            dmax = d[i as usize];
            pmax = i;
        }
    }
    pmax
}

/// C's `%f` for a double: six decimals.
fn fmt_f(v: f64) -> String {
    format!("{:.6}", v)
}

/// `topolorder_lessargs()`: in-order list of the leaves under `pos`.
fn topolorder_lessargs(order: &mut Vec<i32>, topol: &[[i32; 2]], dep: &[Treedep], pos: usize) {
    if dep[pos].child0 == -1 {
        order.push(topol[pos][0]);
    } else {
        topolorder_lessargs(order, topol, dep, dep[pos].child0 as usize);
    }
    if dep[pos].child1 == -1 {
        order.push(topol[pos][1]);
    } else {
        topolorder_lessargs(order, topol, dep, dep[pos].child1 as usize);
    }
}

/// `topolorderz()`: leaves under `pos`, restricted to one child when asked.
///
/// `nchild` is 0, 1 or 2 exactly as in the C.
pub fn topolorderz(topol: &[[i32; 2]], dep: &[Treedep], pos: usize, nchild: i32) -> Vec<i32> {
    let mut order = Vec::new();
    if nchild == 0 || nchild == 2 {
        if dep[pos].child0 == -1 {
            order.push(topol[pos][0]);
        } else {
            topolorder_lessargs(&mut order, topol, dep, dep[pos].child0 as usize);
        }
    }
    if nchild == 1 || nchild == 2 {
        if dep[pos].child1 == -1 {
            order.push(topol[pos][1]);
        } else {
            topolorder_lessargs(&mut order, topol, dep, dep[pos].child1 as usize);
        }
    }
    order
}

const GETA3: f64 = 0.001;

/// `counteff_simple_double_nostatic_memsave()`: per-sequence weights.
pub fn counteff_simple(nseq: usize, tree: &mut Tree) -> Vec<f64> {
    // The C clamps negative branch lengths over 0..nseq, i.e. one past the last
    // real step; len[nseq-1] is calloc'ed zero, so the extra visit is a no-op.
    for i in 0..nseq {
        if tree.len[i][0] < 0.0 {
            eprintln!(
                "WARNING: negative branch length {}, step {}-0",
                tree.len[i][0], i
            );
            tree.len[i][0] = 0.0;
        }
        if tree.len[i][1] < 0.0 {
            eprintln!(
                "WARNING: negative branch length {}, step {}-1",
                tree.len[i][1], i
            );
            tree.len[i][1] = 0.0;
        }
    }

    let mut rootnode = vec![0.0f64; nseq];
    let mut eff = vec![1.0f64; nseq];
    let mut memhist: Vec<Option<Vec<i32>>> = vec![None; nseq.saturating_sub(1)];

    for i in 0..nseq - 1 {
        let mem0: Vec<i32> = if tree.dep[i].child0 == -1 {
            vec![tree.topol[i][0]]
        } else {
            memhist[tree.dep[i].child0 as usize].clone().unwrap()
        };
        let mem1: Vec<i32> = if tree.dep[i].child1 == -1 {
            vec![tree.topol[i][1]]
        } else {
            memhist[tree.dep[i].child1 as usize].clone().unwrap()
        };
        let mut merged = mem0.clone();
        merged.extend_from_slice(&mem1);
        memhist[i] = Some(merged);

        for &s1 in &mem0 {
            let s = s1 as usize;
            rootnode[s] += tree.len[i][0] * eff[s];
            eff[s] *= 0.5;
        }
        for &s2 in &mem1 {
            let s = s2 as usize;
            rootnode[s] += tree.len[i][1] * eff[s];
            eff[s] *= 0.5;
        }
    }

    for v in rootnode.iter_mut() {
        *v += GETA3;
    }
    let mut total = 0.0f64;
    for v in rootnode.iter() {
        total += *v;
    }
    rootnode.iter().map(|v| v / total).collect()
}

/// `increaseintergroupdistanceshalfmtx()`: push every cross-group distance up
/// to at least the widest *within*-group distance, so the UPGMA has a reason
/// to keep each `--merge` group together before it is forced to.
fn increase_intergroup_distances(eff: &mut [Vec<f64>], groups: &[Vec<i32>], nseq: usize) {
    let mut nwarned = 0usize;
    eprintln!(); // "Hitsuyou desu."
    for (gi, g) in groups.iter().enumerate() {
        if g.len() < 2 {
            continue;
        }
        let mut tft = vec![false; nseq];
        for &s in g.iter() {
            tft[s as usize] = true;
        }
        let others: Vec<i32> = (0..nseq as i32).filter(|&m| !tft[m as usize]).collect();

        let mut maxdist = 0.0f64;
        for m in 1..g.len() {
            let s2 = g[m];
            for k in 0..m {
                let s1 = g[k];
                let (sl, ss) = if s2 > s1 { (s2, s1) } else { (s1, s2) };
                let dtmp = eff[ss as usize][(sl - ss) as usize];
                if dtmp > maxdist {
                    maxdist = dtmp;
                }
            }
        }

        for &s2 in g.iter() {
            for &s1 in others.iter() {
                let (sl, ss) = if s2 > s1 { (s2, s1) } else { (s1, s2) };
                let d = &mut eff[ss as usize][(sl - ss) as usize];
                if *d < maxdist {
                    if *d < 0.5 && nwarned < 100 {
                        nwarned += 1;
                        eprintln!("# Sequences {} and {} seem to be closely related, but are not in the same sub MSA ({}) in your setting.", s2 + 1, s1 + 1, gi + 1);
                    } else if *d < 0.5 {
                        nwarned += 1;
                    }
                    *d = maxdist;
                }
            }
        }
    }
    if nwarned > 100 {
        eprintln!("# Sequenc.... (more pairs)");
    }
}

/// `overlapmember()`: the two lists share at least one member.
fn overlapmember(mem1: &[i32], mem2: &[i32]) -> bool {
    mem1.iter().any(|a| mem2.contains(a))
}

/// `includemember()`: every member of `mem` is in `cand`.
fn includemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.iter().all(|m| cand.contains(m))
}

/// `fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained()`: the UPGMA
/// again, but a merge that would split a `--merge` group is rejected and the
/// next-best pair tried instead.
///
/// Returns the tree in the same memsave shape as [`upgma_treeout`] plus, per
/// node and side, the full member list -- which this builder has to keep
/// anyway (it tests candidate merges against the groups) and which
/// `disttbfast`'s sub-alignment checks then read.
#[allow(clippy::type_complexity)]
pub fn upgma_treeout_constrained(
    nseq: usize,
    eff: &mut Vec<Vec<f64>>,
    sueff_global: f64,
    names: &[Vec<u8>],
    outnumber: bool,
    groups: &[Vec<i32>],
) -> (Tree, Vec<[Vec<i32>; 2]>) {
    increase_intergroup_distances(eff, groups, nseq);

    let sueff1 = 1.0 - sueff_global;
    let sueff05 = sueff_global * 0.5;
    let clusterfunc = |d1: f64, d2: f64| -> f64 {
        let m = if d1 < d2 { d1 } else { d2 };
        m * sueff1 + (d1 + d2) * sueff05
    };

    let mut topol = vec![[-1i32, -1i32]; nseq];
    let mut full: Vec<[Vec<i32>; 2]> = (0..nseq).map(|_| [Vec::new(), Vec::new()]).collect();
    let mut len = vec![[0.0f64, 0.0f64]; nseq];
    let mut dep = vec![Treedep::default(); nseq];

    let mut ac = Chain::new(nseq);
    let mut mindisfrom = vec![0.0f64; nseq];
    let mut nearest = vec![0i32; nseq];
    for i in 0..nseq {
        let (m, n) = setnearest(&ac, eff, i);
        mindisfrom[i] = m;
        nearest[i] = n;
    }

    let mut tmptmplen = vec![0.0f64; nseq];
    let mut hist = vec![-1i32; nseq];
    let mut nmemar = vec![1i32; nseq];
    let mut im_g: i32 = -1;
    let mut jm_g: i32 = -1;

    let mut tree: Vec<String> = (0..nseq)
        .map(|i| format!("\n{}_{}\n", i + 1, newick_name(&names[i], outnumber)))
        .collect();
    let mut treetmp = String::new();

    let mut inconsistent = vec![vec![false; nseq]; nseq];
    let mut inconsistentpairs: Vec<(usize, usize)> = Vec::new();
    let mut warned = vec![false; groups.len()];

    for k in 0..nseq - 1 {
        for &(a, b) in inconsistentpairs.iter() {
            inconsistent[a][b] = false;
        }
        inconsistentpairs.clear();
        let mut firsttime = true;
        // `minscore` is function-scope in the C: whichever branch of the retry
        // loop ran last sets the branch length below.
        let mut minscore;

        loop {
            if firsttime {
                firsttime = false;
                minscore = 999.9f64;
                let mut i = 0i32;
                while ac.next[i as usize] != -1 {
                    if mindisfrom[i as usize] < minscore {
                        im_g = i;
                        minscore = mindisfrom[i as usize];
                    }
                    i = ac.next[i as usize];
                }
                jm_g = nearest[im_g as usize];
                if jm_g < im_g {
                    std::mem::swap(&mut im_g, &mut jm_g);
                }
            } else {
                // No cached nearest here: the whole live half-matrix is
                // rescanned, skipping the pairs already found inconsistent.
                minscore = 999.9f64;
                let mut i = 0i32;
                while ac.next[i as usize] != -1 {
                    let iu = i as usize;
                    let mut j = ac.next[iu];
                    while j != -1 {
                        let ju = j as usize;
                        let tmpdouble = eff[iu][ju - iu];
                        if !inconsistent[iu][ju] && tmpdouble < minscore {
                            minscore = tmpdouble;
                            im_g = i;
                            jm_g = j;
                        }
                        j = ac.next[ju];
                    }
                    let mut j = 0i32;
                    while j != -1 && j != i {
                        let ju = j as usize;
                        let tmpdouble = eff[ju][iu - ju];
                        if !inconsistent[ju][iu] && tmpdouble < minscore {
                            minscore = tmpdouble;
                            im_g = j;
                            jm_g = i;
                        }
                        j = ac.next[ju];
                    }
                    i = ac.next[iu];
                }
            }

            // Every remaining pair inconsistent means the grouping cannot be
            // satisfied at all.
            let mut allinconsistent = true;
            'scan: {
                let mut i = 0i32;
                while ac.next[i as usize] != -1 {
                    let mut j = ac.next[i as usize];
                    while j != -1 {
                        if !inconsistent[i as usize][j as usize] {
                            allinconsistent = false;
                            break 'scan;
                        }
                        j = ac.next[j as usize];
                    }
                    i = ac.next[i as usize];
                }
            }
            if allinconsistent {
                eprintln!("\n\n\nPlease check whether the grouping is possible.\n\n");
                std::process::exit(1);
            }

            let im = im_g as usize;
            let jm = jm_g as usize;
            let mut testtopol: Vec<i32> = Vec::new();
            for &side in [im, jm].iter() {
                let prevnode = hist[side];
                if prevnode == -1 {
                    testtopol.push(side as i32);
                } else {
                    let p = prevnode as usize;
                    testtopol.extend_from_slice(&full[p][0]);
                    testtopol.extend_from_slice(&full[p][1]);
                }
            }

            let mut hit = None;
            for (gi, g) in groups.iter().enumerate() {
                if overlapmember(g, &testtopol)
                    && !includemember(&testtopol, g)
                    && !includemember(g, &testtopol)
                {
                    if !warned[gi] {
                        warned[gi] = true;
                        eprintln!(
                            "\n###################################################################"
                        );
                        eprintln!(
                            "# WARNING: Group {} is forced to be a monophyletic cluster.",
                            gi + 1
                        );
                        eprintln!(
                            "###################################################################"
                        );
                    }
                    inconsistent[im][jm] = true;
                    inconsistentpairs.push((im, jm));
                    hit = Some(gi);
                    break;
                }
            }
            if hit.is_none() {
                break;
            }
        }

        let im = im_g as usize;
        let jm = jm_g as usize;

        // `topol[k][side]` keeps the full member list here: the smaller-headed
        // child's members first, then the other's.
        let build = |side: usize| -> Vec<i32> {
            let prevnode = hist[side];
            if prevnode == -1 {
                vec![side as i32]
            } else {
                let p = prevnode as usize;
                let (a, b) = if full[p][0][0] > full[p][1][0] {
                    (1usize, 0usize)
                } else {
                    (0usize, 1usize)
                };
                let mut v = full[p][a].clone();
                v.extend_from_slice(&full[p][b]);
                v
            }
        };
        dep[k].child0 = hist[im];
        let nmemim = nmemar[im];
        let f0 = build(im);
        dep[k].child1 = hist[jm];
        let nmemjm = nmemar[jm];
        let f1 = build(jm);
        topol[k][0] = f0[0];
        topol[k][1] = f1[0];
        full[k] = [f0, f1];

        let minscore = minscore * 0.5;

        len[k][0] = minscore - tmptmplen[im];
        len[k][1] = minscore - tmptmplen[jm];
        dep[k].distfromtip = minscore;
        tmptmplen[im] = minscore;

        hist[im] = k as i32;
        nmemar[im] = nmemim + nmemjm;

        mindisfrom[im] = 999.9;
        let mut i = 0i32;
        while i != -1 {
            let iu = i as usize;
            if iu != im && iu != jm {
                let (miniim, maxiim, minijm, maxijm) = if iu < im {
                    (iu, im, iu, jm)
                } else if iu < jm {
                    (im, iu, iu, jm)
                } else {
                    (im, iu, jm, iu)
                };
                let eff0 = eff[miniim][maxiim - miniim];
                let eff1 = eff[minijm][maxijm - minijm];
                let tmpdouble = clusterfunc(eff0, eff1);
                eff[miniim][maxiim - miniim] = tmpdouble;
                if tmpdouble < mindisfrom[iu] {
                    mindisfrom[iu] = tmpdouble;
                    nearest[iu] = im as i32;
                }
                if tmpdouble < mindisfrom[im] {
                    mindisfrom[im] = tmpdouble;
                    nearest[im] = iu as i32;
                }
                if nearest[iu] == jm as i32 {
                    nearest[iu] = im as i32;
                }
            }
            i = ac.next[iu];
        }

        treetmp = format!(
            "({}:{:7.5},{}:{:7.5})",
            tree[im], len[k][0], tree[jm], len[k][1]
        );
        tree[im] = treetmp.clone();
        tree[jm] = String::new();

        ac.unlink(jm);

        let mut i = 0i32;
        while i != -1 {
            let iu = i as usize;
            if nearest[iu] == im as i32 {
                let (miniim, maxiim) = if iu < im { (iu, im) } else { (im, iu) };
                if eff[miniim][maxiim - miniim] > mindisfrom[iu] {
                    let (m, n) = setnearest(&ac, eff, iu);
                    mindisfrom[iu] = m;
                    nearest[iu] = n;
                }
            }
            i = ac.next[iu];
        }
    }

    std::fs::write("infile.tree", format!("{};\n", treetmp)).expect("Cannot open infile.tree");

    (Tree { topol, len, dep }, full)
}

#[cfg(test)]
mod constrained_tests {
    use super::*;

    /// A half matrix in the C's layout: `eff[i][j-i]` for j >= i.
    fn halfmtx(n: usize, f: impl Fn(usize, usize) -> f64) -> Vec<Vec<f64>> {
        (0..n)
            .map(|i| (0..n - i).map(|d| f(i, i + d)).collect())
            .collect()
    }

    fn names(n: usize) -> Vec<Vec<u8>> {
        (0..n).map(|i| format!("=s{i}").into_bytes()).collect()
    }

    /// The builder writes `infile.tree`, so each case runs in its own
    /// directory. The working directory is process-wide, hence the lock.
    fn in_dir<T>(tag: &str, f: impl FnOnce() -> T) -> T {
        static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
        let _g = LOCK.lock().unwrap_or_else(|e| e.into_inner());
        let d = std::env::temp_dir().join(format!("mafft-rs-tree-{}-{}", tag, std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let old = std::env::current_dir().unwrap();
        std::env::set_current_dir(&d).unwrap();
        let r = f();
        std::env::set_current_dir(old).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        r
    }

    /// Every group ends up as one node's member list -- that is what
    /// "monophyletic" means here, and it is the whole point of this builder.
    #[test]
    fn every_group_comes_out_monophyletic() {
        let n = 8;
        // Neighbours by index are closest, so the plain UPGMA would certainly
        // split these interleaved groups.
        let mut eff = halfmtx(n, |i, j| (j - i) as f64 * 0.1 + 0.05);
        let groups = vec![vec![0, 2, 4, 6], vec![1, 3, 5, 7]];
        let (tree, full) = in_dir("mono", || {
            upgma_treeout_constrained(n, &mut eff, 0.1, &names(n), false, &groups)
        });

        for g in groups.iter() {
            let mut want = g.clone();
            want.sort_unstable();
            let found = (0..n - 1).any(|k| {
                for side in 0..2 {
                    let mut a = full[k][side].clone();
                    a.sort_unstable();
                    if a == want {
                        return true;
                    }
                }
                false
            });
            assert!(found, "group {g:?} is not a node");
        }
        // The memsave view agrees with the full lists.
        for k in 0..n - 1 {
            assert_eq!(tree.topol[k][0], full[k][0][0]);
            assert_eq!(tree.topol[k][1], full[k][1][0]);
        }
    }

    /// Each node's two sides partition its members, and the root covers
    /// everything.
    #[test]
    fn the_constrained_tree_is_still_a_tree() {
        let n = 7;
        let mut eff = halfmtx(n, |i, j| ((i * 7 + j * 13) % 11) as f64 * 0.1 + 0.1);
        let groups = vec![vec![0, 3, 6]];
        let (_tree, full) = in_dir("tree", || {
            upgma_treeout_constrained(n, &mut eff, 0.1, &names(n), false, &groups)
        });

        for k in 0..n - 1 {
            let mut all = full[k][0].clone();
            all.extend_from_slice(&full[k][1]);
            let before = all.len();
            all.sort_unstable();
            all.dedup();
            assert_eq!(all.len(), before, "node {k} repeats a member");
        }
        let mut root = full[n - 2][0].clone();
        root.extend_from_slice(&full[n - 2][1]);
        root.sort_unstable();
        assert_eq!(root, (0..n as i32).collect::<Vec<_>>());
    }

    #[test]
    fn overlap_and_include_are_what_they_say() {
        assert!(overlapmember(&[1, 2], &[2, 3]));
        assert!(!overlapmember(&[1, 2], &[3, 4]));
        assert!(includemember(&[1, 2], &[1, 2, 3]));
        assert!(!includemember(&[1, 4], &[1, 2, 3]));
        // The empty list is included in anything.
        assert!(includemember(&[], &[1]));
    }

    /// Cross-group distances are raised to the widest within-group distance;
    /// within-group ones are left alone.
    #[test]
    fn intergroup_distances_are_raised_to_the_group_diameter() {
        let n = 4;
        let mut eff = halfmtx(n, |i, j| if i == j { 0.0 } else { 0.1 * (j - i) as f64 });
        // group {0,1} has diameter 0.1; 1-2 is 0.1 already, 0-2 is 0.2.
        increase_intergroup_distances(&mut eff, &[vec![0, 1]], n);
        assert!((eff[0][1] - 0.1).abs() < 1e-12, "within group untouched");
        assert!((eff[0][2] - 0.2).abs() < 1e-12, "already above");
        assert!((eff[1][1] - 0.1).abs() < 1e-12, "1-2 stays at the diameter");
    }
}

/// `createchain()`: `--treein pileup` and `--treein shuffle`.
///
/// Not a tree built from distances at all -- it is a comb. Sequence `order[0]`
/// and `order[1]` are joined, then `order[2]` is added to that, then
/// `order[3]`, and so on; `pileup` keeps the input order and `shuffle`
/// permutes it first with glibc's `rand()`.
///
/// Branch lengths are a fixed `2/nseq` per step, so the depth grows linearly.
/// Writes `_guidetree`, and `infile.tree` when `treeout` is on.
pub fn createchain(
    nseq: usize,
    names: Option<&[Vec<u8>]>,
    outnumber: bool,
    shuffle: bool,
    seed: u32,
) -> Tree {
    let mut topol = vec![[-1i32, -1i32]; nseq];
    let mut len = vec![[0.0f64, 0.0f64]; nseq];
    let mut dep = vec![Treedep::default(); nseq];

    // One Newick fragment per sequence, spliced into a comb below.
    let tree: Vec<String> = match names {
        Some(name) => (0..nseq)
            .map(|i| format!("\n{}_{}\n", i + 1, newick_name(&name[i], outnumber)))
            .collect(),
        None => Vec::new(),
    };
    // The C opens all nseq-1 parentheses up front, then fills in from the left.
    let mut instanttree = String::new();
    if names.is_some() {
        for _ in 0..nseq - 1 {
            instanttree.push('(');
        }
    }

    let mut order: Vec<i32> = (0..nseq as i32).collect();
    let mut rng = crate::cerand::CRand::new(seed);
    if shuffle {
        crate::cerand::stringshuffle(&mut order, &mut rng);
    }

    let l = 2.0 / nseq as f64;
    let mut ll = l;

    let im = order[0];
    let jm = order[1];
    // The smaller index goes on side 0, which is what `topolorderz` and the
    // `order` file below assume.
    let mut mm = if im < jm {
        topol[0][0] = im;
        topol[0][1] = jm;
        im
    } else {
        topol[0][0] = jm;
        topol[0][1] = im;
        jm
    };
    len[0][0] = l;
    len[0][1] = l;
    dep[0].child0 = -1;
    dep[0].child1 = -1;
    dep[0].distfromtip = l;
    ll += l;

    if names.is_some() {
        instanttree.push_str(&format!("{}:{:7.5},", tree[im as usize], len[0][0]));
    }

    for i in 1..nseq - 1 {
        let im = order[i];
        let jm = order[i + 1];
        if mm < jm {
            topol[i][0] = mm;
            topol[i][1] = jm;
            len[i][0] = l;
            len[i][1] = ll;
            dep[i].child0 = i as i32 - 1;
            dep[i].child1 = -1;
            dep[i].distfromtip = ll;
        } else {
            topol[i][1] = mm;
            topol[i][0] = jm;
            mm = jm;
            len[i][1] = l;
            len[i][0] = ll;
            dep[i].child1 = i as i32 - 1;
            dep[i].child0 = -1;
            dep[i].distfromtip = ll;
        }
        if names.is_some() {
            instanttree.push_str(&format!("{}:{:7.5}):{:7.5},", tree[im as usize], ll - l, l));
        }
        ll += l;
    }

    if names.is_some() {
        let last = order[nseq - 1];
        instanttree.push_str(&format!("{}:{:7.5})", tree[last as usize], ll - l));
        std::fs::write("infile.tree", format!("{};\n", instanttree))
            .expect("Cannot open infile.tree");
    }

    let mut gt = String::new();
    for i in 0..nseq - 1 {
        gt.push_str(&format!(
            "{} {} {:.6} {:.6}\n",
            topol[i][0] + 1,
            topol[i][1] + 1,
            len[i][0],
            len[i][1]
        ));
    }
    std::fs::write("_guidetree", gt).expect("cannot open _guidetree");

    Tree { topol, len, dep }
}
