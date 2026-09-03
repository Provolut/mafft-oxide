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
    let mut im_g: i32 = -1;

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

    Tree { topol, len, dep }
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
