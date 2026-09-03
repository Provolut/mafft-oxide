//! The progressive alignment driver, transliterated from `treebase()` in
//! `core/disttbfast.c` and `core/tbfast.c` (single-threaded, non-memsave).
//!
//! The two C copies differ in three places, all of them parameters here: the
//! minimum sequence weight, whether the FFT branch may be taken, and whether
//! position-specific constraints from the pairwise stage are applied.

use crate::align::{AlignState, CpmxHist};
use crate::constants::{makedynamicmtx, Constants};
use crate::falign::{falign, falign_givenanchors, FalignState};
use crate::localhom::{fastshrinklocalhom_half, fillimp, LocalHomTable};
use crate::tree::Tree;

/// What varies between the `disttbfast` and `tbfast` call sites.
pub struct TreebaseOpts<'a> {
    /// Terminal gaps are free when this is false (`-O` on the command line).
    pub outgap: bool,
    pub use_fft: bool,
    pub force_fft: bool,
    /// `-N`: suppress the C's automatic switch to the memsave aligner past
    /// 30000 columns.
    pub nevermemsave: bool,
    /// The C's `alg` global: `'A'` for `A__align`, `'M'` for the linear-space
    /// `MSalignmm`. It is a *global* there, and the 30000-column switch below
    /// writes to it, so once any node crosses the threshold every later node
    /// uses the memsave aligner too. Kept as loop state for that reason.
    pub alg: u8,
    /// `-W`: sequence weights are clamped up to this.
    pub minimumweight: f64,
    /// The constraints from the pairwise stage (L-INS-i / G-INS-i).
    pub constraint: Option<&'a LocalHomTable>,
    /// `-l`: how strongly the constraints pull.
    pub fastathreshold: f64,
    /// `--add`: what to do at each tree node -- `'a'` align, `'n'` skip,
    /// `'1'`/`'2'` merge into the group on that side, `'w'` align but both
    /// sides already contain added sequences.
    pub mergeoralign: Option<&'a [u8]>,
    /// The number of sequences appended by `--add`; the first
    /// `njob - nadd` are the existing alignment.
    pub nadd: usize,
    /// `-l <file>`: the external anchors, and the per-pair index into them.
    /// When present every node is aligned with `Falign_givenanchors` instead
    /// of the FFT scan, whatever `use_fft` says.
    pub extanch: Option<(&'a [crate::anchors::ExtAnch], &'a [Vec<i32>])>,
    /// `-x`: how much longer a terminal overhang may be than its partner.
    pub terminalmargin: i32,
    /// Rebuild the score matrix per tree node from its distance from the tips
    /// (`makedynamicmtx`). A no-op unless `--allowshift` is on, but `tbfast`
    /// calls it unconditionally and `disttbfast` only when it matters.
    pub use_dynamic_mtx: bool,
}

/// `eq2dashmatometehayaku()`: the `'='` markers become ordinary gaps once the
/// step that inserted them is finished with.
fn eq2dash(s: &mut [u8]) {
    for c in s.iter_mut() {
        if *c == b'=' {
            *c = b'-';
        }
    }
}

/// `fastconjuction_noname()`: gather a group's members and normalise their
/// weights, returning the pre-normalisation total (`oritotal`).
fn fastconjuction_noname(memlist: &[i32], eff: &[f64], mineff: f64) -> (Vec<f64>, f64) {
    let mut peff: Vec<f64> = memlist
        .iter()
        .map(|&m| {
            let e = eff[m as usize];
            if e < mineff {
                mineff
            } else {
                e
            }
        })
        .collect();
    let mut total = 0.0f64;
    for v in peff.iter() {
        total += *v;
    }
    for v in peff.iter_mut() {
        *v /= total;
    }
    (peff, total)
}

/// `treebase()`: walk the guide tree bottom-up, aligning group to group.
///
/// `aseq` holds the (initially ungapped) sequences and is aligned in place.
#[allow(clippy::too_many_arguments)]
pub fn treebase(
    cn: &Constants,
    opts: &TreebaseOpts,
    njob: usize,
    aseq: &mut [Vec<u8>],
    nlen: &mut [i32],
    tree: &Tree,
    effarr: &[f64],
    st: &mut AlignState,
    fs: &mut FalignState,
) {
    let outgap = opts.outgap;

    let mut alreadyaligned: Vec<bool> = (0..njob).map(|i| i < njob - opts.nadd).collect();
    let mut gapmap: Vec<i32> = Vec::new();

    let mut memhist: Vec<Option<Vec<i32>>> = vec![None; njob.saturating_sub(1)];
    let mut cpmxhist: Vec<Option<CpmxHist>> = (0..njob.saturating_sub(1)).map(|_| None).collect();
    // The C's `alg` global. `treebase` may flip it to 'M' at the 30000-column
    // threshold, and it stays there for every later node.
    let mut alg = opts.alg;
    let mut fftlog = vec![1i32; njob];
    let mut impmtx: Vec<Vec<f64>> = Vec::new();
    let mut impalloclen = 0usize;
    let mut dynmtx: Vec<Vec<f64>> = vec![vec![0.0; cn.nalphabets]; cn.nalphabets];

    for l in 0..njob - 1 {
        let moa = opts.mergeoralign.map(|m| m[l]).unwrap_or(b'a');
        let m1 = tree.topol[l][0] as usize;
        let m2 = tree.topol[l][1] as usize;

        let child0 = tree.dep[l].child0;
        let child1 = tree.dep[l].child1;
        let cpmxchild0 = if child0 == -1 {
            None
        } else {
            cpmxhist[child0 as usize].take()
        };
        let cpmxchild1 = if child1 == -1 {
            None
        } else {
            cpmxhist[child1 as usize].take()
        };

        let mem0: Vec<i32> = if child0 == -1 {
            vec![m1 as i32]
        } else {
            memhist[child0 as usize].clone().unwrap()
        };
        let mem1: Vec<i32> = if child1 == -1 {
            vec![m2 as i32]
        } else {
            memhist[child1 as usize].clone().unwrap()
        };
        if l != njob - 2 {
            let mut merged = mem0.clone();
            merged.extend_from_slice(&mem1);
            memhist[l] = Some(merged);
        }

        // `'n'`: this node is entirely inside a group that is already aligned,
        // so there is nothing to do.
        if moa == b'n' {
            continue;
        }

        if opts.use_dynamic_mtx {
            makedynamicmtx(
                &mut dynmtx,
                &cn.n_dis_consweight_multi,
                tree.dep[l].distfromtip,
                &cn.amino,
                cn.specificityconsideration,
            );
        } else {
            for (i, row) in dynmtx.iter_mut().enumerate() {
                row.copy_from_slice(&cn.n_dis_consweight_multi[i]);
            }
        }

        // `len1 = strlen( aseq[m1] ); len2 = strlen( aseq[m2] );` -- the full
        // alignment rows of the two representatives, not the groups.
        let len1 = aseq[m1].len();
        let len2 = aseq[m2].len();

        let (effarr1, orieff1) = fastconjuction_noname(&mem0, effarr, opts.minimumweight);
        let (effarr2, orieff2) = fastconjuction_noname(&mem1, effarr, opts.minimumweight);
        let clus1 = mem0.len();
        let clus2 = mem1.len();

        // Gather the group members (mseq1 / mseq2 alias rows of aseq in C).
        let mut group1: Vec<Vec<u8>> = mem0.iter().map(|&m| aseq[m as usize].clone()).collect();
        let mut group2: Vec<Vec<u8>> = mem1.iter().map(|&m| aseq[m as usize].clone()).collect();

        // A merge marks the columns it inserts with `'='` so they can be told
        // from the alignment's own gaps and propagated afterwards.
        let newgapstr = if moa == b'1' || moa == b'2' { b'=' } else { b'-' };

        // Collapse the all-gap columns of the side being merged *into*, so the
        // DP does not have to reproduce them, and remember where they were.
        if moa == b'1' {
            gapmap = crate::addfunctions::findcommongaps(&group2);
            crate::addfunctions::commongappick(&mut group2);
        } else if moa == b'2' {
            gapmap = crate::addfunctions::findcommongaps(&group1);
            crate::addfunctions::commongappick(&mut group1);
        }

        // `if( !nevermemsave && ( alg != 'M' && ( len1 > 30000 || len2 > 30000 ) ) )`
        // -- disttbfast.c:2836. Past 30000 columns the C stops using
        // `A__align` and switches, for the rest of the run, to the linear
        // -space aligner. `alg` is a global there and is never switched back,
        // which is why it is loop state here.
        //
        // This is not a memory optimisation with the same answer: forcing the
        // C down the memsave path with `-M -B` on inputs it would otherwise
        // align conventionally changes the alignment on 17 of the 86 corpus
        // inputs, because `MSalignmm.c` compiles with `USE_PENALTY_EX 0` where
        // `Salignmm.c` uses 1.
        if !opts.nevermemsave && alg != b'M' && (len1 > 30000 || len2 > 30000) {
            eprintln!(
                "\nlen1={}, len2={}, Switching to the memsave mode",
                len1, len2
            );
            alg = b'M';
        }

        if let Some(table) = opts.constraint {
            // imp_match_init_strict() -> fillimp(): paint the constraint
            // weights onto the match-bonus matrix for this pair of groups.
            // With constraints the C never takes the FFT branch.
            let lgth1 = group1[0].len();
            let lgth2 = group2[0].len();
            if impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2 {
                impalloclen = lgth1.max(lgth2) + 2;
                impmtx = vec![vec![0.0f64; impalloclen]; impalloclen];
            }
            let shrink = fastshrinklocalhom_half(&mem0, &mem1, table);
            {
                let g1: Vec<&[u8]> = group1.iter().map(|s| s.as_slice()).collect();
                let g2: Vec<&[u8]> = group2.iter().map(|s| s.as_slice()).collect();
                fillimp(
                    &mut impmtx,
                    clus1,
                    clus2,
                    lgth1,
                    lgth2,
                    &g1,
                    &g2,
                    &effarr1,
                    &effarr2,
                    &shrink,
                    &mem0,
                    &mem1,
                    opts.fastathreshold,
                );
            }
            let (_pscore, hist) = st.a_align(
                cn,
                &dynmtx,
                cn.penalty,
                cn.penalty_ex,
                &mut group1,
                &mut group2,
                &effarr1,
                &effarr2,
                outgap,
                outgap,
                None,
                Some(&impmtx),
                mem0[0],
                1,
                cpmxchild0,
                cpmxchild1,
                orieff1,
                orieff2,
                true,
                newgapstr,
            );
            cpmxhist[l] = hist;
        } else {
            // `ffttry`: the FFT scan needs sequences longer than the group is wide.
            let ffttry = nlen[m1] > clus1 as i32 && nlen[m2] > clus2 as i32;
            if let Some((extanch, anchindex)) = opts.extanch {
                let pairanch = crate::anchors::pick_pair_anch(
                    extanch, anchindex, &mem0, &mem1, &group1, &group2,
                );
                falign_givenanchors(
                    &pairanch,
                    st,
                    cn,
                    &dynmtx,
                    &mut group1,
                    &mut group2,
                    &effarr1,
                    &effarr2,
                    outgap,
                    opts.terminalmargin,
                    &mut fftlog[m1],
                    newgapstr,
                );
            } else if opts.force_fft || (opts.use_fft && ffttry) {
                // Under memsave the C anchors with `Falign_udpari_long`
                // rather than `Falign`; `falign` takes that as a mode.
                falign(
                    fs,
                    st,
                    cn,
                    &dynmtx,
                    &mut group1,
                    &mut group2,
                    &effarr1,
                    &effarr2,
                    outgap,
                    &mut fftlog[m1],
                    newgapstr,
                    alg,
                );
            } else if alg == b'M' {
                // `case( 'M' )` has no `clus1 == 1 && clus2 == 1` shortcut --
                // unlike `'A'` and `'d'`, a two-sequence node still goes
                // through the linear-space aligner.
                fftlog[m1] = 0;
                let (_pscore, hist) = crate::msalign::msalignmm(
                    cn,
                    &dynmtx,
                    &mut group1,
                    &mut group2,
                    &effarr1,
                    &effarr2,
                    outgap,
                    outgap,
                    None,
                    cpmxchild0,
                    cpmxchild1,
                    orieff1,
                    orieff2,
                    true,
                    newgapstr,
                );
                cpmxhist[l] = hist;
            } else if alg == b'a' {
                // `case( 'a' )`: like `'M'`, no 1x1 shortcut.
                fftlog[m1] = 0;
                crate::saalign::aalign(
                    cn, &mut group1, &mut group2, &effarr1, &effarr2, outgap, newgapstr,
                );
            } else if clus1 == 1 && clus2 == 1 {
                fftlog[m1] = 0;
                let mut s1 = std::mem::take(&mut group1[0]);
                let mut s2 = std::mem::take(&mut group2[0]);
                st.g_align11(
                    cn,
                    &dynmtx,
                    cn.penalty,
                    cn.penalty_ex,
                    &mut s1,
                    &mut s2,
                    outgap,
                    outgap,
                    newgapstr,
                );
                group1[0] = s1;
                group2[0] = s2;
            } else {
                fftlog[m1] = 0;
                let (_pscore, hist) = st.a_align(
                    cn,
                    &dynmtx,
                    cn.penalty,
                    cn.penalty_ex,
                    &mut group1,
                    &mut group2,
                    &effarr1,
                    &effarr2,
                    outgap,
                    outgap,
                    None,
                    None,
                    mem0[0],
                    1,
                    cpmxchild0,
                    cpmxchild1,
                    orieff1,
                    orieff2,
                    true,
                    newgapstr,
                );
                cpmxhist[l] = hist;
            }
        }

        for (k, &m) in mem0.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group1[k]);
        }
        for (k, &m) in mem1.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group2[k]);
        }

        if moa == b'1' {
            // "jissainiha nai": the added sequences are always last, so the
            // merged side is never group 1.
            eprintln!("Check source!!!");
            std::process::exit(1);
        }
        if moa == b'2' {
            let newlen = aseq[mem0[0] as usize].len();
            crate::addfunctions::adjustgapmap(newlen, &mut gapmap, &aseq[mem0[0] as usize].clone());
            // `restorecommongaps()`'s first act is `if( n0 == 0 ) return;`
            // with `n0 = njob - (clus1+clus2)`: at a node covering every
            // sequence there is nothing outside the two groups to keep in
            // step, so neither the sequences nor `gapmap` are touched.
            if njob - (clus1 + clus2) != 0 {
                let mut allmem = mem0.clone();
                allmem.extend_from_slice(&mem1);
                crate::addfunctions::restorecommongaps(aseq, &allmem, &mut gapmap, b'-');
            }
            let g1: Vec<Vec<u8>> = mem0.iter().map(|&m| aseq[m as usize].clone()).collect();
            let gaplen = crate::addfunctions::findnewgaps(&g1, 0);
            crate::addfunctions::insertnewgaps(
                cn,
                st,
                njob,
                &alreadyaligned,
                aseq,
                &mem0,
                &mem1,
                &gaplen,
                &gapmap,
                b'-',
            );
            for m in mem0.iter().chain(mem1.iter()) {
                eq2dash(&mut aseq[*m as usize]);
            }
            for &m in mem1.iter() {
                alreadyaligned[m as usize] = true;
            }
        }

        // `nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] )` -- assigned to an int.
        nlen[m1] = (0.5 * (nlen[m1] + nlen[m2]) as f64) as i32;
    }
}

/// `dooneiteration()`: `-r` (`--oneiteration`), run after `treebase`.
///
/// Not the tree again -- this is leave-one-out. For each sequence in turn
/// (twice over, `ITERATIVECYCLE` is 2), it is pulled out, realigned against
/// all the others, and the result kept only if it scores at least as well as
/// what was there. So the alignment can only improve, by `intergroup_score`'s
/// reckoning.
#[allow(clippy::too_many_arguments)]
pub fn dooneiteration(
    cn: &Constants,
    opts: &TreebaseOpts,
    njob: usize,
    aseq: &mut [Vec<u8>],
    nlen: &mut [i32],
    tree: &Tree,
    effarr: &[f64],
    st: &mut AlignState,
    fs: &mut FalignState,
) {
    const ITERATIVECYCLE: usize = 2;
    let outgap = opts.outgap;
    let mut fftlog = vec![1i32; njob];
    let mut dynmtx: Vec<Vec<f64>> = vec![vec![0.0; cn.nalphabets]; cn.nalphabets];

    for ll in 0..njob * ITERATIVECYCLE {
        let l = ll % njob;
        let mem0: Vec<i32> = vec![l as i32];
        let mem1: Vec<i32> = (0..njob).filter(|&i| i != l).map(|i| i as i32).collect();
        let m1 = l;
        let m2 = mem1[0] as usize;

        if opts.use_dynamic_mtx {
            makedynamicmtx(
                &mut dynmtx,
                &cn.n_dis_consweight_multi,
                tree.dep[l].distfromtip,
                &cn.amino,
                cn.specificityconsideration,
            );
        } else {
            for (i, row) in dynmtx.iter_mut().enumerate() {
                row.copy_from_slice(&cn.n_dis_consweight_multi[i]);
            }
        }

        // Weights are *not* clamped here: the C passes 0.0 for `minimumweight`.
        let (effarr1, _orieff1) = fastconjuction_noname(&mem0, effarr, 0.0);
        let (effarr2, _orieff2) = fastconjuction_noname(&mem1, effarr, 0.0);
        let clus1 = mem0.len();
        let clus2 = mem1.len();

        let mut group1: Vec<Vec<u8>> = mem0.iter().map(|&m| aseq[m as usize].clone()).collect();
        let mut group2: Vec<Vec<u8>> = mem1.iter().map(|&m| aseq[m as usize].clone()).collect();

        let oscore = crate::score::intergroup_score(
            cn,
            &group1,
            &group2,
            &effarr1,
            &effarr2,
            clus1,
            clus2,
            group1[0].len(),
        );
        let oseq1 = group1.clone();
        let oseq2 = group2.clone();

        crate::addfunctions::commongappick(&mut group2);
        crate::addfunctions::commongappick(&mut group1);

        // `ffttry` has no `fftlog` guard here, unlike `treebase`'s.
        let ffttry = nlen[m1] > clus1 as i32 && nlen[m2] > clus2 as i32;
        if opts.force_fft || (opts.use_fft && ffttry) {
            falign(
                fs,
                st,
                cn,
                &dynmtx,
                &mut group1,
                &mut group2,
                &effarr1,
                &effarr2,
                outgap,
                &mut fftlog[m1],
                b'-',
                // `dooneiteration` uses `Falign_udpari_long` under memsave
                // too (disttbfast.c:2399), but it never reaches the memsave
                // switch: the refinement pass refuses at `l != 0` first, and
                // `alg` cannot have been flipped by `treebase` because the
                // two do not share it here.
                b'A',
            );
        } else {
            fftlog[m1] = 0;
            if clus1 == 1 && clus2 == 1 {
                let mut s1 = std::mem::take(&mut group1[0]);
                let mut s2 = std::mem::take(&mut group2[0]);
                st.g_align11(
                    cn, &dynmtx, cn.penalty, cn.penalty_ex, &mut s1, &mut s2, outgap, outgap, b'-',
                );
                group1[0] = s1;
                group2[0] = s2;
            } else {
                // The C segfaults here for l != 0, and this is the only
                // branch that reaches it: `disttbfast()` frees `cpmxhist` and
                // sets it to NULL at the end of the guide-tree loop
                // (disttbfast.c:4589), then hands that NULL to
                // `dooneiteration`, which passes `cpmxhist + l` to `A__align`
                // as `cpmxresult`. `A__align` guards on `if( cpmxresult )`, so
                // l == 0 is skipped -- but for l != 0 the pointer is a small
                // non-null offset from NULL and writing through it dies
                // (Salignmm.c:2142). The FFT branch and the 1-vs-1 branch pass
                // no `cpmxresult`, which is why most inputs survive `-r`.
                //
                // There is no output to be byte-identical to, so refuse rather
                // than quietly producing an alignment the C cannot.
                if l != 0 {
                    eprintln!(
                        "disttbfast-rs: -r reaches A__align at node {}, where the C \n                           dereferences a freed cpmxhist and segfaults (disttbfast.c:2444).",
                        l
                    );
                    std::process::exit(1);
                }
                st.a_align(
                    cn,
                    &dynmtx,
                    cn.penalty,
                    cn.penalty_ex,
                    &mut group1,
                    &mut group2,
                    &effarr1,
                    &effarr2,
                    outgap,
                    outgap,
                    None,
                    None,
                    mem0[0],
                    1,
                    None,
                    None,
                    0.0,
                    0.0,
                    true,
                    b'-',
                );
            }
        }

        let nscore = crate::score::intergroup_score(
            cn,
            &group1,
            &group2,
            &effarr1,
            &effarr2,
            clus1,
            clus2,
            group1[0].len(),
        );
        if nscore < oscore {
            group1 = oseq1;
            group2 = oseq2;
        }

        for (k, &m) in mem0.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group1[k]);
        }
        for (k, &m) in mem1.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group2[k]);
        }

        nlen[m1] = (0.5 * (nlen[m1] + nlen[m2]) as f64) as i32;
    }
}
