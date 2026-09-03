//! The progressive alignment driver, transliterated from `treebase()` in
//! `core/disttbfast.c` (the single-threaded, unconstrained, non-FFT path).

use crate::align::{AlignState, CpmxHist};
use crate::constants::Constants;
use crate::params::Args;
use crate::tree::Tree;

/// `fastconjuction_noname()`: gather a group's members and normalise their
/// weights, returning the pre-normalisation total (`oritotal`).
fn fastconjuction_noname(memlist: &[i32], eff: &[f64]) -> (Vec<f64>, f64) {
    let mut peff: Vec<f64> = memlist.iter().map(|&m| eff[m as usize]).collect();
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
pub fn treebase(
    cn: &Constants,
    args: &Args,
    njob: usize,
    aseq: &mut [Vec<u8>],
    nlen: &mut [i32],
    tree: &Tree,
    effarr: &[f64],
    st: &mut AlignState,
) {
    let newgapstr = b'-';
    let outgap = args.outgap;

    let mut memhist: Vec<Option<Vec<i32>>> = vec![None; njob.saturating_sub(1)];
    let mut cpmxhist: Vec<Option<CpmxHist>> = (0..njob.saturating_sub(1)).map(|_| None).collect();

    for l in 0..njob - 1 {
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

        let (effarr1, orieff1) = fastconjuction_noname(&mem0, effarr);
        let (effarr2, orieff2) = fastconjuction_noname(&mem1, effarr);
        let clus1 = mem0.len();
        let clus2 = mem1.len();

        // Gather the group members (mseq1 / mseq2 alias rows of aseq in C).
        let mut group1: Vec<Vec<u8>> = mem0.iter().map(|&m| aseq[m as usize].clone()).collect();
        let mut group2: Vec<Vec<u8>> = mem1.iter().map(|&m| aseq[m as usize].clone()).collect();

        if clus1 == 1 && clus2 == 1 {
            let mut s1 = std::mem::take(&mut group1[0]);
            let mut s2 = std::mem::take(&mut group2[0]);
            st.g_align11(
                cn,
                &cn.n_dis_consweight_multi,
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
            let (_pscore, hist) = st.a_align(
                cn,
                &cn.n_dis_consweight_multi,
                cn.penalty,
                cn.penalty_ex,
                &mut group1,
                &mut group2,
                &effarr1,
                &effarr2,
                outgap,
                outgap,
                mem0[0],
                cpmxchild0,
                cpmxchild1,
                orieff1,
                orieff2,
                newgapstr,
            );
            cpmxhist[l] = hist;
        }

        for (k, &m) in mem0.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group1[k]);
        }
        for (k, &m) in mem1.iter().enumerate() {
            aseq[m as usize] = std::mem::take(&mut group2[k]);
        }

        // `nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] )` -- assigned to an int.
        nlen[m1] = (0.5 * (nlen[m1] + nlen[m2]) as f64) as i32;
    }
}
