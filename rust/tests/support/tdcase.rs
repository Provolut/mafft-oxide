//! One staged guide tree through the branch-weighting family, shared by
//! `tests/td.rs` and `examples/tdref.rs` so the golden replay and the live
//! differential drive the port through the same code.
//!
//! The argument vector is `gen_tdcases.py`'s, and the printed lines are
//! `gen_tdharness.c`'s, byte for byte.
use mafft_rs::treeop::{
    calc_branch_weight, dist_from_a_branch, one_cluster_and_the_other_fast, tree_cnv,
    weight_from_a_branch, Topol,
};

/// Parse `"0,3,7"` into the C's `-1`-terminated member list.
fn members(s: &str) -> Vec<i32> {
    let mut v: Vec<i32> = s.split(',').map(|x| x.parse().unwrap()).collect();
    v.push(-1);
    v
}

/// Run one case and return what the C harness prints for it.
///
/// Args:
///     a: The case's argument vector, as `gen_tdcases.py` emits it.
///
/// Returns:
///     Branch weights, the clamped branch lengths, and the distance and weight
///     vectors of every branch, all as f64 bit patterns.
pub fn run(a: &[String]) -> String {
    let mut p = 0usize;
    let nseq: usize = a[p].parse().unwrap();
    p += 1;

    let mut topol: Topol = Vec::with_capacity(nseq);
    let mut len: Vec<[f64; 2]> = Vec::with_capacity(nseq);
    for _ in 0..nseq - 1 {
        let m0 = members(&a[p]);
        let m1 = members(&a[p + 1]);
        let l0: f64 = a[p + 2].parse().unwrap();
        let l1: f64 = a[p + 3].parse().unwrap();
        topol.push([m0, m1]);
        len.push([l0, l1]);
        p += 4;
    }
    // The C allocates `nseq` rows and uses `nseq - 1`; the last row is read by
    // nothing here, but `bw[locnseq-2]` is written, so the row count must match.
    topol.push([vec![-1], vec![-1]]);
    len.push([0.0, 0.0]);

    let stopol = tree_cnv(nseq, &topol, &mut len);
    let mut bw = vec![[0.0f64; 2]; nseq];
    calc_branch_weight(&mut bw, nseq, &stopol, &topol);

    let mut out = String::new();
    for i in 0..nseq - 1 {
        for j in 0..2 {
            out.push_str(&format!("bw {} {} {:016x}\n", i, j, bw[i][j].to_bits()));
        }
    }
    for i in 0..nseq - 1 {
        for j in 0..2 {
            out.push_str(&format!("len {} {} {:016x}\n", i, j, len[i][j].to_bits()));
        }
    }

    let mut result = vec![0.0f64; nseq];
    for i in 0..nseq - 1 {
        for j in 0..2 {
            dist_from_a_branch(nseq, &mut result, &stopol, &topol, &len, i, j);
            for (k, v) in result.iter().enumerate() {
                out.push_str(&format!("dist {} {} {} {:016x}\n", i, j, k, v.to_bits()));
            }
            weight_from_a_branch(nseq, &mut result, &stopol, &bw, &topol, i, j);
            for (k, v) in result.iter().enumerate() {
                out.push_str(&format!("wt {} {} {} {:016x}\n", i, j, k, v.to_bits()));
            }
        }
    }
    for i in 0..nseq - 1 {
        for j in 0..2 {
            let (m1, m2) = one_cluster_and_the_other_fast(nseq, &topol, i, j);
            out.push_str(&format!("split {} {} s {} {}\n", i, j, m1[0], m2[0]));
            for (k, &v) in m1.iter().take_while(|&&v| v != -1).enumerate() {
                out.push_str(&format!("split {} {} m1 {} {}\n", i, j, k, v));
            }
            for (k, &v) in m2.iter().take_while(|&&v| v != -1).enumerate() {
                out.push_str(&format!("split {} {} m2 {} {}\n", i, j, k, v));
            }
        }
    }
    out
}
