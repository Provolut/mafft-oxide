//! One staged run of the refinement loop, shared by `tests/tdi.rs` and
//! `examples/tdiref.rs`.
//!
//! The argument vector is `gen_tdicases.py`'s, and the printed lines are
//! `gen_tdiharness.c`'s, byte for byte.
use mafft_rs::constants::Constants;
use mafft_rs::localhom::{new_table, LocalHom, LocalHomTable};
use mafft_rs::params::Args;
use mafft_rs::tditeration::{tree_dependent_iteration, IterParams};
use mafft_rs::treeop::Topol;

fn members(s: &str) -> Vec<i32> {
    let mut v: Vec<i32> = s.split(',').map(|x| x.parse().unwrap()).collect();
    v.push(-1);
    v
}

/// Run one refinement and return what the C harness prints for it.
///
/// Args:
///     a: The case's argument vector, as `gen_tdicases.py` emits it.
///
/// Returns:
///     The outcome code and every row of the refined alignment.
pub fn run(a: &[String]) -> String {
    let mut p = 0usize;
    let nseq: usize = a[p].parse().unwrap();
    let niter: usize = a[p + 1].parse().unwrap();
    let cut: f64 = a[p + 2].parse().unwrap();
    let constraint: i32 = a[p + 3].parse().unwrap();
    let fastathreshold: f64 = a[p + 4].parse().unwrap();
    p += 5;

    let mut aseq: Vec<Vec<u8>> = (0..nseq).map(|k| a[p + k].as_bytes().to_vec()).collect();
    p += nseq;

    let mut topol: Topol = Vec::with_capacity(nseq);
    let mut len: Vec<[f64; 2]> = Vec::with_capacity(nseq);
    for _ in 0..nseq - 1 {
        topol.push([members(&a[p]), members(&a[p + 1])]);
        len.push([a[p + 2].parse().unwrap(), a[p + 3].parse().unwrap()]);
        p += 4;
    }
    // The C allocates `nseq` rows of each and uses `nseq-1`; `tree_cnv` reads
    // `len[nseq-2]` when it unroots.
    topol.push([vec![-1], vec![-1]]);
    len.push([0.0, 0.0]);

    let mut skip = vec![[false; 2]; nseq];
    for row in skip.iter_mut().take(nseq - 1) {
        row[0] = a[p] == "1";
        row[1] = a[p + 1] == "1";
        p += 2;
    }

    let nlh: usize = a[p].parse().unwrap();
    p += 1;
    let mut table: Option<LocalHomTable> = if constraint != 0 {
        Some(new_table(nseq))
    } else {
        None
    };
    for _ in 0..nlh {
        let i: usize = a[p].parse().unwrap();
        let j: usize = a[p + 1].parse().unwrap();
        let e = LocalHom {
            start1: a[p + 2].parse().unwrap(),
            end1: a[p + 3].parse().unwrap(),
            start2: a[p + 4].parse().unwrap(),
            end2: a[p + 5].parse().unwrap(),
            opt: a[p + 6].parse().unwrap(),
            korh: b'h',
            ..LocalHom::default()
        };
        p += 7;
        if let Some(t) = table.as_mut() {
            let chain = &mut t[i][j - i];
            // The head slot is the one `new_table` made; a second segment
            // between the same pair is appended, as the hat3 reader does.
            let mut e = e;
            e.overlapaa = e.end1 - e.start1 + 1;
            if chain[0].opt == -1.0 {
                chain[0] = e;
            } else {
                chain.push(e);
            }
        }
    }

    let args = Args {
        nblosum: 62,
        ppenalty: -1530,
        ppenalty_ex: 0,
        poffset: 0,
        penalty_shift_factor: 100.0,
        ..Default::default()
    };
    let cn = Constants::new(&args);

    let params = IterParams {
        niter,
        cut,
        cooling: true,
        minimumweight: 0.00001,
        constraint,
        weight: 4,
        score_check: 1,
        use_fft: true,
        utree: true,
        alg: b'A',
        nkozo: 0,
        outgap: true,
        legacygapcost: false,
        newgapstr: b'-',
        fastathreshold,
        // `dvtditr.c`: `alloclen = nlenmax * 9`.
        alloclen: aseq[0].len() * 9,
    };

    let mut bseq: Vec<Vec<u8>> = aseq.clone();
    let outcome = tree_dependent_iteration(
        &cn,
        &params,
        nseq,
        &mut aseq,
        &mut bseq,
        &topol,
        &mut len,
        &skip,
        table.as_mut(),
        &mut |_| {},
    );

    let mut out = String::new();
    out.push_str(&format!("outcome {}\n", outcome));
    for (i, s) in bseq.iter().enumerate() {
        out.push_str(&format!("{} {}\n", i, String::from_utf8_lossy(s)));
    }
    out
}
