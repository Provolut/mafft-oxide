//! One staged `falign_localhom` call, shared by `tests/flh.rs` and
//! `examples/flhref.rs` so the golden replay and the live differential drive
//! the port through exactly the same code.
//!
//! The argument vector is `gen_flhcases.py`'s, and the printed lines are
//! `gen_flhharness.c`'s, byte for byte.
use mafft_rs::constants::Constants;
use mafft_rs::falign::{falign_localhom, FalignState};
use mafft_rs::localhom::LocalHom;
use mafft_rs::params::Args;
use mafft_rs::partalign::{PartAlignState, PartImp};

/// Run one case and return what the C harness prints for it.
///
/// Args:
///     a: The case's argument vector, as `gen_flhcases.py` emits it.
///
/// Returns:
///     The score and constraint score as f64 bit patterns, then every row of
///     both groups.
pub fn run(a: &[String]) -> String {
    let mut p = 0usize;
    let n1: usize = a[p].parse().unwrap();
    p += 1;
    let mut g1: Vec<Vec<u8>> = (0..n1).map(|k| a[p + k].as_bytes().to_vec()).collect();
    p += n1;
    let e1: Vec<f64> = (0..n1).map(|k| a[p + k].parse().unwrap()).collect();
    p += n1;
    let o1: Vec<i32> = (0..n1).map(|k| a[p + k].parse().unwrap()).collect();
    p += n1;
    let n2: usize = a[p].parse().unwrap();
    p += 1;
    let mut g2: Vec<Vec<u8>> = (0..n2).map(|k| a[p + k].as_bytes().to_vec()).collect();
    p += n2;
    let e2: Vec<f64> = (0..n2).map(|k| a[p + k].parse().unwrap()).collect();
    p += n2;
    let o2: Vec<i32> = (0..n2).map(|k| a[p + k].parse().unwrap()).collect();
    p += n2;
    let outgap = a[p] == "1";
    p += 1;
    let legacygapcost = a[p] == "1";
    p += 1;
    let fastathreshold: f64 = a[p].parse().unwrap();
    p += 1;
    let fft_win_size: i32 = a[p].parse().unwrap();
    p += 1;
    let fft_threshold: i32 = a[p].parse().unwrap();
    p += 1;
    let nlh: usize = a[p].parse().unwrap();
    p += 1;

    // One chain per member pair, in the shape `fillimp` walks. `korh` is 'h'
    // and `rimportance` mirrors `importance`, as `msshrinklocalhom_fast_half`
    // leaves them for a run with no structural input.
    let mut chains: Vec<Vec<Vec<LocalHom>>> = vec![vec![Vec::new(); n2]; n1];
    for _ in 0..nlh {
        let i: usize = a[p].parse().unwrap();
        let j: usize = a[p + 1].parse().unwrap();
        let mut e = LocalHom {
            start1: a[p + 2].parse().unwrap(),
            end1: a[p + 3].parse().unwrap(),
            start2: a[p + 4].parse().unwrap(),
            end2: a[p + 5].parse().unwrap(),
            opt: a[p + 6].parse().unwrap(),
            importance: a[p + 7].parse().unwrap(),
            ..LocalHom::default()
        };
        e.rimportance = e.importance;
        e.overlapaa = e.end1 - e.start1 + 1;
        e.korh = b'h';
        chains[i][j].push(e);
        p += 8;
    }
    let shrink: Vec<Vec<Option<&Vec<LocalHom>>>> = chains
        .iter()
        .map(|row| {
            row.iter()
                .map(|c| if c.is_empty() { None } else { Some(c) })
                .collect()
        })
        .collect();

    let args = Args {
        nblosum: 62,
        ppenalty: -1530,
        ppenalty_ex: 0,
        poffset: 0,
        penalty_shift_factor: 100.0,
        fft_win_size,
        fft_threshold,
        ..Default::default()
    };
    let cn = Constants::new(&args);

    let lgth = g1[0].len();
    let mut imp = PartImp::new();
    {
        let r1: Vec<&[u8]> = g1.iter().map(std::vec::Vec::as_slice).collect();
        let r2: Vec<&[u8]> = g2.iter().map(std::vec::Vec::as_slice).collect();
        imp.init_strict(
            n1, n2, lgth, lgth, &r1, &r2, &e1, &e2, &shrink, &o1, &o2, fastathreshold,
        );
    }

    let mut fs = FalignState::new();
    let mut past = PartAlignState::default();
    let (score, totalimp) = falign_localhom(
        &mut fs,
        &mut past,
        &imp,
        &cn,
        &mut g1,
        &mut g2,
        &e1,
        &e2,
        outgap,
        legacygapcost,
        b'-',
    );

    let mut out = String::new();
    out.push_str(&format!("score {:016x}\n", score.to_bits()));
    out.push_str(&format!("imp {:016x}\n", totalimp.to_bits()));
    for (i, s) in g1.iter().enumerate() {
        out.push_str(&format!("1 {} {}\n", i, String::from_utf8_lossy(s)));
    }
    for (i, s) in g2.iter().enumerate() {
        out.push_str(&format!("2 {} {}\n", i, String::from_utf8_lossy(s)));
    }
    out
}
