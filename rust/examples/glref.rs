//! Mirror of `/tmp/glharness.c`: run `gaplencount` and the `calcpfac` family
//! on two groups given on the command line and print the same lines the C
//! harness does, so the two can be diffed.
use mafft_rs::dalign::*;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let mut p = 0usize;
    let n1: usize = a[p].parse().unwrap();
    p += 1;
    let s1: Vec<Vec<u8>> = (0..n1).map(|k| a[p + k].as_bytes().to_vec()).collect();
    p += n1;
    let e1: Vec<f64> = (0..n1).map(|k| a[p + k].parse().unwrap()).collect();
    p += n1;
    let n2: usize = a[p].parse().unwrap();
    p += 1;
    let s2: Vec<Vec<u8>> = (0..n2).map(|k| a[p + k].as_bytes().to_vec()).collect();
    p += n2;
    let e2: Vec<f64> = (0..n2).map(|k| a[p + k].parse().unwrap()).collect();

    let l1 = s1[0].len();
    let l2 = s2[0].len();
    let r1: Vec<&[u8]> = s1.iter().map(|s| s.as_slice()).collect();
    let r2: Vec<&[u8]> = s2.iter().map(|s| s.as_slice()).collect();
    let mut g1 = vec![Col::Null; l1 + 2];
    let mut g2 = vec![Col::Null; l2 + 2];
    gaplencount(n1, l1, &mut g1, &r1, &e1);
    gaplencount(n2, l2, &mut g2, &r2, &e2);

    let dump = |tag: u8, g: &[Col], l: usize| {
        for (j, c) in g.iter().enumerate().take(l + 1) {
            match c.get() {
                None => println!("{} {} NULL", tag, j),
                Some(v) => {
                    for (k, e) in v.iter().enumerate().take_while(|(_, e)| e.idatnext != -1) {
                        println!(
                            "{} {} {} len={} relend={} freq={:016x} next={} prev={} npat={}",
                            tag, j, k, e.len, e.relend, e.freq.to_bits(), e.idatnext, e.idatprev, e.npat
                        );
                    }
                }
            }
        }
    };
    dump(1, &g1, l1);
    dump(2, &g2, l2);
    for i in 0..=l1 {
        for j in 0..=l2 {
            println!(
                "pf {} {} {:016x} gap={:016x} gapex={:016x}",
                i, j,
                calcpfacnoidatend(&g1[i], &g2[j]).to_bits(),
                calcpfac_gap_noidatend(&g1[i], &g2[j], 2).to_bits(),
                calcpfac_gapex_noidatend(&g1[i], &g2[j], 2).to_bits()
            );
        }
    }
}
