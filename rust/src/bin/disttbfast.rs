//! `disttbfast-rs` — a bit-exact Rust port of MAFFT's `disttbfast` for the
//! progressive, non-FFT path (NW-NS-1 / NW-NS-2).
//!
//! It is a drop-in replacement for `libexec/mafft/disttbfast`: same argument
//! vector, FASTA on stdin, FASTA on stdout, and the same `order` file written
//! into the working directory.  Configurations whose C code path is not ported
//! are refused rather than approximated — see `Args::check_supported`.

use std::io::{Read, Write};

use mafft_rs::align::AlignState;
use mafft_rs::constants::{self, Constants};
use mafft_rs::params::Args;
use mafft_rs::treebase::TreebaseOpts;
use mafft_rs::{dist, falign, io, score, tree, treebase};

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match Args::parse(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    if let Err(u) = args.check_supported() {
        eprintln!("disttbfast-rs: unsupported configuration:\n  {}", u.0);
        std::process::exit(1);
    }

    let mut bytes = Vec::new();
    match &args.inputfile {
        Some(path) => {
            let mut f = std::fs::File::open(path).unwrap_or_else(|e| {
                eprintln!("Cannot open {}: {}", path, e);
                std::process::exit(1);
            });
            f.read_to_end(&mut bytes).expect("read input");
        }
        None => {
            std::io::stdin()
                .read_to_end(&mut bytes)
                .expect("read stdin");
        }
    }

    let input = io::read_data(bytes, args.dorp, args.outnumber);
    let njob = input.njob;
    let mut out = Vec::new();

    if njob < 2 {
        let ungapped: Vec<Vec<u8>> = input.seq.iter().map(|s| io::gappick0(s)).collect();
        io::write_data(&mut out, &input.name, &ungapped);
        std::io::stdout().write_all(&out).unwrap();
        eprintln!("Warning: Only {} sequence found.", njob);
        return;
    }

    // `constants( njob, seq )` sees the raw input, gaps and all; only
    // `--fmodel` reads it.
    let cn = Constants::new_seq(&args, input.dorp, &input.seq);
    if cn.trywarp {
        eprintln!("disttbfast-rs: the shift (warp) DP path is not ported; use -Q >= 10");
        std::process::exit(1);
    }

    // seqcheck()
    for (i, s) in input.seq.iter().enumerate() {
        for (p, &c) in s.iter().enumerate() {
            if (c as usize) >= constants::CHARSIZE || cn.amino_n[c as usize] == -1 {
                eprintln!("Illegal character {}", c as char);
                eprintln!("=== Please check site {} in sequence {}.", p + 1, i + 1);
                std::process::exit(1);
            }
        }
    }

    if input.dorp == b'p' && args.tuplesize != 6 {
        eprintln!("tuplesize must be 6 for aa sequence");
        std::process::exit(1);
    }
    if input.dorp == b'd' && args.tuplesize != 6 && args.tuplesize != 10 {
        eprintln!("tuplesize must be 6 or 10 for dna sequence");
        std::process::exit(1);
    }
    if input.dorp == b'd' && args.tuplesize != 6 {
        eprintln!("disttbfast-rs: only -F 6 (6-mer) is ported for nucleotides");
        std::process::exit(1);
    }

    let (mut mtx, _nogaplen) = dist::distance_matrix(&cn, &input.seq, njob);

    let mut nlen: Vec<i32> = input.seq.iter().map(|s| s.len() as i32).collect();
    let mut st = AlignState::new();
    let mut fs = falign::FalignState::new();
    let mut bseq: Vec<Vec<u8>> = Vec::new();

    let nguidetree = args.nguidetree.max(1) as usize;
    for iguidetree in 0..nguidetree {
        let calcpairdists = iguidetree != nguidetree - 1;

        let mut tr = tree::upgma(njob, &mut mtx, args.sueff_global);

        // The `order` file the mafft wrapper uses for --reorder / clustal out.
        {
            let order = tree::topolorderz(&tr.topol, &tr.dep, njob - 2, 2);
            let mut s = String::new();
            for v in order.iter().take(njob) {
                s.push_str(&format!("{}\n", v));
            }
            if let Err(e) = std::fs::write("order", s) {
                eprintln!("Cannot open 'order': {}", e);
                std::process::exit(1);
            }
        }

        let effarr: Vec<f64> = if args.tbrweight != 0 {
            tree::counteff_simple(njob, &mut tr)
        } else {
            vec![1.0; njob]
        };

        bseq = input.seq.iter().map(|s| io::gappick0(s)).collect();

        let selfscore: Vec<i32> = if calcpairdists {
            input
                .seq
                .iter()
                .map(|s| score::naivepairscore11(&cn, s, s, cn.penalty_dist) as i32)
                .collect()
        } else {
            Vec::new()
        };

        eprintln!(
            "Progressive alignment {}/{}... ",
            iguidetree + 1,
            nguidetree
        );
        let opts = TreebaseOpts {
            outgap: args.outgap,
            use_fft: args.use_fft,
            force_fft: args.force_fft,
            minimumweight: 0.0,
            constraint: None,
            fastathreshold: 0.0,
            use_dynamic_mtx: cn.specificityconsideration != 0.0,
        };
        treebase::treebase(
            &cn, &opts, njob, &mut bseq, &mut nlen, &tr, &effarr, &mut st, &mut fs,
        );

        if calcpairdists {
            let skiptable: Vec<Vec<i32>> = bseq.iter().map(|s| score::makeskiptable(s)).collect();
            for i in 0..njob - 1 {
                let ssi = selfscore[i] as f64;
                for j in i + 1..njob {
                    let ssj = selfscore[j] as f64;
                    let bunbo = if ssi < ssj { ssi } else { ssj };
                    let mut v = if bunbo == 0.0 {
                        2.0
                    } else {
                        let s = score::naivepairscorefast(
                            &cn,
                            &bseq[i],
                            &bseq[j],
                            &skiptable[i],
                            &skiptable[j],
                            cn.penalty_dist,
                        );
                        let mut t = (1.0 - s / bunbo) * 2.0;
                        if t > 10.0 {
                            t = 10.0;
                        }
                        t
                    };
                    if v < 0.0 {
                        eprintln!("WARNING: negative distance, iscoretmp = {}", v);
                        v = 0.0;
                    }
                    mtx[i][j - i] = v;
                }
            }
        }
    }

    io::write_data(&mut out, &input.name, &bseq);
    std::io::stdout().write_all(&out).unwrap();
}
