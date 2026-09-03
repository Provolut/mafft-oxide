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
use mafft_rs::{addfunctions, dist, falign, io, score, tree, treebase};

/// `pickup()`: keep the longest sequence on each side of the last few tree
/// nodes and write the rest to `notused`.
///
/// The C's comment says it is not adapted to memsave, and it is not: with the
/// memsave `topol` each side holds a single member, so the "longest" search
/// runs over one candidate. Two more quirks come with it -- `m` counts loop
/// iterations while each iteration selects two sequences, and `m` is then used
/// as the record count for stdout, so the file can be shorter than the
/// selection.
fn pickup(
    n: usize,
    sparsepickup: i32,
    seqlen: &[usize],
    tr: &tree::Tree,
    name: &[Vec<u8>],
    seq: &[Vec<u8>],
) {
    let mut longestseq = vec![[-1i32, -1i32]; n.saturating_sub(1)];
    let mut select = vec![false; n];
    for i in 0..n - 1 {
        for side in 0..2 {
            // memsave: one member, then the -1 terminator, so the C's search
            // for the longest compares a single candidate against -1.
            let m = tr.topol[i][side];
            if m != -1 && seqlen[m as usize] as i64 > -1 {
                longestseq[i][side] = m;
            }
        }
    }
    let mut m = 1i32;
    for i in (0..n - 1).rev() {
        for side in 0..2 {
            let p = longestseq[i][side];
            if p >= 0 {
                select[p as usize] = true;
            }
        }
        m += 1;
        if m >= sparsepickup {
            break;
        }
    }

    let mut nameout = Vec::new();
    let mut seqout = Vec::new();
    let mut namenot = Vec::new();
    let mut seqnot = Vec::new();
    for i in 0..n {
        if select[i] {
            nameout.push(name[i].clone());
            seqout.push(seq[i].clone());
        } else {
            namenot.push(name[i].clone());
            seqnot.push(seq[i].clone());
        }
    }
    let mut out = Vec::new();
    let keep = (m as usize).min(nameout.len());
    io::write_data(&mut out, &nameout[..keep], &seqout[..keep]);
    std::io::stdout().write_all(&out).unwrap();
    let mut nu = Vec::new();
    let keepn = (n - m as usize).min(namenot.len());
    io::write_data(&mut nu, &namenot[..keepn], &seqnot[..keepn]);
    std::fs::write("notused", nu).expect("write notused");
}

/// `includemember()`: is every member of `mem` also in `cand`?
fn includemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.iter().all(|m| cand.contains(m))
}

/// `samemembern()`: `mem` has exactly `cand.len()` members and they are the
/// same set.
fn samemembern(mem: &[i32], cand: &[i32]) -> bool {
    mem.len() == cand.len() && mem.iter().all(|m| cand.contains(m))
}

/// `samemember()`: the same members, in any order. (The C has two definitions
/// of this; the one under `#if 0` compares position by position, the live one
/// compares as sets, and `topolorderz` does not hand them out in `addmem`'s
/// order.)
fn samemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.len() == cand.len() && mem.iter().all(|m| cand.contains(m))
}

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
    // `check_guidetreefile()`: the first line of `_guidetree` selects between a
    // real tree and a family of generated ones. Only the real tree is ported.
    if args.treein {
        let first = std::fs::read_to_string("_guidetree")
            .map(|t| t.lines().next().unwrap_or("").to_string())
            .unwrap_or_else(|_| {
                eprintln!("cannot open _guidetree");
                std::process::exit(1);
            });
        for kw in [
            "shuffle",
            "pileup",
            "auto",
            "test",
            "compact",
            "very compact",
            "stepadd",
            "youngestlinkage",
            "nodepair",
        ] {
            if first.starts_with(kw) {
                eprintln!("disttbfast-rs: '{}' guide trees are not ported", kw);
                std::process::exit(1);
            }
        }
    }

    // With a user tree there is no distance matrix at all.
    let (mut mtx, nogaplen): (Vec<Vec<f64>>, Vec<usize>) = if args.treein {
        (Vec::new(), input.seq.iter().map(|s| s.len()).collect())
    } else {
        dist::distance_matrix(&cn, &input.seq, njob, args.tuplesize as usize)
    };

    let mut nlen: Vec<i32> = input.seq.iter().map(|s| s.len() as i32).collect();
    let mut st = AlignState::new();
    let mut fs = falign::FalignState::new();
    let mut bseq: Vec<Vec<u8>> = Vec::new();

    let nguidetree = if args.treein {
        1
    } else {
        args.nguidetree.max(1) as usize
    };
    for iguidetree in 0..nguidetree {
        let calcpairdists = iguidetree != nguidetree - 1;

        // `-y`: emit the k-tuple distance matrix. This happens before the tree
        // is built, and with `-T` (noalign) and no `-t` the run stops here --
        // before `order` is written.
        if args.distout != 0 && !args.treein {
            io::write_hat2(njob, &input.name, &mtx);
            if args.treeout == 0 && args.noalign {
                io::write_data(&mut out, &input.name, &input.seq);
                std::io::stdout().write_all(&out).unwrap();
                eprintln!();
                return;
            }
        }

        let mut tr = if args.treein {
            match tree::load_tree(
                njob,
                if args.treeout != 0 {
                    Some(&input.name[..])
                } else {
                    None
                },
                args.outnumber,
            ) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
        } else if args.treeout != 0 {
            tree::upgma_treeout(
                njob,
                &mut mtx,
                args.sueff_global,
                Some(&input.name),
                args.outnumber,
                args.treeout == 2,
            )
        } else {
            tree::upgma(njob, &mut mtx, args.sueff_global)
        };

        // `-q N`: keep a sparse subset and stop, before `order` is written.
        if args.sparsepickup != 0 && iguidetree == nguidetree - 1 {
            eprintln!("Sparsepickup! ");
            pickup(
                njob,
                args.sparsepickup,
                &nogaplen,
                &tr,
                &input.name,
                &input.seq,
            );
            eprintln!("done. ");
            return;
        }

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

        if (args.treeout != 0 || args.distout != 0) && args.noalign {
            io::write_data(&mut out, &input.name, &input.seq);
            std::io::stdout().write_all(&out).unwrap();
            eprintln!();
            return;
        }

        let effarr: Vec<f64> = if args.tbrweight != 0 {
            tree::counteff_simple(njob, &mut tr)
        } else {
            vec![1.0; njob]
        };

        // `--add`: the first `njob - nadd` sequences are an existing alignment
        // and keep their columns; only the added ones are ungapped. Everything
        // else takes the plain path.
        let mut mergeoralign: Vec<u8> = vec![b'a'; njob.saturating_sub(1)];
        if args.nadd != 0 {
            let nadd = args.nadd as usize;
            let alignmentlength = input.seq[0].len();
            for i in 0..njob - nadd {
                if input.seq[i].len() != alignmentlength {
                    eprintln!("# The original {} sequences must be aligned", njob - nadd);
                    std::process::exit(1);
                }
            }
            let addmem: Vec<i32> = ((njob - nadd)..njob).map(|i| i as i32).collect();
            if args.addprofile {
                let alignmentlength = input.seq[njob - nadd].len();
                for i in (njob - nadd)..njob {
                    if input.seq[i].len() != alignmentlength {
                        eprintln!("# The {} additional sequences must be aligned", nadd);
                        std::process::exit(1);
                    }
                }
                let mut found = false;
                for i in 0..njob - 1 {
                    let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                    let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                    if samemember(&m0, &addmem) {
                        mergeoralign[i] = b'1';
                        found = true;
                    } else if samemember(&m1, &addmem) {
                        mergeoralign[i] = b'2';
                        found = true;
                    } else {
                        mergeoralign[i] = b'n';
                    }
                }
                if !found {
                    eprintln!("Not found the accurate branch.");
                    std::process::exit(1);
                }
            } else {
                // Plain `--add` starts every node at 'n' -- nodes inside the
                // existing alignment are already done -- and then two passes
                // mark the ones that touch an added sequence.
                for v in mergeoralign.iter_mut() {
                    *v = b'n';
                }
                for j in (njob - nadd)..njob {
                    let one = [j as i32];
                    for i in 0..njob - 1 {
                        let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                        let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                        if samemembern(&m0, &one) {
                            mergeoralign[i] = if mergeoralign[i] != b'n' { b'w' } else { b'1' };
                        } else if samemembern(&m1, &one) {
                            mergeoralign[i] = if mergeoralign[i] != b'n' { b'w' } else { b'2' };
                        }
                    }
                }
                for i in 0..njob - 1 {
                    let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                    let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                    let in0 = includemember(&m0, &addmem);
                    let in1 = includemember(&m1, &addmem);
                    if in0 && in1 {
                        mergeoralign[i] = b'w';
                    } else if in0 {
                        mergeoralign[i] = b'1';
                    } else if in1 {
                        mergeoralign[i] = b'2';
                    }
                }
            }
        }

        bseq = input.seq.iter().map(|s| io::gappick0(s)).collect();
        if args.nadd != 0 {
            let nadd = args.nadd as usize;
            // `--addprofile` keeps the added block's own columns: the C
            // collapses only its all-gap ones and copies it across verbatim,
            // where plain `--add` ungaps each added sequence (done above).
            if args.addprofile {
                let mut tail: Vec<Vec<u8>> = input.seq[njob - nadd..].to_vec();
                addfunctions::commongappick(&mut tail);
                for (i, s) in tail.into_iter().enumerate() {
                    bseq[njob - nadd + i] = s;
                }
            }
            let mut head: Vec<Vec<u8>> = input.seq[..njob - nadd].to_vec();
            addfunctions::commongappick(&mut head);
            for (i, s) in head.into_iter().enumerate() {
                bseq[i] = s;
            }
        }

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
            mergeoralign: if args.nadd != 0 {
                Some(&mergeoralign)
            } else {
                None
            },
            nadd: args.nadd as usize,
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
