//! `dndpre-rs` — a bit-exact Rust port of MAFFT's `dndpre`.
//!
//! Turns an *already aligned* input into the `hat2` distance matrix the second
//! progressive cycle builds its guide tree from. `--auto` reaches it in the
//! FFT-NS-i and FFT-NS-2 brackets, where the script runs a first alignment,
//! re-derives distances from it with this binary, and aligns again.
//!
//! The distance is `maxdist * (1 - score(i,j) / min(self(i), self(j)))`, with
//! the scores coming from `naivescore`. Two clamps follow, and both are the
//! C's: a negative distance becomes 0.0 with a warning, and anything above 9.0
//! becomes **9.9** — not 9.0 — also with a warning.
#![forbid(unsafe_code)]

use std::io::Read;

use mafft_rs::cargs::{self, Opt};
use mafft_rs::constants::Constants;
use mafft_rs::naivescore::{make_skiptable, naive_pairscore11, naive_pairscore_fast};
use mafft_rs::params::{Args, NOTSPECIFIED};
use mafft_rs::io;

/// What the option loop leaves behind.
struct DndArgs {
    inputfile: Option<String>,
    nadd: i32,
    nblosum: i32,
    maxdist: i32,
    dorp: Option<u8>,
    ppenalty: i32,
    ppenalty_ex: i32,
    poffset: i32,
    fmodel: i32,
    /// `-Z`: score without the gap penalty, which the script passes on one of
    /// its two call sites.
    usenaivescore: bool,
}

/// `atof( s ) * 1000 - 0.5`, truncated — the C's way of turning a penalty
/// given in score units into its internal integer scale.
fn scaled(s: Option<&str>) -> i32 {
    let v: f64 = s.unwrap_or("0").trim().parse().unwrap_or_else(|_| {
        // `atof` yields 0.0 for unparseable text rather than complaining.
        let t = s.unwrap_or("");
        let end = t
            .find(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == '+'))
            .unwrap_or(t.len());
        t[..end].parse().unwrap_or(0.0)
    });
    (v * 1000.0 - 0.5) as i32
}

/// Walk `arguments()`.
fn arguments(argv: &[String]) -> DndArgs {
    let mut out = DndArgs {
        inputfile: None,
        nadd: 0,
        nblosum: 62,
        maxdist: 2,
        dorp: None,
        ppenalty: NOTSPECIFIED,
        ppenalty_ex: NOTSPECIFIED,
        poffset: NOTSPECIFIED,
        fmodel: 0,
        usenaivescore: false,
    };
    // Unlike every other binary here, `dndpre`'s switch has **no `default:`**.
    // An unrecognised option character falls through it, the scan carries on,
    // and `argc` is untouched -- so an unknown option is silently ignored
    // rather than fatal. Returning `Flag` for anything unclassified is exactly
    // that behaviour.
    let opts = cargs::parse(argv, |c| match c {
        'y' | 'I' | 'f' | 'g' | 'h' | 'k' | 'b' | 'j' | 'm' | 'i' | 'M' | 'C' => {
            Some(Opt::Value)
        }
        _ => Some(Opt::Flag),
    });
    for o in &opts {
        let v = o.value.as_deref();
        match o.c {
            'y' => {
                let distout = v.and_then(|s| s.bytes().next()).unwrap_or(0);
                eprintln!("distout={}", distout as char);
                if distout != b'c' && distout != b'h' {
                    eprintln!("Set -y c or -y h in v>=7.521.");
                    std::process::exit(1);
                }
            }
            'Z' => out.usenaivescore = true,
            'D' => out.dorp = Some(b'd'),
            'P' => out.dorp = Some(b'p'),
            'a' => out.fmodel = 1,
            'I' => {
                out.nadd = cargs::myatoi(v) as i32;
                eprintln!("nadd = {}", out.nadd);
            }
            'f' => out.ppenalty = scaled(v),
            'g' => out.ppenalty_ex = scaled(v),
            'h' => out.poffset = scaled(v),
            'b' => out.nblosum = cargs::myatoi(v) as i32,
            'M' => out.maxdist = cargs::myatoi(v) as i32,
            'i' => out.inputfile.clone_from(&o.value),
            // `-t`, `-:`, `-K`, `-C`, `-k`, `-j`, `-m` reach this port only in
            // forms the `--auto` script never emits; `-C` in particular selects
            // a thread count, and this port is single-threaded by construction
            // (the C's own `nthread == 0` branch, which is what `-C 0` asks
            // for). The sequential loop is the one that defines the output.
            _ => {}
        }
    }
    out
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = arguments(&argv);

    let bytes = match &args.inputfile {
        Some(path) => match std::fs::read(path) {
            Ok(b) => b,
            Err(_) => {
                eprintln!("Cannot open {path}");
                std::process::exit(1);
            }
        },
        None => {
            let mut buf = Vec::new();
            let _ = std::io::stdin().read_to_end(&mut buf);
            buf
        }
    };

    let input = io::read_data(bytes, args.dorp, false, false);
    if input.njob < 2 {
        eprintln!("At least 2 sequences should be input!");
        eprintln!("Only {} sequence found.", input.njob);
        std::process::exit(1);
    }
    // `njob -= nadd`: sequences being added are not part of the distance
    // matrix, and the C simply shortens the loop bounds rather than filtering.
    let njob = (input.njob as i32 - args.nadd).max(0) as usize;
    let seq = &input.seq[..njob];
    let name = &input.name[..njob];

    for row in seq.iter().skip(1) {
        if row.len() != seq[0].len() {
            eprintln!("Not aligned!");
            std::process::exit(1);
        }
    }

    let cargs_for_constants = Args {
        nblosum: args.nblosum,
        scoremtx: 1,
        fmodel: args.fmodel,
        ppenalty: args.ppenalty,
        ppenalty_ex: args.ppenalty_ex,
        poffset: args.poffset,
        dorp: args.dorp,
        ..Default::default()
    };
    let cn = Constants::new_seq(&cargs_for_constants, input.dorp, seq);

    // `seqcheck()`. Note dndpre prints only the character, where the aligners
    // also name the site -- kept as it is, because the message is the contract.
    for row in seq {
        for &c in row {
            if (c as usize) >= cn.charsize || cn.amino_n[c as usize] == -1 {
                eprintln!("Illegal character {}", c as char);
                std::process::exit(1);
            }
        }
    }

    let penal = if args.usenaivescore { 0 } else { cn.penalty };
    let selfscore: Vec<f64> = seq
        .iter()
        .map(|s| naive_pairscore11(s, s, cn.penalty, &cn.amino_dis))
        .collect();
    let skiptable = make_skiptable(seq);

    let maxdist = f64::from(args.maxdist);
    let mut mtx = vec![vec![0.0f64; njob]; njob];
    for i in 0..njob.saturating_sub(1) {
        let ssi = selfscore[i];
        for j in (i + 1)..njob {
            let bunbo = ssi.min(selfscore[j]);
            let mut mtxv = if bunbo == 0.0 {
                maxdist
            } else {
                maxdist
                    * (1.0
                        - naive_pairscore_fast(
                            &seq[i],
                            &seq[j],
                            &skiptable[i],
                            &skiptable[j],
                            penal,
                            &cn.amino_dis,
                        ) / bunbo)
            };
            if mtxv < 0.0 {
                eprintln!("WARNING: negative distance, mtxv = {mtxv:.6}");
                mtxv = 0.0;
            }
            if mtxv > 9.0 {
                eprintln!("WARNING: Distance {i}-{j} is strange, {mtxv:.6}.");
                mtxv = 9.9;
            }
            // `io::write_hat2` indexes `mtx[i][j - i]`, the offset layout its
            // other caller stores; the C's `WriteHat2_pointer` indexes a
            // square matrix. Same file either way -- what matters is that the
            // writer and the filler agree, so the distance goes where this
            // writer will look for it.
            mtx[i][j - i] = mtxv;
        }
    }

    io::write_hat2(njob, name, &mtx);
}
