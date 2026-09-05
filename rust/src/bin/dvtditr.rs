//! `dvtditr-rs` — a bit-exact Rust port of MAFFT's `dvtditr`.
//!
//! The last C binary a `mafft --auto` run executes. It takes the alignment the
//! progressive stage produced, the distance matrix it was built from, and (for
//! L-INS-i and E-INS-i) the `hat3` constraints, and refines the alignment one
//! branch of the guide tree at a time until nothing improves.
//!
//! Its files are the C's, because the `mafft` script reads them: input and
//! output are both `pre` — read at the start, rewritten at the end — with
//! `hat2` and `hat3` read alongside, and `order`, `trace` and `infile.tree`
//! written.
//!
//! **The implemented configuration is the one the script produces**, and every
//! other one is refused by name rather than approximated. `GOAL-C-FREE.md`
//! records why that is the rule here and what a traced run reports.
#![forbid(unsafe_code)]

use std::io::Read;

use mafft_rs::cargs::{self, Opt};
use mafft_rs::constants::Constants;
use mafft_rs::io;
use mafft_rs::localhom::new_table;
use mafft_rs::params::{Args, NOTSPECIFIED};
use mafft_rs::tditeration::{tree_dependent_iteration, IterParams};
use mafft_rs::tree::upgma_treeout;
use mafft_rs::treeop::{expand_topol, search_anchors};

/// What `dvtditr`'s option loop leaves behind.
struct DvtArgs {
    inputfile: Option<String>,
    niter: usize,
    nadd: i32,
    cut: f64,
    minimumweight: f64,
    fastathreshold: f64,
    constraint: i32,
    bunkatsu: bool,
    use_fft: bool,
    nthread: i32,
    sueff_global: f64,
    nblosum: i32,
    ppenalty: i32,
    ppenalty_ex: i32,
    poffset: i32,
    penalty_shift_factor: f64,
    fft_threshold: i32,
    fft_win_size: i32,
    specificityconsideration: f64,
    dorp: Option<u8>,
    fmodel: i32,
    /// `-:` (`--anysymbol`), which rewrites the `n` row and column of the
    /// nucleotide matrix.
    nwildcard: bool,
}

impl Default for DvtArgs {
    fn default() -> Self {
        // dvtditr.c's own initialisers, in its `arguments()`.
        DvtArgs {
            inputfile: None,
            niter: 2,
            nadd: 0,
            cut: 0.0,
            minimumweight: 0.0006,
            fastathreshold: 2.7,
            constraint: 0,
            bunkatsu: true,
            use_fft: false,
            nthread: 0,
            sueff_global: 0.1,
            nblosum: 62,
            ppenalty: NOTSPECIFIED,
            ppenalty_ex: NOTSPECIFIED,
            poffset: NOTSPECIFIED,
            penalty_shift_factor: 100.0,
            fft_threshold: NOTSPECIFIED,
            fft_win_size: NOTSPECIFIED,
            specificityconsideration: 0.0,
            dorp: None,
            fmodel: 0,
            nwildcard: false,
        }
    }
}

/// `atof( s ) * 1000 - 0.5`, truncated — a penalty in score units on the C's
/// internal integer scale.
fn scaled(s: Option<&str>) -> i32 {
    let v: f64 = s.and_then(|x| x.parse().ok()).unwrap_or(0.0);
    (v * 1000.0 - 0.5) as i32
}

fn atof(s: Option<&str>) -> f64 {
    s.and_then(|x| x.parse().ok()).unwrap_or(0.0)
}

/// `arguments()` from `dvtditr.c`.
///
/// Only the options `scripts/mafft` puts on the `dvtditr` line are accepted.
/// Anything else stops the run naming the flag: this port covers one
/// configuration, and silently ignoring a flag that changes the answer is the
/// one failure mode a refusal cannot have.
fn parse(argv: &[String]) -> DvtArgs {
    let mut a = DvtArgs::default();
    let takes_value = |c: char| -> Option<Opt> {
        match c {
            'i' | 'K' | 'I' | 'f' | 'Q' | 'g' | 'h' | 'b' | 'l' | 'C' | 't' | 'p' | 's' | 'X'
            | 'z' | 'w' | 'E' | 'W' => Some(Opt::Value),
            'F' | 'B' | 'a' | 'D' | 'P' | 'N' | ':' => Some(Opt::Flag),
            _ => None,
        }
    };
    for p in cargs::parse(argv, takes_value) {
        match p.c {
            'i' => a.inputfile = p.value.clone(),
            'K' => a.nadd = cargs::myatoi(p.value.as_deref()) as i32,
            'I' => a.niter = cargs::myatoi(p.value.as_deref()) as usize,
            'f' => a.ppenalty = scaled(p.value.as_deref()),
            'g' => a.ppenalty_ex = scaled(p.value.as_deref()),
            'h' => a.poffset = scaled(p.value.as_deref()),
            'Q' => a.penalty_shift_factor = atof(p.value.as_deref()),
            'b' => a.nblosum = cargs::myatoi(p.value.as_deref()) as i32,
            'l' => {
                // `-l` does three things at once: it is the constraint weight,
                // it turns off segmentation, and it is what switches the
                // constrained path on at all.
                a.bunkatsu = false;
                a.fastathreshold = atof(p.value.as_deref());
                a.constraint = 2;
            }
            'C' => a.nthread = cargs::myatoi(p.value.as_deref()) as i32,
            't' => { /* randomseed: read, and used only by a shuffle the C has disabled */ }
            'p' => {
                let v = p.value.clone().unwrap_or_default();
                assert_eq!(
                    v, "BAATARI2",
                    "only the BAATARI2 parallelization strategy is implemented"
                );
            }
            's' => a.specificityconsideration = atof(p.value.as_deref()),
            'X' => a.sueff_global = atof(p.value.as_deref()),
            'z' => a.fft_threshold = cargs::myatoi(p.value.as_deref()) as i32,
            'w' => a.fft_win_size = cargs::myatoi(p.value.as_deref()) as i32,
            'E' => { /* autosubalignment; the script passes 0.0 */ }
            'W' => a.minimumweight = atof(p.value.as_deref()),
            'F' => a.use_fft = true,
            'B' => a.bunkatsu = false,
            'a' => a.fmodel = 1,
            'D' => a.dorp = Some(b'd'),
            'P' => a.dorp = Some(b'p'),
            'N' => { /* nevermemsave; only reached above 50000 residues */ }
            ':' => a.nwildcard = true,
            c => {
                eprintln!("dvtditr-rs: unsupported option -{c}");
                std::process::exit(1);
            }
        }
    }
    a
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mut opt = parse(&argv);

    // `if( fastathreshold < 0.0001 ) constraint = 0;`
    if opt.fastathreshold < 0.0001 {
        opt.constraint = 0;
    }
    assert!(
        opt.nthread == 0,
        "dvtditr-rs is single-threaded; the script passes -C 0"
    );
    if opt.specificityconsideration != 0.0 {
        // `--allowshift` is the strategy that reaches this: the script passes
        // `-s 0.8`, which switches `tditeration` onto `scoringmatrices`,
        // `classifypairs`, `smalldistmtx` and `partA__align_variousdist` --
        // a second scoring path this port does not implement. It is not on the
        // `--auto` route. See GOAL-C-FREE.md.
        eprintln!(
            "dvtditr-rs: -s {} (specificityconsideration, --allowshift) selects a \
             scoring path this port does not implement",
            opt.specificityconsideration
        );
        std::process::exit(1);
    }
    assert!(
        opt.use_fft,
        "the refinement loop refuses !use_fft in the C too"
    );

    let mut bytes = Vec::new();
    match &opt.inputfile {
        Some(f) => bytes = std::fs::read(f).unwrap_or_else(|e| {
            eprintln!("Cannot open {f}: {e}");
            std::process::exit(1);
        }),
        None => {
            std::io::stdin().read_to_end(&mut bytes).expect("read stdin");
        }
    }
    let input = io::read_data(bytes, opt.dorp, false, false);
    let njob = input.njob;
    let name = input.name;
    let seq_g = input.seq;

    // "Input pre-aligned data": every row must already be the same length.
    if seq_g.iter().any(|s| s.len() != seq_g[0].len()) {
        eprintln!("Input pre-aligned data");
        std::process::exit(1);
    }
    let nlenmax = seq_g[0].len();
    let alloclen = nlenmax * 9;

    // `dvtditr.c`: above 50000 residues an unconstrained run silently switches
    // to `alg = 'M'`, the memsave aligner -- a *different alignment*, not a
    // cheaper route to the same one (`rust/README.md`, "Notes for anyone
    // extending this", item 12). Refuse rather than keep running `'A'` and
    // quietly disagree with the C on the one input class where it matters.
    if nlenmax > 50000 && opt.constraint == 0 {
        eprintln!(
            "dvtditr-rs: nlenmax={nlenmax} > 50000 unconstrained, where the C switches \
             to the memsave aligner (alg='M'); that path is not implemented"
        );
        std::process::exit(1);
    }

    let cargs_for_constants = Args {
        nwildcard: opt.nwildcard,
        nblosum: opt.nblosum,
        ppenalty: opt.ppenalty,
        ppenalty_ex: opt.ppenalty_ex,
        poffset: opt.poffset,
        penalty_shift_factor: opt.penalty_shift_factor,
        fft_threshold: opt.fft_threshold,
        fft_win_size: opt.fft_win_size,
        fmodel: opt.fmodel,
        ..Default::default()
    };
    let cn = Constants::new_seq(&cargs_for_constants, input.dorp, &seq_g);

    // `initFiles()`: the C opens `pre` for writing here and `trace` alongside.
    // `trace` carries progress only, and nothing reads it back.
    std::fs::write("trace", format!("PID = {}\n", std::process::id())).expect("write trace");

    // `pre` is written **once, at the end**. `writePre` takes a `force` flag,
    // and with `USE_XCED` 0 -- which is how this is built -- its whole body is
    // `if( force ) { rewind( prep_g ); writeData_pointer( ... ); }`. Every call
    // from inside the refinement loop passes 0, so it does nothing; only
    // `dvtditr`'s own final call passes 1.
    //
    // That matters more than it looks. A port that wrote on every improvement
    // would leave a *longer* intermediate alignment's tail behind the final,
    // shorter one, because the C rewinds without truncating -- 324 stray bytes
    // on `test/sample`, found by exactly this comparison.
    let write_pre_final = |rows: &[Vec<u8>]| {
        let mut out = Vec::new();
        io::write_data(&mut out, &name, rows);
        std::fs::write("pre", out).expect("write pre");
    };

    // With two sequences the C does not stop early -- that branch is `#if 0`'d
    // out. It drops to equal weights and a single pass, and still builds the
    // tree and writes `infile.tree`.
    let mut weight = 4;
    let mut niter = opt.niter;
    if njob == 2 {
        weight = 0;
        niter = 1;
    }

    // --- the guide tree, from the distance matrix the previous stage wrote ---
    let hat2 = std::fs::read("hat2").unwrap_or_else(|_| {
        eprintln!("Make hat2.");
        std::process::exit(1);
    });
    let full = io::read_hat2(hat2, njob);
    // `upgma_treeout` reads the half layout the memsave C variant uses.
    let mut eff: Vec<Vec<f64>> = vec![vec![0.0; njob]; njob];
    for i in 0..njob {
        for j in i + 1..njob {
            eff[i][j - i] = full[i][j];
        }
    }
    let tree = upgma_treeout(njob, &mut eff, opt.sueff_global, Some(&name), false, false);
    let topol = expand_topol(njob, &tree);
    let mut len: Vec<[f64; 2]> = tree.len.clone();
    len.resize(njob, [0.0, 0.0]);

    // --- constraints ---
    let mut table = if opt.constraint != 0 {
        let hat3 = std::fs::read("hat3").unwrap_or_else(|_| {
            eprintln!("Cannot open hat3.");
            std::process::exit(1);
        });
        let mut t = new_table(njob);
        let kozoarivec = io::read_localhomtable2_half(hat3, njob, &mut t);
        assert!(
            !kozoarivec.iter().any(|&k| k),
            "structural (korh='k') constraints are not implemented"
        );
        Some(t)
    } else {
        None
    };

    // --- `order`: the leaves, in the order the last merge step joined them ---
    {
        let mut order = String::new();
        for side in 0..2 {
            for &m in &topol[njob - 2][side] {
                if m == -1 {
                    break;
                }
                order.push_str(&format!("{m}\n"));
            }
        }
        std::fs::write("order", order).expect("write order");
    }

    // `-K n` adds n sequences to an existing alignment, and then every branch
    // that would break the reference apart is left alone.
    assert_eq!(
        opt.nadd, 0,
        "--add (-K) needs the skipthisbranch computation, which is not implemented"
    );
    let skipthisbranch = vec![[false; 2]; njob];

    // The alignment is refined whole when it is constrained or segmentation is
    // off (`-l` turns it off, and every constrained call passes `-l`), and in
    // segments otherwise -- which is the FFT-NS-i case.
    let mut anchors: Vec<i32> = vec![0];
    if opt.constraint == 0 && opt.bunkatsu {
        anchors.extend(search_anchors(&cn, &seq_g, 20, 65));
    }
    anchors.push(seq_g[0].len() as i32);

    let params = IterParams {
        niter,
        cut: opt.cut,
        cooling: true,
        minimumweight: opt.minimumweight,
        constraint: opt.constraint,
        weight,
        score_check: 1,
        use_fft: opt.use_fft,
        utree: true,
        alg: b'A',
        nkozo: 0,
        outgap: true,
        legacygapcost: false,
        newgapstr: b'-',
        fastathreshold: opt.fastathreshold,
        alloclen,
    };

    let mut len_for_iter = len;
    // `res_g`: the refined segments, concatenated back into whole rows.
    let mut res_g: Vec<Vec<u8>> = vec![Vec::new(); njob];
    for iseg in 0..anchors.len() - 1 {
        let from = anchors[iseg] as usize;
        let to = anchors[iseg + 1] as usize;
        let mut aseq: Vec<Vec<u8>> = seq_g.iter().map(|s| s[from..to].to_vec()).collect();
        let mut bseq = aseq.clone();
        // The C measures this from the *output* built so far, not from the
        // input offsets, so the printed ranges overlap once a segment has been
        // refined shorter than it started.
        let pos = res_g[0].len();
        eprintln!(
            "Segment {:3}/{:3} {:4}-{:4}",
            iseg + 1,
            anchors.len() - 1,
            pos + 1,
            pos + 1 + (to - from)
        );
        tree_dependent_iteration(
            &cn,
            &params,
            njob,
            &mut aseq,
            &mut bseq,
            &topol,
            &mut len_for_iter,
            &skipthisbranch,
            table.as_mut(),
            // The C's `writePre( ..., 0 )` from inside the loop: a no-op here,
            // for the reason above.
            &mut |_| {},
        );
        for (i, row) in bseq.iter().enumerate() {
            res_g[i].extend_from_slice(row);
        }
    }

    eprintln!("done");
    write_pre_final(&res_g);
}
