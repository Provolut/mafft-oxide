//! `tbfast-rs` — a bit-exact Rust port of MAFFT's `tbfast` for the
//! consistency-based strategies L-INS-i and G-INS-i.
//!
//! The C `tbfast` takes two argument groups separated by bare `_` markers: the
//! first configures an internal all-pairs pairwise stage (what used to be a
//! separate `pairlocalalign` binary), the second the progressive stage. Both
//! groups set the same globals, so scoring parameters are parsed, used, and
//! then re-parsed — which is why this port keeps two independent `Constants`.
//!
//! Pipeline: read FASTA, align every pair (locally for `-L`, globally for
//! `-A`), turn each pairwise alignment's ungapped blocks into constraints,
//! build a UPGMA tree from the pairwise distances, then align progressively
//! with those constraints biasing the profile DP.

use std::io::Read;

use mafft_rs::align::AlignState;
use mafft_rs::constants::{self, Constants};
use mafft_rs::falign::FalignState;
use mafft_rs::lalign::LalignState;
use mafft_rs::localhom;
use mafft_rs::pairalign::{self, PairAlg};
use mafft_rs::params::{atof, atoi, NOTSPECIFIED};
use mafft_rs::treebase::TreebaseOpts;
use mafft_rs::{io, tree, treebase};

/// The subset of `tbfast`'s options this port implements. Anything else is
/// refused rather than ignored.
#[derive(Debug, Clone)]
struct TbArgs {
    // pair phase
    pair_alg: Option<PairAlg>,
    pair_ppenalty: i32,
    pair_ppenalty_ex: i32,
    pair_poffset: i32,
    pair_shift_factor: f64,
    pair_nblosum: i32,
    pair_nthread: i32,
    pair_unalignlevel: f64,
    pair_fastathreshold: f64,
    // tree phase
    minimumweight: f64,
    ppenalty_dist: i32,
    ppenalty: i32,
    ppenalty_ex: i32,
    poffset: i32,
    shift_factor: f64,
    nblosum: i32,
    nthread: i32,
    nthreadtb: i32,
    outgap: bool,
    outnumber: bool,
    fastathreshold: f64,
    sueff_global: f64,
    use_fft: bool,
    iterate: i32,
    dorp: Option<u8>,
    unalignlevel: f64,
    unsupported: Vec<String>,
}

impl Default for TbArgs {
    fn default() -> Self {
        TbArgs {
            pair_alg: None,
            pair_ppenalty: NOTSPECIFIED,
            pair_ppenalty_ex: NOTSPECIFIED,
            pair_poffset: NOTSPECIFIED,
            pair_shift_factor: 1000.0,
            pair_nblosum: 62,
            pair_nthread: 1,
            pair_unalignlevel: 0.0,
            pair_fastathreshold: 2.7,
            minimumweight: 0.0,
            ppenalty_dist: NOTSPECIFIED,
            ppenalty: NOTSPECIFIED,
            ppenalty_ex: NOTSPECIFIED,
            poffset: NOTSPECIFIED,
            shift_factor: 1000.0,
            nblosum: 62,
            nthread: 1,
            nthreadtb: 1,
            outgap: true,
            outnumber: false,
            fastathreshold: 0.0,
            sueff_global: 0.1,
            use_fft: false,
            iterate: 0,
            dorp: None,
            unalignlevel: 0.0,
            unsupported: Vec::new(),
        }
    }
}

/// Split `argv` at the bare `_` markers, as `arguments()` in tbfast.c does.
fn split_args(argv: &[String]) -> (Vec<String>, Vec<String>) {
    let mut pav = Vec::new();
    let mut tav = Vec::new();
    let mut i = 0usize;
    let mut in_pair = false;
    let mut seen_pair = false;
    while i < argv.len() {
        if argv[i] == "_" {
            if !seen_pair {
                in_pair = true;
                seen_pair = true;
            } else if in_pair {
                in_pair = false;
            }
            i += 1;
            continue;
        }
        if in_pair {
            pav.push(argv[i].clone());
        } else {
            tav.push(argv[i].clone());
        }
        i += 1;
    }
    (pav, tav)
}

/// Walk clustered single-letter options the way the C parsers do.
fn scan<F>(argv: &[String], mut handle: F) -> Result<(), String>
where
    F: FnMut(char, &mut dyn FnMut() -> Option<String>) -> Result<(), String>,
{
    let mut i = 0usize;
    'outer: while i < argv.len() && argv[i].starts_with('-') {
        let word: Vec<char> = argv[i].chars().collect();
        let mut k = 1usize;
        while k < word.len() {
            let c = word[k];
            k += 1;
            let mut consumed = false;
            {
                let mut take = || -> Option<String> {
                    consumed = true;
                    argv.get(i + 1).cloned()
                };
                handle(c, &mut take)?;
            }
            if consumed {
                i += 2;
                continue 'outer;
            }
        }
        i += 1;
    }
    Ok(())
}

fn parse(argv: &[String]) -> Result<TbArgs, String> {
    let (pav, tav) = split_args(argv);
    let mut a = TbArgs::default();

    // --- pairwise phase (pairlocalalign.c's arguments()) ---
    scan(&pav, |c, take| {
        match c {
            'L' => a.pair_alg = Some(PairAlg::Local),
            'A' => a.pair_alg = Some(PairAlg::Global),
            'u' => a.pair_unalignlevel = atof(&take().ok_or("-u needs a value")?),
            'l' => a.pair_fastathreshold = atof(&take().ok_or("-l needs a value")?),
            'C' => {
                let v = take().ok_or("-C needs a value")?;
                a.pair_nthread = atoi(&v);
            }
            'b' => a.pair_nblosum = atoi(&take().ok_or("-b needs a value")?),
            'g' => a.pair_ppenalty_ex = (atof(&take().ok_or("-g")?) * 1000.0 - 0.5) as i32,
            'f' => a.pair_ppenalty = (atof(&take().ok_or("-f")?) * 1000.0 - 0.5) as i32,
            'Q' => a.pair_shift_factor = atof(&take().ok_or("-Q")?),
            'h' => a.pair_poffset = (atof(&take().ok_or("-h")?) * 1000.0 - 0.5) as i32,
            'D' => a.dorp = Some(b'd'),
            'P' => a.dorp = Some(b'p'),
            // `--anysymbol`; only affects the nucleotide score matrix, so it
            // is a no-op for protein input.
            ':' => {}
            other => a.unsupported.push(format!("pair option -{}", other)),
        }
        Ok(())
    })?;

    // --- progressive phase (tbfast.c's arguments()) ---
    scan(&tav, |c, take| {
        match c {
            '+' => a.iterate = atoi(&take().ok_or("-+ needs a value")?),
            'W' => a.minimumweight = atof(&take().ok_or("-W needs a value")?),
            'V' => a.ppenalty_dist = (atof(&take().ok_or("-V")?) * 1000.0 - 0.5) as i32,
            's' => a.unalignlevel = atof(&take().ok_or("-s")?),
            'C' => {
                let v = take().ok_or("-C needs a value")?;
                a.nthread = atoi(&v);
                a.nthreadtb = match v.find('-') {
                    Some(p) => atoi(&v[p + 1..]),
                    None => a.nthread,
                };
            }
            'b' => a.nblosum = atoi(&take().ok_or("-b")?),
            'f' => a.ppenalty = (atof(&take().ok_or("-f")?) * 1000.0 - 0.5) as i32,
            'g' => a.ppenalty_ex = (atof(&take().ok_or("-g")?) * 1000.0 - 0.5) as i32,
            'Q' => a.shift_factor = atof(&take().ok_or("-Q")?),
            'h' => a.poffset = (atof(&take().ok_or("-h")?) * 1000.0 - 0.5) as i32,
            'l' => a.fastathreshold = atof(&take().ok_or("-l")?),
            'X' => a.sueff_global = atof(&take().ok_or("-X")?),
            'F' => a.use_fft = true,
            'O' => a.outgap = false,
            'n' => a.outnumber = true,
            'D' => a.dorp = Some(b'd'),
            'P' => a.dorp = Some(b'p'),
            'v' | 'B' | 'N' | ':' => {}
            other => a.unsupported.push(format!("option -{}", other)),
        }
        Ok(())
    })?;
    Ok(a)
}

fn refuse(a: &TbArgs) -> Option<String> {
    let mut bad: Vec<String> = a.unsupported.clone();
    if a.pair_alg.is_none() {
        bad.push("only the -L (local) and -A (global) pairwise stages are ported".into());
    }
    if a.fastathreshold < 0.0001 {
        bad.push("this port only covers the constrained (-l) progressive stage".into());
    }
    if a.nthread > 0 || a.nthreadtb > 0 || a.pair_nthread > 0 {
        bad.push("threads change the result; run with -C 0".into());
    }
    if a.nblosum != 62 || a.pair_nblosum != 62 {
        bad.push("only the BLOSUM62 protein matrix is ported".into());
    }
    if a.dorp == Some(b'd') {
        bad.push("nucleotide input is not ported".into());
    }
    if bad.is_empty() {
        None
    } else {
        Some(bad.join("\n  "))
    }
}

/// Build a `Constants` from one phase's penalty settings by reusing the
/// `disttbfast` argument struct, which `constants()` reads.
fn constants_for(
    nblosum: i32,
    ppenalty: i32,
    ppenalty_ex: i32,
    ppenalty_dist: i32,
    poffset: i32,
    shift_factor: f64,
    specificityconsideration: f64,
) -> Constants {
    let mut ca = mafft_rs::params::Args {
        nblosum,
        ppenalty,
        ppenalty_ex,
        ppenalty_dist,
        poffset,
        penalty_shift_factor: shift_factor,
        specificityconsideration,
        ..Default::default()
    };
    ca.scoremtx = 1;
    Constants::new(&ca)
}

/// `hat3`: one line per constraint block, read back by `dvtditr`.
fn write_hat3(njob: usize, table: &localhom::LocalHomTable) {
    use std::fmt::Write as _;
    let mut s = String::new();
    for i in 0..njob - 1 {
        for j in i..njob {
            if j == i {
                continue;
            }
            for e in table[i][j - i].iter() {
                if e.opt == -1.0 {
                    continue;
                }
                let _ = writeln!(
                    s,
                    "{} {} {} {:7.5} {} {} {} {} {}",
                    i,
                    j,
                    e.overlapaa,
                    e.opt / 600.0 * 5.8,
                    e.start1,
                    e.end1,
                    e.start2,
                    e.end2,
                    e.korh as char
                );
            }
        }
    }
    std::fs::write("hat3", s).expect("write hat3");
}

/// `WriteFloatHat2_pointer_halfmtx()`: the distance matrix, for `dvtditr`.
fn write_hat2(njob: usize, name: &[Vec<u8>], mtx: &[Vec<f64>]) {
    use std::fmt::Write as _;
    let mut max = 0.0f64;
    for i in 0..njob - 1 {
        for j in 1..njob - i {
            if mtx[i][j] > max {
                max = mtx[i][j];
            }
        }
    }
    let mut s = String::new();
    let _ = writeln!(s, "{:5}", 1);
    let _ = writeln!(s, "{:5}", njob);
    let _ = writeln!(s, " {:6.3}", max * 2.5);
    for (i, nm) in name.iter().enumerate().take(njob) {
        let _ = writeln!(s, "{:4}. {}", i + 1, String::from_utf8_lossy(nm));
    }
    for i in 0..njob {
        for j in i + 1..njob {
            let _ = write!(s, "{:6.3}", mtx[i][j - i]);
            let ijsa = j - i;
            if ijsa % 12 == 0 || ijsa == njob - i - 1 {
                s.push('\n');
            }
        }
    }
    std::fs::write("hat2", s).expect("write hat2");
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let a = match parse(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    if let Some(msg) = refuse(&a) {
        eprintln!("tbfast-rs: unsupported configuration:\n  {}", msg);
        std::process::exit(1);
    }

    let mut bytes = Vec::new();
    std::io::stdin()
        .read_to_end(&mut bytes)
        .expect("read stdin");
    let input = io::read_data(bytes, a.dorp, a.outnumber);
    let njob = input.njob;
    if njob < 2 {
        eprintln!("At least 2 sequences should be input!");
        std::process::exit(1);
    }
    if input.dorp != b'p' {
        eprintln!("tbfast-rs: only amino-acid input is ported");
        std::process::exit(1);
    }

    // --- pairwise phase ---
    let pair_cn = constants_for(
        a.pair_nblosum,
        a.pair_ppenalty,
        a.pair_ppenalty_ex,
        NOTSPECIFIED,
        a.pair_poffset,
        a.pair_shift_factor,
        a.pair_unalignlevel,
    );
    // `initFiles()` runs here in the C, before `seqcheck()`, so an aborted run
    // still leaves an empty `pre` behind rather than no file at all.
    std::fs::write("pre", b"").expect("Cannot open pre");

    for (i, s) in input.seq.iter().enumerate() {
        for (p, &c) in s.iter().enumerate() {
            if (c as usize) >= constants::CHARSIZE || pair_cn.amino_n[c as usize] == -1 {
                eprintln!("Illegal character {}", c as char);
                eprintln!("=== Please check site {} in sequence {}.", p + 1, i + 1);
                std::process::exit(1);
            }
        }
    }

    let bseq: Vec<Vec<u8>> = input.seq.iter().map(|s| io::gappick0(s)).collect();
    let mut iscore: Vec<Vec<f64>> = (0..njob).map(|i| vec![0.0; njob - i]).collect();
    let mut table = localhom::new_table(njob);
    let mut lst = LalignState::new();
    let mut ast = AlignState::new();

    eprintln!("Pairwise alignment ... ");
    pairalign::pairalign(
        &pair_cn,
        a.pair_alg.unwrap(),
        njob,
        &bseq,
        &mut iscore,
        &mut table,
        true,
        a.outgap,
        &mut lst,
        &mut ast,
    );

    // tbfast rescales `opt` back out of the 5.8/600 units pairalign stores it in.
    // `-+ n` (outputhat23) then asks for the constraints and distances to be
    // written out for the iterative-refinement stage, which is a separate
    // binary (`dvtditr`) that consumes them.
    for row in table.iter_mut() {
        for chain in row.iter_mut() {
            for e in chain.iter_mut() {
                if e.opt == -1.0 {
                    continue;
                }
                e.opt = e.opt / 5.8 * 600.0;
            }
        }
    }

    // `hat2` is written whenever the pairwise stage ran: its argument parser
    // sets `distout = 'h'` unconditionally and tbfast's never resets it, so
    // the `else if( distout )` branch always fires. `hat3` only when asked.
    if a.iterate != 0 {
        write_hat3(njob, &table);
    }
    write_hat2(njob, &input.name, &iscore);

    // --- progressive phase (its own scoring parameters) ---
    let cn = constants_for(
        a.nblosum,
        a.ppenalty,
        a.ppenalty_ex,
        a.ppenalty_dist,
        a.poffset,
        a.shift_factor,
        a.unalignlevel,
    );

    let mut tr = tree::upgma(njob, &mut iscore, a.sueff_global);
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
    let eff = tree::counteff_simple(njob, &mut tr);

    let mut aseq: Vec<Vec<u8>> = bseq.clone();
    localhom::calcimportance_half(njob, &eff, &aseq, &mut table);

    let mut nlen: Vec<i32> = input.seq.iter().map(|s| s.len() as i32).collect();
    let mut fs = FalignState::new();
    let opts = TreebaseOpts {
        outgap: a.outgap,
        use_fft: a.use_fft,
        force_fft: false,
        minimumweight: a.minimumweight,
        constraint: Some(&table),
        fastathreshold: a.fastathreshold,
        use_dynamic_mtx: true,
    };
    eprintln!("Progressive alignment ... ");
    treebase::treebase(
        &cn, &opts, njob, &mut aseq, &mut nlen, &tr, &eff, &mut ast, &mut fs,
    );

    // The C writes the alignment to the file `pre` (prep_g), not to stdout,
    // and leaves an empty `hat3` behind; the mafft wrapper reads `pre`.
    let mut out = Vec::new();
    io::write_data(&mut out, &input.name, &aseq);
    std::fs::write("pre", &out).expect("write pre");
    if a.iterate == 0 {
        // The C's pairwise stage creates an empty hat3 even when nothing
        // asked for one.
        let _ = std::fs::write("hat3", b"");
    }
}
