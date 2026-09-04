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
#![forbid(unsafe_code)]

use std::io::Read;

use mafft_rs::align::AlignState;
use mafft_rs::constants::Constants;
use mafft_rs::falign::FalignState;
use mafft_rs::lalign::LalignState;
use mafft_rs::localhom;
use mafft_rs::pairalign::{self, PairAlg};
use mafft_rs::params::{atof, atoi, NOTSPECIFIED};
use mafft_rs::treebase::TreebaseOpts;
use mafft_rs::{addfunctions, io, keeplength, tree, treebase};

/// The subset of `tbfast`'s options this port implements. Anything else is
/// refused rather than ignored.
#[derive(Debug, Clone)]
struct TbArgs {
    /// `-N`: suppress the C's automatic switch to the memsave aligner past
    /// 30000 columns (tbfast.c:1096). `-B` is accepted and ignored, as in the C.
    nevermemsave: bool,
    // pair phase
    pair_alg: Option<PairAlg>,
    pair_ppenalty: i32,
    pair_ppenalty_ex: i32,
    pair_ppenalty_op: i32,
    /// `-E` (penalty_EX); parsed and unused, see the note at the call site.
    #[allow(dead_code)]
    pair_ppenalty_ex2: i32,
    /// `-Z`: overwrite the pairwise score with `naivepairscore11`.
    usenaivescore: bool,
    /// `-a` (`--fmodel`): take residue frequencies from the input.
    fmodel: i32,
    /// `-:` (`--anysymbol`): score the nucleotide wildcards. Both parsers set
    /// it, and it reaches `constants()` in both phases.
    pair_nwildcard: bool,
    nwildcard: bool,
    /// `-L` (`--legacygapcost`) in the progressive stage.
    legacygapcost: bool,
    /// `-t` / `-^`: write the guide tree.
    treeout: i32,
    /// `-U` (`--treein`): read the guide tree from `_guidetree`.
    treein: bool,
    /// `-T`: stop before aligning.
    noalign: bool,
    /// `-y` in *tbfast's own* parser. Distinct from the `distout` the pairwise
    /// parser sets to 'h' unconditionally: that one is set later, so the exit
    /// below sees only this.
    distout: u8,
    pair_poffset: i32,
    pair_shift_factor: f64,
    pair_nblosum: i32,
    /// `-j` / `-m`: JTT or transmembrane instead of BLOSUM, at `pam_n` PAM.
    pair_scoremtx: i32,
    pair_pam_n: i32,
    pair_tm_or_jtt: i32,
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
    scoremtx: i32,
    pam_n: i32,
    tm_or_jtt: i32,
    nthread: i32,
    nthreadtb: i32,
    outgap: bool,
    outnumber: bool,
    /// `-I <n>` (`--add`): the last `nadd` sequences are appended to the
    /// alignment formed by the first `njob - nadd`.
    nadd: usize,
    /// `-K` clears `addprofile`, which is on by default.
    addprofile: bool,
    /// `-Y` (`--keeplength`) and `-Z` / `-z` (`--mapout`).
    keeplength: bool,
    mapout: i32,
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
            nevermemsave: false,
            pair_alg: None,
            pair_ppenalty: NOTSPECIFIED,
            pair_ppenalty_ex: NOTSPECIFIED,
            pair_ppenalty_op: NOTSPECIFIED,
            pair_ppenalty_ex2: NOTSPECIFIED,
            usenaivescore: false,
            fmodel: 0,
            pair_nwildcard: false,
            nwildcard: false,
            legacygapcost: false,
            treeout: 0,
            treein: false,
            noalign: false,
            distout: 0,
            pair_poffset: NOTSPECIFIED,
            pair_shift_factor: 1000.0,
            pair_nblosum: 62,
            pair_scoremtx: 1,
            pair_pam_n: NOTSPECIFIED,
            pair_tm_or_jtt: mafft_rs::params::JTT,
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
            scoremtx: 1,
            pam_n: NOTSPECIFIED,
            tm_or_jtt: mafft_rs::params::JTT,
            nthread: 1,
            nthreadtb: 1,
            outgap: true,
            outnumber: false,
            nadd: 0,
            addprofile: true,
            keeplength: false,
            mapout: 0,
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
            'j' => {
                a.pair_pam_n = atoi(&take().ok_or("-j needs a value")?);
                a.pair_scoremtx = 0;
                a.pair_tm_or_jtt = mafft_rs::params::JTT;
            }
            'm' => {
                a.pair_pam_n = atoi(&take().ok_or("-m needs a value")?);
                a.pair_scoremtx = 0;
                a.pair_tm_or_jtt = mafft_rs::params::TM;
            }
            'g' => a.pair_ppenalty_ex = (atof(&take().ok_or("-g")?) * 1000.0 - 0.5) as i32,
            'f' => a.pair_ppenalty = (atof(&take().ok_or("-f")?) * 1000.0 - 0.5) as i32,
            'Q' => a.pair_shift_factor = atof(&take().ok_or("-Q")?),
            'h' => a.pair_poffset = (atof(&take().ok_or("-h")?) * 1000.0 - 0.5) as i32,
            'O' => a.pair_ppenalty_op = (atof(&take().ok_or("-O")?) * 1000.0 - 0.5) as i32,
            // `-E` sets penalty_EX, which `genL__align11` does not read (its
            // `+= fpenalty_EX` lines are commented out); parsed so the value
            // is consumed rather than mistaken for a flag.
            'E' => a.pair_ppenalty_ex2 = (atof(&take().ok_or("-E")?) * 1000.0 - 0.5) as i32,
            'N' => a.pair_alg = Some(PairAlg::Generalized),
            'Z' => a.usenaivescore = true,
            'D' => a.dorp = Some(b'd'),
            'P' => a.dorp = Some(b'p'),
            // `--anysymbol`; only affects the nucleotide score matrix, so it
            // is a no-op for protein input.
            'a' => a.fmodel = 1,
            ':' => a.pair_nwildcard = true,
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
            'j' => {
                a.pam_n = atoi(&take().ok_or("-j")?);
                a.scoremtx = 0;
                a.tm_or_jtt = mafft_rs::params::JTT;
            }
            'm' => {
                a.pam_n = atoi(&take().ok_or("-m")?);
                a.scoremtx = 0;
                a.tm_or_jtt = mafft_rs::params::TM;
            }
            'f' => a.ppenalty = (atof(&take().ok_or("-f")?) * 1000.0 - 0.5) as i32,
            'g' => a.ppenalty_ex = (atof(&take().ok_or("-g")?) * 1000.0 - 0.5) as i32,
            'Q' => a.shift_factor = atof(&take().ok_or("-Q")?),
            'h' => a.poffset = (atof(&take().ok_or("-h")?) * 1000.0 - 0.5) as i32,
            'l' => a.fastathreshold = atof(&take().ok_or("-l")?),
            'X' => a.sueff_global = atof(&take().ok_or("-X")?),
            'F' => a.use_fft = true,
            'O' => a.outgap = false,
            'L' => a.legacygapcost = true,
            't' => a.treeout = 1,
            'U' => a.treein = true,
            '^' => a.treeout = 2,
            'T' => a.noalign = true,
            'y' => {
                let v = take().ok_or("-y needs a value")?;
                let c = v.as_bytes().first().copied().unwrap_or(0);
                if c != b'c' && c != b'h' {
                    eprintln!("Set -y c or -y h in v>=7.521.");
                    std::process::exit(1);
                }
                a.distout = c;
            }
            'n' => a.outnumber = true,
            'I' => a.nadd = atoi(&take().ok_or("-I needs a value")?) as usize,
            'K' => a.addprofile = false,
            'Y' => a.keeplength = true,
            'Z' => a.mapout = 1,
            'z' => a.mapout = 2,
            'D' => a.dorp = Some(b'd'),
            'P' => a.dorp = Some(b'p'),
            'a' => a.fmodel = 1,
            ':' => a.nwildcard = true,
            'N' => a.nevermemsave = true,
            'v' | 'B' => {}
            other => a.unsupported.push(format!("option -{}", other)),
        }
        Ok(())
    })?;
    Ok(a)
}

fn refuse(a: &TbArgs) -> Option<String> {
    let mut bad: Vec<String> = a.unsupported.clone();
    if a.pair_alg.is_none() {
        bad.push(
            "only the -L (local), -A (global) and -N (generalised) pairwise stages are ported"
                .into(),
        );
    }
    if a.fastathreshold < 0.0001 {
        bad.push("this port only covers the constrained (-l) progressive stage".into());
    }
    if a.keeplength && a.nadd == 0 {
        // The C reaches `restoreoriginalgaps( njob, bseq, NULL )` and
        // segfaults; there is no output to be identical to.
        bad.push("-Y (--keeplength) without -I <n> crashes the C too".into());
    }
    if a.treein && a.unalignlevel > 0.0 {
        // The C derives `expdist` from the user tree in that case.
        bad.push("-U with --allowshift is not ported".into());
    }
    if a.nthread > 0 || a.nthreadtb > 0 || a.pair_nthread > 0 {
        bad.push("threads change the result; run with -C 0".into());
    }
    // -b is ignored for nucleotides, whose matrix comes from Kimura's model.
    // -1 is --aamatrix and -2 is --text; both read `_aamtx`.
    let known = |n: i32| matches!(n, 30 | 45 | 50 | 62 | 80 | 90 | 100 | -1 | -2);
    if a.dorp != Some(b'd')
        && ((a.scoremtx == 1 && !known(a.nblosum))
            || (a.pair_scoremtx == 1 && !known(a.pair_nblosum)))
    {
        bad.push("only the built-in BLOSUM tables are ported".into());
    }
    if bad.is_empty() {
        None
    } else {
        Some(bad.join("\n  "))
    }
}

/// `includemember()`: every member of `mem` is in `cand`.
fn includemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.iter().all(|m| cand.contains(m))
}

/// `samemember()`: the same members, in any order (the live C definition).
fn samemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.len() == cand.len() && mem.iter().all(|m| cand.contains(m))
}

/// Build a `Constants` from one phase's penalty settings by reusing the
/// `disttbfast` argument struct, which `constants()` reads.
#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn constants_for(
    seq: &[Vec<u8>],
    dorp: u8,
    nblosum: i32,
    ppenalty: i32,
    ppenalty_ex: i32,
    ppenalty_dist: i32,
    poffset: i32,
    ppenalty_op: i32,
    shift_factor: f64,
    specificityconsideration: f64,
    fmodel: i32,
    legacygapcost: bool,
    scoremtx: i32,
    pam_n: i32,
    tm_or_jtt: i32,
    nwildcard: bool,
) -> Constants {
    let mut ca = mafft_rs::params::Args {
        nblosum,
        ppenalty,
        ppenalty_ex,
        ppenalty_dist,
        poffset,
        ppenalty_op,
        penalty_shift_factor: shift_factor,
        specificityconsideration,
        fmodel,
        legacygapcost,
        pam_n,
        tm_or_jtt,
        nwildcard,
        ..Default::default()
    };
    ca.scoremtx = scoremtx;
    Constants::new_seq(&ca, dorp, seq)
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
    let input = io::read_data(bytes, a.dorp, a.outnumber, a.nblosum == -2);
    let njob = input.njob;
    if njob < 2 {
        eprintln!("At least 2 sequences should be input!");
        std::process::exit(1);
    }
    // `-b` is only meaningful for protein: for nucleotides `constants()`
    // ignores it and derives the matrix from Kimura's model, so any value is
    // accepted there. Checked here rather than at parse time because the
    // sequence type is not known until the input has been read -- see
    // `Args::check_blosum`, which does the same for `disttbfast`. Both stages'
    // `-b` values are checked, since they are parsed separately.
    for (stage, nblosum, scoremtx) in [
        ("pairwise", a.pair_nblosum, a.pair_scoremtx),
        ("progressive", a.nblosum, a.scoremtx),
    ] {
        if input.dorp != b'd'
            && scoremtx == 1
            && !matches!(nblosum, 30 | 45 | 50 | 62 | 80 | 90 | 100 | -1 | -2)
        {
            eprintln!(
                "tbfast-rs: unsupported configuration:\n  unknown -b value in the {} stage",
                stage
            );
            std::process::exit(1);
        }
    }
    // --- pairwise phase ---
    let pair_cn = constants_for(
        &input.seq,
        input.dorp,
        a.pair_nblosum,
        a.pair_ppenalty,
        a.pair_ppenalty_ex,
        NOTSPECIFIED,
        a.pair_poffset,
        a.pair_ppenalty_op,
        a.pair_shift_factor,
        a.pair_unalignlevel,
        a.fmodel,
        false,
        a.pair_scoremtx,
        a.pair_pam_n,
        a.pair_tm_or_jtt,
        a.pair_nwildcard,
    );
    // `initFiles()` runs here in the C, before `seqcheck()`, so an aborted run
    // still leaves an empty `pre` behind rather than no file at all.
    std::fs::write("pre", b"").expect("Cannot open pre");

    for (i, s) in input.seq.iter().enumerate() {
        for (p, &c) in s.iter().enumerate() {
            if (c as usize) >= pair_cn.charsize || pair_cn.amino_n[c as usize] == -1 {
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
    let mut gst = mafft_rs::genalign::GenalignState::new();

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
        &mut gst,
        a.usenaivescore,
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
    // The C opens `hat3` unconditionally in the pairwise stage and only fills
    // it when constraints were stored, so an empty file is the normal result
    // of a `-+ 0` run.
    if a.iterate != 0 {
        write_hat3(njob, &table);
    } else {
        std::fs::write("hat3", b"").expect("Cannot open hat3.");
    }
    io::write_hat2(njob, &input.name, &iscore);

    // `pairlocalalign()` runs before this test in the C and sets `distout` to
    // 'h' unconditionally, so with `-T` and no `-t` the run stops here -- after
    // `hat2` and `hat3`, but before the tree and `order`.
    if a.treeout == 0 && a.noalign {
        let mut out = Vec::new();
        io::write_data(&mut out, &input.name, &input.seq);
        std::fs::write("pre", &out).expect("write pre");
        eprintln!();
        return;
    }

    // --- progressive phase (its own scoring parameters) ---
    let cn = constants_for(
        &input.seq,
        input.dorp,
        a.nblosum,
        a.ppenalty,
        a.ppenalty_ex,
        a.ppenalty_dist,
        a.poffset,
        mafft_rs::params::NOTSPECIFIED,
        a.shift_factor,
        a.unalignlevel,
        a.fmodel,
        a.legacygapcost,
        a.scoremtx,
        a.pam_n,
        a.tm_or_jtt,
        a.nwildcard,
    );

    let mut tr = if a.treein {
        match tree::load_tree(
            njob,
            if a.treeout != 0 {
                Some(&input.name[..])
            } else {
                None
            },
            a.outnumber,
        ) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("{}", e);
                std::process::exit(1);
            }
        }
    } else if a.treeout != 0 {
        tree::upgma_treeout(
            njob,
            &mut iscore,
            a.sueff_global,
            Some(&input.name),
            a.outnumber,
            a.treeout == 2,
        )
    } else {
        tree::upgma(njob, &mut iscore, a.sueff_global)
    };
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

    if a.treeout != 0 && a.noalign {
        let mut out = Vec::new();
        io::write_data(&mut out, &input.name, &input.seq);
        std::fs::write("pre", &out).expect("write pre");
        eprintln!();
        return;
    }
    let mut eff = tree::counteff_simple(njob, &mut tr);
    // `--add`: tbfast (and only tbfast) discounts the added sequences by two
    // orders of magnitude, so the reference alignment dominates every node
    // they join.
    for e in eff.iter_mut().skip(njob - a.nadd) {
        *e /= 100.0;
    }

    // `--add`: what to do at each tree node, and how `aseq` starts out. This
    // is tbfast's copy of the block, which differs from disttbfast's: there
    // the per-added-sequence `samemembern` pass is live, here both variants
    // of it sit inside one `#if 0 ... #elseif 0 ... #endif` -- `#elseif` is
    // not a directive, so it does not reopen the block -- and only the
    // `includemember` pass runs.
    let mut mergeoralign: Vec<u8> = vec![b'a'; njob.saturating_sub(1)];
    let mut originalgaps: Option<Vec<u8>> = None;
    let mut addbk: Option<Vec<Vec<u8>>> = None;
    let mut mapseq: Option<Vec<Vec<u8>>> = None;
    let mut aseq: Vec<Vec<u8>> = bseq.clone();
    if a.nadd != 0 {
        let nadd = a.nadd;
        let alignmentlength = input.seq[0].len();
        for i in 0..njob - nadd {
            if input.seq[i].len() != alignmentlength {
                eprintln!("# The original {:4} sequences must be aligned", njob - nadd);
                std::process::exit(1);
            }
        }
        if a.keeplength {
            originalgaps = Some(keeplength::record_original_gaps(njob - nadd, &input.seq));
            if a.mapout != 0 {
                addbk = Some(
                    input.seq[njob - nadd..]
                        .iter()
                        .map(|s| io::gappick0(s))
                        .collect(),
                );
            }
        }
        let addmem: Vec<i32> = ((njob - nadd)..njob).map(|i| i as i32).collect();
        if a.addprofile {
            let alignmentlength = input.seq[njob - nadd].len();
            for i in (njob - nadd)..njob {
                if input.seq[i].len() != alignmentlength {
                    eprintln!("# The {:4} additional sequences must be aligned", nadd);
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
                eprintln!(
                    "# There is no appropriate position to add the {} sequences in the guide tree.",
                    nadd
                );
                std::process::exit(1);
            }
            let mut tail: Vec<Vec<u8>> = input.seq[njob - nadd..].to_vec();
            addfunctions::commongappick(&mut tail);
            mapseq = Some(tail.clone());
            for (i, s) in tail.into_iter().enumerate() {
                aseq[njob - nadd + i] = s;
            }
        } else {
            for v in mergeoralign.iter_mut() {
                *v = b'n';
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
            // Already `gappick0`ed into bseq above.
        }
        let mut head: Vec<Vec<u8>> = input.seq[..njob - nadd].to_vec();
        addfunctions::commongappick(&mut head);
        for (i, s) in head.into_iter().enumerate() {
            aseq[i] = s;
        }
    }

    localhom::calcimportance_half(njob, &eff, &aseq, &mut table);

    let mut nlen: Vec<i32> = input.seq.iter().map(|s| s.len() as i32).collect();
    let mut fs = FalignState::new();
    let opts = TreebaseOpts {
        extanch: None,
        terminalmargin: 0,
        mergeoralign: if a.nadd != 0 {
            Some(&mergeoralign)
        } else {
            None
        },
        nadd: a.nadd,
        outgap: a.outgap,
        use_fft: a.use_fft,
        force_fft: false,
        nevermemsave: a.nevermemsave,
        alg: b'A',
        minimumweight: a.minimumweight,
        constraint: Some(&table),
        fastathreshold: a.fastathreshold,
        use_dynamic_mtx: true,
    };
    eprintln!("Progressive alignment ... ");
    treebase::treebase(
        &cn, &opts, njob, &mut aseq, &mut nlen, &tr, &eff, &mut ast, &mut fs,
    );

    // `--keeplength`: as in disttbfast, and in the same place -- before the
    // alignment is written, since the compact --mapout renames sequences.
    let mut names = input.name.clone();
    if a.keeplength {
        let nadd = a.nadd;
        let (_ndeleted, deletelist) =
            keeplength::delete_new_insertions_whole(njob - nadd, nadd, &mut aseq);
        let mut dlf = String::new();
        for (i, list) in deletelist.iter().enumerate() {
            for e in list.iter() {
                dlf.push_str(&format!("{} {} {}\n", njob - nadd + i, e.pos, e.len));
            }
        }
        std::fs::write("_deletelist", dlf).expect("Cannot open _deletelist");
        keeplength::restore_original_gaps(njob, &mut aseq, originalgaps.as_ref().unwrap());
        if a.mapout != 0 {
            let addbk = addbk.as_ref().unwrap();
            let map = if a.mapout == 1 {
                keeplength::reconstruct_delete_map(
                    addbk,
                    &deletelist,
                    &aseq[njob - nadd..],
                    &input.name[njob - nadd..],
                    a.outnumber,
                )
            } else {
                let mut tail: Vec<Vec<u8>> = names[njob - nadd..].to_vec();
                let realn: &[Vec<u8>] = match mapseq.as_ref() {
                    Some(v) => v,
                    None => &input.seq[njob - nadd..],
                };
                let m = keeplength::reconstruct_delete_map_compact(
                    addbk,
                    &deletelist,
                    realn,
                    &mut tail,
                    a.outnumber,
                );
                names[njob - nadd..].clone_from_slice(&tail);
                m
            };
            std::fs::write("_deletemap", map).expect("Cannot open _deletemap");
        }
    }

    // The C writes the alignment to the file `pre` (prep_g), not to stdout,
    // and leaves an empty `hat3` behind; the mafft wrapper reads `pre`.
    let mut out = Vec::new();
    io::write_data(&mut out, &names, &aseq);
    std::fs::write("pre", &out).expect("write pre");
    if a.iterate == 0 {
        // The C's pairwise stage creates an empty hat3 even when nothing
        // asked for one.
        let _ = std::fs::write("hat3", b"");
    }
}
