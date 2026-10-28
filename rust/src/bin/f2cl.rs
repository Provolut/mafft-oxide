//! `f2cl-rs` — a Rust port of MAFFT's `f2cl`, for the FASTA output path.
//!
//! `f2cl` is the last binary a `mafft` run executes: it reads the finished
//! alignment and writes it in the requested output format. On the `--auto`
//! path the script passes `-n -1 -f -l 60`, which is "FASTA, wrapped at 60
//! columns" -- that is the path ported here, in full.
//!
//! The clustal (`--clustalout`) and phylip (`--phylipout`) writers are here
//! too, with the conservation marking in `setmark_clustal` that clustal needs.
//! They were refused rather than approximated until they had been compared
//! against the C byte for byte, which is the only reason to ship a writer at
//! all: a plausible-looking but different alignment file is the one outcome
//! this port must never produce. `--dash`'s `-d` is still refused -- it
//! re-picks common gaps across a filtered record set, which is a change to the
//! alignment and not an output format.
#![forbid(unsafe_code)]

use std::io::{Read, Write};

use mafft_oxide::cargs::{self, Opt};
use mafft_oxide::constants::Constants;
use mafft_oxide::io;

/// What the option loop leaves behind.
struct Args {
    inputfile: Option<String>,
    orderfile: Option<String>,
    /// `-l`: FASTA line width. Negative means "one line per record", which is
    /// the C's `LineLengthInFASTA < 0` branch, not merely a large width.
    line_length: i64,
    /// `'f'`, `'c'` or `'y'`. The C defaults to clustal.
    format: char,
    /// `-d`: drop DASH records and re-pick common gaps.
    excludedashseq: bool,
    /// `-n`: the name column width. `-1` means "the format's default", which
    /// is 15 for clustal and 10 for phylip -- not the same number, so it
    /// cannot be resolved here.
    namelen: i64,
    /// `-c`: the strategy string the clustal header carries.
    comment: Option<String>,
    /// `-E` / `-N` (`f2clext`). `-E` also sets `nblosum = -2`, which is
    /// `--text`: it forces `dorp` to protein and widens the alphabet, so it
    /// changes `constants()` and not just the marking branch.
    extendedalphabet: bool,
    nblosum: i32,
}

/// Walk `arguments()`, echoing what the C echoes.
///
/// Returns:
///     The parsed options.
fn arguments(argv: &[String]) -> Args {
    let mut out = Args {
        inputfile: None,
        orderfile: None,
        line_length: -1,
        format: 'c',
        excludedashseq: false,
        namelen: -1,
        comment: None,
        extendedalphabet: false,
        nblosum: 62,
    };
    let opts = cargs::parse(argv, |c| match c {
        'i' | 'c' | 'r' | 'n' | 'l' => Some(Opt::Value),
        'f' | 'd' | 'y' | 'E' | 'N' => Some(Opt::Flag),
        _ => None,
    });
    for o in &opts {
        match o.c {
            'i' => {
                out.inputfile.clone_from(&o.value);
                eprintln!("inputfile = {}", o.value.as_deref().unwrap_or("(null)"));
            }
            'c' => {
                out.comment.clone_from(&o.value);
                eprintln!("comment = {}", o.value.as_deref().unwrap_or("(null)"));
            }
            'r' => {
                out.orderfile.clone_from(&o.value);
                eprintln!("orderfile = {}", o.value.as_deref().unwrap_or("(null)"));
            }
            'n' => {
                out.namelen = cargs::myatoi(o.value.as_deref());
                eprintln!("namelen = {}", out.namelen);
            }
            'l' => {
                out.line_length = cargs::myatoi(o.value.as_deref());
                if out.line_length == 0 {
                    eprintln!("LineLengthInFASTA = 0 ??");
                    std::process::exit(1);
                }
            }
            'f' => out.format = 'f',
            'y' => out.format = 'y',
            'd' => out.excludedashseq = true,
            // `-E` widens the alphabet and takes the marking's simple branch;
            // `-N` is the ordinary one. On the FASTA path the marking is
            // computed and thrown away, so neither showed before.
            'E' => {
                out.extendedalphabet = true;
                out.nblosum = -2;
            }
            'N' => out.extendedalphabet = false,
            _ => unreachable!("cargs only yields options this binary classified"),
        }
    }
    out
}

/// `mltaln.h`'s `VERSION`, as the clustal header prints it. Duplicated from
/// `src/bin/version.rs` rather than shared, because that file deliberately
/// keeps the string beside the code that prints it; the two are pinned
/// together by `f2cl_clustal_header_version_matches_version_rs`.
const VERSION: &str = "7.526";

/// `fillspace()` from `core/f2cl.c`: pad to `lenmax` with **spaces**, not gaps.
///
/// Both non-FASTA writers index every record at the same offset and stop at the
/// record's end, so a short record would silently shorten its own line. The C
/// pads first for that reason, and pads with a space rather than a `-` because
/// the padding is beyond the alignment, not part of it.
fn fillspace(seq: &mut Vec<u8>, lenmax: usize) {
    if seq.len() < lenmax {
        seq.resize(lenmax, b' ');
    }
}

/// `setmark_clustal()` from `core/f2cl.c`: the conservation line under each
/// clustal block -- `'*'` identical, `':'` strongly similar, `'.'` weakly so,
/// space otherwise.
///
/// The two group tables are Clustal's own and are transcribed verbatim; their
/// order matters only in that the first match wins, and `'*'` is tested before
/// either. A column containing a gap or a space in any record is left blank
/// without consulting them at all.
fn setmark_clustal(nlen: usize, nseq: usize, seq: &[Vec<u8>], cn: &Constants, extended: bool) -> Vec<u8> {
    let (strong, weaker, nalpha): (&[&[u8]], &[&[u8]], i32) = if cn.dorp == b'd' {
        // The nucleotide weaker list is `"AG"`, `"CT"` -- **two** entries, not
        // three. `core/f2cl.c:35-38` assigns `weaker[0..2]` including
        // `weaker[2] = "CU"` and then sets `nweaker = 2`, so the loop never
        // reaches `"CU"` and a `C`/`U` column is marked by nothing. Including
        // it, which is what the declarations look like they intend, changes the
        // conservation line on nucleotide input: it was 4 of 688 differential
        // cases, all of them `d_rna.fa` under `-c`, and all of them in the mark
        // row only. Transcribed as the C executes it, not as it reads.
        (&[b"TU"], &[b"AG", b"CT"], 10)
    } else {
        (
            &[
                b"STA", b"NEQK", b"NHQK", b"NDEQ", b"QHRK", b"MILV", b"MILF", b"HY", b"FYW",
            ],
            &[
                b"CSA", b"ATV", b"SAG", b"STNK", b"STPA", b"SGND", b"SNDEQK", b"NDEQHK",
                b"NEQHRK", b"FVLIM", b"HFY",
            ],
            20,
        )
    };

    let mut mark = vec![b' '; nlen];
    let aan = |c: u8| -> i32 {
        // `amino_n` is indexed by the raw byte in the C; anything past the
        // table is not in the alphabet.
        if (c as usize) < cn.amino_n.len() {
            cn.amino_n[c as usize]
        } else {
            -1
        }
    };
    for i in 0..nlen {
        // A gap or a space anywhere in the column: no mark, and the tables are
        // never consulted. `fillspace`'s padding lands here.
        if (0..nseq).any(|j| {
            let s = seq[j][i];
            s == b'-' || s == b' '
        }) {
            continue;
        }
        if extended {
            let first = seq[0][i];
            if aan(first) < 0 {
                continue;
            }
            if (0..nseq).all(|j| seq[j][i] == first) {
                mark[i] = b'*';
            }
            continue;
        }
        let first = seq[0][i].to_ascii_uppercase();
        if aan(first) >= nalpha || aan(first) < 0 {
            continue;
        }
        if (0..nseq).all(|j| seq[j][i].to_ascii_uppercase() == first) {
            mark[i] = b'*';
            continue;
        }
        if strong
            .iter()
            .any(|g| (0..nseq).all(|j| g.contains(&seq[j][i].to_ascii_uppercase())))
        {
            mark[i] = b':';
            continue;
        }
        if weaker
            .iter()
            .any(|g| (0..nseq).all(|j| g.contains(&seq[j][i].to_ascii_uppercase())))
        {
            mark[i] = b'.';
        }
    }
    mark
}

/// Refuse a format this port does not implement, naming what to run instead.
fn refuse(what: &str) -> ! {
    eprintln!("f2cl-rs: {what} is not implemented in this port.");
    eprintln!("f2cl-rs: the `--auto` path uses `-f` (FASTA) only; run the C `f2cl` for this.");
    std::process::exit(1);
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = arguments(&argv);

    if args.excludedashseq {
        refuse("-d (--dash record exclusion)");
    }
    match args.format {
        'f' | 'c' | 'y' => {}
        other => {
            eprintln!("Unknown format");
            let _ = other;
            std::process::exit(1);
        }
    }

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

    // Case-preserving, like every other stage after `replaceu`: this is the
    // finished alignment on its way out, and folding case here would undo
    // exactly what `restoreu` just put back.
    let input = io::read_data(bytes, None, false, true);

    let order: Vec<usize> = match &args.orderfile {
        Some(path) => {
            let Ok(text) = std::fs::read_to_string(path) else {
                eprintln!("Cannot open {path}");
                std::process::exit(1);
            };
            let mut ord: Vec<usize> = text
                .lines()
                .take(input.njob)
                .map(|l| cargs::myatoi(Some(l)).max(0) as usize)
                .collect();
            // The C reads exactly njob lines and leaves the rest of the array
            // at its calloc'd zero when the file is short.
            ord.resize(input.njob, 0);
            ord
        }
        None => (0..input.njob).collect(),
    };

    // `if( format == 'c' || format == 'y' ) for(...) fillspace( seq[i], nlenmax )`.
    // The FASTA writer indexes each record on its own length and must not see
    // the padding, so this is conditional in the C and conditional here.
    let mut input = input;
    if args.format == 'c' || args.format == 'y' {
        for s in input.seq.iter_mut() {
            fillspace(s, input.nlenmax);
        }
    }

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());

    if args.format == 'c' || args.format == 'y' {
        // `constants( njob, seq )`, for `amino_n` and `dorp`. The marking needs
        // both, and `-E` reaches them through `nblosum = -2`.
        let cargs_ = mafft_oxide::params::Args {
            nblosum: args.nblosum,
            ..Default::default()
        };
        let cn = Constants::new_seq(&cargs_, input.dorp, &input.seq);
        let mut buf = Vec::new();
        if args.format == 'c' {
            // The C calls `setmark_clustal` unconditionally and only clustal
            // uses the result; computing it only where it is used is the same
            // bytes out.
            let mark = setmark_clustal(
                input.nlenmax,
                input.njob,
                &input.seq,
                &cn,
                args.extendedalphabet,
            );
            io::clustalout_pointer(
                &mut buf,
                input.njob,
                input.nlenmax,
                &input.seq,
                &input.name,
                Some(&mark),
                args.comment.as_deref(),
                &order,
                args.namelen,
                VERSION,
            );
        } else {
            io::phylipout_pointer(
                &mut buf,
                input.njob,
                input.nlenmax,
                &input.seq,
                &input.name,
                &order,
                args.namelen,
            );
        }
        let _ = out.write_all(&buf);
        let _ = out.flush();
        return;
    }

    for &k in &order {
        let (name, seq) = (&input.name[k], &input.seq[k]);
        let _ = out.write_all(b">");
        let _ = out.write_all(&name[1..]);
        let _ = out.write_all(b"\n");
        if args.line_length < 0 {
            let _ = out.write_all(seq);
            let _ = out.write_all(b"\n");
        } else {
            // `for( j=0; j<nalen; j += LineLengthInFASTA )` -- an empty record
            // writes no sequence line at all, where the unwrapped branch above
            // writes a blank one. The difference is visible, so it is kept.
            let width = args.line_length as usize;
            for chunk in seq.chunks(width) {
                let _ = out.write_all(chunk);
                let _ = out.write_all(b"\n");
            }
        }
    }
    let _ = out.flush();
}
