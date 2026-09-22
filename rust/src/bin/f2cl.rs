//! `f2cl-rs` — a Rust port of MAFFT's `f2cl`, for the FASTA output path.
//!
//! `f2cl` is the last binary a `mafft` run executes: it reads the finished
//! alignment and writes it in the requested output format. On the `--auto`
//! path the script passes `-n -1 -f -l 60`, which is "FASTA, wrapped at 60
//! columns" -- that is the path ported here, in full.
//!
//! **The clustal and phylip formats are refused, not approximated.** They are
//! reachable only through `--clustalout` / `--phylipout`, which no `--auto`
//! decision selects, and each needs its own writer plus the conservation
//! marking in `setmark_clustal`. A run that asks for one gets a message
//! naming the C binary rather than a plausible-looking file, because a
//! silently different alignment file is the one outcome this port must never
//! produce. `--dash`'s `-d` is refused for the same reason.
#![forbid(unsafe_code)]

use std::io::{Read, Write};

use mafft_oxide::cargs::{self, Opt};
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
            'c' => eprintln!("comment = {}", o.value.as_deref().unwrap_or("(null)")),
            'r' => {
                out.orderfile.clone_from(&o.value);
                eprintln!("orderfile = {}", o.value.as_deref().unwrap_or("(null)"));
            }
            'n' => eprintln!("namelen = {}", cargs::myatoi(o.value.as_deref())),
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
            // `-E` / `-N` only select the alphabet the clustal conservation
            // marking uses, and that marking is computed and thrown away on
            // the FASTA path. Accepted and ignored, as the C effectively does.
            'E' | 'N' => {}
            _ => unreachable!("cargs only yields options this binary classified"),
        }
    }
    out
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
        'f' => {}
        'c' => refuse("clustal output (the default format)"),
        'y' => refuse("phylip output (-y)"),
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

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
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
