//! `replaceu-rs` — a bit-exact Rust port of MAFFT's `replaceu`.
//!
//! Runs first in every `mafft` invocation with `--anysymbol` / `--preservecase`
//! (which the GEMME/ESCOTT trace always passes). It does two things the rest of
//! the pipeline depends on: it maps every character outside a per-type "usual"
//! set to a single unknown symbol, so the aligners only ever see an alphabet
//! they have a matrix row for; and it rewrites each header as
//! `_os_<index>_oe_<name>`, which is the marker `restoreu` later reads to put
//! the original characters back. That marker is why `tbfast`'s staged input
//! carries names like `>_os_0_oe_P53_HUMAN/1-393`.
//!
//! The reader is the case-preserving one -- `charfilter`, which drops only
//! whitespace -- because preserving case is the whole point: the original
//! bytes have to survive to `restoreu`.
#![forbid(unsafe_code)]

use std::io::{Read, Write};

use mafft_rs::cargs::{self, Opt};
use mafft_rs::io;

/// Protein characters `replace_unusual` leaves alone, upper-cased.
const USUAL_PROTEIN: &[u8] = b"ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-.";

/// Nucleotide characters it leaves alone, lower-cased.
const USUAL_NUCLEOTIDE: &[u8] = b"ATGCUatgcuBDHKMNRSVWYXbdhkmnrsvwyx-";

/// `replace_unusual()`: fold every character to the usual set or to `unknown`.
///
/// Note the case change applies only to characters that *are* usual; an
/// unusual one becomes `unknown` verbatim, without being case-folded first.
/// The two branches are not interchangeable and the C's order is kept.
///
/// Args:
///     seq: The record bodies, edited in place.
///     usual: Characters to keep.
///     unknown: What everything else becomes.
///     upper: Upper-case the kept characters (protein) rather than lower-case
///         them (nucleotide).
fn replace_unusual(seq: &mut [Vec<u8>], usual: &[u8], unknown: u8, upper: bool) {
    for row in seq.iter_mut() {
        for ch in row.iter_mut() {
            if usual.contains(ch) {
                *ch = if upper {
                    ch.to_ascii_uppercase()
                } else {
                    ch.to_ascii_lowercase()
                };
            } else {
                *ch = unknown;
            }
        }
    }
}

/// What the option loop leaves behind.
struct Args {
    /// `-i`: read from this file rather than stdin.
    inputfile: Option<String>,
    /// `-o`: added to each record's index in the `_os_` marker, so a seed
    /// alignment's records keep distinct indices from the main input's.
    seedoffset: i64,
    /// `-D` / `-P`: force nucleotide or protein instead of inferring.
    dorp: Option<u8>,
}

/// Walk `arguments()` and report the options the C reports as it goes.
///
/// Returns:
///     The parsed options.
fn arguments(argv: &[String]) -> Args {
    let mut out = Args {
        inputfile: None,
        seedoffset: 0,
        dorp: None,
    };
    let opts = cargs::parse(argv, |c| match c {
        'o' | 'i' => Some(Opt::Value),
        'D' | 'P' => Some(Opt::Flag),
        _ => None,
    });
    for o in &opts {
        match o.c {
            'o' => {
                out.seedoffset = cargs::myatoi(o.value.as_deref());
                eprintln!("seedoffset = {}", out.seedoffset);
            }
            'i' => {
                out.inputfile.clone_from(&o.value);
                eprintln!("inputfile = {}", o.value.as_deref().unwrap_or("(null)"));
            }
            'D' => out.dorp = Some(b'd'),
            'P' => out.dorp = Some(b'p'),
            _ => unreachable!("cargs only yields options this binary classified"),
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

    // `text = true` is `getnumlen_casepreserve` + `readData_pointer_casepreserve`:
    // the same `charfilter` body, the same 1000000-site sampling cap, and
    // `countATGC` rather than `countATGCandN`.
    let mut input = io::read_data(bytes, args.dorp, false, true);

    // `getnumlen_casepreserve` tracks the shortest record alongside the
    // longest; `read_data` returns only the longest, and the minimum over the
    // records it read is the same number. The sentinel is the C's, and shows
    // through on an input with no records at all.
    let nlenmin = input
        .seq
        .iter()
        .map(std::vec::Vec::len)
        .min()
        .unwrap_or(99_999_999);
    eprintln!(
        "{} x {} - {} {}",
        input.njob, input.nlenmax, nlenmin, input.dorp as char
    );

    if input.dorp == b'p' {
        replace_unusual(&mut input.seq, USUAL_PROTEIN, b'X', true);
    } else {
        replace_unusual(&mut input.seq, USUAL_NUCLEOTIDE, b'n', false);
    }

    let stdout = std::io::stdout();
    let mut out = std::io::BufWriter::new(stdout.lock());
    for (i, (name, seq)) in input.name.iter().zip(&input.seq).enumerate() {
        // `name[i]+1` skips the `'='` the reader puts in front of every name.
        let _ = out.write_all(b">_os_");
        let _ = write!(out, "{}", i as i64 + args.seedoffset);
        let _ = out.write_all(b"_oe_");
        let _ = out.write_all(&name[1..]);
        let _ = out.write_all(b"\n");
        let _ = out.write_all(seq);
        let _ = out.write_all(b"\n");
    }
    let _ = out.flush();
}
