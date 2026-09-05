//! `countlen-rs` — a bit-exact Rust port of MAFFT's `countlen`.
//!
//! Measures an input without aligning it and prints one line the `mafft` shell
//! script parses to pick a strategy: record count, longest and shortest
//! ungapped lengths, inferred sequence type, and the N frequency. The script
//! reads the fields positionally, so the line's shape is a contract.
#![forbid(unsafe_code)]

use std::io::Read;

use mafft_rs::cargs::{self, Opt};
use mafft_rs::io;

/// Format a double the way glibc's `%f` does.
///
/// Rust prints `NaN` and `inf` where C prints `nan` and `inf`, and the script
/// that reads this line does not care -- but a differential test against the C
/// does, and an input with no nucleotide-ish letters at all reaches the
/// `0.0 / 0.0` that produces one.
///
/// Args:
///     value: The number to render.
///
/// Returns:
///     The C's six-decimal rendering.
fn c_format_f64(value: f64) -> String {
    if value.is_nan() {
        return if value.is_sign_negative() {
            "-nan".to_string()
        } else {
            "nan".to_string()
        };
    }
    if value.is_infinite() {
        return if value.is_sign_negative() {
            "-inf".to_string()
        } else {
            "inf".to_string()
        };
    }
    format!("{value:.6}")
}

/// `arguments()`: the C's option loop, which accepts only `-i`.
///
/// The loop itself lives in `cargs` -- its three easy-to-miss behaviours are
/// documented there, and three binaries now depend on getting them right.
///
/// Args:
///     argv: The arguments, program name excluded.
///
/// Returns:
///     The `-i` input file, when one was given.
fn arguments(argv: &[String]) -> Option<String> {
    cargs::parse(argv, |c| match c {
        'i' => Some(Opt::Value),
        _ => None,
    })
    .into_iter()
    .find(|o| o.c == 'i')
    .and_then(|o| o.value)
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let inputfile = arguments(&argv);

    let bytes = match inputfile {
        Some(path) => match std::fs::read(&path) {
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

    // `dorp` is NOTSPECIFIED here in the C, so the type is always inferred.
    let counts = io::getnumlen_nogap_countn(bytes, None);
    println!(
        "{} x {} - {} {} nfreq={}",
        counts.njob,
        counts.nlenmax,
        counts.nlenmin,
        counts.dorp as char,
        c_format_f64(counts.nfreq)
    );
}
