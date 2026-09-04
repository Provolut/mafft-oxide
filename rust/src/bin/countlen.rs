//! `countlen-rs` — a bit-exact Rust port of MAFFT's `countlen`.
//!
//! Measures an input without aligning it and prints one line the `mafft` shell
//! script parses to pick a strategy: record count, longest and shortest
//! ungapped lengths, inferred sequence type, and the N frequency. The script
//! reads the fields positionally, so the line's shape is a contract.
#![forbid(unsafe_code)]

use std::io::Read;

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
/// Transliterated including its failure mode: an unknown option prints its
/// complaint, zeroes the counter so the loop ends, and then trips the
/// `argc != 0` check into an exit -- so an unknown option is fatal, but only
/// after the message that names it.
///
/// Args:
///     argv: The arguments, program name excluded.
///
/// Returns:
///     The `-i` input file, when one was given.
fn arguments(args: &[String]) -> Option<String> {
    // `args` excludes the program name; the C's `argc`/`argv` include it, and
    // both the decrement and the pointer bump happen inside the loop condition,
    // so `argv` still advances on the iteration whose `'-'` test fails. That is
    // what makes a positional argument fatal rather than ignored.
    let mut argc = args.len() as i64 + 1;
    let mut p: i64 = -1;
    let mut inputfile: Option<String> = None;
    loop {
        argc -= 1;
        if argc <= 0 {
            break;
        }
        p += 1;
        let Some(arg) = args.get(p as usize) else {
            break;
        };
        if !arg.starts_with('-') {
            break;
        }
        let bytes = arg.as_bytes().to_vec();
        let mut k = 1usize;
        while k < bytes.len() {
            let c = bytes[k] as char;
            if c == 'i' {
                p += 1;
                inputfile = args.get(p as usize).cloned();
                argc -= 1;
                break;
            }
            // The C's `default:` ends in `break`, which leaves the *switch* --
            // so the scan carries on through the rest of the cluster and
            // complains once per character. `argc = 0` here is what the
            // condition's own decrement then turns into -1, tripping the check
            // below: an unknown option is fatal, but only after it is named.
            eprintln!("illegal option {c}");
            argc = 0;
            k += 1;
        }
    }
    if argc != 0 {
        eprintln!("options: Check source file !");
        std::process::exit(1);
    }
    inputfile
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
