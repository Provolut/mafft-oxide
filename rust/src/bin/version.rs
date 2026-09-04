//! `version-rs` — a bit-exact Rust port of MAFFT's `version`.
//!
//! The whole C program is `fprintf( stdout, VERSION )`: the version string with
//! no trailing newline, because it is read by the `mafft` shell script rather
//! than by a person. The missing newline is therefore load-bearing, not an
//! oversight to correct.
#![forbid(unsafe_code)]

use std::io::Write;

/// `mltaln.h`'s `VERSION`. Kept beside the code that prints it rather than in a
/// shared module: nothing else in the port reads a version, and a second copy
/// somewhere would be a second thing to forget when upstream moves.
const VERSION: &str = "7.526";

fn main() {
    print!("{VERSION}");
    let _ = std::io::stdout().flush();
}
