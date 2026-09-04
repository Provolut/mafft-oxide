//! Golden tests for the small utility binaries, needing no `core/` build.
//!
//! `countlen` and `version` are read by the `mafft` shell script, which parses
//! their output positionally, so their exact bytes are the contract. The
//! expected values in `tests/golden/utils/` came from the real C
//! (`core/countlen`, `core/version` at MAFFT 7.526) and are the ground truth
//! here, exactly as `tests/matrix.rs` does for the aligners: this file is what
//! lets the suite run on a machine with no C toolchain, while `verify_utils.sh`
//! is the wide differential against a live C for anyone extending coverage.
//!
//! Regenerate with `./verify_utils.sh --regenerate ../core` when the C changes
//! or a new case is added -- never to accommodate a diff in the port, which
//! would be the port changing its answer.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Directory holding the inputs and the expected outputs.
const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/utils");

/// Run one utility and return its stdout as bytes.
///
/// Args:
///     bin: Path to the binary under test.
///     args: Its argument vector.
///
/// Returns:
///     Raw stdout. Compared as bytes, because `version`'s contract is partly
///     that it does *not* end in a newline.
fn run(bin: &str, args: &[&str]) -> Vec<u8> {
    let out = Command::new(bin)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("could not run {bin}: {e}"));
    out.stdout
}

/// Read the `<input>\t<expected line>` table `countlen.expected` holds.
fn expected_countlen() -> Vec<(String, String)> {
    let path = PathBuf::from(GOLDEN).join("countlen.expected");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("missing {}: regenerate with verify_utils.sh", path.display()));
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (name, want) = l
                .split_once('\t')
                .unwrap_or_else(|| panic!("malformed fixture line: {l}"));
            (name.to_string(), want.to_string())
        })
        .collect()
}

#[test]
fn countlen_matches_the_c_on_every_fixture() {
    let bin = env!("CARGO_BIN_EXE_countlen-rs");
    let cases = expected_countlen();
    assert!(cases.len() >= 10, "fixture table looks truncated");
    for (name, want) in cases {
        let input = Path::new(GOLDEN).join("in").join(&name);
        let got = run(bin, &["-i", input.to_str().unwrap()]);
        let got = String::from_utf8_lossy(&got).trim_end_matches('\n').to_string();
        assert_eq!(got, want, "countlen disagrees with the C on {name}");
    }
}

#[test]
fn countlen_reads_stdin_when_no_file_is_named() {
    // The `mafft` script uses both forms -- `countlen < infile` at two call
    // sites and `countlen -i infile` at a third -- so both are contracts.
    let bin = env!("CARGO_BIN_EXE_countlen-rs");
    let input = Path::new(GOLDEN).join("in").join("dna.fa");
    let piped = Command::new(bin)
        .stdin(std::fs::File::open(&input).expect("fixture"))
        .output()
        .expect("run countlen-rs");
    let named = run(bin, &["-i", input.to_str().unwrap()]);
    assert_eq!(piped.stdout, named, "the -i and stdin paths disagree");
}

#[test]
fn an_unknown_option_names_it_and_exits_one() {
    // The C's `default:` breaks the switch rather than the scan, then leaves
    // `argc` at -1 so the check below the loop fires. Both halves matter: the
    // message identifies the option, and the exit status is what a caller sees.
    let bin = env!("CARGO_BIN_EXE_countlen-rs");
    let out = Command::new(bin)
        .arg("-Z")
        .stdin(std::process::Stdio::null())
        .output()
        .expect("run countlen-rs");
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("illegal option Z"), "stderr was: {err}");
    assert!(err.contains("options: Check source file !"), "stderr was: {err}");
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stdout.is_empty(), "a refused run must not print a count");
}

#[test]
fn a_positional_argument_is_refused_rather_than_treated_as_the_input() {
    // The C bumps `argv` inside the loop condition even on the iteration whose
    // `'-'` test fails, which leaves `argc` non-zero and makes this fatal. A
    // port that stopped cleanly here would silently read stdin instead.
    let bin = env!("CARGO_BIN_EXE_countlen-rs");
    let input = Path::new(GOLDEN).join("in").join("dna.fa");
    let out = Command::new(bin)
        .arg(input.to_str().unwrap())
        .stdin(std::process::Stdio::null())
        .output()
        .expect("run countlen-rs");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("options: Check source file !"));
}

#[test]
fn a_missing_input_file_is_reported_with_its_name() {
    let bin = env!("CARGO_BIN_EXE_countlen-rs");
    let out = Command::new(bin)
        .args(["-i", "/nonexistent/countlen/input"])
        .stdin(std::process::Stdio::null())
        .output()
        .expect("run countlen-rs");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("Cannot open /nonexistent/countlen/input"));
}

#[test]
fn version_prints_the_string_with_no_trailing_newline() {
    // The missing newline is the C's behaviour and the script depends on it.
    let bin = env!("CARGO_BIN_EXE_version-rs");
    let want = std::fs::read(Path::new(GOLDEN).join("version.expected")).expect("fixture");
    assert_eq!(run(bin, &[]), want);
    assert!(!want.ends_with(b"\n"), "the fixture itself must have no newline");
}
