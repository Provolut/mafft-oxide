//! Golden tests for the `dvtditr` binary, needing no `core/` build.
//!
//! The expected outputs in `tests/golden/dvtditr/` came from the real C
//! (`core/dvtditr` at MAFFT 7.526) and are the ground truth here, exactly as
//! `tests/matrix.rs` is for the aligners: this file is what lets the suite
//! check iterative refinement end to end on a machine with no C toolchain,
//! while `verify_dvtditr.sh` is the wide differential against a live C.
//!
//! Regenerate with `./verify_dvtditr.sh --regenerate ../binaries
//! target/release/dvtditr-rs` when the C changes or a case is added -- never to
//! accommodate a diff in the port, which would be the port changing its answer.
//!
//! **`dvtditr` is not a filter.** It reads the alignment from stdin, reads
//! `hat2` and `hat3` from the working directory the `mafft` script built, and
//! writes `pre`, `order` and `infile.tree` back into that directory. So each
//! case stages a directory of its own and is compared on all three outputs and
//! the exit status, which is what the shell harness compares too.
//!
//! The recorded set is a subset of that harness's 32 cases, chosen for breadth
//! per kilobyte rather than breadth alone -- the selection rule and its
//! reasoning are in `verify_dvtditr.sh` next to the code that applies it.

use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// Directory holding the staged inputs and the C's answers.
const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/dvtditr");

/// The three files the C writes that the harness compares.
const OUTPUTS: [&str; 3] = ["pre", "order", "infile.tree"];

/// A working directory of this case's own, since `dvtditr` reads and writes the
/// current directory and the cases run concurrently.
///
/// Args:
///     tag: The case name, which makes the path unique within one process.
///
/// Returns:
///     An empty directory that did not exist a moment ago.
fn tempdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("dvtditr-rs-test-{}-{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).expect("create the case's working directory");
    d
}

/// Copy a case's staged inputs into a fresh working directory.
///
/// Args:
///     case: The fixture directory, holding `in/`, `argv` and `stdin`.
///     tag: The case name, for the directory's path.
///
/// Returns:
///     The working directory, ready for the binary to run in.
fn stage(case: &Path, tag: &str) -> PathBuf {
    let d = tempdir(tag);
    for e in std::fs::read_dir(case.join("in")).expect("the case has an in/ directory") {
        let p = e.expect("read the staged input").path();
        let name = p.file_name().expect("staged inputs are files");
        std::fs::copy(&p, d.join(name)).expect("copy the staged input");
    }
    d
}

/// Run `dvtditr-rs` in a staged directory.
///
/// Args:
///     dir: The working directory to run in.
///     argv: The argument vector, already split.
///     stdin: The alignment to feed it.
///
/// Returns:
///     The exit status and stderr; the outputs are files, read by the caller.
fn run(dir: &Path, argv: &[String], stdin: &Path) -> (std::process::ExitStatus, String) {
    let f = std::fs::File::open(stdin).expect("the case has a stdin file");
    let out = Command::new(env!("CARGO_BIN_EXE_dvtditr-rs"))
        .args(argv)
        .current_dir(dir)
        .stdin(Stdio::from(f))
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .output()
        .expect("run dvtditr-rs");
    (out.status, String::from_utf8_lossy(&out.stderr).into_owned())
}

/// Read a case's argument vector.
///
/// Args:
///     case: The fixture directory.
///
/// Returns:
///     The vector as the `mafft` script passed it.
fn argv_of(case: &Path) -> Vec<String> {
    std::fs::read_to_string(case.join("argv"))
        .expect("the case has an argv file")
        .split_whitespace()
        .map(str::to_owned)
        .collect()
}

/// Every recorded case, in a stable order.
///
/// Returns:
///     The fixture directories.
fn cases() -> Vec<PathBuf> {
    let mut v: Vec<PathBuf> = std::fs::read_dir(GOLDEN)
        .unwrap_or_else(|e| panic!("no fixtures in {GOLDEN}: {e}"))
        .map(|e| e.expect("read the fixture directory").path())
        .filter(|p| p.is_dir())
        .collect();
    v.sort();
    v
}

#[test]
fn dvtditr_matches_the_c_on_every_fixture() {
    let cases = cases();
    // A fixture directory that emptied itself would otherwise make this test
    // pass while comparing nothing, which is the failure `verify_golden.sh`
    // guards against on the shell side.
    assert!(
        cases.len() >= 10,
        "expected at least 10 recorded cases, found {} -- has the store been \
         truncated? Regenerate with verify_dvtditr.sh --regenerate",
        cases.len()
    );

    for case in &cases {
        let tag = case
            .file_name()
            .and_then(|s| s.to_str())
            .expect("case directories are named");
        let label = std::fs::read_to_string(case.join("label")).unwrap_or_default();
        let label = label.trim();

        let d = stage(case, tag);
        let (status, stderr) = run(&d, &argv_of(case), &case.join("stdin"));

        let want: i32 = std::fs::read_to_string(case.join("exit"))
            .expect("the case has an exit file")
            .trim()
            .parse()
            .expect("the exit file holds a status");
        assert_eq!(
            status.code(),
            Some(want),
            "{label}: exit status\n{stderr}"
        );

        for name in OUTPUTS {
            let got = std::fs::read(d.join(name))
                .unwrap_or_else(|e| panic!("{label}: dvtditr-rs wrote no {name}: {e}"));
            let expected =
                std::fs::read(case.join("out").join(name)).expect("the case records this output");
            // Bytes, not text: these are alignments, and a difference of one
            // residue is the whole thing this comparison exists to catch.
            assert!(
                got == expected,
                "{label}: {name} differs from the C ({} bytes against {})",
                got.len(),
                expected.len()
            );
        }
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// Replace one flag's value in an argument vector.
///
/// Args:
///     argv: The vector to change.
///     flag: The flag whose value moves.
///     value: What it becomes.
///
/// Returns:
///     The changed vector.
fn with_flag(argv: &[String], flag: &str, value: &str) -> Vec<String> {
    let mut v = argv.to_vec();
    let i = v
        .iter()
        .position(|a| a == flag)
        .unwrap_or_else(|| panic!("the recorded argv should carry {flag}"));
    v[i + 1] = value.to_owned();
    v
}

/// The three configurations the port does not implement are refused rather than
/// approximated, on an argument vector that is otherwise a real one.
///
/// This pins the contract that matters -- an unsupported run fails instead of
/// quietly returning a different alignment. It deliberately does not pin *how*
/// each one fails: two of them abort where the third prints a named message,
/// and making that consistent should not have to change this test.
#[test]
fn unported_configurations_are_refused() {
    let case = cases().into_iter().next().expect("at least one case");
    let argv = argv_of(&case);
    let stdin = case.join("stdin");

    for (flag, value, what) in [
        ("-C", "4", "threads"),
        ("-s", "0.8", "the specificity scoring path (--allowshift)"),
        ("-K", "1", "--add"),
    ] {
        let d = stage(&case, &format!("refuse{flag}"));
        let (status, stderr) = run(&d, &with_flag(&argv, flag, value), &stdin);
        assert!(
            !status.success(),
            "{flag} {value} ({what}) should be refused, but the run succeeded\n{stderr}"
        );
        let _ = std::fs::remove_dir_all(&d);
    }
}

/// The one refusal that names what it rejected says so on stderr.
///
/// The other two abort instead, which is why only this one is asserted here.
#[test]
fn the_specificity_refusal_names_the_flag() {
    let case = cases().into_iter().next().expect("at least one case");
    let d = stage(&case, "refusename");
    let (status, stderr) = run(
        &d,
        &with_flag(&argv_of(&case), "-s", "0.8"),
        &case.join("stdin"),
    );
    assert!(!status.success(), "the run should have been refused");
    assert!(
        stderr.contains("--allowshift") && stderr.contains("does not implement"),
        "the refusal should name the flag and say it is not implemented, got:\n{stderr}"
    );
    let _ = std::fs::remove_dir_all(&d);
}
