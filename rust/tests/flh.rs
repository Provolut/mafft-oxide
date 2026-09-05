//! Golden tests for `falign_localhom`, needing no `core/` build.
//!
//! The expected outputs in `tests/golden/flh/` came from the real C
//! (`Falign_localhom` in MAFFT 7.526, via `gen_flhharness.c`) and are the
//! ground truth here, exactly as `tests/matrix.rs` is for the aligners: this
//! file is what lets the suite check the constrained block aligner on a
//! machine with no C toolchain, while `verify_flh.sh` is the wide differential
//! against a live C.
//!
//! Regenerate with `./verify_flh.sh --regenerate ../core` when the C changes or
//! a case is added -- never to accommodate a diff in the port, which would be
//! the port changing its answer.
#[path = "support/flhcase.rs"]
mod flhcase;

use std::path::PathBuf;

/// Directory holding the argument vectors and the expected outputs.
const GOLDEN: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/flh");

/// The C's answers for a run of calls made in one process, replayed here
/// through one state object.
///
/// Separate from the per-case fixtures because it pins something they cannot:
/// what each call leaves in the scratch buffers for the next one to read.
#[test]
fn falign_localhom_matches_the_c_across_a_call_sequence() {
    let dir = PathBuf::from(GOLDEN);
    let raw = std::fs::read_to_string(dir.join("sequence.cases")).unwrap();
    let expected = std::fs::read_to_string(dir.join("sequence.expected")).unwrap();
    let cases: Vec<Vec<String>> = raw
        .lines()
        .filter(|l| !l.is_empty())
        .map(|l| l.split('\t').map(String::from).collect())
        .collect();
    assert!(cases.len() >= 20, "only {} calls in the sequence", cases.len());
    assert_eq!(flhcase::run_sequence(&cases), expected);
}

#[test]
fn falign_localhom_matches_the_c_on_every_fixture() {
    let dir = PathBuf::from(GOLDEN);
    let mut cases: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap_or_else(|e| panic!("no fixtures in {}: {e}", dir.display()))
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "args"))
        .collect();
    cases.sort();
    // A fixture directory that has emptied itself would make this test pass
    // while checking nothing, which is the failure mode `verify_flh.sh`'s case
    // floor guards against on the other side.
    assert!(
        cases.len() >= 50,
        "only {} fixtures; expected the full corpus",
        cases.len()
    );

    for args_path in cases {
        let raw = std::fs::read_to_string(&args_path).unwrap();
        let args: Vec<String> = raw.trim_end_matches('\n').split('\t').map(String::from).collect();
        let expected = std::fs::read_to_string(args_path.with_extension("expected")).unwrap();
        let got = flhcase::run(&args);
        assert_eq!(
            got,
            expected,
            "{} differs from the C",
            args_path.file_name().unwrap().to_string_lossy()
        );
    }
}
