//! Bit-identity matrix: replay every case in `tests/golden/matrix` and require
//! the port's bytes to equal the C's.
//!
//! `MANIFEST` is a tab-separated list of
//! `<key>\t<binary>\t<argv>\t<extra outputs>`, written by `gen_golden.py` from
//! a real `core/` build. `<key>` is `<input>.<tag>`; the expected stdout lives
//! in `<key>.fa` (`disttbfast`) and the expected side files in `<key>.<name>`.
//!
//! This is the differential harness's coverage without the differential
//! harness's requirement: `verify.sh` needs a compiled `core/` on the machine
//! running it, and these fixtures carry the C's answer instead. The inputs are
//! tiny on purpose -- the whole matrix replays in a couple of seconds, so it
//! belongs in `cargo test` rather than in a nightly job.
//!
//! Regenerate with `./gen/gen_golden.py ../core` when a feature lands. A fixture
//! that changes is the port changing its answer, which is the thing this
//! repository exists to prevent.

use std::path::{Path, PathBuf};
use std::process::Command;

const MATRIX: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/golden/matrix");

struct Case {
    key: String,
    binary: String,
    argv: Vec<String>,
    /// Side files the C wrote and the port must reproduce byte for byte.
    extras: Vec<String>,
    /// Files the case needs in the working directory before it runs.
    staged: Vec<String>,
    /// `ok` -- the C succeeded and the bytes must match.
    /// `exit` -- the C exited non-zero, so the port has to refuse too. That is
    /// where the refusal paths and the upstream crashes in issue #4 get their
    /// coverage.
    /// `refuse` -- the C succeeded but the port declines on purpose (threads).
    /// Declared in `gen_golden.py`, so implementing one of these later forces
    /// the declaration to be updated rather than silently passing.
    expect: String,
}

fn manifest() -> Vec<Case> {
    let text = std::fs::read_to_string(Path::new(MATRIX).join("MANIFEST"))
        .expect("tests/golden/matrix/MANIFEST -- regenerate with ./gen/gen_golden.py ../core");
    text.lines()
        .filter(|l| !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert!(f.len() >= 3, "malformed manifest line: {}", l);
            let list = |i: usize| -> Vec<String> {
                f.get(i)
                    .map(|e| {
                        e.split(',')
                            .filter(|s| !s.is_empty())
                            .map(String::from)
                            .collect()
                    })
                    .unwrap_or_default()
            };
            Case {
                key: f[0].to_string(),
                binary: f[1].to_string(),
                argv: f[2].split_whitespace().map(String::from).collect(),
                extras: list(3),
                staged: list(4),
                expect: f.get(5).unwrap_or(&"ok").to_string(),
            }
        })
        .collect()
}

fn workdir(key: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("mafft-oxide-matrix-{}-{}", key, std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// Compare as bytes, but report as text -- these are all FASTA or tables, so a
/// textual diff is what a failure needs to show.
fn same(got: &[u8], want: &[u8], what: &str, key: &str) -> Option<String> {
    if got == want {
        return None;
    }
    let g = String::from_utf8_lossy(got);
    let w = String::from_utf8_lossy(want);
    let mut line = 0;
    for (a, b) in g.lines().zip(w.lines()) {
        line += 1;
        if a != b {
            return Some(format!(
                "{} [{}] differs at line {}\n  port: {}\n  C   : {}",
                key,
                what,
                line,
                a.chars().take(90).collect::<String>(),
                b.chars().take(90).collect::<String>()
            ));
        }
    }
    Some(format!(
        "{} [{}] differs in length: port {} lines / {} bytes, C {} lines / {} bytes",
        key,
        what,
        g.lines().count(),
        got.len(),
        w.lines().count(),
        want.len()
    ))
}

fn replay(c: &Case) -> Vec<String> {
    let input = Path::new(MATRIX).join(format!("{}.fa", c.key.split('.').next().unwrap()));
    let bin = match c.binary.as_str() {
        "disttbfast" => env!("CARGO_BIN_EXE_disttbfast-rs"),
        "tbfast" => env!("CARGO_BIN_EXE_tbfast-rs"),
        other => panic!("unknown binary {}", other),
    };
    let dir = workdir(&c.key);
    for name in &c.staged {
        let from = Path::new(MATRIX).join(format!("{}.stage.{}", c.key, name));
        std::fs::copy(&from, dir.join(name))
            .unwrap_or_else(|e| panic!("staging {:?}: {}", from, e));
    }
    let f = std::fs::File::open(&input).unwrap_or_else(|e| panic!("{:?}: {}", input, e));
    let out = Command::new(bin)
        .args(&c.argv)
        .current_dir(&dir)
        .stdin(std::process::Stdio::from(f))
        .stderr(std::process::Stdio::null())
        .output()
        .expect("run port");

    let mut bad = Vec::new();
    if c.expect == "exit" || c.expect == "refuse" {
        // Either the C failed, so there is nothing to match, or the port
        // declines on purpose. Both mean it must not produce an answer.
        if out.status.success() {
            bad.push(match c.expect.as_str() {
                "exit" => format!("{}: port succeeded where the C exited non-zero", c.key),
                _ => format!(
                    "{}: port produced output, but it is declared as a deliberate \
                     refusal -- if that is no longer true, update gen_golden.py",
                    c.key
                ),
            });
        }
        let _ = std::fs::remove_dir_all(&dir);
        return bad;
    }
    if !out.status.success() {
        bad.push(format!(
            "{}: port exited {:?} where the C succeeded",
            c.key, out.status
        ));
        let _ = std::fs::remove_dir_all(&dir);
        return bad;
    }

    // `disttbfast` writes the alignment to stdout; `tbfast` writes it to `pre`.
    if c.binary == "disttbfast" {
        let want = std::fs::read(Path::new(MATRIX).join(format!("{}.fa", c.key)))
            .unwrap_or_else(|e| panic!("{}.fa: {}", c.key, e));
        if let Some(m) = same(&out.stdout, &want, "stdout", &c.key) {
            bad.push(m);
        }
    }
    for e in &c.extras {
        let want = std::fs::read(Path::new(MATRIX).join(format!("{}.{}", c.key, e)))
            .unwrap_or_else(|err| panic!("{}.{}: {}", c.key, e, err));
        match std::fs::read(dir.join(e)) {
            Ok(got) => {
                if let Some(m) = same(&got, &want, e, &c.key) {
                    bad.push(m);
                }
            }
            Err(err) => bad.push(format!("{}: port wrote no {} ({})", c.key, e, err)),
        }
    }
    let _ = std::fs::remove_dir_all(&dir);
    bad
}

/// One test over the whole matrix rather than one per case: a single run of
/// each binary is milliseconds, and reporting every mismatch at once is more
/// useful than stopping at the first.
#[test]
fn every_case_is_byte_identical_to_the_c() {
    let cases = manifest();
    assert!(
        cases.len() > 100,
        "matrix has only {} cases -- did gen_golden.py run?",
        cases.len()
    );
    let mut bad = Vec::new();
    for c in &cases {
        bad.extend(replay(c));
    }
    if !bad.is_empty() {
        panic!(
            "{} of {} cases differ from the C:\n{}",
            bad.len(),
            cases.len(),
            bad.join("\n")
        );
    }
}

/// The comparison has to actually compare. Verified once by hand -- flipping
/// a single residue in `aa.fftns2.fa` made the run report
/// "1 of 553 cases differ ... differs at line 4" with both versions printed --
/// and pinned here so a refactor of `same()` cannot quietly turn the matrix
/// into an equal number of assertions that nothing is wrong. (It was 553 cases
/// when that was done by hand and is 567 now; the count is deliberately not
/// asserted anywhere, because a feature landing should not have to touch a
/// number in a comment to stay honest. `cases.len() > 100` above is the
/// non-vacuity guard, and the feature list below is what keeps coverage from
/// shrinking.)
#[test]
fn the_comparison_detects_a_one_byte_difference() {
    let a = b">s0\nMNGTEG\n>s1\nMNGTEG\n";
    let b = b">s0\nMNGTEG\n>s1\nMNGTEW\n";
    assert!(same(a, a, "stdout", "k").is_none(), "identical must pass");
    let m = same(a, b, "stdout", "k").expect("a differing byte must be reported");
    assert!(m.contains("line 4"), "should name the line: {}", m);
    assert!(m.contains("MNGTEG") && m.contains("MNGTEW"), "should show both: {}", m);

    // A truncated output is a difference too, and the byte counts say so.
    let m = same(&a[..10], a, "stdout", "k").expect("a short output must be reported");
    assert!(m.contains("differs in length"), "{}", m);

    // An empty output against a real one -- what a silently crashing port
    // would produce -- must not read as a pass.
    assert!(same(b"", a, "stdout", "k").is_some());
}

/// The matrix is only worth what it covers. If a feature is added to
/// `gen_golden.py` but the fixtures are not regenerated, or vice versa, this
/// catches it before the comparison silently tests less than it looks like.
#[test]
fn the_matrix_covers_the_features_it_claims_to() {
    let cases = manifest();
    let keys: Vec<&str> = cases.iter().map(|c| c.key.as_str()).collect();
    for feature in [
        "nwns1",
        "nwns2",
        "fftns1",
        "fftns2",
        "forcefft",
        "memsave",
        "memsavefft",
        "aalign",
        "aalignfft",
        "legacygap",
        "distout",
        "blosum30",
        "blosum80",
        "jtt",
        "fmodel",
        "dynmtx",
        "anysymbol",
        "allowshift",
        "unweighted",
        "linsi",
        "ginsi",
        "einsi",
        // Threads. These were `refuse` fixtures until the port honoured `-C`;
        // they are byte-identity cases now, and naming all three here means a
        // regeneration that quietly dropped one cannot pass. `threadspair` is
        // the `-C 2-0` shape, the only one that tells the two halves of
        // `-C a-b` apart from a transposed parse.
        "threads",
        "threads3",
        "threadspair",
    ] {
        assert!(
            keys.iter().any(|k| k.ends_with(feature)),
            "no matrix case exercises `{}` -- regenerate with ./gen/gen_golden.py",
            feature
        );
    }
    for input in ["aa", "nuc", "ragged", "dup", "tandem", "gapped"] {
        assert!(
            keys.iter().any(|k| k.starts_with(&format!("{}.", input))),
            "no matrix case uses input `{}`",
            input
        );
    }
}
