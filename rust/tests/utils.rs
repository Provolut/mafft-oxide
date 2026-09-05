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
//!
//! `replaceu` and `f2cl` are pinned by output digest rather than by literal
//! text: their output is a whole alignment, and a table of digests stays
//! reviewable where ten embedded FASTA files would not. `restoreu` needs two
//! input files that agree, so it gets a directory of its own.

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


/// Read a `<input>\t<sha256>` table.
fn expected_digests(name: &str) -> Vec<(String, String)> {
    let path = PathBuf::from(GOLDEN).join(name);
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|_| panic!("missing {}: regenerate with verify_utils.sh", path.display()));
    text.lines()
        .filter(|l| !l.is_empty())
        .map(|l| {
            let (n, d) = l.split_once('\t').unwrap_or_else(|| panic!("malformed: {l}"));
            (n.to_string(), d.to_string())
        })
        .collect()
}

/// SHA-256 of some bytes, hex, using the same digest the fixtures were made with.
///
/// Implemented here rather than pulled in as a dependency: the port is
/// deliberately dependency-free, and a test is not a reason to change that.
fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let mut msg = data.to_vec();
    let bitlen = (data.len() as u64) * 8;
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bitlen.to_be_bytes());
    for block in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                block[i * 4],
                block[i * 4 + 1],
                block[i * 4 + 2],
                block[i * 4 + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
            (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        for (dst, v) in h.iter_mut().zip([a, b, c, d, e, f, g, hh]) {
            *dst = dst.wrapping_add(v);
        }
    }
    h.iter().map(|v| format!("{v:08x}")).collect()
}

#[test]
fn the_digest_helper_agrees_with_a_known_vector() {
    // Guards the fixtures: a broken hash here would make every digest
    // comparison below pass or fail for the wrong reason.
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
}

#[test]
fn replaceu_matches_the_c_on_every_fixture() {
    let bin = env!("CARGO_BIN_EXE_replaceu-rs");
    let cases = expected_digests("replaceu.expected");
    assert!(cases.len() >= 10, "fixture table looks truncated");
    for (name, want) in cases {
        let input = Path::new(GOLDEN).join("in").join(&name);
        let got = run(bin, &["-i", input.to_str().unwrap()]);
        assert_eq!(sha256_hex(&got), want, "replaceu disagrees with the C on {name}");
    }
}

#[test]
fn f2cl_matches_the_c_on_the_fasta_path() {
    let bin = env!("CARGO_BIN_EXE_f2cl-rs");
    let cases = expected_digests("f2cl.expected");
    assert!(cases.len() >= 10, "fixture table looks truncated");
    for (name, want) in cases {
        let input = Path::new(GOLDEN).join("in").join(&name);
        let got = run(
            bin,
            &["-n", "-1", "-f", "-l", "60", "-i", input.to_str().unwrap()],
        );
        assert_eq!(sha256_hex(&got), want, "f2cl disagrees with the C on {name}");
    }
}

#[test]
fn f2cl_refuses_the_formats_it_does_not_implement() {
    // Refusing is the contract: a plausible-looking clustal file from a writer
    // that was never compared against the C is the one outcome to avoid.
    let bin = env!("CARGO_BIN_EXE_f2cl-rs");
    let input = Path::new(GOLDEN).join("in").join("sample.fa");
    for args in [
        vec!["-n", "-1", "-i", input.to_str().unwrap()],
        vec!["-n", "-1", "-y", "-i", input.to_str().unwrap()],
    ] {
        let out = Command::new(bin)
            .args(&args)
            .stdin(std::process::Stdio::null())
            .output()
            .expect("run f2cl-rs");
        assert_eq!(out.status.code(), Some(1), "args were {args:?}");
        assert!(out.stdout.is_empty(), "a refused format must write nothing");
        assert!(String::from_utf8_lossy(&out.stderr).contains("not implemented"));
    }
}

#[test]
fn restoreu_puts_the_original_characters_back() {
    let bin = env!("CARGO_BIN_EXE_restoreu-rs");
    let dir = PathBuf::from(GOLDEN).join("restoreu");
    let want = std::fs::read(dir.join("expected.out")).expect("fixture");
    let out = Command::new(bin)
        .args([
            "-a",
            dir.join("pre").to_str().unwrap(),
            "-i",
            dir.join("orig").to_str().unwrap(),
        ])
        .output()
        .expect("run restoreu-rs");
    assert_eq!(out.stdout, want);
}

#[test]
fn restoreu_reports_a_record_count_mismatch_instead_of_crashing() {
    // The C segfaults here with an empty stderr: its own check sits after the
    // loops that overrun (upstream issue #10). The port checks first.
    let bin = env!("CARGO_BIN_EXE_restoreu-rs");
    let dir = std::env::temp_dir().join("restoreu_mismatch_fixture");
    let _ = std::fs::create_dir_all(&dir);
    std::fs::write(dir.join("orig"), ">a\nMKVL\n>b\nMKVI\n").expect("write");
    std::fs::write(
        dir.join("pre"),
        ">_os_0_oe_a\nMKVL\n>_os_1_oe_b\nMKVI\n>_os_0_oe_a\nMKVL\n",
    )
    .expect("write");
    let out = Command::new(bin)
        .args([
            "-a",
            dir.join("pre").to_str().unwrap(),
            "-i",
            dir.join("orig").to_str().unwrap(),
        ])
        .output()
        .expect("run restoreu-rs");
    assert_eq!(out.status.code(), Some(1));
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("njob = 2"), "stderr was: {err}");
    assert!(err.contains("njob = 3"), "stderr was: {err}");
}

#[test]
fn restoreu_refuses_when_no_alignment_is_given() {
    let bin = env!("CARGO_BIN_EXE_restoreu-rs");
    let dir = PathBuf::from(GOLDEN).join("restoreu");
    let out = Command::new(bin)
        .args(["-i", dir.join("orig").to_str().unwrap()])
        .output()
        .expect("run restoreu-rs");
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("No alignment is given."));
}


/// Run a binary and hand back (stdout, stderr, exit code).
fn run_full(bin: &str, args: &[&str]) -> (Vec<u8>, String, Option<i32>) {
    let out = Command::new(bin)
        .args(args)
        .stdin(std::process::Stdio::null())
        .output()
        .unwrap_or_else(|e| panic!("could not run {bin}: {e}"));
    let err = String::from_utf8_lossy(&out.stderr).to_string();
    (out.stdout, err, out.status.code())
}

#[test]
fn f2cl_wraps_at_the_requested_width_and_not_at_all_when_negative() {
    // Two different branches in the C, and they disagree on an empty record:
    // the wrapped loop writes no sequence line, the unwrapped one writes a
    // blank. Both are pinned because both are reachable from the script.
    let bin = env!("CARGO_BIN_EXE_f2cl-rs");
    let dir = std::env::temp_dir().join("f2cl_wrap_fixture");
    let _ = std::fs::create_dir_all(&dir);
    let input = dir.join("in.fa");
    std::fs::write(&input, ">a\n\n>b\nMKVLTA\n").expect("write");
    let path = input.to_str().unwrap();

    let (wrapped, _, code) = run_full(bin, &["-n", "-1", "-f", "-l", "2", "-i", path]);
    assert_eq!(code, Some(0));
    assert_eq!(wrapped, b">a\n>b\nMK\nVL\nTA\n".to_vec());

    let (flat, _, code) = run_full(bin, &["-n", "-1", "-f", "-i", path]);
    assert_eq!(code, Some(0));
    assert_eq!(flat, b">a\n\n>b\nMKVLTA\n".to_vec());
}

#[test]
fn f2cl_refuses_a_zero_line_width() {
    // `-l 0` would be an infinite loop in the C, which is why it checks.
    let bin = env!("CARGO_BIN_EXE_f2cl-rs");
    let input = Path::new(GOLDEN).join("in").join("sample.fa");
    let (out, err, code) = run_full(
        bin,
        &["-n", "-1", "-f", "-l", "0", "-i", input.to_str().unwrap()],
    );
    assert_eq!(code, Some(1));
    assert!(out.is_empty());
    assert!(err.contains("LineLengthInFASTA = 0"), "stderr was: {err}");
}

#[test]
fn f2cl_honours_an_order_file() {
    let bin = env!("CARGO_BIN_EXE_f2cl-rs");
    let dir = std::env::temp_dir().join("f2cl_order_fixture");
    let _ = std::fs::create_dir_all(&dir);
    let input = dir.join("in.fa");
    std::fs::write(&input, ">a\nAA\n>b\nBB\n>c\nCC\n").expect("write");
    let order = dir.join("order");
    std::fs::write(&order, "2\n0\n1\n").expect("write");
    let (out, _, code) = run_full(
        bin,
        &[
            "-n", "-1", "-f", "-r",
            order.to_str().unwrap(),
            "-i",
            input.to_str().unwrap(),
        ],
    );
    assert_eq!(code, Some(0));
    assert_eq!(out, b">c\nCC\n>a\nAA\n>b\nBB\n".to_vec());
}

#[test]
fn restoreu_reports_a_missing_marker_as_a_format_error() {
    let bin = env!("CARGO_BIN_EXE_restoreu-rs");
    let dir = std::env::temp_dir().join("restoreu_format_fixture");
    let _ = std::fs::create_dir_all(&dir);
    std::fs::write(dir.join("orig"), ">a\nMKVL\n").expect("write");
    std::fs::write(dir.join("pre"), ">a\nMKVL\n").expect("write");
    let (out, err, code) = run_full(
        bin,
        &[
            "-a",
            dir.join("pre").to_str().unwrap(),
            "-i",
            dir.join("orig").to_str().unwrap(),
        ],
    );
    assert_eq!(code, Some(1));
    assert!(out.is_empty());
    assert!(err.contains("Format error!"), "stderr was: {err}");
}

#[test]
fn restoreu_reports_an_original_that_does_not_fit_the_alignment() {
    // Both directions of `fillorichar`'s consistency check: the original
    // running out mid-row, and having residues left over at the end.
    let bin = env!("CARGO_BIN_EXE_restoreu-rs");
    let dir = std::env::temp_dir().join("restoreu_inconsistent_fixture");
    let _ = std::fs::create_dir_all(&dir);

    std::fs::write(dir.join("orig"), ">a\nMK\n").expect("write");
    std::fs::write(dir.join("pre"), ">_os_0_oe_a\nMKVL\n").expect("write");
    let (_, err, code) = run_full(
        bin,
        &[
            "-a",
            dir.join("pre").to_str().unwrap(),
            "-i",
            dir.join("orig").to_str().unwrap(),
        ],
    );
    assert_eq!(code, Some(1));
    assert!(err.contains("(1)"), "stderr was: {err}");

    std::fs::write(dir.join("orig2"), ">a\nMKVLTA\n").expect("write");
    std::fs::write(dir.join("pre2"), ">_os_0_oe_a\nMKVL\n").expect("write");
    let (_, err, code) = run_full(
        bin,
        &[
            "-a",
            dir.join("pre2").to_str().unwrap(),
            "-i",
            dir.join("orig2").to_str().unwrap(),
        ],
    );
    assert_eq!(code, Some(1));
    assert!(err.contains("(2)"), "stderr was: {err}");
}


#[test]
fn dndpre_matches_the_c_on_every_aligned_fixture() {
    // `dndpre` writes `hat2` into the working directory rather than to stdout,
    // so each case runs in a directory of its own -- and that is also what
    // stops one case reading the previous one's file and passing on it.
    let bin = env!("CARGO_BIN_EXE_dndpre-rs");
    let cases = expected_digests("dndpre.expected");
    assert!(!cases.is_empty(), "fixture table is empty");
    for (name, want) in cases {
        let input = Path::new(GOLDEN).join("aligned").join(&name);
        let dir = std::env::temp_dir().join(format!("dndpre_case_{name}"));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("case dir");
        let status = Command::new(bin)
            .args(["-y", "hat2", "-b", "62", "-M", "2", "-C", "0", "-i"])
            .arg(&input)
            .current_dir(&dir)
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .expect("run dndpre-rs");
        assert!(status.success(), "dndpre-rs failed on {name}");
        let got = std::fs::read(dir.join("hat2")).expect("dndpre-rs wrote no hat2");
        assert_eq!(sha256_hex(&got), want, "dndpre disagrees with the C on {name}");
    }
}

#[test]
fn dndpre_refuses_fewer_than_two_records() {
    let bin = env!("CARGO_BIN_EXE_dndpre-rs");
    let input = Path::new(GOLDEN).join("in").join("one.fa");
    let (_, err, code) = run_full(bin, &["-y", "hat2", "-i", input.to_str().unwrap()]);
    assert_eq!(code, Some(1));
    assert!(err.contains("At least 2 sequences"), "stderr was: {err}");
}

#[test]
fn dndpre_refuses_an_unaligned_input() {
    // Its whole job is re-deriving distances from an alignment, so ragged
    // input is a caller error rather than something to measure.
    let bin = env!("CARGO_BIN_EXE_dndpre-rs");
    let dir = std::env::temp_dir().join("dndpre_ragged_fixture");
    let _ = std::fs::create_dir_all(&dir);
    let input = dir.join("ragged.fa");
    std::fs::write(&input, ">a\nMKVL\n>b\nMK\n").expect("write");
    let (_, err, code) = run_full(bin, &["-y", "hat2", "-i", input.to_str().unwrap()]);
    assert_eq!(code, Some(1));
    assert!(err.contains("Not aligned!"), "stderr was: {err}");
}

#[test]
fn dndpre_rejects_a_distance_format_it_does_not_know() {
    let bin = env!("CARGO_BIN_EXE_dndpre-rs");
    let input = Path::new(GOLDEN).join("aligned").join("sample.fa");
    let (_, err, code) = run_full(bin, &["-y", "zzz", "-i", input.to_str().unwrap()]);
    assert_eq!(code, Some(1));
    assert!(err.contains("Set -y c or -y h"), "stderr was: {err}");
}
