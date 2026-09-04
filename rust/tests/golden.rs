//! Golden-output tests: run the built binary and compare against output
//! captured from the C `disttbfast` (MAFFT 7.526) on `test/sample`.

use std::process::Command;

const ARGS_NWNS2: &[&str] = &[
    "-q", "0", "-E", "2", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O", "-C", "0-0", "-b", "62",
    "-g", "0", "-f", "-1.53", "-Q", "100.0", "-h", "0", "-X", "0.1", "-x", "100",
];
const ARGS_NWNS1: &[&str] = &[
    "-q", "0", "-E", "1", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O", "-C", "0-0", "-b", "62",
    "-g", "0", "-f", "-1.53", "-Q", "100.0", "-h", "0", "-X", "0.1", "-x", "100",
];
const ARGS_FFTNS2: &[&str] = &[
    "-q", "0", "-E", "2", "-V", "-1.53", "-s", "0.0", "-W", "6", "-O", "-C", "0-0", "-b", "62",
    "-g", "0", "-f", "-1.53", "-Q", "100.0", "-h", "0", "-F", "-X", "0.1", "-x", "100",
];

fn run(args: &[&str], dir: &std::path::Path) -> String {
    let input = concat!(env!("CARGO_MANIFEST_DIR"), "/../test/sample");
    let f = std::fs::File::open(input).expect("test/sample");
    let out = Command::new(env!("CARGO_BIN_EXE_disttbfast-rs"))
        .args(args)
        .current_dir(dir)
        .stdin(std::process::Stdio::from(f))
        .stderr(std::process::Stdio::null())
        .output()
        .expect("run disttbfast-rs");
    assert!(out.status.success(), "exit {:?}", out.status);
    String::from_utf8(out.stdout).expect("utf-8 output")
}

fn tempdir(tag: &str) -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("disttbfast-rs-test-{}-{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

#[test]
fn nwns2_matches_golden() {
    let d = tempdir("nwns2");
    let got = run(ARGS_NWNS2, &d);
    assert_eq!(got, include_str!("golden/sample.nwns2.fa"));
    let order = std::fs::read_to_string(d.join("order")).expect("order file");
    assert_eq!(order, include_str!("golden/sample.nwns2.order"));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn nwns1_matches_golden() {
    let d = tempdir("nwns1");
    let got = run(ARGS_NWNS1, &d);
    assert_eq!(got, include_str!("golden/sample.nwns1.fa"));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn fftns2_matches_golden() {
    let d = tempdir("fftns2");
    let got = run(ARGS_FFTNS2, &d);
    assert_eq!(got, include_str!("golden/sample.fftns2.fa"));
    let _ = std::fs::remove_dir_all(&d);
}

#[test]
fn unsupported_configuration_is_refused() {
    // -Y (--keeplength) is not ported: it must fail loudly rather than
    // quietly returning a different alignment.
    let d = tempdir("unsupported");
    let input = concat!(env!("CARGO_MANIFEST_DIR"), "/../test/sample");
    let f = std::fs::File::open(input).unwrap();
    let mut args: Vec<&str> = ARGS_NWNS2.to_vec();
    args.push("-Y");
    let out = Command::new(env!("CARGO_BIN_EXE_disttbfast-rs"))
        .args(&args)
        .current_dir(&d)
        .stdin(std::process::Stdio::from(f))
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("keeplength"));
    let _ = std::fs::remove_dir_all(&d);
}

// ---------------------------------------------------------------------------
// tbfast: the consistency-based strategies. Unlike `disttbfast`, `tbfast`
// writes its alignment to a file called `pre` rather than to stdout, which is
// why these read the working directory instead of capturing output.
//
// The argument vectors are the ones the `mafft` wrapper builds, split at the
// bare `_` markers into a pairwise stage and a progressive stage. They match
// the L-INS-i, G-INS-i and E-INS-i entries in `verify_tbfast.sh`.

const P_LOCAL: &[&str] = &[
    "-u", "0.0", "-l", "2.7", "-C", "0", "-b", "62", "-g", "-0.100", "-f", "-2.00", "-Q", "100.0",
    "-h", "0.1",
];
const P_GLOBAL: &[&str] = &[
    "-u", "0.0", "-l", "2.7", "-C", "0", "-b", "62", "-g", "-0.10", "-f", "-2.00", "-Q", "100.0",
    "-h", "0.1",
];
const P_GENAF: &[&str] = &[
    "-u", "0.0", "-l", "2.7", "-C", "0", "-b", "62", "-g", "0.0", "-f", "-2.00", "-Q", "100.0",
    "-h", "0.0", "-O", "-6.00", "-E", "-0.000", "-N", "-Z",
];
const T_STAGE: &[&str] = &[
    "-W", "0.00001", "-V", "-1.53", "-s", "0.0", "-C", "0", "-b", "62", "-f", "-1.53", "-Q",
    "100.0", "-h", "0", "-F", "-l", "2.7", "-X", "0.1",
];
const T_GENAF: &[&str] = &[
    "-W", "0.00001", "-V", "-1.53", "-s", "0.0", "-O", "-C", "0", "-b", "62", "-f", "-1.53", "-Q",
    "100.0", "-h", "0.000", "-F", "-l", "2.7", "-X", "0.1",
];

fn run_tbfast(pair: &[&str], pair_alg: &[&str], prog: &[&str], dir: &std::path::Path) {
    let input = concat!(env!("CARGO_MANIFEST_DIR"), "/../test/sample");
    let f = std::fs::File::open(input).expect("test/sample");
    let mut args: Vec<&str> = vec!["_"];
    args.extend_from_slice(pair);
    args.extend_from_slice(pair_alg);
    args.push("_");
    args.extend_from_slice(prog);
    let out = Command::new(env!("CARGO_BIN_EXE_tbfast-rs"))
        .args(&args)
        .current_dir(dir)
        .stdin(std::process::Stdio::from(f))
        .stderr(std::process::Stdio::null())
        .output()
        .expect("run tbfast-rs");
    assert!(out.status.success(), "exit {:?}", out.status);
}

fn read(dir: &std::path::Path, name: &str) -> String {
    std::fs::read_to_string(dir.join(name)).unwrap_or_else(|e| panic!("{}: {}", name, e))
}

#[test]
fn linsi_matches_golden() {
    let d = tempdir("linsi");
    run_tbfast(P_LOCAL, &["-L"], &[&["-+", "2"][..], T_STAGE, &["-O"]].concat(), &d);
    assert_eq!(read(&d, "pre"), include_str!("golden/sample.linsi.pre"));
    assert_eq!(read(&d, "order"), include_str!("golden/sample.linsi.order"));
}

#[test]
fn ginsi_matches_golden() {
    let d = tempdir("ginsi");
    run_tbfast(P_GLOBAL, &["-A"], &[&["-+", "2"][..], T_STAGE].concat(), &d);
    assert_eq!(read(&d, "pre"), include_str!("golden/sample.ginsi.pre"));
    assert_eq!(read(&d, "order"), include_str!("golden/sample.ginsi.order"));
}

#[test]
fn einsi_matches_golden() {
    let d = tempdir("einsi");
    run_tbfast(P_GENAF, &[], &[&["-+", "2"][..], T_GENAF].concat(), &d);
    assert_eq!(read(&d, "pre"), include_str!("golden/sample.einsi.pre"));
    assert_eq!(read(&d, "order"), include_str!("golden/sample.einsi.order"));
}
