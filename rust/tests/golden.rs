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
