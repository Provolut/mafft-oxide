//! Unpack the golden fixtures before anything reads them.
//!
//! `tests/golden/<suite>.pack` is tracked; `tests/golden/<suite>/` is not. The
//! fixtures were 1914 loose files, 95% of everything this repository tracks,
//! and git pays for each of them on every checkout and every `status` -- while
//! being 2 MB of text in total, so the file count was the whole cost. Packed,
//! they are five files that still diff line by line, which is the property they
//! exist for: a fixture that changes is the port changing its answer, and a
//! reviewer has to be able to see it.
//!
//! The tests are unchanged and still read ordinary files. This is deliberately
//! the smaller of the two ways to do it: the other one teaches six test files
//! and their support modules to read from a container instead of the
//! filesystem, which is a great deal more code to get wrong for the same
//! result.
//!
//! Written in Rust rather than calling `gen/pack_golden.py` because this runs
//! on every `cargo build`, including the release builds that produce the
//! shipped binaries. Those must not acquire a Python dependency.

use std::fs;
use std::path::{Path, PathBuf};

const HEADER: &[u8] = b"# mafft-rs golden pack v1\n";

fn main() {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let Ok(entries) = fs::read_dir(&root) else {
        return; // No fixtures here at all; nothing to do.
    };
    for e in entries.flatten() {
        let pack = e.path();
        if pack.extension().and_then(|s| s.to_str()) != Some("pack") {
            continue;
        }
        println!("cargo:rerun-if-changed={}", pack.display());
        let dst = pack.with_extension("");
        if is_current(&pack, &dst) {
            continue;
        }
        if let Err(e) = unpack(&pack, &dst) {
            // A fixture directory that silently failed to appear would make
            // every test using it fail with something that looks like a port
            // bug. Say what actually happened.
            panic!("cannot unpack {}: {e}", pack.display());
        }
    }
}

/// Has `dst` already been unpacked from a `pack` no newer than itself?
fn is_current(pack: &Path, dst: &Path) -> bool {
    let (Ok(p), Ok(d)) = (fs::metadata(pack), fs::metadata(dst)) else {
        return false;
    };
    if !d.is_dir() {
        return false;
    }
    match (p.modified(), d.modified()) {
        (Ok(pt), Ok(dt)) => dt >= pt,
        _ => false,
    }
}

/// Recreate the directory the pack was made from.
///
/// Records are length-prefixed rather than delimited, so no fixture content can
/// imitate a boundary -- these hold arbitrary alignment output, and any
/// sentinel line would turn up inside one of them eventually.
fn unpack(pack: &Path, dst: &Path) -> std::io::Result<()> {
    let raw = fs::read(pack)?;
    if !raw.starts_with(HEADER) {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "not a golden pack",
        ));
    }
    // Rebuilt from scratch: a file dropped from the pack must disappear here
    // too, or a stale fixture would go on being tested against.
    if dst.exists() {
        fs::remove_dir_all(dst)?;
    }
    fs::create_dir_all(dst)?;

    let mut i = HEADER.len();
    while i < raw.len() {
        let nl = raw[i..]
            .iter()
            .position(|&b| b == b'\n')
            .map(|n| i + n)
            .ok_or_else(|| bad("truncated record header"))?;
        let head = std::str::from_utf8(&raw[i..nl]).map_err(|_| bad("non-utf8 record header"))?;
        let (rel, count) = head.rsplit_once(' ').ok_or_else(|| bad("malformed header"))?;
        let count: usize = count.parse().map_err(|_| bad("bad byte count"))?;
        i = nl + 1;
        let end = i + count;
        if end > raw.len() {
            return Err(bad("record runs past the end of the pack"));
        }
        let target = dst.join(rel);
        if let Some(parent) = target.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(&target, &raw[i..end])?;
        i = end;
        // The packer adds a newline after a payload that lacks one, purely so
        // the next header starts a line; the count above is what delimits.
        if raw.get(i) == Some(&b'\n') {
            i += 1;
        }
    }
    Ok(())
}

fn bad(msg: &str) -> std::io::Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, msg.to_string())
}
