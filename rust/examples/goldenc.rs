//! Stand in for the C binary in a differential harness: record what it does
//! once, replay it thereafter.
//!
//! ```text
//! goldenc record <store> <real-c> [args...]
//! goldenc replay <store> [args...]
//! ```
//!
//! The harnesses take the C binary as an argument and otherwise care nothing
//! about it, so a program that reproduces its bytes is indistinguishable from
//! it. That is the whole idea here: `verify.sh` and `verify_tbfast.sh` are not
//! modified at all, keep making exactly the comparisons they made before, and
//! keep reporting the same tallies -- they simply stop waiting for an aligner
//! that already answered this question once.
//!
//! Why a compiled binary rather than a script: this is executed once per case,
//! 13708 times. A Python stand-in costs ~25ms of interpreter startup, which is
//! more than the ~21ms C run it is replacing -- it would be slower than doing
//! the work. This starts in about a millisecond.
//!
//! # What a case is keyed on
//!
//! Everything that can change what the C writes: the argument vector, the bytes
//! on stdin, and the working directory it starts in -- the harnesses stage
//! `_guidetree`, `_aamtx`, `_externalanchors` and friends there before running.
//! Nothing else reaches these binaries.
//!
//! The key is a 128-bit FNV-1a over that material rather than a cryptographic
//! hash, which needs no dependency and is safe here for a reason worth stating:
//! a collision replays some *other* case's output, the port's real output then
//! differs from it, and the harness reports a failure. Collisions fail loudly
//! and towards refusing, never towards a silent pass.
//!
//! # What it does not do
//!
//! It does not verify that the recording is current. A store made from one
//! `core/` and replayed against a port built from a newer upstream import would
//! be comparing against history. `verify_golden.sh` records the C's identity
//! alongside the store and checks it; this program only does the mechanics.

use std::collections::BTreeMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// FNV-1a, two independent 64-bit lanes making a 128-bit key.
struct Key {
    a: u64,
    b: u64,
}

impl Key {
    fn new() -> Self {
        Self {
            a: 0xcbf2_9ce4_8422_2325,
            b: 0x9e37_79b9_7f4a_7c15,
        }
    }

    fn eat(&mut self, bytes: &[u8]) {
        for &x in bytes {
            self.a = (self.a ^ u64::from(x)).wrapping_mul(0x0000_0100_0000_01b3);
            self.b = (self.b ^ u64::from(x)).wrapping_mul(0x8864_3f65_0000_002f);
        }
        // Length-separated, so that ("ab", "c") and ("a", "bc") differ.
        let n = bytes.len() as u64;
        self.a ^= n;
        self.b = self.b.rotate_left(7) ^ n;
    }

    fn hex(&self) -> String {
        format!("{:016x}{:016x}", self.a, self.b)
    }
}

/// Every file under `dir`, by path relative to it.
fn walk(dir: &Path, prefix: &Path, out: &mut BTreeMap<PathBuf, Vec<u8>>) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for e in entries.flatten() {
        let p = e.path();
        let rel = prefix.join(e.file_name());
        if p.is_dir() {
            walk(&p, &rel, out);
        } else if let Ok(b) = std::fs::read(&p) {
            out.insert(rel, b);
        }
    }
}

fn snapshot() -> BTreeMap<PathBuf, Vec<u8>> {
    let mut m = BTreeMap::new();
    walk(Path::new("."), Path::new(""), &mut m);
    m
}

fn key_of(args: &[String], stdin: &[u8], before: &BTreeMap<PathBuf, Vec<u8>>) -> String {
    let mut k = Key::new();
    for a in args {
        k.eat(a.as_bytes());
    }
    k.eat(stdin);
    for (p, b) in before {
        k.eat(p.to_string_lossy().as_bytes());
        k.eat(b);
    }
    k.hex()
}

fn die(msg: &str) -> ! {
    eprintln!("goldenc: {msg}");
    // Distinct from anything the C returns, so a harness sees a status
    // mismatch and reports a failure rather than mistaking it for agreement.
    std::process::exit(97)
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 3 {
        die("usage: goldenc record <store> <real-c> [args...] | replay <store> [args...]");
    }
    let mode = argv[1].clone();
    let store = PathBuf::from(&argv[2]);

    let mut stdin = Vec::new();
    let _ = std::io::stdin().read_to_end(&mut stdin);
    let before = snapshot();

    match mode.as_str() {
        "record" => {
            if argv.len() < 4 {
                die("record needs the real binary");
            }
            let real = &argv[3];
            let args: Vec<String> = argv[4..].to_vec();
            let key = key_of(&args, &stdin, &before);

            // A real file, not a pipe. The harness hands these binaries a
            // redirected file, and MAFFT's readers seek on it -- give the child
            // a pipe instead and it reads a prefix, or nothing, and writes a
            // different alignment. That is not a subtle difference: it fails
            // every case, which is how it was found.
            let tmp_in = std::env::temp_dir().join(format!("goldenc.{}.in", std::process::id()));
            std::fs::write(&tmp_in, &stdin)
                .unwrap_or_else(|e| die(&format!("cannot stage stdin: {e}")));
            let fh = std::fs::File::open(&tmp_in)
                .unwrap_or_else(|e| die(&format!("cannot reopen stdin: {e}")));
            let child = Command::new(real)
                .args(&args)
                .stdin(Stdio::from(fh))
                .stdout(Stdio::piped())
                .spawn()
                .unwrap_or_else(|e| die(&format!("cannot run {real}: {e}")));
            let out = child.wait_with_output().unwrap_or_else(|e| die(&format!("{real}: {e}")));
            let _ = std::fs::remove_file(&tmp_in);
            // A killed child has no exit code, and 16 of `verify.sh`'s cases
            // exist precisely because the C crashes there -- the harness scores
            // "C died, port refused cleanly" as a pass, keyed on a status at or
            // above 128. Recording those as a plain 1 would turn every one of
            // them into a failure on replay.
            let code = match out.status.code() {
                Some(c) => c,
                None => {
                    use std::os::unix::process::ExitStatusExt;
                    128 + out.status.signal().unwrap_or(0)
                }
            };

            // Only what this run produced: the harness stages input files into
            // the same directory, and storing those back would make the replay
            // rewrite its own inputs.
            let after = snapshot();
            let mut produced: Vec<(PathBuf, Vec<u8>)> = Vec::new();
            for (p, b) in after {
                if before.get(&p) != Some(&b) {
                    produced.push((p, b));
                }
            }

            // Written aside and renamed into place, because the recording runs
            // under `verify_parallel.sh` with eight workers sharing one store.
            // The key is derived from content, so two workers racing on one
            // case are writing the same bytes; last one wins, harmlessly.
            let dst = store.join(&key);
            if !dst.exists() {
                let tmp = store.join(format!(".tmp.{}.{key}", std::process::id()));
                let _ = std::fs::create_dir_all(&tmp);
                let _ = std::fs::write(tmp.join("__exit"), code.to_string());
                let _ = std::fs::write(tmp.join("__stdout"), &out.stdout);
                for (p, b) in &produced {
                    let f = tmp.join("f").join(p);
                    if let Some(d) = f.parent() {
                        let _ = std::fs::create_dir_all(d);
                    }
                    let _ = std::fs::write(&f, b);
                }
                if std::fs::rename(&tmp, &dst).is_err() {
                    let _ = std::fs::remove_dir_all(&tmp);
                }
            }

            std::io::stdout().write_all(&out.stdout).ok();
            std::process::exit(code);
        }
        "replay" => {
            let args: Vec<String> = argv[3..].to_vec();
            let key = key_of(&args, &stdin, &before);
            let dir = store.join(&key);
            if !dir.is_dir() {
                // Never fall back to running anything: a store that has lost a
                // case must stop the run, not quietly become a different test.
                die(&format!("no recording for {key} ({})", args.join(" ")));
            }
            let mut files = BTreeMap::new();
            walk(&dir.join("f"), Path::new(""), &mut files);
            for (p, b) in files {
                if let Some(d) = p.parent() {
                    if !d.as_os_str().is_empty() {
                        let _ = std::fs::create_dir_all(d);
                    }
                }
                if let Err(e) = std::fs::write(&p, &b) {
                    die(&format!("cannot restore {}: {e}", p.display()));
                }
            }
            let stdout = std::fs::read(dir.join("__stdout")).unwrap_or_default();
            std::io::stdout().write_all(&stdout).ok();
            let code: i32 = std::fs::read_to_string(dir.join("__exit"))
                .ok()
                .and_then(|s| s.trim().parse().ok())
                .unwrap_or_else(|| die("recording has no exit status"));
            std::process::exit(code);
        }
        other => die(&format!("unknown mode {other}")),
    }
}
