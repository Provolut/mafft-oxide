//! Mirror of `gen_flhharness.c`: run `falign_localhom` on staged groups and
//! print the same lines the C harness does, so `verify_flh.sh` can diff them.
//!
//! With `--sequence`, cases are read from stdin (one tab-separated argument
//! vector per line) and run in order through one state, which is what the C
//! harness does with its process-lifetime statics.
#[path = "../tests/support/flhcase.rs"]
mod flhcase;

use std::io::Read;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    if a.first().map(String::as_str) == Some("--sequence") {
        let mut buf = String::new();
        std::io::stdin().read_to_string(&mut buf).unwrap();
        let cases: Vec<Vec<String>> = buf
            .lines()
            .filter(|l| !l.is_empty())
            .map(|l| l.split('\t').map(String::from).collect())
            .collect();
        print!("{}", flhcase::run_sequence(&cases));
        return;
    }
    print!("{}", flhcase::run(&a));
}
