//! Mirror of `gen_tdiharness.c`: run the refinement loop on one staged
//! alignment and print the same lines the C harness does, so `verify_tdi.sh`
//! can diff them.
#[path = "../tests/support/tdicase.rs"]
mod tdicase;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    print!("{}", tdicase::run(&a));
}
