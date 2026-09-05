//! Mirror of `gen_tdharness.c`: run the branch-weighting family on one staged
//! guide tree and print the same lines the C harness does, so `verify_td.sh`
//! can diff them.
#[path = "../tests/support/tdcase.rs"]
mod tdcase;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    print!("{}", tdcase::run(&a));
}
