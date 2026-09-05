//! Mirror of `gen_flhharness.c`: run `falign_localhom` on one staged pair of
//! groups and print the same lines the C harness does, so `verify_flh.sh` can
//! diff them.
#[path = "../tests/support/flhcase.rs"]
mod flhcase;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    print!("{}", flhcase::run(&a));
}
