//! Bit-exact Rust ports of MAFFT 7.526's alignment algorithms.
//!
//! Each module names the C file it was transliterated from. The rule
//! throughout is that the C's evaluation order, its persistent scratch
//! buffers and its edge-case behaviour are reproduced exactly, because the
//! output has to match byte for byte -- not because any of it is good style.
//!
//! No `unsafe`, and the lint below keeps it that way. The C's pointer walks
//! are reproduced with indices and slices, its `static` scratch buffers with
//! owned state passed explicitly, and its `char **` matrices with `Vec<Vec<_>>`
//! -- so the transliteration is faithful in arithmetic without inheriting the
//! memory model. That is what lets the port stop at upstream's out-of-bounds
//! writes (issue #4) instead of reproducing them.
// `deny`, not `forbid`, for exactly two modules: `lalign_avx512` and
// `partalign_avx512` allow themselves,
// and nothing else in the crate can. `forbid` cannot be lifted locally, which
// is the point of the two lints -- so this is the whole of the licence the
// kernel gets, and it is visible here rather than buried.
#![deny(unsafe_code)]

pub mod addfunctions;
pub mod align;
pub mod anchors;
pub mod blosum;
pub mod cargs;
pub mod cerand;
pub mod compacttree;
pub mod constants;
pub mod dalign;
pub mod dist;
pub mod falign;
pub mod fft;
pub mod genalign;
pub mod io;
pub mod jtt;
pub mod keeplength;
pub mod lalign;
pub mod lalign_avx512;
pub mod localhom;
pub mod msalign;
pub mod naivescore;
pub mod pairalign;
pub mod partalign;
pub mod partalign_avx512;
pub mod params;
pub mod profile;
pub mod saalign;
pub mod score;
pub mod tditeration;
pub mod tree;
pub mod treeop;
pub mod treebase;
pub mod warp;

#[cfg(test)]
mod tests_unit;

/// One lock for every test that swaps the working directory.
///
/// The working directory is process-wide, and four helpers across three
/// modules change it to run a case in a scratch directory. Each of them had a
/// `static LOCK` of its own, which serialises a helper against itself and
/// against nothing else -- so two of them would interleave, and a test would
/// find itself running in another helper's directory. It failed the way that
/// implies: the file the helper had just written was "not found", and the
/// restore afterwards hit `ENOENT` on a directory the other thread had already
/// removed. Roughly one run in eight.
///
/// Every one of those helpers takes this lock instead.
#[cfg(test)]
pub(crate) fn cwd_lock() -> std::sync::MutexGuard<'static, ()> {
    static LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    LOCK.lock().unwrap_or_else(|e| e.into_inner())
}
