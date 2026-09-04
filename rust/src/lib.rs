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
#![forbid(unsafe_code)]

pub mod addfunctions;
pub mod align;
pub mod anchors;
pub mod blosum;
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
pub mod localhom;
pub mod msalign;
pub mod pairalign;
pub mod params;
pub mod profile;
pub mod saalign;
pub mod score;
pub mod tree;
pub mod treebase;
pub mod warp;

#[cfg(test)]
mod tests_unit;
