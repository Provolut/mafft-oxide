//! Scoring-matrix and gap-penalty setup, transliterated from `constants()` in
//! `core/constants.c` (the `dorp == 'p' && scoremtx == 1` branch).
//!
//! Every arithmetic step is kept in the original order so the resulting integer
//! matrix is identical to the C one.  The matrix is derived, not tabulated,
//! precisely so that it stays identical if the inputs ever change.

use crate::blosum;
use crate::params::Args;

pub const NALPHABETS: usize = 26;
pub const CHARSIZE: usize = 0x80;

/// `shishagonyuu()`: round half away from zero, truncating like the C cast.
fn shishagonyuu(v: f64) -> i32 {
    if v > 0.0 {
        (v + 0.5) as i32
    } else if v == 0.0 {
        0
    } else {
        (v - 0.5) as i32
    }
}

/// Everything `constants()` publishes as a global.
pub struct Constants {
    /// Residue alphabet, index -> character.
    pub amino: [u8; NALPHABETS],
    /// Character -> alphabet index, or -1.
    pub amino_n: [i32; CHARSIZE],
    /// Character -> residue group (used by the 6-mer distance).
    pub amino_grp: [i32; 0x100],
    /// Integer score matrix over alphabet indices.  The two matrices below are
    /// derived from it; it is kept because it is what the golden test compares.
    #[allow(dead_code)]
    pub n_dis: [[i32; NALPHABETS]; NALPHABETS],
    /// `n_dis` scaled by `consweight_multi` (1.0), used by the profile DP.
    pub n_dis_consweight_multi: Vec<Vec<f64>>,
    /// Score matrix indexed by character, used by the 1-vs-1 DP and by
    /// `naivepairscore*`.
    pub amino_dis: Vec<Vec<i32>>,

    pub penalty: i32,
    pub penalty_ex: i32,
    pub penalty_dist: i32,
    /// Only read by the shift ("warp") DP, which this port refuses to run;
    /// kept because `constants()` computes it and it gates `trywarp`.
    #[allow(dead_code)]
    pub penalty_shift: i32,
    /// Folded into the score matrix above; kept for parity with the C global.
    #[allow(dead_code)]
    pub offset: i32,
    pub trywarp: bool,
}

// Defaults from blosum.h / constants.c for the protein BLOSUM branch.
const DEFAULTGOP_B: i32 = -1530;
const DEFAULTGEP_B: i32 = 0;
const DEFAULTOFS_B: i32 = -123;
const CONSWEIGHT_MULTI: f64 = 1.0;

impl Constants {
    pub fn new(args: &Args) -> Constants {
        let trywarp = !(args.penalty_shift_factor >= 10.0);

        let mut ppenalty = args.ppenalty;
        let mut ppenalty_dist = args.ppenalty_dist;
        let mut ppenalty_ex = args.ppenalty_ex;
        let mut poffset = args.poffset;
        if ppenalty == crate::params::NOTSPECIFIED {
            ppenalty = DEFAULTGOP_B;
        }
        if ppenalty_dist == crate::params::NOTSPECIFIED {
            ppenalty_dist = ppenalty;
        }
        if ppenalty_ex == crate::params::NOTSPECIFIED {
            ppenalty_ex = DEFAULTGEP_B;
        }
        if poffset == crate::params::NOTSPECIFIED {
            poffset = DEFAULTOFS_B;
        }

        let penalty = (600.0 / 1000.0 * ppenalty as f64 + 0.5) as i32;
        let penalty_dist = (600.0 / 1000.0 * ppenalty_dist as f64 + 0.5) as i32;
        let penalty_shift = (args.penalty_shift_factor * penalty as f64) as i32;
        let penalty_ex = (600.0 / 1000.0 * ppenalty_ex as f64 + 0.5) as i32;
        let offset = (600.0 / 1000.0 * poffset as f64 + 0.5) as i32;

        let mut n_distmp = [[0.0f64; 20]; 20];
        let mut freq = [0.0f64; 20];
        let mut amino = [0u8; NALPHABETS];
        let mut amino_grp = [0i32; 0x100];
        blosum::blosum_mtx(
            args.nblosum,
            &mut n_distmp,
            &mut freq,
            &mut amino,
            &mut amino_grp,
        );

        // freq1 = freq, since fmodel != 1 is enforced by check_supported().
        let freq1 = &freq;

        let mut average = 0.0f64;
        for i in 0..20 {
            for j in 0..20 {
                average += n_distmp[i][j] * freq1[i] * freq1[j];
            }
        }
        // rescale == 1 for the built-in tables.
        for i in 0..20 {
            for j in 0..20 {
                n_distmp[i][j] -= average;
            }
        }

        let mut average = 0.0f64;
        for i in 0..20 {
            average += n_distmp[i][i] * freq1[i];
        }
        for i in 0..20 {
            for j in 0..20 {
                n_distmp[i][j] *= 600.0 / average;
            }
        }

        for i in 0..20 {
            for j in 0..20 {
                n_distmp[i][j] -= offset as f64;
            }
        }
        for i in 0..20 {
            for j in 0..20 {
                n_distmp[i][j] = shishagonyuu(n_distmp[i][j]) as f64;
            }
        }

        let mut n_dis = [[0i32; NALPHABETS]; NALPHABETS];
        for i in 0..20 {
            for j in 0..20 {
                n_dis[i][j] = n_distmp[i][j] as i32;
            }
        }

        // --- common tail of constants() ---
        let mut amino_n = [-1i32; CHARSIZE];
        for i in 0..NALPHABETS {
            amino_n[amino[i] as usize] = i as i32;
        }

        let mut amino_dis = vec![vec![0i32; CHARSIZE]; CHARSIZE];
        let mut n_dis_consweight_multi = vec![vec![0.0f64; NALPHABETS]; NALPHABETS];
        for i in 0..NALPHABETS {
            for j in 0..NALPHABETS {
                amino_dis[amino[i] as usize][amino[j] as usize] = n_dis[i][j];
                n_dis_consweight_multi[i][j] = n_dis[i][j] as f64 * CONSWEIGHT_MULTI;
            }
        }

        Constants {
            amino,
            amino_n,
            amino_grp,
            n_dis,
            n_dis_consweight_multi,
            amino_dis,
            penalty,
            penalty_ex,
            penalty_dist,
            penalty_shift,
            offset,
            trywarp,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Args;

    /// Values dumped from the C `constants()` (see `tests/golden/constants.txt`,
    /// produced by linking against `core/constants.c`).  The matrix is derived
    /// by the same floating-point steps in both languages, so any divergence in
    /// rounding shows up here rather than as a mysterious alignment difference.
    #[test]
    fn matches_c_constants() {
        let golden = include_str!("../tests/golden/constants.txt");
        let args = Args {
            nblosum: 62,
            scoremtx: 1,
            ppenalty: -1530,
            ppenalty_ex: 0,
            ppenalty_dist: -1530,
            poffset: 0,
            penalty_shift_factor: 100.0,
            ..Args::default()
        };
        let cn = Constants::new(&args);

        let mut lines = golden.lines();
        let mut kv = || -> (String, i64) {
            let l = lines.next().expect("golden line");
            let mut it = l.split_whitespace();
            let k = it.next().unwrap().to_string();
            let v: i64 = it.next().unwrap().parse().unwrap();
            (k, v)
        };
        assert_eq!(kv(), ("penalty".into(), cn.penalty as i64));
        assert_eq!(kv(), ("penalty_ex".into(), cn.penalty_ex as i64));
        assert_eq!(kv(), ("penalty_dist".into(), cn.penalty_dist as i64));
        assert_eq!(kv(), ("offset".into(), cn.offset as i64));
        assert_eq!(kv(), ("trywarp".into(), cn.trywarp as i64));

        for i in 0..NALPHABETS {
            let row: Vec<i32> = lines
                .next()
                .expect("golden matrix row")
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            assert_eq!(row.len(), NALPHABETS, "row {}", i);
            for j in 0..NALPHABETS {
                assert_eq!(row[j], cn.n_dis[i][j], "n_dis[{}][{}]", i, j);
            }
        }

        let aminoline = lines.next().unwrap();
        assert_eq!(
            aminoline.split_whitespace().nth(1).unwrap().as_bytes(),
            &cn.amino[..]
        );
    }
}
