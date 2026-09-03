//! Scoring-matrix and gap-penalty setup, transliterated from `constants()` in
//! `core/constants.c` (the `dorp == 'p' && scoremtx == 1` branch).
//!
//! Every arithmetic step is kept in the original order so the resulting integer
//! matrix is identical to the C one.  The matrix is derived, not tabulated,
//! precisely so that it stays identical if the inputs ever change.

use crate::blosum;
use crate::params::Args;

/// Miyata's polarity scale (`polarity_` in `core/miyata.h`), in `LOCAMINOD`
/// order.
const POLARITY_: [f64; 20] = [
    8.1, 10.5, 11.6, 13.0, 5.5, 10.5, 12.3, 9.0, 10.4, 5.2, 4.9, 11.3, 5.7, 5.2, 8.0, 9.2, 8.6,
    5.4, 6.2, 5.9,
];

/// Miyata's volume scale (`volume_` in `core/miyata.h`).
const VOLUME_: [f64; 20] = [
    31.0, 124.0, 56.0, 54.0, 55.0, 85.0, 83.0, 3.0, 96.0, 111.0, 111.0, 119.0, 105.0, 132.0, 32.5,
    32.0, 61.0, 170.0, 136.0, 84.0,
];

/// `FFT_THRESHOLD` and `FFT_WINSIZE_P` from `mltaln.h`.
pub const FFT_THRESHOLD: i32 = 80;
pub const FFT_WINSIZE_P: i32 = 20;

/// Centre and scale a 20-element property vector, as the `fftscore` block of
/// `constants()` does.
fn normalise(raw: &[f64; 20]) -> [f64; 0x100] {
    let mut out = [0.0f64; 0x100];
    let mut v = *raw;
    let mut av = 0.0f64;
    for x in v.iter() {
        av += *x;
    }
    av /= 20.0;
    let mut sd = 0.0f64;
    for x in v.iter() {
        sd += (*x - av) * (*x - av);
    }
    sd /= 20.0;
    sd = sd.sqrt();
    for x in v.iter_mut() {
        *x -= av;
    }
    for x in v.iter_mut() {
        *x /= sd;
    }
    out[..20].copy_from_slice(&v);
    out
}

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

    /// Score matrix used by the FFT segment scan: `n_dis + offset - offsetFFT`.
    pub n_dis_fft: Vec<Vec<i32>>,
    /// Normalised polarity and volume, the two components of the sequence
    /// vector the FFT correlates.
    /// 256 entries like the C globals; only the 20 scored residues are set, and
    /// `seq_vec_5` does read index 20 (`B`), which is therefore 0.0.
    pub polarity: [f64; 0x100],
    pub volume: [f64; 0x100],
    pub fft_win_size: i32,
    pub fft_threshold: i32,
    /// `-s` / `-u`: how far a shift is allowed to pay for itself
    /// (`--allowshift` / `--unalignlevel`). Zero disables the whole mechanism.
    pub specificityconsideration: f64,
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

        // offsetFFT = (int)( 600.0 / 1000.0 * (-0) + 0.5 ) == 0
        let offset_fft = (600.0 / 1000.0 * (-0.0f64) + 0.5) as i32;
        let mut n_dis_fft = vec![vec![0i32; NALPHABETS]; NALPHABETS];
        for i in 0..20 {
            for j in 0..20 {
                n_dis_fft[i][j] = n_dis[i][j] + offset - offset_fft;
            }
        }

        let fft_win_size = if args.fft_win_size == crate::params::NOTSPECIFIED {
            FFT_WINSIZE_P
        } else {
            args.fft_win_size
        };
        let fft_threshold = if args.fft_threshold == crate::params::NOTSPECIFIED {
            FFT_THRESHOLD
        } else {
            args.fft_threshold
        };

        Constants {
            n_dis_fft,
            polarity: normalise(&POLARITY_),
            volume: normalise(&VOLUME_),
            fft_win_size,
            fft_threshold,
            specificityconsideration: args.specificityconsideration,
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

/// `dist2offset()`: how much of a bonus a shifted match may earn.
pub fn dist2offset(dist: f64, specificityconsideration: f64) -> f64 {
    let val = dist * 0.5 - specificityconsideration;
    if val > 0.0 {
        0.0
    } else {
        val
    }
}

/// `makedynamicmtx()`: raise every residue-residue score by an offset that
/// depends on how diverged the node being aligned is.
///
/// Under `--allowshift` the score matrix is no longer a constant: closely
/// related groups get a flatter matrix, so a shift costs relatively less.
pub fn makedynamicmtx(
    out: &mut [Vec<f64>],
    input: &[Vec<f64>],
    offset: f64,
    amino: &[u8; NALPHABETS],
    specificityconsideration: f64,
) {
    let offset = dist2offset(offset * 2.0, specificityconsideration);
    for i in 0..NALPHABETS {
        for j in 0..NALPHABETS {
            out[i][j] = input[i][j];
        }
    }
    if offset == 0.0 {
        return;
    }
    for i in 0..NALPHABETS {
        if amino[i] == b'-' {
            continue;
        }
        for j in 0..NALPHABETS {
            if amino[j] == b'-' {
                continue;
            }
            out[i][j] = input[i][j] + offset * 600.0;
        }
    }
}
