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

/// `locaminon` from `core/DNA.h`: the nucleotide alphabet, index -> character.
/// Note `'n'` appears twice (10 and 17); `amino_n['n']` ends up 17.
pub const LOCAMINON: &[u8; 26] = b"agctuAGCTUnNbdhkmnrsvwyx-O";

/// `locgrpn` from `core/DNA.h`: character -> group for the 6-mer distance.
const LOCGRPN: [i32; 26] = [
    0, 1, 2, 3, 3, 0, 1, 2, 3, 3, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
];

/// `FFT_WINSIZE_D` from `mltaln.h`.
pub const FFT_WINSIZE_D: i32 = 100;

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

/// `generatenuc1pam()`: the 1-PAM transition matrix of Kimura's two-parameter
/// model, with transition/transversion ratio `kimura_r`.
fn generatenuc1pam(kimura_r: i32, freq: &[f64; 4]) -> [[f64; 11]; 11] {
    let k = kimura_r as f64;
    let r = [
        [0.0, k, 1.0, 1.0],
        [k, 0.0, 1.0, 1.0],
        [1.0, 1.0, 0.0, k],
        [1.0, 1.0, k, 0.0],
    ];
    let mut mutv = [0.0f64; 4];
    let mut total = 0.0f64;
    for i in 0..4 {
        let mut tmp = 0.0f64;
        for j in 0..4 {
            tmp += r[i][j] * freq[j];
        }
        mutv[i] = tmp;
        total += tmp * freq[i];
    }
    let mut pam1 = [[0.0f64; 11]; 11];
    for i in 0..4 {
        for j in 0..4 {
            if i != j {
                pam1[i][j] = 0.01 / total * r[i][j] * freq[j];
            } else {
                pam1[i][j] = 1.0 - 0.01 / total * mutv[i];
            }
        }
    }
    pam1
}

/// `MtxmltDouble( mtx1, mtx2, 4 )`: `mtx1 <- mtx1 * mtx2`, top-left 4x4 only.
fn mtxmlt_double(mtx1: &mut [[f64; 11]; 11], mtx2: &[[f64; 11]; 11]) {
    let mut tmp = [0.0f64; 4];
    for i in 0..4 {
        tmp[..4].copy_from_slice(&mtx1[i][..4]);
        for j in 0..4 {
            let mut sum = 0.0f64;
            for k in 0..4 {
                sum += tmp[k] * mtx2[k][j];
            }
            mtx1[i][j] = sum;
        }
    }
}

/// `ambiguousscore()`: give each IUPAC ambiguity code the mean score of the
/// bases it stands for. Row `i` is filled first and then mirrored, so the
/// diagonal entries below are computed from the *unambiguous* diagonal.
fn ambiguousscore(amino_n: &[i32; CHARSIZE], n_dis: &mut [[i32; NALPHABETS]; NALPHABETS]) {
    let ix = |c: u8| amino_n[c as usize] as usize;
    let (a, g, c, t) = (ix(b'a'), ix(b'g'), ix(b'c'), ix(b't'));
    let two: [(u8, usize, usize); 6] = [
        (b'r', a, g),
        (b'y', c, t),
        (b'k', g, t),
        (b'm', a, c),
        (b's', g, c),
        (b'w', a, t),
    ];
    let three: [(u8, usize, usize, usize); 4] = [
        (b'b', c, g, t),
        (b'd', a, g, t),
        (b'h', a, c, t),
        (b'v', a, c, g),
    ];
    for i in 0..NALPHABETS {
        for &(code, p, q) in two.iter() {
            let k = ix(code);
            n_dis[i][k] = shishagonyuu(0.5 * (n_dis[p][i] + n_dis[q][i]) as f64);
        }
        for &(code, p, q, r) in three.iter() {
            let k = ix(code);
            n_dis[i][k] =
                shishagonyuu((1.0 / 3.0) * (n_dis[p][i] + n_dis[q][i] + n_dis[r][i]) as f64);
        }
        for &(code, _, _) in two.iter() {
            let k = ix(code);
            n_dis[k][i] = n_dis[i][k];
        }
        for &(code, _, _, _) in three.iter() {
            let k = ix(code);
            n_dis[k][i] = n_dis[i][k];
        }
    }
    for &(code, p, q) in two.iter() {
        let k = ix(code);
        n_dis[k][k] = shishagonyuu(0.5 * (n_dis[p][p] + n_dis[q][q]) as f64);
    }
    for &(code, p, q, r) in three.iter() {
        let k = ix(code);
        n_dis[k][k] = shishagonyuu((1.0 / 3.0) * (n_dis[p][p] + n_dis[q][q] + n_dis[r][r]) as f64);
    }
}

/// `nscore()`: `--anysymbol` for nucleotides, scoring `n` against everything.
fn nscore(amino_n: &[i32; CHARSIZE], n_dis: &mut [[i32; NALPHABETS]; NALPHABETS]) {
    let n = amino_n[b'n' as usize] as usize;
    for i in 0..NALPHABETS {
        n_dis[i][n] = shishagonyuu(0.25 * n_dis[i][i] as f64);
        n_dis[n][i] = n_dis[i][n];
    }
    n_dis[n][n] =
        shishagonyuu(0.25 * (n_dis[0][0] + n_dis[1][1] + n_dis[2][2] + n_dis[3][3]) as f64);
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
    /// `'p'` or `'d'`: which branch of `constants()` built this.
    pub dorp: u8,
    /// `-s` / `-u`: how far a shift is allowed to pay for itself
    /// (`--allowshift` / `--unalignlevel`). Zero disables the whole mechanism.
    pub specificityconsideration: f64,
}

// Defaults from blosum.h / constants.c for the protein BLOSUM branch.
const DEFAULTGOP_B: i32 = -1530;
const DEFAULTGEP_B: i32 = 0;
const DEFAULTOFS_B: i32 = -123;
// Defaults from DNA.h for the nucleotide branch.
const DEFAULTGOP_N: i32 = -1530;
const DEFAULTGEP_N: i32 = 0;
const DEFAULTOFS_N: i32 = -369;
const DEFAULTPAMN: i32 = 200;
const CONSWEIGHT_MULTI: f64 = 1.0;

impl Constants {
    /// The protein branch (`dorp == 'p'`), which is what every caller that does
    /// not read its sequences first still wants.
    pub fn new(args: &Args) -> Constants {
        Constants::new_dorp(args, b'p')
    }

    /// `constants()`, for either sequence type.
    pub fn new_dorp(args: &Args, dorp: u8) -> Constants {
        let trywarp = !(args.penalty_shift_factor >= 10.0);

        let mut amino = [0u8; NALPHABETS];
        let mut amino_grp = [0i32; 0x100];
        let mut n_dis = [[0i32; NALPHABETS]; NALPHABETS];
        let penalty;
        let penalty_ex;
        let penalty_dist;
        let penalty_shift;
        let offset;
        // `nscoredalphabets`: 20 for protein, 10 for DNA. Only the FFT score
        // matrix below is built over it.
        let nscored;

        if dorp == b'd' {
            nscored = 10;
            let mut ppenalty = args.ppenalty;
            let mut ppenalty_dist = args.ppenalty_dist;
            let mut ppenalty_ex = args.ppenalty_ex;
            let mut poffset = args.poffset;
            if ppenalty == crate::params::NOTSPECIFIED {
                ppenalty = DEFAULTGOP_N;
            }
            if ppenalty_dist == crate::params::NOTSPECIFIED {
                ppenalty_dist = ppenalty;
            }
            if ppenalty_ex == crate::params::NOTSPECIFIED {
                ppenalty_ex = DEFAULTGEP_N;
            }
            if poffset == crate::params::NOTSPECIFIED {
                poffset = DEFAULTOFS_N;
            }
            let pam_n = if args.pam_n == crate::params::NOTSPECIFIED {
                DEFAULTPAMN
            } else {
                args.pam_n
            };
            let kimura_r = if args.kimura_r == crate::params::NOTSPECIFIED {
                2
            } else {
                args.kimura_r
            };

            penalty = (3.0 * 600.0 / 1000.0 * ppenalty as f64 + 0.5) as i32;
            penalty_dist = (3.0 * 600.0 / 1000.0 * ppenalty_dist as f64 + 0.5) as i32;
            penalty_shift = (args.penalty_shift_factor * penalty as f64) as i32;
            penalty_ex = (3.0 * 600.0 / 1000.0 * ppenalty_ex as f64 + 0.5) as i32;
            offset = (1.0 * 600.0 / 1000.0 * poffset as f64 + 0.5) as i32;

            amino.copy_from_slice(LOCAMINON);
            let mut amino_n = [-1i32; CHARSIZE];
            for i in 0..NALPHABETS {
                amino_n[amino[i] as usize] = i as i32;
            }
            for i in 0..NALPHABETS {
                amino_grp[amino[i] as usize] = LOCGRPN[i];
            }

            // fmodel != 1 is refused, so the frequencies are flat.
            let freq = [0.25f64; 4];

            // `pamx` is an 11x11 calloc'd matrix; only [0,10)x[0,10) is read
            // back, and the halves that neither branch below writes stay zero.
            let mut pamx = [[0.0f64; 11]; 11];
            let pam1 = generatenuc1pam(kimura_r, &freq);
            // MtxuntDouble( pamx, 4 )
            for i in 0..4 {
                pamx[i][i] = 1.0;
            }
            for _ in 0..pam_n {
                mtxmlt_double(&mut pamx, &pam1);
            }
            for i in 0..4 {
                for j in 0..4 {
                    pamx[i][j] /= freq[j];
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    if pamx[i][j] == 0.0 {
                        pamx[i][j] = 0.00001;
                    }
                    pamx[i][j] = pamx[i][j].log10() * 1000.0;
                }
            }
            let mut average = 0.0f64;
            for i in 0..4 {
                for j in 0..4 {
                    average += pamx[i][j] * freq[i] * freq[j];
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    pamx[i][j] -= average;
                }
            }
            let mut average = 0.0f64;
            for i in 0..4 {
                average += pamx[i][i] / 4.0;
            }
            for i in 0..4 {
                for j in 0..4 {
                    pamx[i][j] *= 600.0 / average;
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    pamx[i][j] -= offset as f64;
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    pamx[i][j] = shishagonyuu(pamx[i][j]) as f64;
                }
            }

            // 'u' scores as 't', and the upper-case half repeats the lower.
            // The order matters: at i == 4 both assignments read cells this
            // same loop has already overwritten.
            for i in 0..5 {
                pamx[4][i] = pamx[3][i];
                pamx[i][4] = pamx[i][3];
            }
            for i in 5..10 {
                for j in 5..10 {
                    pamx[i][j] = pamx[i - 5][j - 5];
                }
            }

            for i in 0..10 {
                for j in 0..10 {
                    n_dis[i][j] = shishagonyuu(pamx[i][j]);
                }
            }
            ambiguousscore(&amino_n, &mut n_dis);
            if args.nwildcard {
                nscore(&amino_n, &mut n_dis);
            }
        } else {
            nscored = 20;
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

            penalty = (600.0 / 1000.0 * ppenalty as f64 + 0.5) as i32;
            penalty_dist = (600.0 / 1000.0 * ppenalty_dist as f64 + 0.5) as i32;
            penalty_shift = (args.penalty_shift_factor * penalty as f64) as i32;
            penalty_ex = (600.0 / 1000.0 * ppenalty_ex as f64 + 0.5) as i32;
            offset = (600.0 / 1000.0 * poffset as f64 + 0.5) as i32;

            let mut n_distmp = [[0.0f64; 20]; 20];
            let mut freq = [0.0f64; 20];
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

            for i in 0..20 {
                for j in 0..20 {
                    n_dis[i][j] = n_distmp[i][j] as i32;
                }
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
        for i in 0..nscored {
            for j in 0..nscored {
                n_dis_fft[i][j] = n_dis[i][j] + offset - offset_fft;
            }
        }

        let fft_win_size = if args.fft_win_size == crate::params::NOTSPECIFIED {
            if dorp == b'd' {
                FFT_WINSIZE_D
            } else {
                FFT_WINSIZE_P
            }
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
            dorp,
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

    /// The same, for the nucleotide branch: Kimura's model raised to the 200th
    /// power, taken through `log10`, rescaled and rounded. Nothing about that
    /// chain is obvious enough to trust without pinning it.
    #[test]
    fn matches_c_constants_dna() {
        let golden = include_str!("../tests/golden/constants_dna.txt");
        let args = Args {
            nblosum: 62,
            scoremtx: 1,
            penalty_shift_factor: 100.0,
            ..Args::default()
        };
        let cn = Constants::new_dorp(&args, b'd');

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
        let mut kv2 = || -> (String, i64) {
            let l = lines.next().expect("golden line");
            let mut it = l.split_whitespace();
            (
                it.next().unwrap().to_string(),
                it.next().unwrap().parse().unwrap(),
            )
        };
        assert_eq!(kv2(), ("fftWinSize".into(), cn.fft_win_size as i64));
        assert_eq!(kv2(), ("fftThreshold".into(), cn.fft_threshold as i64));
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
