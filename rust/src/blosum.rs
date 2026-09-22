//! BLOSUM matrices, transliterated from `core/blosum.c`.
//!
//! Every built-in table `BLOSUMmtx()` can select is reproduced; the values and
//! their order are copied verbatim so that the floating-point arithmetic in
//! `constants()` reproduces the C result bit for bit. `--aamatrix` (`n == -1`,
//! a matrix read from a file) is refused rather than ported.

/// Amino-acid ordering used by the BLOSUM tables (`locaminod` in blosum.c).
pub const LOCAMINOD: &[u8; 26] = b"ARNDCQEGHILKMFPSTWYVBZX.-J";

/// Residue group membership (`locgrpd` in blosum.c), parallel to `LOCAMINOD`.
pub const LOCGRPD: [i32; 26] = [
    0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2, 6, 6, 6, 1,
];

/// Background amino-acid frequencies (`freqd` in blosum.c).
pub const FREQD: [f64; 20] = [
    0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062, 0.074, 0.023, 0.052, 0.091, 0.059, 0.024,
    0.040, 0.051, 0.069, 0.059, 0.014, 0.032, 0.066,
];

/// BLOSUM62, lower triangle in row-major order (`tmpmtx62` in blosum.c).
pub const TMPMTX62: [f64; 210] = [
    5.893685, -2.120252, 8.210189, -2.296072, -0.659672, 8.479856, -2.630151, -2.408668, 1.907550,
    8.661363, -0.612761, -5.083814, -3.989626, -5.189966, 12.873172, -1.206025, 1.474162, 0.002529,
    -0.470069, -4.352838, 7.927704, -1.295821, -0.173087, -0.402015, 2.265459, -5.418729, 2.781955,
    7.354247, 0.239392, -3.456163, -0.634136, -1.970281, -3.750621, -2.677743, -3.165266, 8.344902,
    -2.437724, -0.374792, 0.867735, -1.678363, -4.481724, 0.672051, -0.176497, -3.061315,
    11.266586, -1.982718, -4.485360, -4.825558, -4.681732, -1.841495, -4.154454, -4.791538,
    -5.587336, -4.847345, 5.997760, -2.196882, -3.231860, -5.068375, -5.408471, -1.916207,
    -3.200863, -4.269723, -5.440437, -4.180099, 2.282412, 5.774148, -1.101017, 3.163105, -0.268534,
    -1.052724, -4.554510, 1.908859, 1.163010, -2.291924, -1.081539, -4.005209, -3.670219, 6.756827,
    -1.402897, -2.050705, -3.226290, -4.587785, -2.129758, -0.631437, -2.997038, -4.014898,
    -2.326896, 1.690191, 2.987638, -2.032119, 8.088951, -3.315080, -4.179521, -4.491005, -5.225795,
    -3.563219, -4.746598, -4.788639, -4.661029, -1.851231, -0.241317, 0.622170, -4.618016,
    0.018880, 9.069126, -1.221394, -3.162863, -3.000581, -2.220163, -4.192770, -1.922917,
    -1.674258, -3.200320, -3.241363, -4.135001, -4.290107, -1.520445, -3.714633, -5.395930,
    11.046892, 1.673639, -1.147170, 0.901353, -0.391548, -1.312485, -0.151708, -0.220375,
    -0.438748, -1.322366, -3.522266, -3.663923, -0.305170, -2.221304, -3.553533, -1.213470,
    5.826527, -0.068042, -1.683495, -0.069138, -1.576054, -1.299983, -1.012997, -1.294878,
    -2.363065, -2.528844, -1.076382, -1.796229, -1.004336, -0.999449, -3.161436, -1.612919,
    2.071710, 6.817956, -3.790328, -4.019108, -5.543911, -6.321502, -3.456164, -2.919725,
    -4.253197, -3.737232, -3.513238, -3.870811, -2.447829, -4.434676, -2.137255, 1.376341,
    -5.481260, -4.127804, -3.643382, 15.756041, -2.646022, -2.540799, -3.122641, -4.597428,
    -3.610671, -2.131601, -3.030688, -4.559647, 2.538948, -1.997058, -1.593097, -2.730047,
    -1.492308, 4.408690, -4.379667, -2.528713, -2.408996, 3.231335, 9.892544, -0.284140, -3.753871,
    -4.314525, -4.713963, -1.211518, -3.297575, -3.663425, -4.708118, -4.676220, 3.820569,
    1.182672, -3.393535, 1.030861, -1.273542, -3.523054, -2.469318, -0.083276, -4.251392,
    -1.811267, 5.653391,
];

/// Fills `matrix` (20x20) and `freq` from the requested BLOSUM table, and sets
/// `amino` / `amino_grp` exactly as `BLOSUMmtx()` does.
///
/// Only `n == 62` is supported here; the other tables in blosum.c are selected
/// by options this port does not accept, and are rejected loudly rather than
/// silently approximated.
/// `tmpmtx30` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX30: [f64; 210] = [
    4.0, -1.0, 8.0, 0.0, -2.0, 8.0, 0.0, -1.0, 1.0, 9.0, -3.0, -2.0, -1.0, -3.0, 17.0, 1.0, 3.0,
    -1.0, -1.0, -2.0, 8.0, 0.0, -1.0, -1.0, 1.0, 1.0, 2.0, 6.0, 0.0, -2.0, 0.0, -1.0, -4.0, -2.0,
    -2.0, 8.0, -2.0, -1.0, -1.0, -2.0, -5.0, 0.0, 0.0, -3.0, 14.0, 0.0, -3.0, 0.0, -4.0, -2.0,
    -2.0, -3.0, -1.0, -2.0, 6.0, -1.0, -2.0, -2.0, -1.0, 0.0, -2.0, -1.0, -2.0, -1.0, 2.0, 4.0,
    0.0, 1.0, 0.0, 0.0, -3.0, 0.0, 2.0, -1.0, -2.0, -2.0, -2.0, 4.0, 1.0, 0.0, 0.0, -3.0, -2.0,
    -1.0, -1.0, -2.0, 2.0, 1.0, 2.0, 2.0, 6.0, -2.0, -1.0, -1.0, -5.0, -3.0, -3.0, -4.0, -3.0,
    -3.0, 0.0, 2.0, -1.0, -2.0, 10.0, -1.0, -1.0, -3.0, -1.0, -3.0, 0.0, 1.0, -1.0, 1.0, -3.0,
    -3.0, 1.0, -4.0, -4.0, 11.0, 1.0, -1.0, 0.0, 0.0, -2.0, -1.0, 0.0, 0.0, -1.0, -1.0, -2.0, 0.0,
    -2.0, -1.0, -1.0, 4.0, 1.0, -3.0, 1.0, -1.0, -2.0, 0.0, -2.0, -2.0, -2.0, 0.0, 0.0, -1.0, 0.0,
    -2.0, 0.0, 2.0, 5.0, -5.0, 0.0, -7.0, -4.0, -2.0, -1.0, -1.0, 1.0, -5.0, -3.0, -2.0, -2.0,
    -3.0, 1.0, -3.0, -3.0, -5.0, 20.0, -4.0, 0.0, -4.0, -1.0, -6.0, -1.0, -2.0, -3.0, 0.0, -1.0,
    3.0, -1.0, -1.0, 3.0, -2.0, -2.0, -1.0, 5.0, 9.0, 1.0, -1.0, -2.0, -2.0, -2.0, -3.0, -3.0,
    -3.0, -3.0, 4.0, 1.0, -2.0, 0.0, 1.0, -4.0, -1.0, 1.0, -3.0, 1.0, 5.0,
];

/// `tmpmtx45` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX45: [f64; 210] = [
    5.0, -2.0, 7.0, -1.0, 0.0, 6.0, -2.0, -1.0, 2.0, 7.0, -1.0, -3.0, -2.0, -3.0, 12.0, -1.0, 1.0,
    0.0, 0.0, -3.0, 6.0, -1.0, 0.0, 0.0, 2.0, -3.0, 2.0, 6.0, 0.0, -2.0, 0.0, -1.0, -3.0, -2.0,
    -2.0, 7.0, -2.0, 0.0, 1.0, 0.0, -3.0, 1.0, 0.0, -2.0, 10.0, -1.0, -3.0, -2.0, -4.0, -3.0, -2.0,
    -3.0, -4.0, -3.0, 5.0, -1.0, -2.0, -3.0, -3.0, -2.0, -2.0, -2.0, -3.0, -2.0, 2.0, 5.0, -1.0,
    3.0, 0.0, 0.0, -3.0, 1.0, 1.0, -2.0, -1.0, -3.0, -3.0, 5.0, -1.0, -1.0, -2.0, -3.0, -2.0, 0.0,
    -2.0, -2.0, 0.0, 2.0, 2.0, -1.0, 6.0, -2.0, -2.0, -2.0, -4.0, -2.0, -4.0, -3.0, -3.0, -2.0,
    0.0, 1.0, -3.0, 0.0, 8.0, -1.0, -2.0, -2.0, -1.0, -4.0, -1.0, 0.0, -2.0, -2.0, -2.0, -3.0,
    -1.0, -2.0, -3.0, 9.0, 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, -2.0, -3.0, -1.0, -2.0,
    -2.0, -1.0, 4.0, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0, -1.0,
    -1.0, -1.0, 2.0, 5.0, -2.0, -2.0, -4.0, -4.0, -5.0, -2.0, -3.0, -2.0, -3.0, -2.0, -2.0, -2.0,
    -2.0, 1.0, -3.0, -4.0, -3.0, 15.0, -2.0, -1.0, -2.0, -2.0, -3.0, -1.0, -2.0, -3.0, 2.0, 0.0,
    0.0, -1.0, 0.0, 3.0, -3.0, -2.0, -1.0, 3.0, 8.0, 0.0, -2.0, -3.0, -3.0, -1.0, -3.0, -3.0, -3.0,
    -3.0, 3.0, 1.0, -2.0, 1.0, 0.0, -3.0, -1.0, 0.0, -3.0, -1.0, 5.0,
];

/// `tmpmtx50` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX50: [f64; 210] = [
    5.0, -2.0, 7.0, -1.0, -1.0, 7.0, -2.0, -2.0, 2.0, 8.0, -1.0, -4.0, -2.0, -4.0, 13.0, -1.0, 1.0,
    0.0, 0.0, -3.0, 7.0, -1.0, 0.0, 0.0, 2.0, -3.0, 2.0, 6.0, 0.0, -3.0, 0.0, -1.0, -3.0, -2.0,
    -3.0, 8.0, -2.0, 0.0, 1.0, -1.0, -3.0, 1.0, 0.0, -2.0, 10.0, -1.0, -4.0, -3.0, -4.0, -2.0,
    -3.0, -4.0, -4.0, -4.0, 5.0, -2.0, -3.0, -4.0, -4.0, -2.0, -2.0, -3.0, -4.0, -3.0, 2.0, 5.0,
    -1.0, 3.0, 0.0, -1.0, -3.0, 2.0, 1.0, -2.0, 0.0, -3.0, -3.0, 6.0, -1.0, -2.0, -2.0, -4.0, -2.0,
    0.0, -2.0, -3.0, -1.0, 2.0, 3.0, -2.0, 7.0, -3.0, -3.0, -4.0, -5.0, -2.0, -4.0, -3.0, -4.0,
    -1.0, 0.0, 1.0, -4.0, 0.0, 8.0, -1.0, -3.0, -2.0, -1.0, -4.0, -1.0, -1.0, -2.0, -2.0, -3.0,
    -4.0, -1.0, -3.0, -4.0, 10.0, 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, -1.0, 0.0, -1.0, -3.0, -3.0, 0.0,
    -2.0, -3.0, -1.0, 5.0, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, -1.0, -2.0, -2.0, -1.0, -1.0, -1.0,
    -1.0, -2.0, -1.0, 2.0, 5.0, -3.0, -3.0, -4.0, -5.0, -5.0, -1.0, -3.0, -3.0, -3.0, -3.0, -2.0,
    -3.0, -1.0, 1.0, -4.0, -4.0, -3.0, 15.0, -2.0, -1.0, -2.0, -3.0, -3.0, -1.0, -2.0, -3.0, 2.0,
    -1.0, -1.0, -2.0, 0.0, 4.0, -3.0, -2.0, -2.0, 2.0, 8.0, 0.0, -3.0, -3.0, -4.0, -1.0, -3.0,
    -3.0, -4.0, -4.0, 4.0, 1.0, -3.0, 1.0, -1.0, -3.0, -2.0, 0.0, -3.0, -1.0, 5.0,
];

/// `tmpmtx80` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX80: [f64; 210] = [
    7.0, -3.0, 9.0, -3.0, -1.0, 9.0, -3.0, -3.0, 2.0, 10.0, -1.0, -6.0, -5.0, -7.0, 13.0, -2.0,
    1.0, 0.0, -1.0, -5.0, 9.0, -2.0, -1.0, -1.0, 2.0, -7.0, 3.0, 8.0, 0.0, -4.0, -1.0, -3.0, -6.0,
    -4.0, -4.0, 9.0, -3.0, 0.0, 1.0, -2.0, -7.0, 1.0, 0.0, -4.0, 12.0, -3.0, -5.0, -6.0, -7.0,
    -2.0, -5.0, -6.0, -7.0, -6.0, 7.0, -3.0, -4.0, -6.0, -7.0, -3.0, -4.0, -6.0, -7.0, -5.0, 2.0,
    6.0, -1.0, 3.0, 0.0, -2.0, -6.0, 2.0, 1.0, -3.0, -1.0, -5.0, -4.0, 8.0, -2.0, -3.0, -4.0, -6.0,
    -3.0, -1.0, -4.0, -5.0, -4.0, 2.0, 3.0, -3.0, 9.0, -4.0, -5.0, -6.0, -6.0, -4.0, -5.0, -6.0,
    -6.0, -2.0, -1.0, 0.0, -5.0, 0.0, 10.0, -1.0, -3.0, -4.0, -3.0, -6.0, -3.0, -2.0, -5.0, -4.0,
    -5.0, -5.0, -2.0, -4.0, -6.0, 12.0, 2.0, -2.0, 1.0, -1.0, -2.0, -1.0, -1.0, -1.0, -2.0, -4.0,
    -4.0, -1.0, -3.0, -4.0, -2.0, 7.0, 0.0, -2.0, 0.0, -2.0, -2.0, -1.0, -2.0, -3.0, -3.0, -2.0,
    -3.0, -1.0, -1.0, -4.0, -3.0, 2.0, 8.0, -5.0, -5.0, -7.0, -8.0, -5.0, -4.0, -6.0, -6.0, -4.0,
    -5.0, -4.0, -6.0, -3.0, 0.0, -7.0, -6.0, -5.0, 16.0, -4.0, -4.0, -4.0, -6.0, -5.0, -3.0, -5.0,
    -6.0, 3.0, -3.0, -2.0, -4.0, -3.0, 4.0, -6.0, -3.0, -3.0, 3.0, 11.0, -1.0, -4.0, -5.0, -6.0,
    -2.0, -4.0, -4.0, -6.0, -5.0, 4.0, 1.0, -4.0, 1.0, -2.0, -4.0, -3.0, 0.0, -5.0, -3.0, 7.0,
];

/// `tmpmtx90` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX90: [f64; 210] = [
    5.0, -2.0, 6.0, -2.0, -1.0, 7.0, -3.0, -3.0, 1.0, 7.0, -1.0, -5.0, -4.0, -5.0, 9.0, -1.0, 1.0,
    0.0, -1.0, -4.0, 7.0, -1.0, -1.0, -1.0, 1.0, -6.0, 2.0, 6.0, 0.0, -3.0, -1.0, -2.0, -4.0, -3.0,
    -3.0, 6.0, -2.0, 0.0, 0.0, -2.0, -5.0, 1.0, -1.0, -3.0, 8.0, -2.0, -4.0, -4.0, -5.0, -2.0,
    -4.0, -4.0, -5.0, -4.0, 5.0, -2.0, -3.0, -4.0, -5.0, -2.0, -3.0, -4.0, -5.0, -4.0, 1.0, 5.0,
    -1.0, 2.0, 0.0, -1.0, -4.0, 1.0, 0.0, -2.0, -1.0, -4.0, -3.0, 6.0, -2.0, -2.0, -3.0, -4.0,
    -2.0, 0.0, -3.0, -4.0, -3.0, 1.0, 2.0, -2.0, 7.0, -3.0, -4.0, -4.0, -5.0, -3.0, -4.0, -5.0,
    -5.0, -2.0, -1.0, 0.0, -4.0, -1.0, 7.0, -1.0, -3.0, -3.0, -3.0, -4.0, -2.0, -2.0, -3.0, -3.0,
    -4.0, -4.0, -2.0, -3.0, -4.0, 8.0, 1.0, -1.0, 0.0, -1.0, -2.0, -1.0, -1.0, -1.0, -2.0, -3.0,
    -3.0, -1.0, -2.0, -3.0, -2.0, 5.0, 0.0, -2.0, 0.0, -2.0, -2.0, -1.0, -1.0, -3.0, -2.0, -1.0,
    -2.0, -1.0, -1.0, -3.0, -2.0, 1.0, 6.0, -4.0, -4.0, -5.0, -6.0, -4.0, -3.0, -5.0, -4.0, -3.0,
    -4.0, -3.0, -5.0, -2.0, 0.0, -5.0, -4.0, -4.0, 11.0, -3.0, -3.0, -3.0, -4.0, -4.0, -3.0, -4.0,
    -5.0, 1.0, -2.0, -2.0, -3.0, -2.0, 3.0, -4.0, -3.0, -2.0, 2.0, 8.0, -1.0, -3.0, -4.0, -5.0,
    -2.0, -3.0, -3.0, -5.0, -4.0, 3.0, 0.0, -3.0, 0.0, -2.0, -3.0, -2.0, -1.0, -3.0, -3.0, 5.0,
];

/// `tmpmtx100` from `core/blosum.c`, lower triangle in row-major order.
pub const TMPMTX100: [f64; 210] = [
    8.0, -3.0, 10.0, -4.0, -2.0, 11.0, -5.0, -5.0, 1.0, 10.0, -2.0, -8.0, -5.0, -8.0, 14.0, -2.0,
    0.0, -1.0, -2.0, -7.0, 11.0, -3.0, -2.0, -2.0, 2.0, -9.0, 2.0, 10.0, -1.0, -6.0, -2.0, -4.0,
    -7.0, -5.0, -6.0, 9.0, -4.0, -1.0, 0.0, -3.0, -8.0, 1.0, -2.0, -6.0, 13.0, -4.0, -7.0, -7.0,
    -8.0, -3.0, -6.0, -7.0, -9.0, -7.0, 8.0, -4.0, -6.0, -7.0, -8.0, -5.0, -5.0, -7.0, -8.0, -6.0,
    2.0, 8.0, -2.0, 3.0, -1.0, -3.0, -8.0, 2.0, 0.0, -5.0, -3.0, -6.0, -6.0, 10.0, -3.0, -4.0,
    -5.0, -8.0, -4.0, -2.0, -5.0, -7.0, -5.0, 1.0, 3.0, -4.0, 12.0, -5.0, -6.0, -7.0, -8.0, -4.0,
    -6.0, -8.0, -8.0, -4.0, -2.0, 0.0, -6.0, -1.0, 11.0, -2.0, -5.0, -5.0, -5.0, -8.0, -4.0, -4.0,
    -6.0, -5.0, -7.0, -7.0, -3.0, -5.0, -7.0, 12.0, 1.0, -3.0, 0.0, -2.0, -3.0, -2.0, -2.0, -2.0,
    -3.0, -5.0, -6.0, -2.0, -4.0, -5.0, -3.0, 9.0, -1.0, -3.0, -1.0, -4.0, -3.0, -3.0, -3.0, -5.0,
    -4.0, -3.0, -4.0, -3.0, -2.0, -5.0, -4.0, 2.0, 9.0, -6.0, -7.0, -8.0, -10.0, -7.0, -5.0, -8.0,
    -7.0, -5.0, -6.0, -5.0, -8.0, -4.0, 0.0, -8.0, -7.0, -7.0, 17.0, -5.0, -5.0, -5.0, -7.0, -6.0,
    -4.0, -7.0, -8.0, 1.0, -4.0, -4.0, -5.0, -5.0, 4.0, -7.0, -5.0, -5.0, 2.0, 12.0, -2.0, -6.0,
    -7.0, -8.0, -3.0, -5.0, -5.0, -8.0, -7.0, 4.0, 0.0, -5.0, 0.0, -3.0, -6.0, -4.0, -1.0, -5.0,
    -5.0, 8.0,
];

/// `loadaamtx()`: the user matrix `--aamatrix` puts in `_aamtx`.
///
/// Returns the 210-entry lower triangle in the canonical
/// `ARNDCQEGHILKMFPSTWYV` order followed by 20 background frequencies (all
/// `-1.0` unless the file carried a `frequency` line), and whether the caller
/// should rescale -- the file can turn that off with a `norescale` line.
///
/// The file's own residue order comes from its first non-`#` line and may be
/// anything; the matrix is permuted into ours. Only the lower triangle of the
/// file is read, so a value is looked up at (max, min) of the two mapped
/// indices.
pub fn loadaamtx(dorp: u8) -> (Vec<f64>, bool) {
    const AAORDER: &[u8; 20] = b"ARNDCQEGHILKMFPSTWYV";

    if dorp != b'p' {
        eprintln!("User-defined matrix is not supported for DNA");
        std::process::exit(1);
    }
    let text = match std::fs::read("_aamtx") {
        Ok(t) => t,
        Err(_) => {
            eprintln!("Cannot open the _aamtx file");
            std::process::exit(1);
        }
    };
    // `fgets` keeps the newline; splitting on it and dropping the trailing
    // empty piece gives the same sequence of lines.
    let mut lines: Vec<&[u8]> = text.split(|&c| c == b'\n').collect();
    if lines.last().map(|l| l.is_empty()).unwrap_or(false) {
        lines.pop();
    }
    let mut li = 0usize;

    // The header: the first line that does not start with '#'.
    let mut inorder: Vec<u8> = Vec::new();
    while li < lines.len() {
        let l = lines[li];
        li += 1;
        if l.first() != Some(&b'#') {
            inorder = l
                .iter()
                .filter(|c| c.is_ascii_alphabetic())
                .map(|c| c.to_ascii_uppercase())
                .collect();
            break;
        }
    }
    inorder.truncate(20);

    let mut map = [0usize; 20];
    for i in 0..20 {
        match inorder.iter().position(|&c| c == AAORDER[i]) {
            Some(p) => map[i] = p,
            None => {
                eprintln!("{}: not found in the first 20 letters.", AAORDER[i] as char);
                showaamtxexample();
            }
        }
    }

    // The lower triangle, in the file's own order.
    let mut raw = [[0.0f64; 20]; 21];
    let mut i = 0usize;
    while li < lines.len() {
        let line = lines[li];
        li += 1;
        if line.first() == Some(&b'#') {
            continue;
        }
        let mut p = 0usize;
        for j in 0..=i {
            while p < line.len() && !line[p].is_ascii_digit() && line[p] != b'-' && line[p] != b'.'
            {
                p += 1;
            }
            raw[i][j] = atof_at(line, p);
            match line[p..].iter().position(|&c| c == b' ') {
                Some(q) => p += q,
                None => {
                    if j < i {
                        showaamtxexample();
                    }
                    p = line.len();
                }
            }
        }
        i += 1;
        if i > 19 {
            break;
        }
    }

    // Trailing keywords. `norescale` turns the rescaling off; `frequency`
    // overrides the background frequencies. Either ends the scan.
    let mut rescale = true;
    for v in raw[20].iter_mut() {
        *v = -1.0;
    }
    while li < lines.len() {
        let line = lines[li];
        li += 1;
        let key: &[u8] = line
            .split(|c| c.is_ascii_whitespace())
            .find(|t| !t.is_empty())
            .unwrap_or(b"");
        if key == b"norescale" {
            eprintln!("no rescale");
            rescale = false;
            break;
        } else if key == b"frequency" {
            let mut p = 0usize;
            for j in 0..20 {
                while p < line.len()
                    && !line[p].is_ascii_digit()
                    && line[p] != b'-'
                    && line[p] != b'.'
                {
                    p += 1;
                }
                raw[20][j] = atof_at(line, p);
                match line[p..].iter().position(|&c| c == b' ') {
                    Some(q) => p += q,
                    None => {
                        if j < 19 {
                            showaamtxexample();
                        }
                        p = line.len();
                    }
                }
            }
            break;
        }
    }

    let mut val = vec![0.0f64; 420];
    let mut k = 0usize;
    for i in 0..20 {
        for j in 0..=i {
            let (ii, jj) = if i != j {
                (map[i].max(map[j]), map[i].min(map[j]))
            } else {
                (map[i], map[i])
            };
            val[k] = raw[ii][jj];
            k += 1;
        }
    }
    for i in 0..20 {
        val[400 + i] = raw[20][map[i]];
    }

    eprintln!("inorder = {}", String::from_utf8_lossy(&inorder));
    (val, rescale)
}

/// `atof()` starting at `p`, stopping wherever `strtod` would.
fn atof_at(line: &[u8], p: usize) -> f64 {
    if p >= line.len() {
        return 0.0;
    }
    let rest = &line[p..];
    let mut e = 0usize;
    if rest[0] == b'-' || rest[0] == b'+' {
        e = 1;
    }
    while e < rest.len() && rest[e].is_ascii_digit() {
        e += 1;
    }
    if e < rest.len() && rest[e] == b'.' {
        e += 1;
        while e < rest.len() && rest[e].is_ascii_digit() {
            e += 1;
        }
    }
    if e < rest.len() && (rest[e] == b'e' || rest[e] == b'E') {
        let mut t = e + 1;
        if t < rest.len() && (rest[t] == b'-' || rest[t] == b'+') {
            t += 1;
        }
        let ds = t;
        while t < rest.len() && rest[t].is_ascii_digit() {
            t += 1;
        }
        if t > ds {
            e = t;
        }
    }
    std::str::from_utf8(&rest[..e])
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// `showaamtxexample()`: the C prints the expected file format and exits.
fn showaamtxexample() -> ! {
    eprintln!("Format error in aa matrix");
    eprintln!("# Example:");
    eprintln!("# comment");
    eprintln!("   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V");
    eprintln!("A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0");
    eprintln!("R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3");
    eprintln!("...");
    eprintln!("V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4");
    eprintln!("frequency 0.07 0.05 0.04 0.05 0.02 .. ");
    eprintln!("# Example end");
    eprintln!("Only the lower half is loaded");
    eprintln!("The last line (frequency) is optional.");
    std::process::exit(1);
}

pub fn blosum_mtx(
    n: i32,
    dorp: u8,
    matrix: &mut [[f64; 20]; 20],
    freq: &mut [f64; 20],
    amino: &mut [u8],
    amino_grp: &mut [i32; 0x100],
) -> bool {
    // `BLOSUMmtx()`'s dispatch. `n == -1` is `--aamatrix`: the table comes
    // from `_aamtx`, which can also carry the background frequencies and turn
    // the rescaling off.
    let loaded;
    let mut rescale = true;
    let tmpmtx: &[f64] = match n {
        30 => &TMPMTX30,
        45 => &TMPMTX45,
        50 => &TMPMTX50,
        62 => &TMPMTX62,
        80 => &TMPMTX80,
        90 => &TMPMTX90,
        100 => &TMPMTX100,
        -1 => {
            let (v, r) = loadaamtx(dorp);
            rescale = r;
            loaded = v;
            &loaded
        }
        _ => {
            eprintln!("blosum {} ?", n);
            std::process::exit(1);
        }
    };

    let mut count = 0usize;
    for i in 0..20 {
        for j in 0..=i {
            let v = tmpmtx[count];
            count += 1;
            matrix[i][j] = v;
            matrix[j][i] = v;
        }
    }
    if n == -1 && tmpmtx[400] != -1.0 {
        freq.copy_from_slice(&tmpmtx[400..420]);
        let mut av = 0.0f64;
        for v in freq.iter() {
            av += *v;
        }
        for v in freq.iter_mut() {
            *v /= av;
        }
    } else {
        freq.copy_from_slice(&FREQD);
    }

    for i in 0..26 {
        amino[i] = LOCAMINOD[i];
    }
    for i in 0..26 {
        amino_grp[amino[i] as usize] = LOCGRPD[i];
    }
    rescale
}

#[cfg(test)]
mod aamtx_tests {
    use super::*;
    use std::io::Write as _;

    /// `loadaamtx` reads `_aamtx` from the working directory, so each case
    /// runs in its own temporary one.
    fn in_dir<T>(name: &str, contents: &str, f: impl FnOnce() -> T) -> T {
        // Serialise: the working directory is process-wide.
        let _g = crate::cwd_lock();
        let dir =
            std::env::temp_dir().join(format!("mafft-oxide-aamtx-{}-{}", name, std::process::id()));
        let _ = std::fs::create_dir_all(&dir);
        let mut fh = std::fs::File::create(dir.join("_aamtx")).unwrap();
        fh.write_all(contents.as_bytes()).unwrap();
        drop(fh);
        let old = std::env::current_dir().unwrap();
        std::env::set_current_dir(&dir).unwrap();
        let r = f();
        std::env::set_current_dir(old).unwrap();
        let _ = std::fs::remove_dir_all(&dir);
        r
    }

    const AA: &str = "ARNDCQEGHILKMFPSTWYV";

    /// A file in the canonical order whose (i,j) entry is a known function of
    /// i and j, so a permutation is easy to check for.
    fn triangle(order: &str, extra: &str) -> String {
        let idx = |c: char| AA.find(c).unwrap() as i32;
        let mut s = String::from("# test\n   ");
        s.push_str(
            &order
                .chars()
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join("  "),
        );
        s.push('\n');
        for (i, a) in order.chars().enumerate() {
            s.push(a);
            for b in order.chars().take(i + 1) {
                // Symmetric and distinct for every unordered pair.
                s.push_str(&format!(" {:4}", idx(a) * idx(b) + idx(a) + idx(b)));
            }
            s.push('\n');
        }
        s.push_str(extra);
        s
    }

    fn expect(i: usize, j: usize) -> f64 {
        (i * j + i + j) as f64
    }

    #[test]
    fn a_canonical_file_loads_in_order() {
        let (val, rescale) = in_dir("plain", &triangle(AA, ""), || loadaamtx(b'p'));
        assert!(rescale);
        let mut k = 0;
        for i in 0..20 {
            for j in 0..=i {
                assert_eq!(val[k], expect(i, j), "at {i},{j}");
                k += 1;
            }
        }
        // No `frequency` line: the sentinel stays.
        assert!(val[400..420].iter().all(|&v| v == -1.0));
    }

    /// The file's own residue order is arbitrary; the matrix is permuted into
    /// ours, so the values land in the same places as the canonical file's.
    #[test]
    fn a_permuted_file_gives_the_same_matrix() {
        let permuted = "WYVARNDCQEGHILKMFPST";
        let (val, _) = in_dir("perm", &triangle(permuted, ""), || loadaamtx(b'p'));
        let mut k = 0;
        for i in 0..20 {
            for j in 0..=i {
                assert_eq!(val[k], expect(i, j), "at {i},{j}");
                k += 1;
            }
        }
    }

    #[test]
    fn a_frequency_line_is_read_in_the_files_order() {
        let permuted = "WYVARNDCQEGHILKMFPST";
        let freq: Vec<String> = (0..20).map(|i| format!("0.0{:02}", i + 1)).collect();
        let extra = format!("frequency {}\n", freq.join(" "));
        let (val, rescale) = in_dir("freq", &triangle(permuted, &extra), || loadaamtx(b'p'));
        assert!(rescale);
        // freq[k] belongs to permuted[k]; val[400+i] belongs to AA[i].
        for (i, a) in AA.chars().enumerate() {
            let p = permuted.find(a).unwrap();
            let want: f64 = freq[p].parse().unwrap();
            assert!(
                (val[400 + i] - want).abs() < 1e-12,
                "{a}: {} vs {want}",
                val[400 + i]
            );
        }
    }

    #[test]
    fn norescale_is_picked_up() {
        let (_, rescale) = in_dir("nores", &triangle(AA, "norescale\n"), || loadaamtx(b'p'));
        assert!(!rescale);
    }

    #[test]
    fn atof_at_stops_where_strtod_would() {
        assert_eq!(atof_at(b"  -12.5 3", 2), -12.5);
        assert_eq!(atof_at(b"7", 0), 7.0);
        assert_eq!(atof_at(b".5x", 0), 0.5);
        assert_eq!(atof_at(b"1e2 ", 0), 100.0);
        // A stray 'e' with no exponent digits is not part of the number.
        assert_eq!(atof_at(b"3e", 0), 3.0);
        assert_eq!(atof_at(b"abc", 3), 0.0);
    }
}

/// `checkchar()`: which byte values `_aamtx` may not name in `--text` mode.
fn checkchar(i: u32) -> bool {
    i > 0xff
        || i == 0
        || i == 0x3E
        || i == 0x3D
        || i == 0x3C
        || i == 0x2D
        || i == 0x20
        || i == 0x0d
        || i == 0x0a
}

/// `overridematrix()`: `--textmatrix`'s `_aamtx`, one `<hex> <hex> <score>`
/// per line. A missing file is a warning, not an error.
///
/// Returns whether any cell was actually overridden, which is what gates the
/// "alphabet used but unscored" check.
fn overridematrix(matrix: &mut [Vec<f64>]) -> bool {
    let text = match std::fs::read("_aamtx") {
        Ok(t) => t,
        Err(_) => {
            eprintln!("warning: cannot open scorematrix. Use the default one.");
            return false;
        }
    };
    let mut userdefined = false;
    for raw in text.split_inclusive(|&c| c == b'\n') {
        if raw.last() != Some(&b'\n') {
            // `fgets` then `feof`: a last line with no newline is not read.
            break;
        }
        // `#` truncates the line wherever it appears.
        let line: &[u8] = match raw.iter().position(|&c| c == b'#') {
            Some(p) => &raw[..p],
            None => raw,
        };
        let mut it = line
            .split(|c: &u8| c.is_ascii_whitespace())
            .filter(|t| !t.is_empty());
        let (t1, t2, t3) = (it.next(), it.next(), it.next());
        // `sscanf` returning EOF (nothing at all on the line) is skipped;
        // anything else short of three fields is fatal.
        if t1.is_none() {
            continue;
        }
        let parse_hex = |t: Option<&[u8]>| -> Option<u32> {
            let t = t?;
            let s = std::str::from_utf8(t).ok()?;
            let s = s
                .strip_prefix("0x")
                .or_else(|| s.strip_prefix("0X"))
                .unwrap_or(s);
            u32::from_str_radix(s, 16).ok()
        };
        let i1 = parse_hex(t1);
        let i2 = parse_hex(t2);
        let v = t3
            .and_then(|t| std::str::from_utf8(t).ok())
            .and_then(|s| s.parse::<f64>().ok());
        let (i1, i2, v) = match (i1, i2, v) {
            (Some(a), Some(b), Some(c)) => (a, b, c),
            _ => {
                eprintln!("Format error in this line?");
                eprintln!("{}", String::from_utf8_lossy(raw).trim_end());
                std::process::exit(1);
            }
        };
        if checkchar(i1) {
            eprintln!("{}=0x{:x} cannot be used (1)", i1 as u8 as char, i1);
            std::process::exit(1);
        }
        if checkchar(i2) {
            eprintln!("{}=0x{:x} cannot be used (2)", i2 as u8 as char, i2);
            std::process::exit(1);
        }
        matrix[i1 as usize][i2 as usize] = v;
        matrix[i2 as usize][i1 as usize] = v;
        userdefined = true;
    }
    userdefined
}

/// `extendedmtx()`: the `--text` alphabet is the byte values themselves, with
/// a flat 1.0/-1.0 matrix that `_aamtx` then overrides.
pub fn extendedmtx(
    nalphabets: usize,
    matrix: &mut [Vec<f64>],
    freq: &mut [f64],
    amino: &mut [u8],
    amino_grp: &mut [i32; 0x100],
) -> bool {
    for i in 0..nalphabets {
        amino[i] = i as u8;
    }
    for i in 0..nalphabets {
        amino_grp[amino[i] as usize] = (i % 6) as i32;
    }
    for f in freq.iter_mut().take(nalphabets) {
        *f = 1.0 / nalphabets as f64;
    }
    for i in 0..nalphabets {
        for j in 0..=i {
            matrix[i][j] = -1.0;
            matrix[j][i] = -1.0;
        }
    }
    for i in 0..nalphabets {
        matrix[i][i] = 1.0;
    }
    overridematrix(matrix)
}

#[cfg(test)]
mod text_tests {
    use super::*;
    use std::io::Write as _;

    fn in_dir<T>(tag: &str, aamtx: Option<&str>, f: impl FnOnce() -> T) -> T {
        let _g = crate::cwd_lock();
        let d = std::env::temp_dir().join(format!("mafft-oxide-text-{}-{}", tag, std::process::id()));
        let _ = std::fs::remove_dir_all(&d);
        let _ = std::fs::create_dir_all(&d);
        if let Some(c) = aamtx {
            let mut fh = std::fs::File::create(d.join("_aamtx")).unwrap();
            fh.write_all(c.as_bytes()).unwrap();
        }
        let old = std::env::current_dir().unwrap();
        std::env::set_current_dir(&d).unwrap();
        let r = f();
        std::env::set_current_dir(old).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        r
    }

    fn blank() -> (Vec<Vec<f64>>, Vec<f64>, Vec<u8>, [i32; 0x100]) {
        (
            vec![vec![0.0; 0x100]; 0x100],
            vec![0.0; 0x100],
            vec![0u8; 0x100],
            [0i32; 0x100],
        )
    }

    /// Without `_aamtx` the alphabet is the byte values, the matrix is 1.0 on
    /// the diagonal and -1.0 elsewhere, and nothing is "user defined".
    #[test]
    fn the_default_text_matrix_is_identity_ish() {
        let (mut m, mut f, mut a, mut g) = blank();
        let user = in_dir("default", None, || {
            extendedmtx(0x100, &mut m, &mut f, &mut a, &mut g)
        });
        assert!(!user);
        for i in 0..0x100 {
            assert_eq!(a[i], i as u8);
            assert_eq!(g[i], (i % 6) as i32);
            assert!((f[i] - 1.0 / 256.0).abs() < 1e-15);
            assert_eq!(m[i][i], 1.0);
        }
        assert_eq!(m[b'a' as usize][b'b' as usize], -1.0);
    }

    /// `_aamtx` overrides cells symmetrically, and `#` ends a line.
    #[test]
    fn the_matrix_file_overrides_cells_both_ways() {
        let (mut m, mut f, mut a, mut g) = blank();
        let file = "# a comment\n0x41 0x42 7.5\n0x41 0x41 -3 # trailing comment\n";
        let user = in_dir("over", Some(file), || {
            extendedmtx(0x100, &mut m, &mut f, &mut a, &mut g)
        });
        assert!(user);
        assert_eq!(m[0x41][0x42], 7.5);
        assert_eq!(m[0x42][0x41], 7.5);
        assert_eq!(m[0x41][0x41], -3.0);
        // Untouched cells keep the defaults.
        assert_eq!(m[0x43][0x43], 1.0);
        assert_eq!(m[0x43][0x44], -1.0);
    }

    /// A last line with no newline is not read, matching `fgets` + `feof`.
    #[test]
    fn a_final_line_without_a_newline_is_dropped() {
        let (mut m, mut f, mut a, mut g) = blank();
        in_dir("nonl", Some("0x41 0x42 7\n0x43 0x44 9"), || {
            extendedmtx(0x100, &mut m, &mut f, &mut a, &mut g)
        });
        assert_eq!(m[0x41][0x42], 7.0);
        assert_eq!(m[0x43][0x44], -1.0, "the unterminated line is ignored");
    }

    #[test]
    fn checkchar_rejects_the_reserved_bytes() {
        // '=' '<' '>' '-' ' ' CR LF and 0, plus anything past a byte.
        for c in [0x3d, 0x3c, 0x3e, 0x2d, 0x20, 0x0d, 0x0a, 0x00, 0x100] {
            assert!(checkchar(c), "0x{c:x} should be rejected");
        }
        for c in [0x41u32, 0x7e, 0xff, 0x01] {
            assert!(!checkchar(c), "0x{c:x} should be allowed");
        }
    }
}
