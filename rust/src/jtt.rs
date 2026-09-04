//! JTT and transmembrane substitution counts, transliterated from
//! `core/JTT.c`.
//!
//! One table carries both models: the lower triangle is JTT's accepted point
//! mutations, the upper triangle the transmembrane ones. `JTTmtx()` normalises
//! whichever half it wants by the corresponding background frequencies and
//! mirrors it into a symmetric rate matrix, so the two halves never mix.

/// `r` in `JTTmtx()`: lower triangle JTT, upper triangle transmembrane.
pub const R_RAW: [[f64; 20]; 20] = [
    [
        0.0, 21.0, 2.0, 7.0, 13.0, 4.0, 6.0, 160.0, 6.0, 44.0, 43.0, 5.0, 10.0, 21.0, 34.0, 198.0,
        202.0, 0.0, 1.0, 292.0,
    ],
    [
        247.0, 0.0, 0.0, 1.0, 2.0, 21.0, 3.0, 22.0, 21.0, 4.0, 8.0, 53.0, 19.0, 0.0, 1.0, 5.0, 5.0,
        28.0, 0.0, 0.0,
    ],
    [
        216.0, 116.0, 0.0, 14.0, 1.0, 7.0, 0.0, 0.0, 8.0, 4.0, 5.0, 11.0, 3.0, 1.0, 2.0, 32.0,
        19.0, 1.0, 1.0, 2.0,
    ],
    [
        386.0, 48.0, 1433.0, 0.0, 0.0, 0.0, 12.0, 15.0, 4.0, 1.0, 0.0, 2.0, 1.0, 0.0, 1.0, 0.0,
        6.0, 0.0, 1.0, 4.0,
    ],
    [
        106.0, 125.0, 32.0, 13.0, 0.0, 0.0, 0.0, 13.0, 2.0, 4.0, 11.0, 0.0, 1.0, 34.0, 0.0, 48.0,
        13.0, 8.0, 23.0, 47.0,
    ],
    [
        208.0, 750.0, 159.0, 130.0, 9.0, 0.0, 16.0, 1.0, 26.0, 1.0, 16.0, 6.0, 3.0, 0.0, 5.0, 7.0,
        2.0, 0.0, 0.0, 0.0,
    ],
    [
        600.0, 119.0, 180.0, 2914.0, 8.0, 1027.0, 0.0, 21.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        4.0, 2.0, 0.0, 0.0, 7.0,
    ],
    [
        1183.0, 614.0, 291.0, 577.0, 98.0, 84.0, 610.0, 0.0, 1.0, 10.0, 0.0, 0.0, 3.0, 4.0, 7.0,
        64.0, 12.0, 5.0, 0.0, 53.0,
    ],
    [
        46.0, 446.0, 466.0, 144.0, 40.0, 635.0, 41.0, 41.0, 0.0, 3.0, 2.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        4.0, 0.0, 29.0, 2.0,
    ],
    [
        173.0, 76.0, 130.0, 37.0, 19.0, 20.0, 43.0, 25.0, 26.0, 0.0, 273.0, 0.0, 161.0, 66.0, 4.0,
        22.0, 150.0, 1.0, 4.0, 883.0,
    ],
    [
        257.0, 205.0, 63.0, 34.0, 36.0, 314.0, 65.0, 56.0, 134.0, 1324.0, 0.0, 1.0, 153.0, 251.0,
        37.0, 43.0, 26.0, 20.0, 6.0, 255.0,
    ],
    [
        200.0, 2348.0, 758.0, 102.0, 7.0, 858.0, 754.0, 142.0, 85.0, 75.0, 94.0, 0.0, 4.0, 0.0,
        0.0, 1.0, 2.0, 0.0, 5.0, 1.0,
    ],
    [
        100.0, 61.0, 39.0, 27.0, 23.0, 52.0, 30.0, 27.0, 21.0, 704.0, 974.0, 103.0, 0.0, 8.0, 0.0,
        1.0, 32.0, 1.0, 5.0, 89.0,
    ],
    [
        51.0, 16.0, 15.0, 8.0, 66.0, 9.0, 13.0, 18.0, 50.0, 196.0, 1093.0, 7.0, 49.0, 0.0, 0.0,
        32.0, 9.0, 2.0, 54.0, 37.0,
    ],
    [
        901.0, 217.0, 31.0, 39.0, 15.0, 395.0, 71.0, 93.0, 157.0, 31.0, 578.0, 77.0, 23.0, 36.0,
        0.0, 9.0, 10.0, 0.0, 1.0, 1.0,
    ],
    [
        2413.0, 413.0, 1738.0, 244.0, 353.0, 182.0, 156.0, 1131.0, 138.0, 172.0, 436.0, 228.0,
        54.0, 309.0, 1138.0, 0.0, 134.0, 1.0, 22.0, 13.0,
    ],
    [
        2440.0, 230.0, 693.0, 151.0, 66.0, 149.0, 142.0, 164.0, 76.0, 930.0, 172.0, 398.0, 343.0,
        39.0, 412.0, 2258.0, 0.0, 1.0, 3.0, 48.0,
    ],
    [
        11.0, 109.0, 2.0, 5.0, 38.0, 12.0, 12.0, 69.0, 5.0, 12.0, 82.0, 9.0, 8.0, 37.0, 6.0, 36.0,
        8.0, 0.0, 2.0, 18.0,
    ],
    [
        41.0, 46.0, 114.0, 89.0, 164.0, 40.0, 15.0, 15.0, 514.0, 61.0, 84.0, 20.0, 17.0, 850.0,
        22.0, 164.0, 45.0, 41.0, 0.0, 2.0,
    ],
    [
        1766.0, 69.0, 55.0, 127.0, 99.0, 58.0, 226.0, 276.0, 22.0, 3938.0, 1261.0, 58.0, 559.0,
        189.0, 84.0, 219.0, 526.0, 27.0, 42.0, 0.0,
    ],
];

/// `freq0`: JTT background frequencies.
pub const FREQ0: [f64; 20] = [
    0.077, 0.051, 0.043, 0.052, 0.020, 0.041, 0.062, 0.074, 0.023, 0.052, 0.091, 0.059, 0.024,
    0.040, 0.051, 0.069, 0.059, 0.014, 0.032, 0.066,
];

/// `freq0_TM`: transmembrane background frequencies.
pub const FREQ0_TM: [f64; 20] = [
    0.1051, 0.0157, 0.0185, 0.0089, 0.0219, 0.0141, 0.0097, 0.0758, 0.0168, 0.1188, 0.1635, 0.0112,
    0.0333, 0.0777, 0.0260, 0.0568, 0.0523, 0.0223, 0.0324, 0.1195,
];

/// `JTTmtx()`: build the symmetric rate matrix `rsr` and the frequencies.
///
/// The two branches read opposite triangles of `r`, which is why the raw table
/// above holds both models at once.
pub fn jtt_mtx(is_tm: bool) -> ([[f64; 20]; 20], [f64; 20]) {
    let mut r = R_RAW;
    for (i, row) in r.iter_mut().enumerate() {
        row[i] = 0.0;
    }
    let freq = if is_tm { FREQ0_TM } else { FREQ0 };
    if is_tm {
        for i in 1..20 {
            for j in 0..i {
                r[j][i] /= 400.0 * freq[i] * freq[j];
                r[i][j] = r[j][i];
            }
        }
    } else {
        for i in 1..20 {
            for j in 0..i {
                r[i][j] /= 400.0 * freq[i] * freq[j];
                r[j][i] = r[i][j];
            }
        }
    }
    (r, freq)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// XOR of the upper triangle's f64 bit patterns. The matrix is symmetric,
    /// so XOR-ing the whole thing cancels to zero and pins nothing; the upper
    /// triangle alone catches any cell that moves.
    fn upper_xor(r: &[[f64; 20]; 20]) -> u64 {
        let mut x = 0u64;
        for i in 0..20 {
            for j in i + 1..20 {
                x ^= r[i][j].to_bits();
            }
        }
        x
    }

    /// Every expected value came from the C, by linking `core/JTT.c` against a
    /// harness that calls `JTTmtx()` and dumps bit patterns.
    #[test]
    fn jtt_matches_the_c() {
        let (r, freq) = jtt_mtx(false);
        assert_eq!(upper_xor(&r), 0x0019_7346_79cc_c860);
        assert_eq!(r[0][1].to_bits(), 0x4063_a7d4_b6ca_bd71);
        assert_eq!(r[3][7].to_bits(), 0x4077_6deb_c685_0e5f);
        assert_eq!(r[18][19].to_bits(), 0x4048_dba2_e8ba_2e8b);
        assert_eq!(r[0][19].to_bits(), 0x408b_2605_0992_2f5a);
        assert_eq!(freq[0].to_bits(), 0x3fb3_b645_a1ca_c083);
        assert_eq!(freq[19].to_bits(), 0x3fb0_e560_4189_374c);
    }

    #[test]
    fn transmembrane_matches_the_c() {
        let (r, freq) = jtt_mtx(true);
        assert_eq!(upper_xor(&r), 0x3f16_c579_8535_2c94);
        assert_eq!(r[0][1].to_bits(), 0x403f_d11b_e71e_f724);
        assert_eq!(r[3][7].to_bits(), 0x404b_cb1d_dd50_f10a);
        assert_eq!(r[18][19].to_bits(), 0x3ff4_a987_87ed_7c95);
        assert_eq!(r[0][19].to_bits(), 0x404d_0fd0_f4b5_f9d3);
        assert_eq!(freq[0].to_bits(), 0x3fba_e7d5_66cf_41f2);
        assert_eq!(freq[19].to_bits(), 0x3fbe_978d_4fdf_3b64);
    }

    #[test]
    fn the_two_models_read_opposite_triangles_of_one_table() {
        // R_RAW carries JTT below the diagonal and transmembrane above it, so
        // the two calls must not agree anywhere off-diagonal -- if they did,
        // one branch would be reading the wrong half.
        let (jtt, _) = jtt_mtx(false);
        let (tm, _) = jtt_mtx(true);
        let agree = (0..20)
            .flat_map(|i| (i + 1..20).map(move |j| (i, j)))
            .filter(|&(i, j)| jtt[i][j] == tm[i][j])
            .count();
        assert_eq!(agree, 0, "no off-diagonal cell should match between models");
    }

    #[test]
    fn both_models_are_symmetric_with_a_zero_diagonal() {
        for tm in [false, true] {
            let (r, _) = jtt_mtx(tm);
            for i in 0..20 {
                assert_eq!(r[i][i], 0.0, "diagonal {} (tm={})", i, tm);
                for j in 0..20 {
                    assert_eq!(
                        r[i][j].to_bits(),
                        r[j][i].to_bits(),
                        "({},{}) tm={}",
                        i,
                        j,
                        tm
                    );
                }
            }
        }
    }
}
