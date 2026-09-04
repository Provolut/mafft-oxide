//! glibc's `rand()`, because `--treein shuffle` seeds it and shuffles the
//! input order with it. The sequence is part of the output, so it has to be
//! the same sequence, not merely a random one.
//!
//! This is glibc's default TYPE_3 additive-feedback generator: a 34-word
//! state seeded by a Lehmer generator, then `r[i] = r[i-3] + r[i-31]` with
//! the top bit dropped on output. Reproduced exactly, including the 310
//! discarded warm-up draws.

pub struct CRand {
    r: Vec<i32>,
    f: usize,
    rr: usize,
}

impl CRand {
    /// `srand( seed )`.
    pub fn new(seed: u32) -> CRand {
        // glibc: a seed of 0 is treated as 1.
        let seed = if seed == 0 { 1 } else { seed };
        let mut r = vec![0i32; 344];
        r[0] = seed as i32;
        for i in 1..31 {
            // r[i] = (16807 * r[i-1]) % 2147483647, by Schrage's method so it
            // fits in 32 bits -- and glibc really does compute it this way.
            let prev = r[i - 1] as i64;
            let hi = prev / 127_773;
            let lo = prev % 127_773;
            let mut word = 16_807 * lo - 2_836 * hi;
            if word < 0 {
                word += 2_147_483_647;
            }
            r[i] = word as i32;
        }
        for i in 31..34 {
            r[i] = r[i - 31];
        }
        for i in 34..344 {
            r[i] = r[i - 31].wrapping_add(r[i - 3]);
        }
        CRand {
            r,
            f: 344 - 31,
            rr: 344 - 3,
        }
    }

    /// `rand()`.
    pub fn next(&mut self) -> i32 {
        let v = self.r[self.f].wrapping_add(self.r[self.rr]);
        self.r.push(v);
        self.f += 1;
        self.rr += 1;
        ((v as u32) >> 1) as i32
    }
}

/// `stringshuffle()`: Fisher-Yates with a biased `rand() % size`, run over
/// every position rather than downwards, exactly as the C has it.
pub fn stringshuffle(ary: &mut [i32], rng: &mut CRand) {
    let size = ary.len();
    for i in 0..size {
        let j = (rng.next() as usize) % size;
        ary.swap(i, j);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The first twelve draws for a handful of seeds, taken from glibc's own
    /// `rand()` on this machine. If these ever stop matching, `--treein
    /// shuffle` is producing a different alignment from the C's.
    #[test]
    fn it_matches_glibc() {
        let cases: [(u32, [i32; 12]); 4] = [
            (
                // glibc treats seed 0 as seed 1, so these two agree.
                0,
                [
                    1804289383, 846930886, 1681692777, 1714636915, 1957747793, 424238335,
                    719885386, 1649760492, 596516649, 1189641421, 1025202362, 1350490027,
                ],
            ),
            (
                1,
                [
                    1804289383, 846930886, 1681692777, 1714636915, 1957747793, 424238335,
                    719885386, 1649760492, 596516649, 1189641421, 1025202362, 1350490027,
                ],
            ),
            (
                42,
                [
                    71876166, 708592740, 1483128881, 907283241, 442951012, 537146758, 1366999021,
                    1854614940, 647800535, 53523743, 783815874, 1643643143,
                ],
            ),
            (
                12345,
                [
                    383100999, 858300821, 357768173, 455528251, 133005921, 116285904, 591987137,
                    102557902, 689413528, 585691128, 789708827, 477528897,
                ],
            ),
        ];
        for (seed, want) in cases.iter() {
            let mut r = CRand::new(*seed);
            for (k, w) in want.iter().enumerate() {
                assert_eq!(r.next(), *w, "seed {seed}, draw {k}");
            }
        }
    }

    /// The shuffle is a permutation -- it never loses or duplicates an entry.
    #[test]
    fn the_shuffle_is_a_permutation() {
        for seed in [0u32, 1, 7, 99] {
            let mut a: Vec<i32> = (0..20).collect();
            let mut r = CRand::new(seed);
            stringshuffle(&mut a, &mut r);
            let mut sorted = a.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, (0..20).collect::<Vec<i32>>());
        }
    }
}
