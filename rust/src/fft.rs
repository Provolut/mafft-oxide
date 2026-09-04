//! The FFT itself, transliterated from `core/fft.c`.
//!
//! MAFFT locates homologous segments by correlating two sequences that have
//! been converted into vectors of normalised volume and polarity (Katoh,
//! Misawa, Kuma & Miyata, *NAR* 30:3059, 2002) — the correlation is computed
//! through this radix-2 FFT rather than by direct summation.
//!
//! The sine table is built by the recurrence in the original (Okumura's
//! algorithm book), not by calling `sin` per entry: only `sin(PI/n)` and
//! `sqrt` are evaluated. `sqrt` is correctly rounded by IEEE-754, and the
//! single `sin` argument is `PI/n` with `n` a power of two, so it is an exact
//! double. `tests/golden/sintbl.txt` pins the resulting table against the C.

/// `Fukusosuu` — a complex number, laid out as in `mltaln.h`.
#[derive(Clone, Copy, Debug, Default)]
pub struct Cpx {
    pub r: f64,
    pub i: f64,
}

/// `PI` from `fft.h`. The C's literal `3.14159265358979323846` rounds to the
/// same `f64` as `std::f64::consts::PI` (both `0x400921fb54442d18`), so this
/// is the same number, written the way Rust writes it.
pub const PI: f64 = std::f64::consts::PI;

/// `make_sintbl()`.
fn make_sintbl(n: usize, sintbl: &mut [f64]) {
    let n2 = n / 2;
    let n4 = n / 4;
    let n8 = n / 8;
    let mut t = (PI / n as f64).sin();
    let mut dc = 2.0 * t * t;
    let mut ds = (dc * (2.0 - dc)).sqrt();
    t = 2.0 * dc;
    let mut c = 1.0;
    sintbl[n4] = 1.0;
    let mut s = 0.0;
    sintbl[0] = 0.0;
    for i in 1..n8 {
        c -= dc;
        dc += t * c;
        s += ds;
        ds -= t * s;
        sintbl[i] = s;
        sintbl[n4 - i] = c;
    }
    if n8 != 0 {
        sintbl[n8] = 0.5f64.sqrt();
    }
    for i in 0..n4 {
        sintbl[n2 - i] = sintbl[i];
    }
    for i in 0..n2 + n4 {
        sintbl[i + n2] = -sintbl[i];
    }
}

/// `make_bitrev()`.
fn make_bitrev(n: usize, bitrev: &mut [usize]) {
    let n2 = n / 2;
    let mut i = 0usize;
    let mut j = 0usize;
    loop {
        bitrev[i] = j;
        i += 1;
        if i >= n {
            break;
        }
        let mut k = n2;
        while k <= j {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

/// The `static` tables of `fft()`, kept across calls and rebuilt only when the
/// transform length changes.
#[derive(Default)]
pub struct Fft {
    last_n: usize,
    bitrev: Vec<usize>,
    sintbl: Vec<f64>,
}

impl Fft {
    pub fn new() -> Fft {
        Fft::default()
    }

    /// `fft()`. A negative `n` requests the inverse transform, exactly as in C;
    /// the forward transform divides by `n`.
    pub fn run(&mut self, n_signed: isize, x: &mut [Cpx]) {
        let (n, inverse) = if n_signed < 0 {
            ((-n_signed) as usize, true)
        } else {
            (n_signed as usize, false)
        };
        let n4 = n / 4;

        if n != self.last_n || n == 0 {
            self.last_n = n;
            if n == 0 {
                return;
            }
            // realloc(): the C does not preserve or clear the old contents, and
            // both tables are fully rewritten below.
            self.sintbl = vec![0.0; n + n4];
            self.bitrev = vec![0usize; n];
            make_sintbl(n, &mut self.sintbl);
            make_bitrev(n, &mut self.bitrev);
        }

        for i in 0..n {
            let j = self.bitrev[i];
            if i < j {
                x.swap(i, j);
            }
        }

        let mut k = 1usize;
        while k < n {
            let k2 = k + k;
            let mut h = 0usize;
            let d = n / k2;
            for j in 0..k {
                let c = self.sintbl[h + n4];
                let s = if inverse {
                    -self.sintbl[h]
                } else {
                    self.sintbl[h]
                };
                let mut i = j;
                while i < n {
                    let ik = i + k;
                    let d_r = s * x[ik].i + c * x[ik].r;
                    let d_i = c * x[ik].i - s * x[ik].r;
                    x[ik].r = x[i].r - d_r;
                    x[i].r += d_r;
                    x[ik].i = x[i].i - d_i;
                    x[i].i += d_i;
                    i += k2;
                }
                h += d;
            }
            k = k2;
        }

        if !inverse {
            let nf = n as f64;
            for v in x.iter_mut().take(n) {
                v.r /= nf;
                v.i /= nf;
            }
        }
    }
}

/// `calcNaiseki()`: the conjugate product used to build the cross-correlation.
pub fn calc_naiseki(x: &Cpx, y: &Cpx) -> Cpx {
    Cpx {
        r: x.r * y.r + x.i * y.i,
        i: -x.r * y.i + x.i * y.r,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The sine table is the only place the port depends on a libm
    /// transcendental (`sin(PI/n)`, with `n` a power of two). If a platform's
    /// `sin` differs from the reference by even one ulp, alignments can change --
    /// so pin every bit of the table against a dump from `core/fft.c`.
    #[test]
    fn sintbl_matches_c_bit_for_bit() {
        let golden = include_str!("../tests/golden/sintbl.txt");
        let mut lines = golden.lines().peekable();
        let mut checked = 0usize;
        while let Some(header) = lines.next() {
            let mut it = header.split_whitespace();
            assert_eq!(it.next(), Some("n"));
            let n: usize = it.next().unwrap().parse().unwrap();
            let mut tbl = vec![0.0f64; n + n / 4];
            make_sintbl(n, &mut tbl);
            for (i, v) in tbl.iter().enumerate() {
                let want = u64::from_str_radix(lines.next().unwrap(), 16).unwrap();
                assert_eq!(
                    v.to_bits(),
                    want,
                    "sintbl[{}] for n={}: {:016x} != {:016x}",
                    i,
                    n,
                    v.to_bits(),
                    want
                );
            }
            checked += 1;
        }
        assert!(checked >= 5);
    }

    /// The bit-reversal permutation is integer-only, but it decides the whole
    /// transform, so check it against a direct definition.
    #[test]
    fn bitrev_is_a_bit_reversal() {
        for bits in 1..12u32 {
            let n = 1usize << bits;
            let mut br = vec![0usize; n];
            make_bitrev(n, &mut br);
            for i in 0..n {
                let want = i.reverse_bits() >> (usize::BITS - bits);
                assert_eq!(br[i], want, "n={} i={}", n, i);
            }
        }
    }
}
