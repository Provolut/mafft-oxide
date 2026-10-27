//! The generalised affine-gap local pairwise DP, transliterated from
//! `genL__align11` in `core/genalign11.c`.
//!
//! This is what E-INS-i (`--genafpair`) runs on every pair of input sequences.
//! It differs from `L__align11` by a third recurrence: besides the usual
//! "gap in seq1" and "gap in seq2" states it keeps a *pair* of unaligned runs
//! (`tbk`, tracked by `largeM` / `Mi`), which can be entered from either side
//! for a single `penalty_OP` however long the two runs are. That is what lets
//! it leave large, mutually unalignable insertions out of the alignment
//! instead of paying an extension penalty across them.
//!
//! Three details decide the result and are reproduced exactly:
//!
//! * the traceback needs two matrices, not one: a cell records both the row
//!   and the column it came from, because a `tbk` transition moves in both
//!   axes at once;
//! * `tbk` is compared against `wm` *before* it is updated from this cell's
//!   `Mi` / `largeM[j]`, so a run can never be entered and left in one step;
//! * `Mi` and `largeM[j]` are never charged an extension -- the
//!   `+= fpenalty_EX` lines are commented out in the C, which is what makes
//!   the unaligned run free after its opening cost.

use crate::constants::Constants;

fn max_usize(a: usize, b: usize) -> usize {
    if a > b {
        a
    } else {
        b
    }
}

fn fill_amino_dynamicmtx(cn: &Constants, dst: &mut [Vec<f64>], src: &[Vec<f64>]) {
    for i in 0..cn.nalphabets {
        for j in 0..cn.nalphabets {
            dst[cn.amino[i] as usize][cn.amino[j] as usize] = src[i][j];
        }
    }
}

/// One generalised local alignment.
pub struct GenAln {
    pub score: f64,
    pub a1: Vec<u8>,
    pub a2: Vec<u8>,
    pub off1: usize,
    pub off2: usize,
}

/// The scratch buffers `genL__align11` keeps in `static`s.
pub struct GenalignState {
    orlgth1: usize,
    orlgth2: usize,
    w1: Vec<f64>,
    w2: Vec<f64>,
    initverticalw: Vec<f64>,
    lastverticalw: Vec<f64>,
    m: Vec<f64>,
    mp: Vec<i32>,
    large_m: Vec<f64>,
    mp_big: Vec<i32>,
    amino_dynamicmtx: Vec<Vec<f64>>,

    /// `commonIP` and `commonJP`. Every cell the traceback reads is written
    /// either by the DP or by the border pass at the top of `gentracking`, so
    /// unlike the other aligners nothing stale is ever read back -- but the
    /// growth rule is kept the same anyway.
    ijpi: Vec<Vec<i32>>,
    ijpj: Vec<Vec<i32>>,
    common_alloc1: usize,
    common_alloc2: usize,
}

impl Default for GenalignState {
    fn default() -> Self {
        Self::new()
    }
}

impl GenalignState {
    pub fn new() -> GenalignState {
        GenalignState {
            orlgth1: 0,
            orlgth2: 0,
            w1: Vec::new(),
            w2: Vec::new(),
            initverticalw: Vec::new(),
            lastverticalw: Vec::new(),
            m: Vec::new(),
            mp: Vec::new(),
            large_m: Vec::new(),
            mp_big: Vec::new(),
            amino_dynamicmtx: Vec::new(),
            ijpi: Vec::new(),
            ijpj: Vec::new(),
            common_alloc1: 0,
            common_alloc2: 0,
        }
    }

    /// `genL__align11()`.
    pub fn gen_l_align11(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        seq1: &[u8],
        seq2: &[u8],
    ) -> GenAln {
        let lgth1 = seq1.len();
        let lgth2 = seq2.len();

        let fpenalty = cn.penalty as f64;
        let fpenalty_op = cn.penalty_op as f64;
        let fpenalty_ex = cn.penalty_ex as f64;
        let foffset = cn.offset as f64;
        let localthr = -foffset;
        let localthr2 = -foffset;

        if lgth1 > self.orlgth1 || lgth2 > self.orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.orlgth2) + 100;
            self.w1 = vec![0.0; ll2 + 2];
            self.w2 = vec![0.0; ll2 + 2];
            self.initverticalw = vec![0.0; ll1 + 2];
            self.lastverticalw = vec![0.0; ll1 + 2];
            self.m = vec![0.0; ll2 + 2];
            self.mp = vec![0i32; ll2 + 2];
            self.large_m = vec![0.0; ll2 + 2];
            self.mp_big = vec![0i32; ll2 + 2];
            // Freed before the replacement is built: an assignment evaluates its
            // right-hand side first, so the old buffer and the new one are both
            // live for the duration -- and a growth happens at the moment the
            // working set is largest, so the overlap lands on the peak. Worth
            // 2.8 MB where `commonip.rs` had the same shape.
            self.amino_dynamicmtx = Vec::new();
            self.amino_dynamicmtx = vec![vec![0.0; 0x100]; 0x100];
            self.orlgth1 = ll1 - 100;
            self.orlgth2 = ll2 - 100;
        }

        let mut mtx = std::mem::take(&mut self.amino_dynamicmtx);
        fill_amino_dynamicmtx(cn, &mut mtx, n_dynamicmtx);

        if self.orlgth1 > self.common_alloc1 || self.orlgth2 > self.common_alloc2 {
            let ll1 = max_usize(self.orlgth1, self.common_alloc1);
            let ll2 = max_usize(self.orlgth2, self.common_alloc2);
            self.ijpi = Vec::new();
            self.ijpi = vec![vec![0i32; ll2 + 10]; ll1 + 10];
            self.ijpj = Vec::new();
            self.ijpj = vec![vec![0i32; ll2 + 10]; ll1 + 10];
            self.common_alloc1 = ll1;
            self.common_alloc2 = ll2;
        }

        let mut w1 = std::mem::take(&mut self.w1);
        let mut w2 = std::mem::take(&mut self.w2);
        let mut initverticalw = std::mem::take(&mut self.initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.lastverticalw);
        let mut m = std::mem::take(&mut self.m);
        let mut mp = std::mem::take(&mut self.mp);
        let mut large_m = std::mem::take(&mut self.large_m);
        let mut mp_big = std::mem::take(&mut self.mp_big);
        let mut ijpi = std::mem::take(&mut self.ijpi);
        let mut ijpj = std::mem::take(&mut self.ijpj);

        // match_calc_mtx( mtx, initverticalw, seq2, seq1, 0, lgth1 )
        {
            let row = &mtx[seq2[0] as usize];
            for i in 0..lgth1 {
                initverticalw[i] = row[seq1[i] as usize];
            }
        }
        // match_calc_mtx( mtx, currentw, seq1, seq2, 0, lgth2 )
        {
            let row = &mtx[seq1[0] as usize];
            for j in 0..lgth2 {
                w1[j] = row[seq2[j] as usize];
            }
        }

        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1];
            mp[j] = 0;
            large_m[j] = w1[j - 1];
            mp_big[j] = 0;
        }
        lastverticalw[0] = w1[lgth2 - 1];

        let localstop = (lgth1 + lgth2 + 1) as i32;
        let mut maxwm = -999999999.9f64;
        let mut endali = 0usize;
        let mut endalj = 0usize;

        let lastj = lgth2 + 1;
        let mut cur_is_w1 = true;
        for i in 1..lgth1 + 1 {
            cur_is_w1 = !cur_is_w1;
            let (currentw, previousw) = if cur_is_w1 {
                (&mut w1, &mut w2)
            } else {
                (&mut w2, &mut w1)
            };
            previousw[0] = initverticalw[i - 1];
            {
                // `lasti` is lgth1+1, so the last row reads seq1[lgth1] -- the
                // NUL terminator, whose row of the matrix is all zeros.
                let c1 = if i < lgth1 { seq1[i] } else { 0u8 };
                let row = &mtx[c1 as usize];
                for j in 0..lgth2 {
                    currentw[j] = row[seq2[j] as usize];
                }
            }
            currentw[0] = initverticalw[i];

            let mut mi = previousw[0];
            let mut mpi = 0i32;
            let mut big_mi = previousw[0];
            let mut mpi_big = 0i32;
            let mut tbk = -999999.9f64;
            let mut tbki = 0i32;
            let mut tbkj = 0i32;

            let ipi = &mut ijpi[i];
            let ipj = &mut ijpj[i];

            for j in 1..lastj {
                let mut wm = previousw[j - 1];
                ipi[j] = (i - 1) as i32;
                ipj[j] = (j - 1) as i32;

                let g = mi + fpenalty;
                if g > wm {
                    wm = g;
                    ipj[j] = mpi;
                }
                let g = previousw[j - 1];
                if g > mi {
                    mi = g;
                    mpi = (j - 1) as i32;
                }
                mi += fpenalty_ex;

                let g = m[j] + fpenalty;
                if g > wm {
                    wm = g;
                    ipi[j] = mp[j];
                    ipj[j] = (j - 1) as i32;
                }
                let g = previousw[j - 1];
                if g > m[j] {
                    m[j] = g;
                    mp[j] = (i - 1) as i32;
                }
                m[j] += fpenalty_ex;

                // The generalised part: enter or extend a pair of unaligned
                // runs. `tbk` is read before it is updated from this cell.
                let g = tbk + fpenalty_op;
                if g > wm {
                    wm = g;
                    ipi[j] = tbki;
                    ipj[j] = tbkj;
                }
                if big_mi > tbk {
                    tbk = big_mi;
                    tbki = (i - 1) as i32;
                    tbkj = mpi_big;
                }
                if large_m[j] > tbk {
                    tbk = large_m[j];
                    tbki = mp_big[j];
                    tbkj = (j - 1) as i32;
                }
                if previousw[j - 1] > large_m[j] {
                    large_m[j] = previousw[j - 1];
                    mp_big[j] = (i - 1) as i32;
                }
                if previousw[j - 1] > big_mi {
                    big_mi = previousw[j - 1];
                    mpi_big = (j - 1) as i32;
                }

                if maxwm < wm {
                    maxwm = wm;
                    endali = i;
                    endalj = j;
                }
                if wm < localthr {
                    ipi[j] = localstop;
                    wm = localthr2;
                }

                currentw[j] += wm;
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        let stopped = ijpi[endali][endalj] == localstop;
        let (a1, a2, off1, off2) = if stopped {
            (Vec::new(), Vec::new(), 0, 0)
        } else {
            gentracking(
                seq1, seq2, &mut ijpi, &mut ijpj, endali, endalj, localstop, b'-',
            )
        };

        self.amino_dynamicmtx = mtx;
        self.w1 = w1;
        self.w2 = w2;
        self.initverticalw = initverticalw;
        self.lastverticalw = lastverticalw;
        self.m = m;
        self.mp = mp;
        self.large_m = large_m;
        self.mp_big = mp_big;
        self.ijpi = ijpi;
        self.ijpj = ijpj;

        GenAln {
            score: if stopped { 0.0 } else { maxwm },
            a1,
            a2,
            off1,
            off2,
        }
    }
}

/// `gentracking()`: walk back through the two index matrices.
#[allow(clippy::too_many_arguments)]
fn gentracking(
    seq1: &[u8],
    seq2: &[u8],
    ijpi: &mut [Vec<i32>],
    ijpj: &mut [Vec<i32>],
    endi: usize,
    endj: usize,
    localstop: i32,
    gap: u8,
) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let lgth1 = seq1.len();
    let lgth2 = seq2.len();

    for row in ijpi.iter_mut().take(lgth1 + 1) {
        row[0] = localstop;
    }
    for row in ijpj.iter_mut().take(lgth1 + 1) {
        row[0] = localstop;
    }
    for j in 0..lgth2 + 1 {
        ijpi[0][j] = localstop;
        ijpj[0][j] = localstop;
    }

    // The C writes backwards from the end of a scratch buffer; collect in
    // reverse and flip at the end.
    let mut r1: Vec<u8> = Vec::with_capacity(lgth1 + lgth2);
    let mut r2: Vec<u8> = Vec::with_capacity(lgth1 + lgth2);

    let mut iin = endi as i32;
    let mut jin = endj as i32;
    let mut ifi = 0i32;
    let mut jfi = 0i32;
    let limk = (lgth1 + lgth2) as i32;
    let mut k = 0i32;
    while k <= limk {
        ifi = ijpi[iin as usize][jin as usize];
        jfi = ijpj[iin as usize][jin as usize];

        let mut l = iin - ifi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            r1.push(seq1[(ifi + l) as usize]);
            r2.push(gap);
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            r1.push(gap);
            r2.push(seq2[(jfi + l) as usize]);
            k += 1;
        }

        if iin <= 0 || jin <= 0 {
            break;
        }
        r1.push(seq1[ifi as usize]);
        r2.push(seq2[jfi as usize]);

        if ijpi[ifi as usize][jfi as usize] == localstop {
            break;
        }
        if ijpj[ifi as usize][jfi as usize] == localstop {
            break;
        }
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1; // the `for` header's own increment
    }

    r1.reverse();
    r2.reverse();
    let off1 = if ifi == -1 { 0 } else { ifi as usize };
    let off2 = if jfi == -1 { 0 } else { jfi as usize };
    (r1, r2, off1, off2)
}
