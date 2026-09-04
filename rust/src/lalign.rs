//! Local (Smith-Waterman style) pairwise DP, transliterated from
//! `core/Lalign11.c`, plus the score-only global variant from
//! `core/Galign11.c`.
//!
//! `L__align11` is the function L-INS-i spends most of its time in: it is run
//! once for every pair of input sequences to find the best-scoring local
//! alignment, whose ungapped blocks then become the position-specific
//! constraints ("localhom") that steer the progressive stage.
//!
//! Three details decide the result and are reproduced exactly:
//!
//! * The recurrence clamps at `localthr`; a cell below it is marked
//!   `localstop` and its score replaced by `localthr2`, which is what makes
//!   the alignment local rather than global.
//! * `maxwm` is compared with strict `<`, so the *first* cell achieving the
//!   maximum wins, and the traceback starts there.
//! * The `mi` / `mpi` recurrence reads only the previous row, so it is
//!   computed as a prefix scan ahead of the cells that consume it. Same
//!   operations in the same order; every value is an exact integer here.

use crate::constants::Constants;
use crate::warp::Warp;

/// The scratch buffers `L__align11` and friends keep in `static`s.
pub struct LalignState {
    orlgth1: usize,
    orlgth2: usize,
    w1: Vec<f64>,
    w2: Vec<f64>,
    initverticalw: Vec<f64>,
    lastverticalw: Vec<f64>,
    m: Vec<f64>,
    mp: Vec<i32>,
    miv: Vec<f64>,
    mpiv: Vec<i32>,
    amino_dynamicmtx: Vec<Vec<f64>>,

    // L__align11_noalign has its own set.
    n_orlgth1: usize,
    n_orlgth2: usize,
    n_w1: Vec<f64>,
    n_w2: Vec<f64>,
    n_initverticalw: Vec<f64>,
    n_lastverticalw: Vec<f64>,
    n_m: Vec<f64>,
    n_amino_dynamicmtx: Vec<Vec<f64>>,

    // G__align11_noalign likewise.
    g_orlgth1: usize,
    g_orlgth2: usize,
    g_w1: Vec<f64>,
    g_w2: Vec<f64>,
    g_initverticalw: Vec<f64>,
    g_lastverticalw: Vec<f64>,
    g_m: Vec<f64>,
    g_amino_dynamicmtx: Vec<Vec<f64>>,

    /// `commonIP`, shared with the other aligners in the same process.
    ijp: Vec<Vec<i32>>,
    common_alloc1: usize,
    common_alloc2: usize,
}

impl Default for LalignState {
    fn default() -> Self {
        Self::new()
    }
}

fn max_usize(a: usize, b: usize) -> usize {
    if a > b {
        a
    } else {
        b
    }
}

impl LalignState {
    pub fn new() -> LalignState {
        LalignState {
            orlgth1: 0,
            orlgth2: 0,
            w1: Vec::new(),
            w2: Vec::new(),
            initverticalw: Vec::new(),
            lastverticalw: Vec::new(),
            m: Vec::new(),
            mp: Vec::new(),
            miv: Vec::new(),
            mpiv: Vec::new(),
            amino_dynamicmtx: Vec::new(),
            n_orlgth1: 0,
            n_orlgth2: 0,
            n_w1: Vec::new(),
            n_w2: Vec::new(),
            n_initverticalw: Vec::new(),
            n_lastverticalw: Vec::new(),
            n_m: Vec::new(),
            n_amino_dynamicmtx: Vec::new(),
            g_orlgth1: 0,
            g_orlgth2: 0,
            g_w1: Vec::new(),
            g_w2: Vec::new(),
            g_initverticalw: Vec::new(),
            g_lastverticalw: Vec::new(),
            g_m: Vec::new(),
            g_amino_dynamicmtx: Vec::new(),
            ijp: Vec::new(),
            common_alloc1: 0,
            common_alloc2: 0,
        }
    }

    fn ensure_ijp(&mut self, orlgth1: usize, orlgth2: usize) {
        if orlgth1 > self.common_alloc1 || orlgth2 > self.common_alloc2 {
            let ll1 = max_usize(orlgth1, self.common_alloc1);
            let ll2 = max_usize(orlgth2, self.common_alloc2);
            self.ijp = vec![vec![0i32; ll2 + 10]; ll1 + 10];
            self.common_alloc1 = ll1;
            self.common_alloc2 = ll2;
        }
    }
}

/// Fill a 256x256 character-indexed matrix from the alphabet-indexed one.
fn fill_amino_dynamicmtx(cn: &Constants, dst: &mut [Vec<f64>], src: &[Vec<f64>]) {
    for i in 0..cn.nalphabets {
        for j in 0..cn.nalphabets {
            dst[cn.amino[i] as usize][cn.amino[j] as usize] = src[i][j];
        }
    }
}

/// The local alignment of one pair.
pub struct LocalAln {
    pub score: f64,
    /// The two aligned strings (gapped), or empty if nothing scored.
    pub a1: Vec<u8>,
    pub a2: Vec<u8>,
    /// Offsets of the alignment start within the ungapped inputs.
    pub off1: usize,
    pub off2: usize,
}

impl LalignState {
    /// `L__align11()`.
    ///
    /// `scoreoffset` shifts the local threshold; the pair loop passes 0.0
    /// except under `--allowshift`, which this port refuses.
    pub fn l_align11(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        scoreoffset: f64,
        seq1: &[u8],
        seq2: &[u8],
    ) -> LocalAln {
        let lgth1 = seq1.len();
        let lgth2 = seq2.len();

        let localthr = -(cn.offset as f64) + scoreoffset * 600.0;
        let localthr2 = localthr;
        let fpenalty = cn.penalty as f64;
        let fpenalty_ex = cn.penalty_ex as f64;
        let fpenalty_shift = cn.penalty_shift as f64;
        let mut warp = if cn.trywarp {
            Some(Warp::new(lgth1, lgth2))
        } else {
            None
        };

        if lgth1 > self.orlgth1 || lgth2 > self.orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.orlgth2) + 100;
            self.w1 = vec![0.0; ll2 + 2];
            self.w2 = vec![0.0; ll2 + 2];
            self.initverticalw = vec![0.0; ll1 + 2];
            self.lastverticalw = vec![0.0; ll1 + 2];
            self.m = vec![0.0; ll2 + 2];
            self.mp = vec![0; ll2 + 2];
            self.miv = vec![0.0; ll2 + 2];
            self.mpiv = vec![0; ll2 + 2];
            self.amino_dynamicmtx = vec![vec![0.0; 0x100]; 0x100];
            self.orlgth1 = ll1 - 100;
            self.orlgth2 = ll2 - 100;
        }
        {
            let mut mtx = std::mem::take(&mut self.amino_dynamicmtx);
            fill_amino_dynamicmtx(cn, &mut mtx, n_dynamicmtx);
            self.amino_dynamicmtx = mtx;
        }
        self.ensure_ijp(self.orlgth1, self.orlgth2);

        let mtx = std::mem::take(&mut self.amino_dynamicmtx);
        let mut w1 = std::mem::take(&mut self.w1);
        let mut w2 = std::mem::take(&mut self.w2);
        let mut initverticalw = std::mem::take(&mut self.initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.lastverticalw);
        let mut m = std::mem::take(&mut self.m);
        let mut mp = std::mem::take(&mut self.mp);
        let mut miv = std::mem::take(&mut self.miv);
        let mut mpiv = std::mem::take(&mut self.mpiv);
        let mut ijp = std::mem::take(&mut self.ijp);

        {
            let row = &mtx[seq2[0] as usize];
            for i in 0..lgth1 {
                initverticalw[i] = row[seq1[i] as usize];
            }
        }
        {
            let row = &mtx[seq1[0] as usize];
            for j in 0..lgth2 {
                w1[j] = row[seq2[j] as usize];
            }
        }

        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1];
            mp[j] = 0;
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
            // match_calc_mtx( mtx, currentw, seq1, seq2, i, lgth2 ): row i of
            // seq1. i runs up to lgth1, and seq1[lgth1] is the C terminator,
            // whose matrix row is all zeros.
            {
                let c1 = if i < lgth1 { seq1[i] } else { 0u8 };
                let row = &mtx[c1 as usize];
                for j in 0..lgth2 {
                    currentw[j] = row[seq2[j] as usize];
                }
            }
            currentw[0] = initverticalw[i];

            // Prefix scan of the mi/mpi recurrence.
            {
                let mut miscan = previousw[0];
                let mut mpiscan = 0i32;
                for jj in 1..lastj {
                    miv[jj] = miscan;
                    mpiv[jj] = mpiscan;
                    if previousw[jj - 1] > miscan {
                        miscan = previousw[jj - 1];
                        mpiscan = jj as i32 - 1;
                    }
                    miscan += fpenalty_ex;
                }
            }

            let ijprow = &mut ijp[i];
            // The C splits this in two on `trywarp` rather than testing it per
            // cell (Lalign11.c: `if( !trywarp ) { ... } else ...`). Testing it
            // inside the row costs a branch on every one of lgth1*lgth2 cells
            // and stops the body being elementwise in j.
            if warp.is_none() {
                // Straight-line form, matching the C's scalar tail exactly:
                // every conditional is a select and every store is
                // unconditional. Same comparisons on the same operands in the
                // same order, so the values are identical -- and every one of
                // them is an exact integer here (the score matrix and both
                // penalties are ints), so no reassociation question arises.
                //
                // Sliced to `lastj` up front so the bounds checks fold away:
                // the compiler cannot otherwise see that `jj < lastj` is in
                // range for buffers sized from a previous, larger call.
                let pw = &previousw[..lastj];
                let cw = &mut currentw[..lastj];
                let mv = &miv[..lastj];
                let mpv = &mpiv[..lastj];
                let mm = &mut m[..lastj];
                let mpp = &mut mp[..lastj];
                let ijr = &mut ijprow[..lastj];
                for jj in 1..lastj {
                    let diag = pw[jj - 1];
                    let gi = mv[jj] + fpenalty;
                    let gj = mm[jj] + fpenalty;
                    let taki = gi > diag;
                    let w1v = if taki { gi } else { diag };
                    let ip1 = if taki { -(jj as i32 - mpv[jj]) } else { 0 };
                    let takj = gj > w1v;
                    let w2v = if takj { gj } else { w1v };
                    let ip2 = if takj { i as i32 - mpp[jj] } else { ip1 };
                    let takm = diag > mm[jj];

                    mpp[jj] = if takm { i as i32 - 1 } else { mpp[jj] };
                    mm[jj] = if takm { diag } else { mm[jj] } + fpenalty_ex;

                    if maxwm < w2v {
                        maxwm = w2v;
                        endali = i;
                        endalj = jj;
                    }

                    let stop = w2v < localthr;
                    ijr[jj] = if stop { localstop } else { ip2 };
                    cw[jj] += if stop { localthr2 } else { w2v };
                }
            } else {
                for jj in 1..lastj {
                    let diag = previousw[jj - 1];
                    let gi = miv[jj] + fpenalty;
                    let gj = m[jj] + fpenalty;
                    let taki = gi > diag;
                    let w1v = if taki { gi } else { diag };
                    let ip1 = if taki { -(jj as i32 - mpiv[jj]) } else { 0 };
                    let takj = gj > w1v;
                    let w2v = if takj { gj } else { w1v };
                    let ip2 = if takj { i as i32 - mp[jj] } else { ip1 };
                    let takm = diag > m[jj];

                    if takm {
                        mp[jj] = i as i32 - 1;
                    }
                    m[jj] = if takm { diag } else { m[jj] } + fpenalty_ex;

                    if maxwm < w2v {
                        maxwm = w2v;
                        endali = i;
                        endalj = jj;
                    }

                    let mut w2v = w2v;
                    let mut ip2 = ip2;
                    let stop = w2v < localthr;
                    if stop {
                        ip2 = localstop;
                        w2v = localthr2;
                    }
                    if let Some(wp) = warp.as_mut() {
                        let (nw, code) =
                            wp.cell(i, jj, w2v, currentw[jj], fpenalty_shift, fpenalty_ex);
                        w2v = nw;
                        if let Some(c) = code {
                            ip2 = c;
                        }
                    }
                    ijprow[jj] = ip2;
                    currentw[jj] += w2v;
                }
            }
            if let Some(wp) = warp.as_mut() {
                wp.end_row(lastj);
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        let mut result = LocalAln {
            score: maxwm,
            a1: Vec::new(),
            a2: Vec::new(),
            off1: 0,
            off2: 0,
        };

        if ijp[endali][endalj] == localstop {
            result.score = 0.0;
            eprintln!("maxwm <- 0.0 ");
        } else {
            let (a1, a2, off1, off2) = ltracking(
                seq1,
                seq2,
                &ijp,
                endali,
                endalj,
                localstop,
                b'-',
                warp.as_ref(),
            );
            result.a1 = a1;
            result.a2 = a2;
            result.off1 = off1;
            result.off2 = off2;
        }

        self.amino_dynamicmtx = mtx;
        self.w1 = w1;
        self.w2 = w2;
        self.initverticalw = initverticalw;
        self.lastverticalw = lastverticalw;
        self.m = m;
        self.mp = mp;
        self.miv = miv;
        self.mpiv = mpiv;
        self.ijp = ijp;

        result
    }

    /// `L__align11_noalign()`: the same recurrence without traceback.
    pub fn l_align11_noalign(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        seq1: &[u8],
        seq2: &[u8],
    ) -> f64 {
        let lgth1 = seq1.len();
        let lgth2 = seq2.len();
        let localthr = -(cn.offset as f64);
        let localthr2 = localthr;
        let fpenalty = cn.penalty as f64;
        let fpenalty_ex = cn.penalty_ex as f64;

        if lgth1 > self.n_orlgth1 || lgth2 > self.n_orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.n_orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.n_orlgth2) + 100;
            self.n_w1 = vec![0.0; ll2 + 2];
            self.n_w2 = vec![0.0; ll2 + 2];
            self.n_initverticalw = vec![0.0; ll1 + 2];
            self.n_lastverticalw = vec![0.0; ll1 + 2];
            self.n_m = vec![0.0; ll2 + 2];
            self.n_amino_dynamicmtx = vec![vec![0.0; 0x100]; 0x100];
            self.n_orlgth1 = ll1 - 100;
            self.n_orlgth2 = ll2 - 100;
        }
        let mut mtx = std::mem::take(&mut self.n_amino_dynamicmtx);
        fill_amino_dynamicmtx(cn, &mut mtx, n_dynamicmtx);
        let mut w1 = std::mem::take(&mut self.n_w1);
        let mut w2 = std::mem::take(&mut self.n_w2);
        let mut initverticalw = std::mem::take(&mut self.n_initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.n_lastverticalw);
        let mut m = std::mem::take(&mut self.n_m);

        {
            let row = &mtx[seq2[0] as usize];
            for i in 0..lgth1 {
                initverticalw[i] = row[seq1[i] as usize];
            }
        }
        {
            let row = &mtx[seq1[0] as usize];
            for j in 0..lgth2 {
                w1[j] = row[seq2[j] as usize];
            }
        }
        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1];
        }
        lastverticalw[0] = w1[lgth2 - 1];

        let mut maxwm = -999999999.9f64;
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
                let c1 = if i < lgth1 { seq1[i] } else { 0u8 };
                let row = &mtx[c1 as usize];
                for j in 0..lgth2 {
                    currentw[j] = row[seq2[j] as usize];
                }
            }
            currentw[0] = initverticalw[i];

            let mut mi = previousw[0];
            for j in 1..lastj {
                let mut wm = previousw[j - 1];
                let g = mi + fpenalty;
                if g > wm {
                    wm = g;
                }
                if previousw[j - 1] > mi {
                    mi = previousw[j - 1];
                }
                mi += fpenalty_ex;

                let g = m[j] + fpenalty;
                if g > wm {
                    wm = g;
                }
                if previousw[j - 1] > m[j] {
                    m[j] = previousw[j - 1];
                }
                m[j] += fpenalty_ex;

                if maxwm < wm {
                    maxwm = wm;
                }
                if wm < localthr {
                    wm = localthr2;
                }
                currentw[j] += wm;
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        self.n_amino_dynamicmtx = mtx;
        self.n_w1 = w1;
        self.n_w2 = w2;
        self.n_initverticalw = initverticalw;
        self.n_lastverticalw = lastverticalw;
        self.n_m = m;

        if maxwm < localthr {
            eprintln!("maxwm <- 0.0 ");
            return 0.0;
        }
        maxwm
    }

    /// `G__align11_noalign()` from `core/Galign11.c`: global score, no traceback.
    pub fn g_align11_noalign(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        penal: i32,
        penal_ex: i32,
        seq1: &[u8],
        seq2: &[u8],
    ) -> f64 {
        let lgth1 = seq1.len();
        let lgth2 = seq2.len();
        let fpenalty = penal as f64;
        let fpenalty_ex = penal_ex as f64;

        if lgth1 > self.g_orlgth1 || lgth2 > self.g_orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.g_orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.g_orlgth2) + 100;
            self.g_w1 = vec![0.0; ll2 + 2];
            self.g_w2 = vec![0.0; ll2 + 2];
            self.g_initverticalw = vec![0.0; ll1 + 2];
            self.g_lastverticalw = vec![0.0; ll1 + 2];
            self.g_m = vec![0.0; ll2 + 2];
            self.g_amino_dynamicmtx = vec![vec![0.0; 0x100]; 0x100];
            self.g_orlgth1 = ll1 - 100;
            self.g_orlgth2 = ll2 - 100;
        }
        let mut mtx = std::mem::take(&mut self.g_amino_dynamicmtx);
        fill_amino_dynamicmtx(cn, &mut mtx, n_dynamicmtx);
        let mut w1 = std::mem::take(&mut self.g_w1);
        let mut w2 = std::mem::take(&mut self.g_w2);
        let mut initverticalw = std::mem::take(&mut self.g_initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.g_lastverticalw);
        let mut m = std::mem::take(&mut self.g_m);

        {
            let row = &mtx[seq2[0] as usize];
            for i in 0..lgth1 {
                initverticalw[i] = row[seq1[i] as usize];
            }
        }
        {
            let row = &mtx[seq1[0] as usize];
            for j in 0..lgth2 {
                w1[j] = row[seq2[j] as usize];
            }
        }
        for i in 1..lgth1 + 1 {
            initverticalw[i] += fpenalty;
        }
        for j in 1..lgth2 + 1 {
            w1[j] += fpenalty;
        }
        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1];
        }
        lastverticalw[0] = w1[lgth2 - 1];

        let lastj = lgth2 + 1;
        let mut cur_is_w1 = true;
        let mut wm = 0.0f64;
        for i in 1..lgth1 + 1 {
            cur_is_w1 = !cur_is_w1;
            let (currentw, previousw) = if cur_is_w1 {
                (&mut w1, &mut w2)
            } else {
                (&mut w2, &mut w1)
            };
            previousw[0] = initverticalw[i - 1];
            {
                let c1 = if i < lgth1 { seq1[i] } else { 0u8 };
                let row = &mtx[c1 as usize];
                for j in 0..lgth2 {
                    currentw[j] = row[seq2[j] as usize];
                }
            }
            currentw[0] = initverticalw[i];

            let mut mi = previousw[0];
            // The last row and the last column charge no extension (2018/May/11).
            let fpenalty_ex_i = if i < lgth1 { fpenalty_ex } else { 0.0 };
            for j in 1..lastj {
                wm = previousw[j - 1];
                let g = mi + fpenalty;
                if g > wm {
                    wm = g;
                }
                if previousw[j - 1] > mi {
                    mi = previousw[j - 1];
                }
                mi += fpenalty_ex_i;

                let g = m[j] + fpenalty;
                if g > wm {
                    wm = g;
                }
                if previousw[j - 1] > m[j] {
                    m[j] = previousw[j - 1];
                }
                if j < lgth2 {
                    m[j] += fpenalty_ex;
                }

                currentw[j] += wm;
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        self.g_amino_dynamicmtx = mtx;
        self.g_w1 = w1;
        self.g_w2 = w2;
        self.g_initverticalw = initverticalw;
        self.g_lastverticalw = lastverticalw;
        self.g_m = m;
        wm
    }
}

/// `Ltracking()`: walk back from the best cell to the first `localstop`.
fn ltracking(
    seq1: &[u8],
    seq2: &[u8],
    ijp: &[Vec<i32>],
    endi: usize,
    endj: usize,
    localstop: i32,
    gap: u8,
    warp: Option<&Warp>,
) -> (Vec<u8>, Vec<u8>, usize, usize) {
    let lgth1 = seq1.len();
    let lgth2 = seq2.len();

    // ijp[i][0] and ijp[0][j] are set to localstop -- but the caller's matrix
    // is shared, so make a local view of just those two edges.
    let edge_is_stop = |i: usize, j: usize| -> i32 {
        if i == 0 || j == 0 {
            localstop
        } else {
            ijp[i][j]
        }
    };

    let mut r1: Vec<u8> = Vec::new();
    let mut r2: Vec<u8> = Vec::new();
    let mut iin = endi as i32;
    let mut jin = endj as i32;
    let mut ifi = 0i32;
    let mut jfi = 0i32;
    let limk = (lgth1 + lgth2) as i32;
    let mut k = 0i32;
    let warpbase = warp.map(|w| w.warpbase).unwrap_or(i32::MAX);
    while k <= limk {
        let v = edge_is_stop(iin as usize, jin as usize);
        if v >= warpbase && v != localstop {
            let (a, b) = warp.unwrap().target(v);
            ifi = a;
            jfi = b;
        } else if v < 0 {
            ifi = iin - 1;
            jfi = jin + v;
        } else if v > 0 {
            ifi = iin - v;
            jfi = jin - 1;
        } else {
            ifi = iin - 1;
            jfi = jin - 1;
        }
        // v == localstop takes the v > 0 branch and lands outside the matrix;
        // the loop below breaks on the `ijp[ifi][jfi] == localstop` test before
        // that can matter, exactly as in the C.

        #[allow(unused_assignments)]
        if ifi == -warpbase && jfi == -warpbase {
            let mut l = iin;
            while {
                l -= 1;
                l >= 0
            } {
                r1.push(seq1[l as usize]);
                r2.push(gap);
                k += 1;
            }
            let mut l = jin;
            while {
                l -= 1;
                l >= 0
            } {
                r1.push(gap);
                r2.push(seq2[l as usize]);
                k += 1;
            }
            break;
        }

        let mut l = iin - ifi;
        while {
            l -= 1;
            l > 0
        } {
            r1.push(seq1[(ifi + l) as usize]);
            r2.push(gap);
            k += 1;
        }
        let mut l = jin - jfi;
        while {
            l -= 1;
            l > 0
        } {
            r1.push(gap);
            r2.push(seq2[(jfi + l) as usize]);
            k += 1;
        }

        if iin <= 0 || jin <= 0 {
            break;
        }
        r1.push(seq1[ifi as usize]);
        r2.push(seq2[jfi as usize]);
        if edge_is_stop(ifi as usize, jfi as usize) == localstop {
            break;
        }
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1;
    }
    r1.reverse();
    r2.reverse();
    let off1 = if ifi == -1 { 0 } else { ifi as usize };
    let off2 = if jfi == -1 { 0 } else { jfi as usize };
    (r1, r2, off1, off2)
}
