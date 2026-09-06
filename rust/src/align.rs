//! The dynamic-programming core.
//!
//! * [`AlignState::a_align`] is `A__align()` from `core/Salignmm.c` — the
//!   group-to-group affine DP with position-specific gap penalties that does
//!   the actual work of a progressive alignment step.
//! * [`AlignState::g_align11`] is `G__align11()` from `core/Galign11.c`, the
//!   1-vs-1 special case that `treebase()` uses when both groups are single
//!   sequences.
//!
//! Both keep their scratch buffers between calls exactly as the C `static`s do:
//! the buffers grow but are never cleared, so a stale value left beyond the
//! current sequence length is read back the same way it is in C.  That is not
//! an accident of the port — several loops in `A__align` index one element past
//! the current length.

use crate::constants::Constants;
use crate::profile;
use crate::warp::Warp;

// Salignmm.c / Galign11.c
const TERMGAPFAC: f64 = 0.0;
const TERMGAPFAC_EX: f64 = 0.0;

/// A cached profile for one subtree: `nalphabets` composition rows, then gap
/// frequency, opening-gap counts and final-gap counts (`cpmxhist[]` in the C).
pub type CpmxHist = Vec<Vec<f64>>;

fn max_usize(a: usize, b: usize) -> usize {
    if a > b {
        a
    } else {
        b
    }
}

/// All the `static TLS` state of `A__align()`, `G__align11()` and the shared
/// `commonIP` traceback matrix.
pub struct AlignState {
    // --- A__align ---
    orlgth1: usize,
    orlgth2: usize,
    w1: Vec<f64>,
    w2: Vec<f64>,
    initverticalw: Vec<f64>,
    lastverticalw: Vec<f64>,
    m: Vec<f64>,
    mp: Vec<i32>,
    ogcp1: Vec<f64>,
    ogcp1o: Vec<f64>,
    ogcp2: Vec<f64>,
    ogcp2o: Vec<f64>,
    fgcp1: Vec<f64>,
    fgcp1o: Vec<f64>,
    fgcp2: Vec<f64>,
    fgcp2o: Vec<f64>,
    cpmx1: Vec<Vec<f64>>,
    cpmx2: Vec<Vec<f64>>,
    gapfreq1: Vec<f64>,
    gapfreq2: Vec<f64>,
    doublework: Vec<Vec<f64>>,
    intwork: Vec<Vec<i32>>,
    previousfirstlen: i32,
    previousicyc: i32,
    previousfirstmem: i32,
    previouscall: i32,

    // --- commonIP, shared by both aligners ---
    ijp: Vec<Vec<i32>>,
    common_alloc1: usize,
    common_alloc2: usize,

    // --- G__align11 ---
    g_orlgth1: usize,
    g_orlgth2: usize,
    g_w1: Vec<f64>,
    g_w2: Vec<f64>,
    g_initverticalw: Vec<f64>,
    g_lastverticalw: Vec<f64>,
    g_m: Vec<f64>,
    g_mp: Vec<i32>,
    g_amino_dynamicmtx: Vec<Vec<f64>>,
}

impl Default for AlignState {
    fn default() -> Self {
        Self::new()
    }
}

impl AlignState {
    pub fn new() -> AlignState {
        AlignState {
            orlgth1: 0,
            orlgth2: 0,
            w1: Vec::new(),
            w2: Vec::new(),
            initverticalw: Vec::new(),
            lastverticalw: Vec::new(),
            m: Vec::new(),
            mp: Vec::new(),
            ogcp1: Vec::new(),
            ogcp1o: Vec::new(),
            ogcp2: Vec::new(),
            ogcp2o: Vec::new(),
            fgcp1: Vec::new(),
            fgcp1o: Vec::new(),
            fgcp2: Vec::new(),
            fgcp2o: Vec::new(),
            cpmx1: Vec::new(),
            cpmx2: Vec::new(),
            gapfreq1: Vec::new(),
            gapfreq2: Vec::new(),
            doublework: Vec::new(),
            intwork: Vec::new(),
            previousfirstlen: -1,
            previousicyc: -1,
            previousfirstmem: 0,
            previouscall: 0,
            ijp: Vec::new(),
            common_alloc1: 0,
            common_alloc2: 0,
            g_orlgth1: 0,
            g_orlgth2: 0,
            g_w1: Vec::new(),
            g_w2: Vec::new(),
            g_initverticalw: Vec::new(),
            g_lastverticalw: Vec::new(),
            g_m: Vec::new(),
            g_mp: Vec::new(),
            g_amino_dynamicmtx: Vec::new(),
        }
    }

    /// `commonIP` growth, driven by whichever aligner's `orlgth` is larger.
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

/// `match_calc()`: profile-profile score of column `i1` of group 1 against
/// every column of group 2.
#[allow(clippy::too_many_arguments)]
/// Shared with `partalign`: `partSalignmm.c`'s `match_calc` is this same
/// function, differing only in that it reads the global
/// `n_dis_consweight_multi` where `A__align` is handed `n_dynamicmtx`.
/// Exposed rather than copied -- two transliterations of one DP inner
/// loop would be free to drift apart.
pub(crate) fn match_calc(
    n_dynamicmtx: &[Vec<f64>],
    out: &mut [f64],
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    i1: usize,
    lgth2: usize,
    doublework: &mut [Vec<f64>],
    intwork: &mut [Vec<i32>],
    initialize: bool,
) {
    // The alphabet size is however many rows the matrix has: 26 normally,
    // 0x100 under `--text`.
    let nalphabets = n_dynamicmtx.len();
    if initialize {
        // The C marks the end of each column's non-zero entries with a `-1`
        // and re-tests for it on every step of the innermost loop below. Here
        // the length carries that instead: `clear` keeps the capacity `ensure`
        // allocated, so the pushes never reallocate, and the hot loop walks a
        // slice whose end is known rather than checking a sentinel per element.
        // Same entries, same order, one branch per element less.
        for j in 0..lgth2 {
            let vals = &mut doublework[j];
            let codes = &mut intwork[j];
            vals.clear();
            codes.clear();
            for l in 0..nalphabets {
                let v = cpmx2[l][j];
                if v != 0.0 {
                    vals.push(v);
                    codes.push(l as i32);
                }
            }
        }
    }

    // The C writes this l-outer, j-inner:
    //     for( l=0; l<nalphabets; l++ ) { scarr[l] = 0.0;
    //         for( j=0; j<nalphabets; j++ )
    //             scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1]; }
    // Transposed here. Each `scarr[l]` still accumulates over j in increasing
    // j, so every sum is the same sequence of additions in the same order --
    // the one thing that must not move. What changes is the addressing:
    // `n_dynamicmtx[j]` is walked along its row rather than down a column, and
    // `cpmx1[j][i1]` is loaded once per j instead of once per (l, j).
    // On the stack rather than a fresh `vec!` per call: this runs once per DP
    // row, which is ~128 thousand times on a 36-sequence refinement, and the
    // alphabet is 26 (0x100 under `--text`), so the array is sized for the
    // larger case and only `nalphabets` of it is ever touched.
    let mut scarr = [0.0f64; 0x100];
    for j in 0..nalphabets {
        let c = cpmx1[j][i1];
        let row = &n_dynamicmtx[j];
        for (l, sc) in scarr[..nalphabets].iter_mut().enumerate() {
            *sc += row[l] * c;
        }
    }

    // The C fetches each column's two row pointers once and then walks them
    // (`cpmxpdnpt = *cpmxpdnptpt++;` ... `while( *cpmxpdnpt>-1 ) *matchpt +=
    // scarr[*cpmxpdnpt++] * *cpmxpdpt++;`). Indexing `intwork[j][k]` inside
    // the inner loop instead re-walks the outer `Vec` on every step, and does
    // it twice per iteration -- once for the test and once for the value.
    // Same arithmetic in the same order; only the addressing changes.
    for (j, o) in out.iter_mut().enumerate().take(lgth2) {
        let codes = &intwork[j];
        let vals = &doublework[j];
        let mut acc = 0.0f64;
        for (&code, &val) in codes.iter().zip(vals.iter()) {
            // `& 0xff` costs one instruction and lets the compiler drop the
            // bounds check on an array of exactly 0x100: a code is an alphabet
            // index, always below `nalphabets`, so the mask never changes a
            // value. Without it this innermost line carries a branch that can
            // never be taken.
            acc += scarr[(code as usize) & 0xff] * val;
        }
        *o = acc;
    }
}

/// The result of `Atracking()`: the score plus the two gap tables.
struct Track {
    wm: f64,
    gt1: Vec<u8>,
    gt2: Vec<u8>,
    ngap1: bool,
    ngap2: bool,
}

/// `Atracking()` from Salignmm.c.  The warp branch is absent because a
/// `trywarp` configuration is rejected before we get here.
#[allow(clippy::too_many_arguments)]
fn atracking(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    fpenalty: f64,
    fpenalty_ex: f64,
    lgth1: usize,
    lgth2: usize,
    ijp: &mut [Vec<i32>],
    tailgp: bool,
    warp: Option<&Warp>,
) -> Track {
    let mut wm = 0.0f64;

    if !tailgp {
        wm = lasthorizontalw[lgth2 - 1] - 1.0;
        for j in (0..lgth2 - 1).rev() {
            let g = lasthorizontalw[j]
                + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2 - 1 - j) as f64 * TERMGAPFAC_EX);
            if g > wm {
                wm = g;
                ijp[lgth1][lgth2] = -((lgth2 - j) as i32);
            }
        }
        for i in (0..lgth1 - 1).rev() {
            let g = lastverticalw[i]
                + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1 - 1 - i) as f64 * TERMGAPFAC_EX);
            if g > wm {
                wm = g;
                ijp[lgth1][lgth2] = (lgth1 - i) as i32;
            }
        }
        if lasthorizontalw[lgth2 - 1] > wm {
            wm = lasthorizontalw[lgth2 - 1];
            ijp[lgth1][lgth2] = 0;
        }
    }

    for (i, row) in ijp.iter_mut().enumerate().take(lgth1 + 1) {
        row[0] = i as i32 + 1;
    }
    for j in 0..lgth2 + 1 {
        ijp[0][j] = -(j as i32 + 1);
    }

    // The C builds the tables backwards into a buffer; we push and reverse.
    let mut rev1: Vec<u8> = Vec::with_capacity(lgth1 + lgth2 + 1);
    let mut rev2: Vec<u8> = Vec::with_capacity(lgth1 + lgth2 + 1);

    let mut iin = lgth1 as i32;
    let mut jin = lgth2 as i32;
    let limk = (lgth1 + lgth2 + 1) as i32;
    let mut k = 0i32;
    let warpbase = warp.map(|w| w.warpbase).unwrap_or(i32::MAX);
    while k < limk {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi) = if v >= warpbase {
            warp.unwrap().target(v)
        } else if v < 0 {
            (iin - 1, jin + v)
        } else if v > 0 {
            (iin - v, jin - 1)
        } else {
            (iin - 1, jin - 1)
        };

        #[allow(unused_assignments)]
        if ifi == -warpbase && jfi == -warpbase {
            // A jump straight off the front: everything left is terminal gap.
            let mut l = iin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'o');
                rev2.push(b'-');
                k += 1;
            }
            let mut l = jin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'-');
                rev2.push(b'o');
                k += 1;
            }
            break;
        }

        let mut l = iin - ifi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'o');
            rev2.push(b'-');
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'-');
            rev2.push(b'o');
            k += 1;
        }

        if iin <= 0 || jin <= 0 {
            break;
        }
        rev1.push(b'o');
        rev2.push(b'o');
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1;
    }

    rev1.reverse();
    rev2.reverse();
    let ngap1 = rev1.contains(&b'-');
    let ngap2 = rev2.contains(&b'-');
    Track {
        wm,
        gt1: rev1,
        gt2: rev2,
        ngap1,
        ngap2,
    }
}

/// `createcpmxresult()`: merged composition profile of the two aligned groups.
pub(crate) fn create_cpmx_result(
    limk: usize,
    eff1: f64,
    eff2: f64,
    cpmx1: &[Vec<f64>],
    cpmx2: &[Vec<f64>],
    gt1: &[u8],
    gt2: &[u8],
) -> Vec<Vec<f64>> {
    let alen = gt1.len();
    // `nalphabets + 3`: the profile rows plus gapfreq, opening and final gap
    // counts, which the caller appends. Only the profile rows are built here.
    let mut out = vec![vec![0.0f64; limk + 1]; cpmx1.len()];
    for (i, row) in out.iter_mut().enumerate() {
        let mut p = 0usize;
        for j in 0..alen {
            if gt1[j] != b'-' {
                row[j] += cpmx1[i][p] * eff1;
                p += 1;
            }
        }
        let mut p = 0usize;
        for j in 0..alen {
            if gt2[j] != b'-' {
                row[j] += cpmx2[i][p] * eff2;
                p += 1;
            }
        }
    }
    out
}

/// `creategapfreqresult()`.
pub(crate) fn create_gapfreq_result(
    limk: usize,
    eff1: f64,
    eff2: f64,
    gapf1: &[f64],
    gapf2: &[f64],
    gt1: &[u8],
    gt2: &[u8],
) -> Vec<f64> {
    let alen = gt1.len();
    let mut out = vec![0.0f64; limk + 1];
    // The C reads gaptable1[alen], i.e. the NUL terminator, which is not '-',
    // so exactly one extra value of gapf1 is consumed.
    let mut p = 0usize;
    for j in 0..alen + 1 {
        let c = if j < alen { gt1[j] } else { 0 };
        if c != b'-' {
            out[j] += gapf1[p] * eff1;
            p += 1;
        }
    }
    let mut p = 0usize;
    for j in 0..alen {
        if gt2[j] != b'-' {
            out[j] += gapf2[p] * eff2;
            p += 1;
        }
    }
    out[alen] = 1.0;
    out
}

/// `createogresult()`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn create_og_result(
    limk: usize,
    eff1: f64,
    eff2: f64,
    ori1: &[f64],
    ori2: &[f64],
    gf1: &[f64],
    gf2: &[f64],
    gt1: &[u8],
    gt2: &[u8],
) -> Vec<f64> {
    let alen = gt1.len();
    let mut out = vec![0.0f64; limk + 1];
    for (gt, ori, gf, eff) in [(gt1, ori1, gf1, eff1), (gt2, ori2, gf2, eff2)] {
        let mut p = 0usize;
        for j in 0..alen {
            if gt[j] == b'-' {
                if j == 0 {
                    out[j] += eff;
                } else if gt[j - 1] != b'-' {
                    out[j] += gf[p - 1] * eff;
                }
            } else {
                if j == 0 || gt[j - 1] != b'-' {
                    out[j] += ori[p] * eff;
                }
                p += 1;
            }
        }
    }
    out
}

/// `createfgresult()`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn create_fg_result(
    limk: usize,
    eff1: f64,
    eff2: f64,
    ori1: &[f64],
    ori2: &[f64],
    gf1: &[f64],
    gf2: &[f64],
    gt1: &[u8],
    gt2: &[u8],
) -> Vec<f64> {
    let alen = gt1.len();
    let mut out = vec![0.0f64; limk + 1];
    for (gt, ori, gf, eff) in [(gt1, ori1, gf1, eff1), (gt2, ori2, gf2, eff2)] {
        let mut p = 0usize;
        for j in 0..alen {
            // gaptable[alen] is the NUL terminator in the C, hence != '-'.
            let next_is_gap = j + 1 < alen && gt[j + 1] == b'-';
            if gt[j] == b'-' {
                if j == alen - 1 {
                    out[j] += eff;
                } else if !next_is_gap {
                    out[j] += gf[p] * eff;
                }
            } else {
                if !next_is_gap {
                    out[j] += ori[p] * eff;
                }
                p += 1;
            }
        }
    }
    out
}

impl AlignState {
    /// `A__align()`, with no constraints and no warp.
    ///
    /// Two call sites reach it and they take different branches of the C:
    ///
    /// * `treebase()` supplies a profile-history slot (`want_cpmxresult`) and
    ///   no boundary-gap patterns, so a child's cached profile can be reused
    ///   and the parent's is produced.
    /// * `Falign()` aligns one block between two anchors, so it supplies the
    ///   `gaps` patterns describing what lies just outside the block, and no
    ///   history slot.
    ///
    /// `group1` / `group2` hold the aligned members of each group; on return
    /// they carry the new alignment. The returned profile is the parent's
    /// `cpmxhist` entry, or `None` when the C would store a null.
    #[allow(clippy::too_many_arguments)]
    pub fn a_align(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        penalty_l: i32,
        penalty_ex_l: i32,
        group1: &mut [Vec<u8>],
        group2: &mut [Vec<u8>],
        eff1: &[f64],
        eff2: &[f64],
        headgp: bool,
        tailgp: bool,
        gaps: Option<(&[u8], &[u8], &[u8], &[u8])>,
        // `imp` is `impmtx` when the caller runs with constraints: the
        // per-cell match bonus derived from the pairwise local alignments.
        imp: Option<&[Vec<f64>]>,
        firstmem: i32,
        calledbyfulltreebase: i32,
        cpmxchild0: Option<CpmxHist>,
        cpmxchild1: Option<CpmxHist>,
        orieff1: f64,
        orieff2: f64,
        want_cpmxresult: bool,
        newgapstr: u8,
    ) -> (f64, Option<CpmxHist>) {
        let icyc = group1.len();
        let jcyc = group2.len();
        let lgth1 = group1[0].len();
        let lgth2 = group2[0].len();

        if lgth1 == 0 && lgth2 == 0 {
            return (0.0, None);
        }
        if lgth1 == 0 {
            for s in group1.iter_mut() {
                *s = vec![newgapstr; lgth2];
            }
            return (0.0, None);
        }
        if lgth2 == 0 {
            for s in group2.iter_mut() {
                *s = vec![newgapstr; lgth1];
            }
            return (0.0, None);
        }

        let fpenalty = penalty_l as f64;
        let fpenalty_ex = penalty_ex_l as f64;
        let fpenalty_shift = cn.penalty_shift as f64;

        let mut warp = if cn.trywarp {
            if !headgp || !tailgp {
                eprintln!("At present, headgp and tailgp must be 1 to allow shift.");
                std::process::exit(1);
            }
            Some(Warp::new(lgth1, lgth2))
        } else {
            None
        };

        // --- buffer growth, mirroring the C `static` reallocation rule ---
        if lgth1 > self.orlgth1 || lgth2 > self.orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.orlgth2) + 100;
            self.w1 = vec![0.0; ll2 + 2];
            self.w2 = vec![0.0; ll2 + 2];
            self.initverticalw = vec![0.0; ll1 + 2];
            self.lastverticalw = vec![0.0; ll1 + 2];
            self.m = vec![0.0; ll2 + 2];
            self.mp = vec![0; ll2 + 2];
            self.ogcp1 = vec![0.0; ll1 + 2];
            self.ogcp1o = vec![0.0; ll1 + 2];
            self.ogcp2 = vec![0.0; ll2 + 2];
            self.ogcp2o = vec![0.0; ll2 + 2];
            self.fgcp1 = vec![0.0; ll1 + 2];
            self.fgcp1o = vec![0.0; ll1 + 2];
            self.fgcp2 = vec![0.0; ll2 + 2];
            self.fgcp2o = vec![0.0; ll2 + 2];
            self.cpmx1 = vec![vec![0.0; ll1 + 2]; cn.nalphabets];
            self.cpmx2 = vec![vec![0.0; ll2 + 2]; cn.nalphabets];
            self.previousfirstlen = -1;
            self.previousicyc = -1;
            self.gapfreq1 = vec![0.0; ll1 + 2];
            self.gapfreq2 = vec![0.0; ll2 + 2];
            let mw = max_usize(ll1, ll2) + 2;
            self.doublework = vec![vec![0.0; cn.nalphabets]; mw];
            self.intwork = vec![vec![0; cn.nalphabets + 1]; mw];
            self.orlgth1 = ll1 - 100;
            self.orlgth2 = ll2 - 100;
        }
        self.ensure_ijp(self.orlgth1, self.orlgth2);

        // The C's `previouscall` holds the *previous* call's
        // `calledbyfulltreebase`, which Falign passes as -1: truthy in C, so
        // the flag is stored as-is and tested against zero.
        let reuseprofiles = calledbyfulltreebase == 1
            && self.previouscall != 0
            && firstmem >= 0
            && firstmem == self.previousfirstmem
            && lgth1 as i32 == self.previousfirstlen
            && icyc as i32 == self.previousicyc + 1;

        // Move every scratch buffer out of `self` so the borrow checker lets us
        // hold views into them while writing others.
        let mut w1 = std::mem::take(&mut self.w1);
        let mut w2 = std::mem::take(&mut self.w2);
        let mut initverticalw = std::mem::take(&mut self.initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.lastverticalw);
        let mut m = std::mem::take(&mut self.m);
        let mut mp = std::mem::take(&mut self.mp);
        let mut ogcp1 = std::mem::take(&mut self.ogcp1);
        let mut ogcp1o = std::mem::take(&mut self.ogcp1o);
        let mut ogcp2 = std::mem::take(&mut self.ogcp2);
        let mut ogcp2o = std::mem::take(&mut self.ogcp2o);
        let mut fgcp1 = std::mem::take(&mut self.fgcp1);
        let mut fgcp1o = std::mem::take(&mut self.fgcp1o);
        let mut fgcp2 = std::mem::take(&mut self.fgcp2);
        let mut fgcp2o = std::mem::take(&mut self.fgcp2o);
        let mut cpmx1 = std::mem::take(&mut self.cpmx1);
        let mut cpmx2 = std::mem::take(&mut self.cpmx2);
        let mut gapfreq1 = std::mem::take(&mut self.gapfreq1);
        let mut gapfreq2 = std::mem::take(&mut self.gapfreq2);
        let mut doublework = std::mem::take(&mut self.doublework);
        let mut intwork = std::mem::take(&mut self.intwork);
        let mut ijp = std::mem::take(&mut self.ijp);

        // --- profiles ---
        let g1: Vec<&[u8]> = group1.iter().map(|s| s.as_slice()).collect();
        let g2: Vec<&[u8]> = group2.iter().map(|s| s.as_slice()).collect();

        let mut hist1: Option<CpmxHist> = None;
        let mut hist2: Option<CpmxHist> = None;
        let headgapfreq1: f64;
        let headgapfreq2: f64;

        // The C selects this branch on `cpmxresult && specificityconsideration
        // == 0.0`: with --allowshift the score matrix varies per node, so a
        // child's cached profile cannot be reused -- but the parent's is still
        // produced.
        if want_cpmxresult && cn.specificityconsideration == 0.0 {
            assert!(
                gaps.is_none(),
                "The combination of sgap1 and cpmxhist is not supported"
            );
            hist1 = cpmxchild0;
            if hist1.is_none() {
                profile::cpmx_calc_new(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
                profile::gapcountf(&mut gapfreq1, &g1, icyc, eff1, lgth1);
                for v in gapfreq1.iter_mut().take(lgth1) {
                    *v = 1.0 - *v;
                }
                profile::st_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1);
                profile::st_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
            }
            hist2 = cpmxchild1;
            if hist2.is_none() {
                profile::cpmx_calc_new(cn, &g2, &mut cpmx2, eff2, lgth2, jcyc);
                profile::gapcountf(&mut gapfreq2, &g2, jcyc, eff2, lgth2);
                for v in gapfreq2.iter_mut().take(lgth2) {
                    *v = 1.0 - *v;
                }
                profile::st_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2);
                profile::st_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
            }

            // Head and tail gaps are taken to be absent on this branch.
            headgapfreq1 = 1.0;
            headgapfreq2 = 1.0;
            match hist1.as_mut() {
                Some(h) => h[cn.nalphabets][lgth1] = 1.0,
                None => gapfreq1[lgth1] = 1.0,
            }
            match hist2.as_mut() {
                Some(h) => h[cn.nalphabets][lgth2] = 1.0,
                None => gapfreq2[lgth2] = 1.0,
            }
        } else {
            // The Falign block branch, and the --allowshift branch. When the
            // previous call built the same group one member smaller, the C
            // updates the profile incrementally rather than recomputing it --
            // a different summation order, so it has to be reproduced.
            if reuseprofiles {
                profile::cpmx_calc_add(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
            } else {
                profile::cpmx_calc_new(cn, &g1, &mut cpmx1, eff1, lgth1, icyc);
            }
            profile::cpmx_calc_new(cn, &g2, &mut cpmx2, eff2, lgth2, jcyc);

            let mut hgf1 = 0.0f64;
            let mut hgf2 = 0.0f64;
            match gaps {
                Some((sgap1, sgap2, egap1, egap2)) => {
                    profile::new_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1, sgap1);
                    profile::new_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
                    profile::new_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2, sgap2);
                    profile::new_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
                    hgf1 = profile::outgapcount(icyc, sgap1, eff1);
                    hgf2 = profile::outgapcount(jcyc, sgap2, eff2);
                    gapfreq1[lgth1] = profile::outgapcount(icyc, egap1, eff1);
                    gapfreq2[lgth2] = profile::outgapcount(jcyc, egap2, eff2);
                }
                None => {
                    if reuseprofiles {
                        profile::st_opening_gap_add(&mut ogcp1o, icyc, &g1, eff1, lgth1);
                        profile::st_final_gap_add(&mut fgcp1o, icyc, &g1, eff1, lgth1);
                    } else {
                        profile::st_opening_gap_count(&mut ogcp1o, icyc, &g1, eff1, lgth1);
                        profile::st_final_gap_count(&mut fgcp1o, icyc, &g1, eff1, lgth1);
                    }
                    profile::st_opening_gap_count(&mut ogcp2o, jcyc, &g2, eff2, lgth2);
                    profile::st_final_gap_count(&mut fgcp2o, jcyc, &g2, eff2, lgth2);
                    gapfreq1[lgth1] = 0.0;
                    gapfreq2[lgth2] = 0.0;
                }
            }

            // `--legacygapcost` drops the weighting entirely. Note the C only
            // consults it here: the cpmxhist branch above assumes head and
            // tail gaps are absent "regardless of legacygap", so the two
            // branches disagree once the flag is on.
            if !cn.legacygapcost {
                if reuseprofiles {
                    profile::gapcountadd(&mut gapfreq1, &g1, icyc, eff1, lgth1);
                } else {
                    profile::gapcountf(&mut gapfreq1, &g1, icyc, eff1, lgth1);
                }
                profile::gapcountf(&mut gapfreq2, &g2, jcyc, eff2, lgth2);
                for v in gapfreq1.iter_mut().take(lgth1 + 1) {
                    *v = 1.0 - *v;
                }
                for v in gapfreq2.iter_mut().take(lgth2 + 1) {
                    *v = 1.0 - *v;
                }
                headgapfreq1 = 1.0 - hgf1;
                headgapfreq2 = 1.0 - hgf2;
            } else {
                // The counting calls are skipped, not just overwritten: they
                // update a `static` buffer incrementally, so skipping them is
                // what the next call sees.
                for v in gapfreq1.iter_mut().take(lgth1 + 1) {
                    *v = 1.0;
                }
                for v in gapfreq2.iter_mut().take(lgth2 + 1) {
                    *v = 1.0;
                }
                headgapfreq1 = 1.0;
                headgapfreq2 = 1.0;
            }
        }

        let (cpmx1v, gapfreq1v, ogcp1ov, fgcp1ov): (&[Vec<f64>], &[f64], &[f64], &[f64]) =
            match hist1.as_ref() {
                Some(h) => (
                    &h[0..cn.nalphabets],
                    &h[cn.nalphabets],
                    &h[cn.nalphabets + 1],
                    &h[cn.nalphabets + 2],
                ),
                None => (&cpmx1[..], &gapfreq1[..], &ogcp1o[..], &fgcp1o[..]),
            };
        let (cpmx2v, gapfreq2v, ogcp2ov, fgcp2ov): (&[Vec<f64>], &[f64], &[f64], &[f64]) =
            match hist2.as_ref() {
                Some(h) => (
                    &h[0..cn.nalphabets],
                    &h[cn.nalphabets],
                    &h[cn.nalphabets + 1],
                    &h[cn.nalphabets + 2],
                ),
                None => (&cpmx2[..], &gapfreq2[..], &ogcp2o[..], &fgcp2o[..]),
            };

        // --- position-specific gap penalties ---
        for i in 0..lgth1 {
            ogcp1[i] = 0.5 * (1.0 - ogcp1ov[i]) * fpenalty * gapfreq1v[i];
            fgcp1[i] = 0.5 * (1.0 - fgcp1ov[i]) * fpenalty * gapfreq1v[i];
        }
        for i in 0..lgth2 {
            ogcp2[i] = 0.5 * (1.0 - ogcp2ov[i]) * fpenalty * gapfreq2v[i];
            fgcp2[i] = 0.5 * (1.0 - fgcp2ov[i]) * fpenalty * gapfreq2v[i];
        }

        // --- DP: `currentw` starts as w1 and the two swap on every row ---
        let mut cur_is_w1 = true;

        match_calc(
            n_dynamicmtx,
            &mut initverticalw,
            cpmx2v,
            cpmx1v,
            0,
            lgth1,
            &mut doublework,
            &mut intwork,
            true,
        );
        if let Some(im) = imp {
            // imp_match_out_vead_tate( initverticalw, 0, lgth1 )
            for (i, v) in initverticalw.iter_mut().enumerate().take(lgth1) {
                *v += im[i][0];
            }
        }
        match_calc(
            n_dynamicmtx,
            &mut w1,
            cpmx1v,
            cpmx2v,
            0,
            lgth2,
            &mut doublework,
            &mut intwork,
            true,
        );
        if let Some(im) = imp {
            // imp_match_out_vead( currentw, 0, lgth2 )
            for (j, v) in w1.iter_mut().enumerate().take(lgth2) {
                *v += im[0][j];
            }
        }

        if headgp {
            for i in 1..lgth1 + 1 {
                initverticalw[i] += ogcp1[0] * headgapfreq2 + fgcp1[i - 1] * gapfreq2v[0];
                initverticalw[i] += fpenalty_ex * i as f64;
            }
            for j in 1..lgth2 + 1 {
                w1[j] += ogcp2[0] * headgapfreq1 + fgcp2[j - 1] * gapfreq1v[0];
                w1[j] += fpenalty_ex * j as f64;
            }
        } else {
            for i in 1..lgth1 + 1 {
                initverticalw[i] += fpenalty * TERMGAPFAC;
                initverticalw[i] += fpenalty_ex * i as f64 * TERMGAPFAC_EX;
            }
            for j in 1..lgth2 + 1 {
                w1[j] += fpenalty * TERMGAPFAC;
                w1[j] += fpenalty_ex * j as f64 * TERMGAPFAC_EX;
            }
        }

        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1] + ogcp1[1] * gapfreq2v[j - 1];
            mp[j] = 0;
        }
        lastverticalw[0] = w1[lgth2 - 1];

        let lasti = if tailgp { lgth1 + 1 } else { lgth1 };
        let lastj = lgth2 + 1;

        let mut wm = 0.0f64;
        for i in 1..lasti {
            cur_is_w1 = !cur_is_w1;
            let (currentw, previousw) = if cur_is_w1 {
                (&mut w1, &mut w2)
            } else {
                (&mut w2, &mut w1)
            };

            previousw[0] = initverticalw[i - 1];
            match_calc(
                n_dynamicmtx,
                currentw,
                cpmx1v,
                cpmx2v,
                i,
                lgth2,
                &mut doublework,
                &mut intwork,
                false,
            );
            if let Some(im) = imp {
                // imp_match_out_vead( currentw, i, lgth2 )
                for (j, v) in currentw.iter_mut().enumerate().take(lgth2) {
                    *v += im[i][j];
                }
            }
            currentw[0] = initverticalw[i];

            let mut mi = previousw[0] + ogcp2[1] * gapfreq1v[i - 1];
            let mut mpi = 0i32;
            let fgcp1va = fgcp1[i - 1];
            let ogcp1va = ogcp1[i];
            let gf1va = gapfreq1v[i];
            let gf1vapre = gapfreq1v[i - 1];
            let ijprow = &mut ijp[i];

            for j in 1..lastj {
                let mut w = previousw[j - 1];
                ijprow[j] = 0;

                let g = mi + fgcp2[j - 1] * gf1va;
                if g > w {
                    w = g;
                    ijprow[j] = -(j as i32 - mpi);
                }
                let g = previousw[j - 1] + ogcp2[j] * gf1vapre;
                if g >= mi {
                    mi = g;
                    mpi = j as i32 - 1;
                }
                mi += fpenalty_ex;

                let g = m[j] + fgcp1va * gapfreq2v[j];
                if g > w {
                    w = g;
                    ijprow[j] = i as i32 - mp[j];
                }
                let g = previousw[j - 1] + ogcp1va * gapfreq2v[j - 1];
                if g >= m[j] {
                    m[j] = g;
                    mp[j] = i as i32 - 1;
                }
                m[j] += fpenalty_ex;

                if let Some(wp) = warp.as_mut() {
                    let (nw, code) = wp.cell(i, j, w, currentw[j], fpenalty_shift, fpenalty_ex);
                    w = nw;
                    if let Some(c) = code {
                        ijprow[j] = c;
                    }
                }

                currentw[j] += w;
                wm = w;
            }
            if let Some(wp) = warp.as_mut() {
                wp.end_row(lastj);
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        let track = {
            let currentw: &[f64] = if cur_is_w1 { &w1 } else { &w2 };
            match imp {
                // Under constraints the C uses Atracking_localhom, which picks
                // the end cell with `>=` and consults the global `outgap`
                // rather than the `tailgp` argument, and leaves `wm` as the
                // DP's last value.
                Some(im) => atracking_localhom(
                    currentw,
                    &lastverticalw,
                    lgth1,
                    lgth2,
                    &mut ijp,
                    tailgp,
                    im,
                    warp.as_ref(),
                ),
                None => {
                    let t = atracking(
                        currentw,
                        &lastverticalw,
                        fpenalty,
                        fpenalty_ex,
                        lgth1,
                        lgth2,
                        &mut ijp,
                        tailgp,
                        warp.as_ref(),
                    );
                    if !tailgp {
                        wm = t.wm;
                    }
                    t
                }
            }
        };

        // --- profile history for the parent node ---
        let cpmxresult = if want_cpmxresult && icyc + jcyc > 20 {
            let alen = track.gt1.len();
            let limk = ((alen + 1000) / 1000) * 1000;

            let mut totaleff1 = 0.0f64;
            for v in eff1.iter().take(icyc) {
                totaleff1 += *v;
            }
            let mut totaleff2 = 0.0f64;
            for v in eff2.iter().take(jcyc) {
                totaleff2 += *v;
            }
            if (totaleff1 - 1.0).abs() > 0.001 || (totaleff2 - 1.0).abs() > 0.001 {
                panic!("totaleff1 = {}, totaleff2 = {}", totaleff1, totaleff2);
            }
            let totaleff1 = totaleff1 * orieff1 / (orieff1 + orieff2);
            let totaleff2 = totaleff2 * orieff2 / (orieff1 + orieff2);

            let mut res = create_cpmx_result(
                limk, totaleff1, totaleff2, cpmx1v, cpmx2v, &track.gt1, &track.gt2,
            );
            res.push(create_gapfreq_result(
                limk, totaleff1, totaleff2, gapfreq1v, gapfreq2v, &track.gt1, &track.gt2,
            ));
            res.push(create_og_result(
                limk, totaleff1, totaleff2, ogcp1ov, ogcp2ov, gapfreq1v, gapfreq2v, &track.gt1,
                &track.gt2,
            ));
            res.push(create_fg_result(
                limk, totaleff1, totaleff2, fgcp1ov, fgcp2ov, gapfreq1v, gapfreq2v, &track.gt1,
                &track.gt2,
            ));
            Some(res)
        } else {
            None
        };

        // --- write the new alignment back ---
        // The C skips the copy when a group gained no gaps and the profile was
        // reused; `gapireru` with an all-`o` table is a no-op, so a single
        // branch reproduces both cases.
        if track.ngap1 || !reuseprofiles {
            for s in group1.iter_mut() {
                *s = profile::gapireru(s, &track.gt1, newgapstr);
            }
        }
        if track.ngap2 || !reuseprofiles {
            for s in group2.iter_mut() {
                *s = profile::gapireru(s, &track.gt2, newgapstr);
            }
        }

        self.previousfirstlen = lgth1 as i32;
        self.previousfirstmem = firstmem;
        self.previousicyc = icyc as i32;
        self.previouscall = calledbyfulltreebase;

        // Give the scratch buffers back.
        self.w1 = w1;
        self.w2 = w2;
        self.initverticalw = initverticalw;
        self.lastverticalw = lastverticalw;
        self.m = m;
        self.mp = mp;
        self.ogcp1 = ogcp1;
        self.ogcp1o = ogcp1o;
        self.ogcp2 = ogcp2;
        self.ogcp2o = ogcp2o;
        self.fgcp1 = fgcp1;
        self.fgcp1o = fgcp1o;
        self.fgcp2 = fgcp2;
        self.fgcp2o = fgcp2o;
        self.cpmx1 = cpmx1;
        self.cpmx2 = cpmx2;
        self.gapfreq1 = gapfreq1;
        self.gapfreq2 = gapfreq2;
        self.doublework = doublework;
        self.intwork = intwork;
        self.ijp = ijp;

        (wm, cpmxresult)
    }

    /// `G__align11()`: the 1-vs-1 global DP.
    #[allow(clippy::too_many_arguments)]
    pub fn g_align11(
        &mut self,
        cn: &Constants,
        n_dynamicmtx: &[Vec<f64>],
        penalty: i32,
        penalty_ex: i32,
        seq1: &mut Vec<u8>,
        seq2: &mut Vec<u8>,
        headgp: bool,
        tailgp: bool,
        newgapstr: u8,
    ) -> f64 {
        let lgth1 = seq1.len();
        let lgth2 = seq2.len();

        if lgth1 == 0 && lgth2 == 0 {
            return 0.0;
        }
        if lgth1 == 0 {
            *seq1 = vec![newgapstr; lgth2];
            return 0.0;
        }
        if lgth2 == 0 {
            *seq2 = vec![newgapstr; lgth1];
            return 0.0;
        }

        let fpenalty = penalty as f64;
        let fpenalty_ex = penalty_ex as f64;
        let fpenalty_shift = cn.penalty_shift as f64;
        let mut warp = if cn.trywarp {
            if !headgp || !tailgp {
                eprintln!("At present, headgp and tailgp must be 1.");
                std::process::exit(1);
            }
            Some(Warp::new(lgth1, lgth2))
        } else {
            None
        };

        if lgth1 > self.g_orlgth1 || lgth2 > self.g_orlgth2 {
            let ll1 = max_usize((1.3 * lgth1 as f64) as usize, self.g_orlgth1) + 100;
            let ll2 = max_usize((1.3 * lgth2 as f64) as usize, self.g_orlgth2) + 100;
            self.g_w1 = vec![0.0; ll2 + 2];
            self.g_w2 = vec![0.0; ll2 + 2];
            self.g_initverticalw = vec![0.0; ll1 + 2];
            self.g_lastverticalw = vec![0.0; ll1 + 2];
            self.g_m = vec![0.0; ll2 + 2];
            self.g_mp = vec![0; ll2 + 2];
            self.g_amino_dynamicmtx = vec![vec![0.0; 0x100]; 0x100];
            self.g_orlgth1 = ll1 - 100;
            self.g_orlgth2 = ll2 - 100;
        }
        for i in 0..cn.nalphabets {
            for j in 0..cn.nalphabets {
                self.g_amino_dynamicmtx[cn.amino[i] as usize][cn.amino[j] as usize] =
                    n_dynamicmtx[i][j];
            }
        }
        self.ensure_ijp(self.g_orlgth1, self.g_orlgth2);

        let mtx = std::mem::take(&mut self.g_amino_dynamicmtx);
        let mut w1 = std::mem::take(&mut self.g_w1);
        let mut w2 = std::mem::take(&mut self.g_w2);
        let mut initverticalw = std::mem::take(&mut self.g_initverticalw);
        let mut lastverticalw = std::mem::take(&mut self.g_lastverticalw);
        let mut m = std::mem::take(&mut self.g_m);
        let mut mp = std::mem::take(&mut self.g_mp);
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

        if headgp {
            for i in 1..lgth1 + 1 {
                initverticalw[i] += fpenalty;
            }
            for j in 1..lgth2 + 1 {
                w1[j] += fpenalty;
            }
        } else {
            for i in 1..lgth1 + 1 {
                initverticalw[i] += fpenalty * TERMGAPFAC;
                initverticalw[i] += fpenalty_ex * i as f64 * TERMGAPFAC_EX;
            }
            for j in 1..lgth2 + 1 {
                w1[j] += fpenalty * TERMGAPFAC;
                w1[j] += fpenalty_ex * j as f64 * TERMGAPFAC_EX;
            }
        }

        for j in 1..lgth2 + 1 {
            m[j] = w1[j - 1];
            mp[j] = 0;
        }
        lastverticalw[0] = w1[lgth2 - 1];

        let lasti = if tailgp { lgth1 + 1 } else { lgth1 };
        let lastj = lgth2 + 1;
        let mut cur_is_w1 = true;
        let mut wm = 0.0f64;

        for i in 1..lasti {
            cur_is_w1 = !cur_is_w1;
            let (currentw, previousw) = if cur_is_w1 {
                (&mut w1, &mut w2)
            } else {
                (&mut w2, &mut w1)
            };
            previousw[0] = initverticalw[i - 1];
            {
                // With tailgp the last row reads seq1[lgth1], i.e. the C string
                // terminator; row 0 of the matrix is all zeros.
                let c1 = if i < lgth1 { seq1[i] } else { 0u8 };
                let row = &mtx[c1 as usize];
                for j in 0..lgth2 {
                    currentw[j] = row[seq2[j] as usize];
                }
            }
            currentw[0] = initverticalw[i];

            let mut mi = previousw[0];
            let mut mpi = 0i32;
            let fpenalty_ex_i = if i < lgth1 { fpenalty_ex } else { 0.0 };
            let ijprow = &mut ijp[i];

            for j in 1..lastj {
                let mut w = previousw[j - 1];
                ijprow[j] = 0;

                let g = mi + fpenalty;
                if g > w {
                    w = g;
                    ijprow[j] = -(j as i32 - mpi);
                }
                let g = previousw[j - 1];
                if g >= mi {
                    mi = g;
                    mpi = j as i32 - 1;
                }
                mi += fpenalty_ex_i;

                let g = m[j] + fpenalty;
                if g > w {
                    w = g;
                    ijprow[j] = i as i32 - mp[j];
                }
                let g = previousw[j - 1];
                if g >= m[j] {
                    m[j] = g;
                    mp[j] = i as i32 - 1;
                }
                if j < lgth2 {
                    m[j] += fpenalty_ex;
                }

                if let Some(wp) = warp.as_mut() {
                    let (nw, code) = wp.cell(i, j, w, currentw[j], fpenalty_shift, fpenalty_ex);
                    w = nw;
                    if let Some(c) = code {
                        ijprow[j] = c;
                    }
                }

                currentw[j] += w;
                wm = w;
            }
            if let Some(wp) = warp.as_mut() {
                wp.end_row(lastj);
            }
            lastverticalw[i] = currentw[lgth2 - 1];
        }

        let (wmo, mseq1, mseq2) = {
            let currentw: &[f64] = if cur_is_w1 { &w1 } else { &w2 };
            g_atracking(
                currentw,
                &lastverticalw,
                seq1,
                seq2,
                &mut ijp,
                tailgp,
                fpenalty,
                fpenalty_ex,
                newgapstr,
                warp.as_ref(),
            )
        };
        if !tailgp {
            wm = wmo;
        }

        *seq1 = mseq1;
        *seq2 = mseq2;

        self.g_amino_dynamicmtx = mtx;
        self.g_w1 = w1;
        self.g_w2 = w2;
        self.g_initverticalw = initverticalw;
        self.g_lastverticalw = lastverticalw;
        self.g_m = m;
        self.g_mp = mp;
        self.ijp = ijp;

        wm
    }
}

/// `Atracking()` from Galign11.c — emits the two aligned strings directly.
#[allow(clippy::too_many_arguments)]
fn g_atracking(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    seq1: &[u8],
    seq2: &[u8],
    ijp: &mut [Vec<i32>],
    tailgp: bool,
    fpenalty: f64,
    fpenalty_ex: f64,
    gap: u8,
    warp: Option<&Warp>,
) -> (f64, Vec<u8>, Vec<u8>) {
    let lgth1 = seq1.len();
    let lgth2 = seq2.len();
    let mut wm = 0.0f64;

    for (i, row) in ijp.iter_mut().enumerate().take(lgth1 + 1) {
        row[0] = i as i32 + 1;
    }
    for j in 0..lgth2 + 1 {
        ijp[0][j] = -(j as i32 + 1);
    }

    if !tailgp {
        wm = lasthorizontalw[lgth2 - 1] - 1.0;
        for j in (0..lgth2 - 1).rev() {
            let g = lasthorizontalw[j]
                + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth2 - 1 - j) as f64 * TERMGAPFAC_EX);
            if g > wm {
                wm = g;
                ijp[lgth1][lgth2] = -((lgth2 - j) as i32);
            }
        }
        for i in (0..lgth1 - 1).rev() {
            let g = lastverticalw[i]
                + (fpenalty * TERMGAPFAC + fpenalty_ex * (lgth1 - 1 - i) as f64 * TERMGAPFAC_EX);
            if g > wm {
                wm = g;
                ijp[lgth1][lgth2] = (lgth1 - i) as i32;
            }
        }
        if lasthorizontalw[lgth2 - 1] > wm {
            wm = lasthorizontalw[lgth2 - 1];
            ijp[lgth1][lgth2] = 0;
        }
    }

    let mut r1: Vec<u8> = Vec::new();
    let mut r2: Vec<u8> = Vec::new();
    let mut iin = lgth1 as i32;
    let mut jin = lgth2 as i32;
    let limk = (lgth1 + lgth2 + 1) as i32;
    let mut k = 0i32;
    let warpbase = warp.map(|w| w.warpbase).unwrap_or(i32::MAX);
    while k < limk {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi) = if v >= warpbase {
            warp.unwrap().target(v)
        } else if v < 0 {
            (iin - 1, jin + v)
        } else if v > 0 {
            (iin - v, jin - 1)
        } else {
            (iin - 1, jin - 1)
        };

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
        loop {
            l -= 1;
            if l <= 0 {
                break;
            }
            r1.push(seq1[(ifi + l) as usize]);
            r2.push(gap);
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l <= 0 {
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
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1;
    }
    r1.reverse();
    r2.reverse();
    (wm, r1, r2)
}

/// `Atracking_localhom()` from Salignmm.c.
///
/// Differs from [`atracking`] in two ways that change the answer: the end cell
/// is chosen with `>=` (so the *last* equal-scoring cell wins, not the first),
/// and the choice is governed by the global `outgap` — which `treebase` passes
/// as `tailgp`, so the parameter carries it here.
fn atracking_localhom(
    lasthorizontalw: &[f64],
    lastverticalw: &[f64],
    lgth1: usize,
    lgth2: usize,
    ijp: &mut [Vec<i32>],
    outgap: bool,
    imp: &[Vec<f64>],
    warp: Option<&Warp>,
) -> Track {
    if !outgap {
        let mut wm = lastverticalw[0];
        for i in 0..lgth1 {
            if lastverticalw[i] >= wm {
                wm = lastverticalw[i];
                ijp[lgth1][lgth2] = (lgth1 - i) as i32;
            }
        }
        for j in 0..lgth2 {
            if lasthorizontalw[j] >= wm {
                wm = lasthorizontalw[j];
                ijp[lgth1][lgth2] = -((lgth2 - j) as i32);
            }
        }
    }

    for (i, row) in ijp.iter_mut().enumerate().take(lgth1 + 1) {
        row[0] = i as i32 + 1;
    }
    for j in 0..lgth2 + 1 {
        ijp[0][j] = -(j as i32 + 1);
    }

    let mut rev1: Vec<u8> = Vec::with_capacity(lgth1 + lgth2 + 1);
    let mut rev2: Vec<u8> = Vec::with_capacity(lgth1 + lgth2 + 1);

    let mut iin = lgth1 as i32;
    let mut jin = lgth2 as i32;
    let limk = (lgth1 + lgth2 + 1) as i32;
    let mut k = 0i32;
    let mut impwm = 0.0f64;
    let warpbase = warp.map(|w| w.warpbase).unwrap_or(i32::MAX);
    while k < limk {
        let v = ijp[iin as usize][jin as usize];
        let (ifi, jfi) = if v >= warpbase {
            warp.unwrap().target(v)
        } else if v < 0 {
            (iin - 1, jin + v)
        } else if v > 0 {
            (iin - v, jin - 1)
        } else {
            (iin - 1, jin - 1)
        };

        #[allow(unused_assignments)]
        if ifi == -warpbase && jfi == -warpbase {
            let mut l = iin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'o');
                rev2.push(b'-');
                k += 1;
            }
            let mut l = jin;
            while {
                l -= 1;
                l >= 0
            } {
                rev1.push(b'-');
                rev2.push(b'o');
                k += 1;
            }
            break;
        }

        let mut l = iin - ifi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'o');
            rev2.push(b'-');
            k += 1;
        }
        let mut l = jin - jfi;
        loop {
            l -= 1;
            if l == 0 {
                break;
            }
            rev1.push(b'-');
            rev2.push(b'o');
            k += 1;
        }

        if !(iin as usize == lgth1 || jin as usize == lgth2) {
            impwm += imp[iin as usize][jin as usize];
        }

        if iin <= 0 || jin <= 0 {
            break;
        }
        rev1.push(b'o');
        rev2.push(b'o');
        k += 1;
        iin = ifi;
        jin = jfi;
        k += 1;
    }
    let _ = impwm; // the caller discards it (tbfast passes a dummy)

    rev1.reverse();
    rev2.reverse();
    let ngap1 = rev1.contains(&b'-');
    let ngap2 = rev2.contains(&b'-');
    Track {
        wm: 0.0,
        gt1: rev1,
        gt2: rev2,
        ngap1,
        ngap2,
    }
}
