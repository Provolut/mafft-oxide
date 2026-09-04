//! The shift ("warp") extension of the pairwise and profile DPs, used by
//! `--allowshift` / `--unalignlevel > 0`.
//!
//! Transliterated from the `trywarp` blocks in `core/Salignmm.c`,
//! `core/Galign11.c` and `core/Lalign11.c`, which are the same code three
//! times.
//!
//! Ordinarily a DP cell may only come from its three neighbours. With warp
//! enabled a cell may instead jump to the best cell seen anywhere in the
//! previous row, paying `penalty_shift` plus a gap-extension cost proportional
//! to how far it jumped. That is what lets MAFFT align sequences whose
//! homologous blocks appear in a different order. The jump targets are
//! recorded in `warpis`/`warpjs` and referenced from the traceback matrix by
//! indices at or above `warpbase`.

/// Per-row bookkeeping for the warp option.
pub struct Warp {
    pub warpis: Vec<i32>,
    pub warpjs: Vec<i32>,
    warpi: Vec<i32>,
    warpj: Vec<i32>,
    prevwarpi: Vec<i32>,
    prevwarpj: Vec<i32>,
    wmrecords: Vec<f64>,
    prevwmrecords: Vec<f64>,
    pub warpbase: i32,
}

impl Warp {
    /// Allocate for one alignment of the given lengths.
    pub fn new(lgth1: usize, lgth2: usize) -> Warp {
        let warpbase = (lgth1 + lgth2) as i32;
        Warp {
            warpis: Vec::new(),
            warpjs: Vec::new(),
            warpi: vec![-warpbase; lgth2 + 1],
            warpj: vec![-warpbase; lgth2 + 1],
            prevwarpi: vec![-warpbase; lgth2 + 1],
            prevwarpj: vec![-warpbase; lgth2 + 1],
            wmrecords: vec![0.0; lgth2 + 1],
            prevwmrecords: vec![0.0; lgth2 + 1],
            warpbase,
        }
    }

    /// The per-cell part of the `trywarp` block.
    ///
    /// Returns the (possibly improved) score and, when the jump wins, the
    /// traceback code to store. `curw` is the match score of the current cell,
    /// i.e. `*curpt` before `wm` is added to it.
    pub fn cell(
        &mut self,
        i: usize,
        j: usize,
        wm: f64,
        curw: f64,
        fpenalty_shift: f64,
        fpenalty_ex: f64,
    ) -> (f64, Option<i32>) {
        let mut wm = wm;
        let mut code = None;
        let g = self.prevwmrecords[j - 1]
            + fpenalty_shift
            + fpenalty_ex
                * (i as f64 - self.prevwarpi[j - 1] as f64 + j as f64
                    - self.prevwarpj[j - 1] as f64);
        if g > wm {
            let n = self.warpis.len();
            if n > 0
                && self.prevwarpi[j - 1] == self.warpis[n - 1]
                && self.prevwarpj[j - 1] == self.warpjs[n - 1]
            {
                code = Some(self.warpbase + n as i32 - 1);
            } else {
                code = Some(self.warpbase + n as i32);
                self.warpis.push(self.prevwarpi[j - 1]);
                self.warpjs.push(self.prevwarpj[j - 1]);
            }
            wm = g;
        }

        let curm = curw + wm;
        if self.wmrecords[j - 1] > self.wmrecords[j] {
            self.wmrecords[j] = self.wmrecords[j - 1];
            self.warpi[j] = self.warpi[j - 1];
            self.warpj[j] = self.warpj[j - 1];
        }
        if curm > self.wmrecords[j] {
            self.wmrecords[j] = curm;
            self.warpi[j] = i as i32;
            self.warpj[j] = j as i32;
        }
        (wm, code)
    }

    /// End of row: `fltncpy` / `intncpy` of the records into the `prev` copies.
    pub fn end_row(&mut self, lastj: usize) {
        self.prevwmrecords[..lastj].copy_from_slice(&self.wmrecords[..lastj]);
        self.prevwarpi[..lastj].copy_from_slice(&self.warpi[..lastj]);
        self.prevwarpj[..lastj].copy_from_slice(&self.warpj[..lastj]);
    }

    /// Resolve a traceback code at or above `warpbase` to its jump target.
    pub fn target(&self, code: i32) -> (i32, i32) {
        let k = (code - self.warpbase) as usize;
        (self.warpis[k], self.warpjs[k])
    }
}
