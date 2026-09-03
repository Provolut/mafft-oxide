//! `D__align()` (`core/Dalignmm.c`): the generalized-affine DP, reached by
//! `-@` (`mafft --sp`).
//!
//! Where the other aligners carry two running gap scores, this one carries a
//! *distribution of gap lengths* per column and charges each candidate move
//! against it. `Gaplen` is one entry of that distribution.
//!
//! `Dalignmm.c` compiles with `USEGAPLENMTX 0` and `USEGAPLENHALF 0`, so the
//! `idatend` field and the `copygaplenrestricted*` helpers are not in the
//! build at all; what is live is the "compactx" family below. Those two zeros
//! also kill every `if( pfac != pfactmp )` cross-check in the DP, so
//! `PFACERROR` is always 0 -- which matters, because the `WMCHECK` block at
//! the end of `D__align` *is* live and tests it.
//!
//! ## The six cubes, for whoever writes the DP
//!
//! `D__align` threads six `Gaplen ***` through the DP, each `[pos][col]`, and
//! half of them are stored *relative to* `pos` rather than from column 0. The
//! allocation sizes are what say which is which, and the C's `// half`
//! comments agree:
//!
//! | cube                | rows      | columns per row | indexing        |
//! |---------------------|-----------|-----------------|-----------------|
//! | `gaplen1ibestkamo`  | `ll1+2`   | `ll1+2-i`       | half, offset i  |
//! | `gaplen2ibestkamo`  | `ll1+2`   | `ll2+2`         | full            |
//! | `gaplen1icurr`      | `ll1+2`   | `ll1+2-i`       | half, offset i  |
//! | `gaplen2icurr`      | `ll1+2`   | `ll2+2`         | full            |
//! | `gaplen1jbestkamo`  | `ll2+2`   | `ll1+2`         | full            |
//! | `gaplen2jbestkamo`  | `ll2+2`   | `ll2+2-j`       | half, offset j  |
//! | `gaplen1jbest`      | `ll2+2`   | `ll1+2`         | full            |
//! | `gaplen2jbest`      | `ll2+2`   | `ll2+2-j`       | half, offset j  |
//! | `gaplen1jcurr`      | `ll2+2`   | `ll1+2`         | full            |
//! | `gaplen2jcurr`      | `ll2+2`   | `ll2+2-j`       | half, offset j  |
//! | `gaplen1jprev`      | `ll2+2`   | `ll1+2`         | full            |
//! | `gaplen2jprev`      | `ll2+2`   | `ll2+2-j`       | half, offset j  |
//!
//! Every row ends with a `Col::Barrier` one past its last real column, which
//! is what stops `duplicategaplencompactx` and `freegaplenpartly` at the row
//! end without either of them being told the length.

/// One gap-length pattern at a column: `len` residues of gap, seen with
/// weight `freq`, and this column sits `relend` positions before its end.
///
/// A column's list is terminated by an entry with `idatnext == -1`, and entry
/// 0 carries the list's length in `npat`. `idatnext`/`idatprev` chain the same
/// pattern across the columns it spans, which is how a length is updated in
/// one place and read in another.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Gaplen {
    pub relend: i32,
    pub idatnext: i32,
    pub idatprev: i32,
    pub npat: i32,
    pub len: i32,
    pub freq: f64,
}

impl Gaplen {
    /// The C's terminator: `len` and `relend` -1, `freq` 0, `npat` -1.
    pub fn terminator() -> Gaplen {
        Gaplen { relend: -1, idatnext: -1, idatprev: -1, npat: -1, len: -1, freq: 0.0 }
    }
}

/// One column of a `Gaplen **`.
///
/// The C stores three things in that pointer: NULL, an owned array, and the
/// sentinel `(Gaplen *)1`. The sentinel is a *barrier* -- `duplicategaplencompactx`
/// and `freegaplenpartly` stop when they reach it -- not a shared pointer, so
/// nothing here is aliased and an owned `Vec` per column is faithful.
#[derive(Clone, Debug, Default, PartialEq)]
pub enum Col {
    #[default]
    Null,
    /// `(Gaplen *)1`
    Barrier,
    Data(Vec<Gaplen>),
}

impl Col {
    pub fn get(&self) -> Option<&Vec<Gaplen>> {
        match self {
            Col::Data(v) => Some(v),
            _ => None,
        }
    }
    pub fn get_mut(&mut self) -> Option<&mut Vec<Gaplen>> {
        match self {
            Col::Data(v) => Some(v),
            _ => None,
        }
    }
    pub fn is_null(&self) -> bool {
        matches!(self, Col::Null)
    }
}

/// `strralpha()`: how many non-gap characters run backwards from `pos`.
fn strralpha(seq: &[u8], pos: usize) -> i32 {
    let mut n = 0i32;
    let mut i = pos;
    while i > 0 {
        i -= 1;
        if seq[i] != b'-' {
            break;
        }
        n += 1;
    }
    n
}

/// `fillgaplen()`: propagate each column's patterns backwards over the columns
/// they span, chaining them with `idatnext`/`idatprev`.
pub fn fillgaplen(mtx: &mut [Col], l: usize) {
    for i in 0..=l {
        if mtx[i].is_null() {
            continue;
        }
        let entries: Vec<(usize, i32, f64)> = match mtx[i].get() {
            Some(v) => v
                .iter()
                .enumerate()
                .take_while(|(_, g)| g.len != -1)
                .map(|(n, g)| (n, g.len, g.freq))
                .collect(),
            None => continue,
        };
        for (n, len, freq) in entries {
            let mut idatnext = n as i32;
            for j in 0..len {
                let pos = i - 1 - j as usize;
                let k = match mtx[pos].get() {
                    None => {
                        mtx[pos] = Col::Data(vec![Gaplen::terminator()]);
                        0usize
                    }
                    Some(v) => v[0].npat as usize,
                };
                {
                    let v = mtx[pos].get_mut().unwrap();
                    if v.len() < k + 2 {
                        v.resize(k + 2, Gaplen::terminator());
                    }
                    v[k].len = len;
                    v[k].relend = j + 1;
                    v[k].freq = freq;
                    v[k].idatnext = idatnext;
                    v[k].idatprev = -1;
                    v[k].npat = -1;
                    v[k + 1] = Gaplen::terminator();
                    v[0].npat = k as i32 + 1;
                }
                // `mtx[pos+1][idatnext].idatprev = k` -- "kanarazu aru?"
                // ("is it always there?"): the C does not check.
                if let Some(v) = mtx[pos + 1].get_mut() {
                    v[idatnext as usize].idatprev = k as i32;
                }
                idatnext = k as i32;
            }
        }
    }
}

/// `gaplencount()`: the gap-length distribution of one group, per column.
pub fn gaplencount(n: usize, l: usize, mtx: &mut [Col], seq: &[&[u8]], eff: &[f64]) {
    for j in 0..=l {
        if matches!(mtx[j], Col::Barrier) {
            break;
        }
        mtx[j] = Col::Null;
    }
    let mut known = vec![0i32; l + 1];
    for j in 0..=l {
        for v in known.iter_mut().take(j) {
            *v = 0;
        }
        let mut nknown = 0usize;
        for i in 0..n {
            // `j` runs to `l` inclusive, so the last read is the C's NUL
            // terminator -- not a gap, so the column past the end still
            // records the gap that ends there.
            if seq[i].get(j).copied().unwrap_or(0) == b'-' {
                continue;
            }
            let gl = strralpha(seq[i], j);
            if gl > 0 {
                if known[gl as usize] != 0 {
                    let v = mtx[j].get_mut().expect("known implies a list");
                    let mut k = 0usize;
                    while v[k].len != -1 && v[k].len != gl {
                        k += 1;
                    }
                    if v[k].len == -1 {
                        eprintln!("Unexpected error!");
                        std::process::exit(1);
                    }
                    v[k].freq += eff[i];
                } else {
                    if mtx[j].is_null() {
                        mtx[j] = Col::Data(Vec::new());
                    }
                    let v = mtx[j].get_mut().unwrap();
                    v.resize(nknown + 2, Gaplen::terminator());
                    v[nknown].len = gl;
                    v[nknown].relend = 0;
                    v[nknown].freq = eff[i];
                    v[nknown].idatnext = -2;
                    v[nknown + 1] = Gaplen::terminator();
                    known[gl as usize] += 1;
                    nknown += 1;
                    v[0].npat = nknown as i32;
                }
            }
        }
    }
    fillgaplen(mtx, l);
}

/// `duplicategaplencompactx()`: replace `cpy[start..=end]` with a copy of
/// `orig`, stopping early at a barrier.
pub fn duplicategaplencompactx(cpy: &mut [Col], orig: &[Col], start: i32, end: i32) {
    let start = start.max(0) as usize;
    let mut end = end;
    for i in start..=(end.max(start as i32) as usize) {
        if i as i32 > end {
            break;
        }
        if matches!(cpy[i], Col::Barrier) {
            end = i as i32 - 1;
            break;
        }
        cpy[i] = Col::Null;
    }
    if end < start as i32 {
        return;
    }
    for i in start..=end as usize {
        match orig[i].get() {
            None => cpy[i] = Col::Null,
            Some(v) => {
                let mut out: Vec<Gaplen> = Vec::new();
                let mut l = 0usize;
                while v[l].idatnext != -1 {
                    out.push(v[l]);
                    l += 1;
                }
                out.push(v[l]);
                cpy[i] = Col::Data(out);
            }
        }
    }
}

/// `gaplenextendnoidatend()`: lengthen every pattern that starts at
/// `gapstartpos` by `insertionlen`, writing the new length at this column and
/// at the previous column's linked entry.
pub fn gaplenextendnoidatend(cpy: &mut [Col], gapstartpos: usize, insertionlen: i32) {
    let heads: Vec<(usize, i32, i32, i32)> = match cpy[gapstartpos].get() {
        None => return,
        Some(v) => {
            let mut out = Vec::new();
            let mut l = 0usize;
            while v[l].idatnext != -1 {
                out.push((l, v[l].len + insertionlen, v[l].relend, v[l].idatprev));
                l += 1;
            }
            out
        }
    };
    for (l, len, relend, idatprev) in heads {
        // The C's inner loop `break`s unconditionally on its first iteration,
        // so it only ever touches `cpy[gapstartpos][l]`, and only when that
        // entry's `relend` is non-zero.
        if relend != 0 {
            if let Some(v) = cpy[gapstartpos].get_mut() {
                v[l].len = len;
            }
        }
        if gapstartpos != 0 && idatprev != -1 {
            if let Some(v) = cpy[gapstartpos - 1].get_mut() {
                v[idatprev as usize].len = len;
            }
        }
    }
}

/// `copygaplencompactx()`: copy only the `len` fields of one column, then
/// extend. The C walks the *source* list and writes into the destination
/// without checking its length, so the destination must already be at least
/// as long -- every call site prepares it with `duplicategaplencompactx`.
#[allow(clippy::too_many_arguments)]
pub fn copygaplencompactx(
    cpy: &mut [Col],
    orig: &[Col],
    gapstartpos: i32,
    insertionlen: i32,
    posincopy: usize,
    posinori: usize,
) {
    let lens: Vec<i32> = match orig[posinori].get() {
        None => return,
        Some(v) => {
            let mut out = Vec::new();
            let mut l = 0usize;
            while v[l].idatnext != -1 {
                out.push(v[l].len);
                l += 1;
            }
            out
        }
    };
    if let Some(v) = cpy[posincopy].get_mut() {
        for (k, len) in lens.iter().enumerate() {
            if k < v.len() {
                v[k].len = *len;
            }
        }
    }
    let gapstartpos = if gapstartpos == -1 {
        posincopy
    } else {
        gapstartpos as usize
    };
    gaplenextendnoidatend(cpy, gapstartpos, insertionlen);
}

/// `extendgaplencompactx()`: make `cpy[start]` exist if it does not, then
/// carry the previous column's lengths into it through `idatprev`.
pub fn extendgaplencompactx(cpy: &mut [Col], orig: &[Col], start: usize) {
    let ov = match orig[start].get() {
        None => {
            cpy[start] = Col::Null;
            return;
        }
        Some(v) => v.clone(),
    };
    if cpy[start].is_null() {
        let mut out: Vec<Gaplen> = Vec::new();
        let mut l = 0usize;
        while ov[l].idatnext != -1 {
            out.push(ov[l]);
            l += 1;
        }
        out.push(ov[l]);
        cpy[start] = Col::Data(out);
    }
    if start == 0 || cpy[start - 1].is_null() {
        return;
    }
    let prev_lens: Vec<i32> = cpy[start - 1].get().unwrap().iter().map(|g| g.len).collect();
    let v = cpy[start].get_mut().unwrap();
    let mut l = 0usize;
    while ov[l].idatnext != -1 {
        let id = ov[l].idatprev;
        if id != -1 && (id as usize) < prev_lens.len() && l < v.len() {
            v[l].len = prev_lens[id as usize];
        }
        l += 1;
    }
}

/// `calcpfacnoidatend()`: the penalty factor for extending both sides
/// straight through -- the fraction of the two groups for which this column
/// pair is *not* inside a gap that already covers it.
pub fn calcpfacnoidatend(gaplen1: &Col, gaplen2: &Col) -> f64 {
    let mut pfac1 = 0.0f64;
    let mut pfac2 = 0.0f64;
    if let Some(v) = gaplen1.get() {
        for g in v.iter().take_while(|g| g.idatnext != -1) {
            if g.relend == 0 {
                pfac1 += g.freq;
            }
        }
    }
    if let Some(v) = gaplen2.get() {
        for g in v.iter().take_while(|g| g.idatnext != -1) {
            if g.relend == 0 {
                pfac2 += g.freq;
            }
        }
    }
    let mut pfac = pfac1 * pfac2 + pfac1 * (1.0 - pfac2) + pfac2 * (1.0 - pfac1);
    if let (Some(v1), Some(v2)) = (gaplen1.get(), gaplen2.get()) {
        for g1 in v1.iter().take_while(|g| g.idatnext != -1) {
            for g2 in v2.iter().take_while(|g| g.idatnext != -1) {
                let (p1, p2) = (g1.relend, g2.relend);
                let hit = if p1 == 0 && p2 == 0 {
                    g1.len == g2.len
                } else if p1 == 0 && p2 != 0 {
                    g2.len - (p2 - 1) > g1.len
                } else if p1 != 0 && p2 == 0 {
                    g1.len - (p1 - 1) > g2.len
                } else {
                    false
                };
                if hit {
                    pfac -= g1.freq * g2.freq;
                }
            }
        }
    }
    pfac
}

/// `calcpfac_gapex_noidatend()`: the factor for *extending* a gap in group 1
/// against a gap that ends at this column in group 2.
pub fn calcpfac_gapex_noidatend(gaplen1i: &Col, gaplen2j: &Col, newgaplen: i32) -> f64 {
    let mut pfac = 0.0f64;
    let v2 = match gaplen2j.get() {
        None => return pfac,
        Some(v) => v,
    };
    for g2 in v2.iter().take_while(|g| g.idatnext != -1) {
        let pos2 = g2.relend;
        if pos2 != 0 {
            continue;
        }
        let mut pfac1 = 1.0f64;
        let mut pfac10 = 1.0f64;
        if let Some(v1) = gaplen1i.get() {
            for g1 in v1.iter().take_while(|g| g.idatnext != -1) {
                let pos1 = g1.relend;
                pfac10 -= g1.freq;
                if newgaplen + g1.len - pos1 > g2.len - pos2 {
                    pfac1 -= g1.freq;
                }
            }
        }
        pfac += pfac1 * g2.freq;
        // The C marks this block "????" -- kept as written.
        if newgaplen >= g2.len - (pos2 - 1) {
            pfac -= pfac10 * g2.freq;
        }
    }
    pfac
}

/// `calcpfac_gap_noidatend()`: the factor for *opening* a gap of
/// `newgaplen` in group 1 at this column pair.
pub fn calcpfac_gap_noidatend(gaplen1i: &Col, gaplen2j: &Col, newgaplen: i32) -> f64 {
    let mut pfac = 0.0f64;
    let mut pfac1 = 0.0f64;
    let mut pfac10 = 0.0f64;
    if let Some(v1) = gaplen1i.get() {
        for g1 in v1.iter().take_while(|g| g.idatnext != -1) {
            let pos1 = g1.relend;
            if pos1 != 0 {
                let mut pfac2 = 0.0f64;
                if let Some(v2) = gaplen2j.get() {
                    for g2 in v2.iter().take_while(|g| g.idatnext != -1) {
                        if g2.relend == 0 && g2.len > g1.len - pos1 + newgaplen {
                            pfac2 += g2.freq;
                        }
                    }
                }
                pfac += pfac2 * g1.freq;
                pfac1 += g1.freq;
            } else {
                let mut pfac2 = 1.0f64;
                if let Some(v2) = gaplen2j.get() {
                    for g2 in v2.iter().take_while(|g| g.idatnext != -1) {
                        let pos2 = g2.relend;
                        // a gap ending here, then one continuing through
                        if pos2 == 0 && g2.len == g1.len + newgaplen {
                            pfac2 -= g2.freq;
                        }
                        if pos2 != 0 && g2.len - (pos2 - 1) > g1.len + newgaplen {
                            pfac2 -= g2.freq;
                        }
                    }
                }
                pfac += pfac2 * g1.freq;
                pfac10 += g1.freq;
            }
        }
    }
    pfac1 = 1.0 - pfac1 - pfac10;
    let mut pfac2 = 1.0f64;
    if let Some(v2) = gaplen2j.get() {
        for g2 in v2.iter().take_while(|g| g.idatnext != -1) {
            let pos2 = g2.relend;
            if pos2 == 0 && g2.len == newgaplen {
                pfac2 -= g2.freq;
            }
            if pos2 != 0 && g2.len - (pos2 - 1) > newgaplen {
                pfac2 -= g2.freq;
            }
        }
    }
    pfac += pfac1 * pfac2;
    pfac
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build(rows: &[&str], eff: &[f64]) -> Vec<Col> {
        let owned: Vec<Vec<u8>> = rows.iter().map(|s| s.as_bytes().to_vec()).collect();
        let refs: Vec<&[u8]> = owned.iter().map(|s| s.as_slice()).collect();
        let l = owned[0].len();
        let mut g = vec![Col::Null; l + 2];
        gaplencount(owned.len(), l, &mut g, &refs, eff);
        g
    }

    /// Every expected value here was taken from the C, through
    /// `gen_glharness.c`, as an f64 bit pattern. `verify_gaplen.sh` checks the
    /// same code against a built `core/` over 203 fuzzed pairs; these pin a
    /// few of those answers so `cargo test` alone still covers the layer.
    #[test]
    fn gaplencount_spreads_one_gap_back_over_the_columns_it_spans() {
        // "MNGTEG---DNFY" at weight 0.25: a 3-residue gap whose run ends at
        // column 9, recorded at columns 6..9 with relend counting down.
        let g = build(&["MNGTEG---DNFY", "MNGTEGAAADNFY"], &[0.25, 0.75]);
        for c in g.iter().take(6) {
            assert!(c.is_null());
        }
        let expect = [(6, 3, 3, -1), (7, 3, 2, 0), (8, 3, 1, 0), (9, 3, 0, 0)];
        for (col, len, relend, prev) in expect {
            let v = g[col].get().unwrap_or_else(|| panic!("column {} empty", col));
            assert_eq!(v[0].len, len, "col {} len", col);
            assert_eq!(v[0].relend, relend, "col {} relend", col);
            assert_eq!(v[0].freq, 0.25, "col {} freq", col);
            assert_eq!(v[0].idatprev, prev, "col {} idatprev", col);
            assert_eq!(v[0].npat, 1, "col {} npat", col);
        }
        // The run's own column chains forward to nothing: idatnext -2 is the
        // marker `gaplencount` writes before `fillgaplen` links the rest.
        assert_eq!(g[9].get().unwrap()[0].idatnext, -2);
        assert_eq!(g[6].get().unwrap()[0].idatnext, 0);
    }

    #[test]
    fn a_gap_free_group_records_nothing() {
        let g = build(&["MNGTEGDNFY"], &[1.0]);
        assert!(g.iter().all(|c| c.is_null()));
    }

    #[test]
    fn two_different_gap_lengths_at_one_column_are_two_patterns() {
        // Both rows end a gap at column 4, of different lengths, so column 4
        // carries two entries and `npat` is 2.
        let g = build(&["MN--TEG", "M---TEG"], &[0.5, 0.5]);
        let v = g[4].get().expect("column 4 has both runs");
        assert_eq!(v[0].npat, 2);
        let mut lens: Vec<i32> = v.iter().take_while(|e| e.idatnext != -1).map(|e| e.len).collect();
        lens.sort();
        assert_eq!(lens, vec![2, 3]);
    }

    #[test]
    fn calcpfac_matches_the_c_on_a_pinned_pair() {
        // Group 1 "MNGTEG---DNFY"/"MNGTEGAAADNFY" at 0.25/0.75, group 2
        // "MN--TEGDNF" at 1.0. Bit patterns straight from gen_glharness.c.
        let g1 = build(&["MNGTEG---DNFY", "MNGTEGAAADNFY"], &[0.25, 0.75]);
        let g2 = build(&["MN--TEGDNF"], &[1.0]);
        // Column 4 of group 2 ends its 2-residue gap; against a group-1
        // column with no pattern the straight factor is 1 and opening is 0.
        assert_eq!(calcpfacnoidatend(&g1[0], &g2[4]).to_bits(), 0x3ff0000000000000);
        assert_eq!(calcpfac_gap_noidatend(&g1[0], &g2[4], 2).to_bits(), 0);
        assert_eq!(calcpfac_gapex_noidatend(&g1[0], &g2[4], 2).to_bits(), 0x3ff0000000000000);
        // Away from that column the roles invert.
        assert_eq!(calcpfacnoidatend(&g1[0], &g2[0]).to_bits(), 0);
        assert_eq!(calcpfac_gap_noidatend(&g1[0], &g2[0], 2).to_bits(), 0x3ff0000000000000);
        // Group-1 column 6 opens its own run, which drops the opening factor
        // to 1 - 0.25 = 0.75.
        assert_eq!(calcpfac_gap_noidatend(&g1[6], &g2[0], 2).to_bits(), 0x3fe8000000000000);
    }

    #[test]
    fn duplicate_then_extend_leaves_the_source_untouched() {
        let src = build(&["MNGTEG---DNFY", "MNGTEGAAADNFY"], &[0.25, 0.75]);
        let mut dst = vec![Col::Null; src.len()];
        duplicategaplencompactx(&mut dst, &src, 0, src.len() as i32 - 2);
        assert_eq!(dst[6], src[6]);
        assert_eq!(dst[9], src[9]);
        let before = src[6].clone();
        gaplenextendnoidatend(&mut dst, 6, 5);
        assert_eq!(src[6], before, "the source must not move");
        // relend != 0 at column 6, so its own length grows.
        assert_eq!(dst[6].get().unwrap()[0].len, 8);
    }

    #[test]
    fn duplicate_stops_at_a_barrier() {
        let src = build(&["MNGTEG---DNFY", "MNGTEGAAADNFY"], &[0.25, 0.75]);
        let mut dst = vec![Col::Null; src.len()];
        dst[8] = Col::Barrier;
        duplicategaplencompactx(&mut dst, &src, 0, src.len() as i32 - 2);
        assert_eq!(dst[6], src[6]);
        assert_eq!(dst[7], src[7]);
        assert_eq!(dst[8], Col::Barrier, "the barrier survives");
        assert!(dst[9].is_null(), "nothing past the barrier is copied");
    }
}
