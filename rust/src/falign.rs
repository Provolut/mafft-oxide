//! FFT-based anchoring, transliterated from `core/Falign.c` (`Falign`) and
//! `core/fftFunctions.c` (`getKouho`, `zurasu2`, `alignableReagion`,
//! `blockAlign2`).
//!
//! This is the step MAFFT is named for. Both groups are turned into a complex
//! signal whose real and imaginary parts are the weighted, normalised volume
//! and polarity of each column; the cross-correlation of the two signals is
//! computed by FFT, and its highest peaks are candidate relative shifts
//! ("lags") at which homologous regions line up. For each candidate lag a
//! sliding window scores the two groups column by column and keeps the runs
//! that beat a threshold; those become anchors. `blockAlign2` then picks a
//! consistent, monotone subset of anchors by dynamic programming, and the
//! groups are cut at the chosen anchors and aligned block by block with the
//! ordinary profile DP. See Katoh et al., *NAR* 30:3059-3066 (2002).
//!
//! Everything below keeps the C's evaluation order, its persistent scratch
//! buffers, and its quirks — including `blockAlign2` comparing against `maxj`
//! in the loop that fills `maxi`, and reading one element of `result1` left
//! over from a previous call.

use crate::align::AlignState;
use crate::constants::{Constants, NALPHABETS};
use crate::fft::{calc_naiseki, Cpx, Fft};

/// `NKOUHO` from `fft.h`: how many correlation peaks to try.
const NKOUHO: usize = 20;
/// `MAXSEG` from `mltaln.h`.
const MAXSEG: usize = 100000;
/// `SEGMENTSIZE` from `fftFunctions.c`.
const SEGMENTSIZE: i32 = 150;

/// `Segment` from `mltaln.h`, minus the two `skip*` fields the C writes but
/// never reads.
#[derive(Clone, Copy, Debug, Default)]
struct Segment {
    start: i32,
    end: i32,
    center: i32,
    score: f64,
    number: i32,
}

/// `getKouho()`: the `nkouho` highest correlation peaks, as lags. Destroys
/// `soukan`, as the C does.
fn get_kouho(kouho: &mut [i32], nkouho: usize, soukan: &mut [f64], nlen2: usize) {
    let nlen4 = nlen2 / 2;
    for k in kouho.iter_mut().take(nkouho) {
        let mut max = -9999.9f64;
        let mut ikouho = 0usize;
        for i in 0..nlen2 {
            let tmp = soukan[i];
            if tmp > max {
                ikouho = i;
                max = tmp;
            }
        }
        soukan[ikouho] = -9999.9;
        *k = ikouho as i32 - nlen4 as i32;
    }
}

/// The `static` state of `alignableReagion()`.
struct AlignableState {
    stra: Vec<f64>,
    alloclen: usize,
    threshold: f64,
    threshold_set: bool,
}

impl AlignableState {
    fn new() -> AlignableState {
        AlignableState {
            stra: Vec::new(),
            alloclen: 0,
            threshold: 0.0,
            threshold_set: false,
        }
    }
}

/// `alignableReagion()`: sliding-window scan for runs of high profile-profile
/// similarity at a fixed lag. Appends to `seg` and returns how many it added.
#[allow(clippy::too_many_arguments)]
fn alignable_reagion(
    st: &mut AlignableState,
    cn: &Constants,
    clus1: usize,
    clus2: usize,
    seq1: &[&[u8]],
    seq2: &[&[u8]],
    eff1: &[f64],
    eff2: &[f64],
    seg: &mut [Segment],
) -> usize {
    let fft_win_size = cn.fft_win_size;
    let len = seq1[0].len().min(seq2[0].len());
    let maxlen = seq1[0].len().max(seq2[0].len()) + fft_win_size as usize;
    if st.alloclen < maxlen {
        if !st.threshold_set {
            st.threshold = cn.fft_threshold as f64 / 100.0 * 600.0 * fft_win_size as f64;
            st.threshold_set = true;
        }
        // The C frees and callocs, so the fresh buffer is zeroed and any values
        // beyond `len` are zero until a later call leaves something there.
        st.stra = vec![0.0; maxlen];
        st.alloclen = maxlen;
    }
    let threshold = st.threshold;

    let mut totaleff = 0.0f64;
    for e1 in eff1.iter().take(clus1) {
        for e2 in eff2.iter().take(clus2) {
            totaleff += *e1 * *e2;
        }
    }

    let mut hat1 = [0i32; NALPHABETS + 1];
    let mut hat2 = [0i32; NALPHABETS + 1];

    for i in 0..len {
        let mut prf1 = [0.0f64; NALPHABETS];
        let mut prf2 = [0.0f64; NALPHABETS];
        for j in 0..clus1 {
            prf1[cn.amino_n[seq1[j][i] as usize] as usize] += eff1[j];
        }
        for j in 0..clus2 {
            prf2[cn.amino_n[seq2[j][i] as usize] as usize] += eff2[j];
        }

        // Linked lists over the occupied letters, highest index first.
        let mut pre1 = NALPHABETS;
        let mut pre2 = NALPHABETS;
        for j in (0..=25usize).rev() {
            if prf1[j] != 0.0 {
                hat1[pre1] = j as i32;
                pre1 = j;
            }
            if prf2[j] != 0.0 {
                hat2[pre2] = j as i32;
                pre2 = j;
            }
        }
        hat1[pre1] = -1;
        hat2[pre2] = -1;

        let mut s = 0.0f64;
        let mut k = hat1[NALPHABETS];
        while k != -1 {
            let mut j = hat2[NALPHABETS];
            while j != -1 {
                s += cn.n_dis_fft[k as usize][j as usize] as f64
                    * prf1[k as usize]
                    * prf2[j as usize];
                j = hat2[j as usize];
            }
            k = hat1[k as usize];
        }
        st.stra[i] = s / totaleff;
    }

    let mut value = 0usize;
    let mut status = false;
    let mut starttmp = 0i32;
    let mut length = 0i32;
    let mut cumscore = 0.0f64;
    let mut score = 0.0f64;
    for j in 0..fft_win_size as usize {
        score += st.stra[j];
    }

    let lim = len as isize - fft_win_size as isize;
    let mut i: isize = 1;
    while i < lim {
        let iu = i as usize;
        score = score - st.stra[iu - 1] + st.stra[iu + fft_win_size as usize - 1];

        if score > threshold {
            if !status {
                status = true;
                starttmp = i as i32;
                length = 0;
                cumscore = 0.0;
            }
            length += 1;
            cumscore += score;
        }
        if score <= threshold || length > SEGMENTSIZE {
            if status {
                if length > fft_win_size {
                    seg[value].start = starttmp;
                    seg[value].end = i as i32;
                    seg[value].center = (starttmp + i as i32 + fft_win_size) / 2;
                    seg[value].score = cumscore;
                    value += 1;
                }
                length = 0;
                cumscore = 0.0;
                status = false;
                starttmp = i as i32;
                if value > MAXSEG - 3 {
                    eprintln!("TOO MANY SEGMENTS!");
                    std::process::exit(1);
                }
            }
        }
        i += 1;
    }
    // C leaves `i` at max(1, lim) after the loop.
    let iend = if lim > 1 { lim as i32 } else { 1 };
    if status && length > fft_win_size {
        seg[value].end = iend;
        seg[value].start = starttmp;
        seg[value].center = (starttmp + iend + fft_win_size) / 2;
        seg[value].score = cumscore;
        value += 1;
    }
    value
}

/// `mymergesort()`: stable merge sort of segment indices by `center`.
fn mymergesort(first: isize, last: isize, order: &mut [usize], center: &dyn Fn(usize) -> i32) {
    if first < last {
        let middle = (first + last) / 2;
        mymergesort(first, middle, order, center);
        mymergesort(middle + 1, last, order, center);
        let work: Vec<usize> = order[first as usize..=middle as usize].to_vec();
        let mut i = (middle + 1) as usize;
        let mut j = 0usize;
        let mut k = first as usize;
        while i <= last as usize && j < work.len() {
            if center(work[j]) <= center(order[i]) {
                order[k] = work[j];
                j += 1;
            } else {
                order[k] = order[i];
                i += 1;
            }
            k += 1;
        }
        while j < work.len() {
            order[k] = work[j];
            j += 1;
            k += 1;
        }
    }
}

/// The `static` state of `blockAlign2()`.
///
/// `result1` / `result2` persist between calls and one element of each is read
/// back before being written; `maxi` / `maxj` / `pointi` / `pointj` persist too
/// and `maxj` is read in the loop that computes `maxi`. Both are faithful.
pub struct BlockAlignState {
    result1: Vec<i32>,
    result2: Vec<i32>,
    maxi: f64,
    maxj: f64,
    pointi: usize,
    pointj: usize,
}

impl BlockAlignState {
    fn new() -> BlockAlignState {
        BlockAlignState {
            result1: vec![0; MAXSEG],
            result2: vec![0; MAXSEG],
            maxi: 0.0,
            maxj: 0.0,
            pointi: 0,
            pointj: 0,
        }
    }
}

/// `blockAlign2()`: choose a consistent, monotone subset of the anchors by DP,
/// and rewrite `cut1` / `cut2` / `ncut` with the chosen cut positions.
fn block_align2(
    st: &mut BlockAlignState,
    penalty: i32,
    cut1: &mut [i32],
    cut2: &mut [i32],
    ocrossscore: &[Vec<f64>],
    ncut: &mut usize,
) {
    let n = *ncut;
    let mut crossscore = vec![vec![0.0f64; n + 2]; n + 2];
    let mut track = vec![vec![0i32; n + 2]; n + 2];
    for i in 0..n {
        crossscore[i][..n].clone_from_slice(&ocrossscore[i][..n]);
    }
    let ocut1: Vec<i32> = cut1[..n].to_vec();
    let ocut2: Vec<i32> = cut2[..n].to_vec();

    for i in 1..n {
        for j in 1..n {
            // `permit()` returns 0 unconditionally, so the guard reduces to
            // "only k == 0 unless we are on the last row/column".
            st.pointi = 0;
            st.maxi = 0.0;
            let klim = j as isize - 2;
            let mut k: isize = 0;
            while k < klim {
                let ku = k as usize;
                if !(k != 0 && ku < n - 1 && j < n - 1) {
                    // The C compares against maxj here, not maxi.
                    if crossscore[i - 1][ku] > st.maxj {
                        st.pointi = ku;
                        st.maxi = crossscore[i - 1][ku];
                    }
                }
                k += 1;
            }

            st.pointj = 0;
            st.maxj = 0.0;
            let klim = i as isize - 2;
            let mut k: isize = 0;
            while k < klim {
                let ku = k as usize;
                if !(k != 0 && ku < n - 1 && i < n - 1) {
                    if crossscore[ku][j - 1] > st.maxj {
                        st.pointj = ku;
                        st.maxj = crossscore[ku][j - 1];
                    }
                }
                k += 1;
            }

            // `maxi += penalty` / `maxj += penalty` in the C mutate the
            // statics, and the penalised `maxj` is what the *next* iteration's
            // first k-loop compares against.
            st.maxi += penalty as f64;
            st.maxj += penalty as f64;

            let mut maximum = crossscore[i - 1][j - 1];
            track[i][j] = 0;
            if maximum < st.maxi {
                maximum = st.maxi;
                track[i][j] = j as i32 - st.pointi as i32;
            }
            if maximum < st.maxj {
                maximum = st.maxj;
                track[i][j] = st.pointj as i32 - i as i32;
            }
            crossscore[i][j] += maximum;
        }
    }

    st.result1[MAXSEG - 1] = n as i32 - 1;
    st.result2[MAXSEG - 1] = n as i32 - 1;

    let mut i = MAXSEG - 1;
    while i >= 1 {
        let cur1 = st.result1[i];
        let cur2 = st.result2[i];
        if cur1 == 0 || cur2 == 0 {
            break;
        }
        let shift = track[cur1 as usize][cur2 as usize];
        if shift == 0 {
            st.result1[i - 1] = cur1 - 1;
            st.result2[i - 1] = cur2 - 1;
        } else if shift > 0 {
            st.result1[i - 1] = cur1 - 1;
            st.result2[i - 1] = cur2 - shift;
        } else {
            st.result1[i - 1] = cur1 + shift;
            st.result2[i - 1] = cur2 - 1;
        }
        i -= 1;
    }

    let mut count: i32 = 0;
    for j in i..MAXSEG {
        let r1 = st.result1[j] as usize;
        let r2 = st.result2[j] as usize;
        if ocrossscore[r1][r2] == 0.0 {
            continue;
        }
        // Reads result1[j-1]; at j == i that element was last written by a
        // previous call (or is the calloc zero on the first).
        let p1 = st.result1[j - 1] as usize;
        let p2 = st.result2[j - 1] as usize;
        if (st.result1[j] == st.result1[j - 1] || st.result2[j] == st.result2[j - 1])
            && ocrossscore[r1][r2] > ocrossscore[p1][p2]
        {
            count -= 1;
        }
        cut1[count as usize] = ocut1[r1];
        cut2[count as usize] = ocut2[r2];
        count += 1;
    }
    *ncut = count as usize;
}

/// All persistent state `Falign()` keeps in `static`s.
pub struct FalignState {
    fft: Fft,
    alignable: AlignableState,
    block: BlockAlignState,
    n_alloc: usize,
    seq_vector1: Vec<Cpx>,
    seq_vector2: Vec<Cpx>,
    naiseki: Vec<Cpx>,
    naiseki_no_wa: Vec<Cpx>,
    soukan: Vec<f64>,
    kouho: Vec<i32>,
    cut1: Vec<i32>,
    cut2: Vec<i32>,
    segment: Vec<Segment>,
    segment1: Vec<Segment>,
    segment2: Vec<Segment>,
}

impl Default for FalignState {
    fn default() -> Self {
        Self::new()
    }
}

impl FalignState {
    pub fn new() -> FalignState {
        FalignState {
            fft: Fft::new(),
            alignable: AlignableState::new(),
            block: BlockAlignState::new(),
            n_alloc: 0,
            seq_vector1: Vec::new(),
            seq_vector2: Vec::new(),
            naiseki: Vec::new(),
            naiseki_no_wa: Vec::new(),
            soukan: Vec::new(),
            kouho: vec![0; NKOUHO],
            cut1: vec![0; MAXSEG],
            cut2: vec![0; MAXSEG],
            segment: vec![Segment::default(); MAXSEG],
            segment1: vec![Segment::default(); MAXSEG],
            segment2: vec![Segment::default(); MAXSEG],
        }
    }
}

/// `getkyokaigap()`: the column of `seq` at `pos`, one character per member.
fn getkyokaigap(seq: &[Vec<u8>], pos: usize) -> Vec<u8> {
    seq.iter().map(|s| s[pos]).collect()
}

/// `Falign()` for the `treebase()` call site: `alg == 'A'`, no constraints, no
/// codon options, not called from `tditeration`.
///
/// `group1` / `group2` are the aligned members of each group and are replaced
/// by the new alignment. Returns the summed block scores.
#[allow(clippy::too_many_arguments)]
pub fn falign(
    fs: &mut FalignState,
    ast: &mut AlignState,
    cn: &Constants,
    n_dynamicmtx: &[Vec<f64>],
    group1: &mut Vec<Vec<u8>>,
    group2: &mut Vec<Vec<u8>>,
    eff1: &[f64],
    eff2: &[f64],
    outgap: bool,
    fftlog: &mut i32,
    newgapstr: u8,
) -> f64 {
    let clus1 = group1.len();
    let clus2 = group2.len();
    let len1 = group1[0].len();
    let len2 = group2[0].len();

    // nlen = smallest power of two strictly greater than max(len1, len2)
    let nlentmp = len1.max(len2);
    let mut nlen = 1usize;
    while nlentmp >= nlen {
        nlen <<= 1;
    }
    let nlen2 = nlen / 2;

    if fs.n_alloc < nlen {
        fs.naiseki_no_wa = vec![Cpx::default(); nlen];
        fs.naiseki = vec![Cpx::default(); nlen];
        fs.seq_vector1 = vec![Cpx::default(); nlen + 1];
        fs.seq_vector2 = vec![Cpx::default(); nlen + 1];
        fs.soukan = vec![0.0; nlen + 1];
        fs.n_alloc = nlen;
    }

    // --- build the two signals ---
    // `n20or4or2` is 1 in this function for both sequence types, so there is a
    // single vector each side. For protein it carries polarity in the real part
    // and volume in the imaginary one; for nucleotides the branch taken is
    // `seq_vec_3`, which -- because `n20or4or2` is 1, not 4 -- fills only the
    // vector for alphabet index 0, i.e. the indicator of `a`.
    let nuc = cn.dorp == b'd';
    for v in fs.seq_vector1.iter_mut().take(nlen) {
        *v = Cpx::default();
    }
    for i in 0..clus1 {
        if nuc {
            seq_vec_3(&mut fs.seq_vector1, eff1[i], &group1[i], cn);
        } else {
            seq_vec_5(
                &mut fs.seq_vector1,
                &cn.polarity,
                &cn.volume,
                eff1[i],
                &group1[i],
                cn,
            );
        }
    }
    for v in fs.seq_vector2.iter_mut().take(nlen) {
        *v = Cpx::default();
    }
    for i in 0..clus2 {
        if nuc {
            seq_vec_3(&mut fs.seq_vector2, eff2[i], &group2[i], cn);
        } else {
            seq_vec_5(
                &mut fs.seq_vector2,
                &cn.polarity,
                &cn.volume,
                eff2[i],
                &group2[i],
                cn,
            );
        }
    }

    fs.fft.run(nlen as isize, &mut fs.seq_vector2);
    fs.fft.run(nlen as isize, &mut fs.seq_vector1);

    for l in 0..nlen {
        fs.naiseki[l] = calc_naiseki(&fs.seq_vector1[l], &fs.seq_vector2[l]);
    }
    for l in 0..nlen {
        fs.naiseki_no_wa[l].r = fs.naiseki[l].r;
        fs.naiseki_no_wa[l].i = fs.naiseki[l].i;
    }

    fs.fft.run(-(nlen as isize), &mut fs.naiseki_no_wa);

    for m in 0..=nlen2 {
        fs.soukan[m] = fs.naiseki_no_wa[nlen2 - m].r;
    }
    for m in nlen2 + 1..nlen {
        fs.soukan[m] = fs.naiseki_no_wa[nlen + nlen2 - m].r;
    }

    get_kouho(&mut fs.kouho, NKOUHO, &mut fs.soukan, nlen);

    // --- collect anchors at each candidate lag ---
    let mut count = 0usize;
    for k in 0..NKOUHO {
        let lag = fs.kouho[k];
        if lag <= -(len1 as i32) || len2 as i32 <= lag {
            continue;
        }
        // zurasu2(): shift whichever group starts later.
        let (off1, off2) = if lag > 0 {
            (0usize, lag as usize)
        } else {
            ((-lag) as usize, 0usize)
        };
        let s1: Vec<&[u8]> = group1.iter().map(|s| &s[off1..]).collect();
        let s2: Vec<&[u8]> = group2.iter().map(|s| &s[off2..]).collect();

        let tmpint = alignable_reagion(
            &mut fs.alignable,
            cn,
            clus1,
            clus2,
            &s1,
            &s2,
            eff1,
            eff2,
            &mut fs.segment[count..],
        );

        if count + tmpint > MAXSEG - 3 {
            eprintln!("TOO MANY SEGMENTS.");
            std::process::exit(1);
        }
        if tmpint == 0 {
            break;
        }
        for _ in 0..tmpint {
            let seg = fs.segment[count];
            if lag > 0 {
                fs.segment1[count] = Segment {
                    start: seg.start,
                    end: seg.end,
                    center: seg.center,
                    score: seg.score,
                    number: 0,
                };
                fs.segment2[count] = Segment {
                    start: seg.start + lag,
                    end: seg.end + lag,
                    center: seg.center + lag,
                    score: seg.score,
                    number: 0,
                };
            } else {
                fs.segment1[count] = Segment {
                    start: seg.start - lag,
                    end: seg.end - lag,
                    center: seg.center - lag,
                    score: seg.score,
                    number: 0,
                };
                fs.segment2[count] = Segment {
                    start: seg.start,
                    end: seg.end,
                    center: seg.center,
                    score: seg.score,
                    number: 0,
                };
            }
            count += 1;
        }
    }

    // --- sort anchors along each sequence and index them ---
    let mut sorted1: Vec<usize> = (0..count).collect();
    let mut sorted2: Vec<usize> = (0..count).collect();
    {
        let c1: Vec<i32> = (0..count).map(|i| fs.segment1[i].center).collect();
        let c2: Vec<i32> = (0..count).map(|i| fs.segment2[i].center).collect();
        mymergesort(0, count as isize - 1, &mut sorted1, &|i| c1[i]);
        mymergesort(0, count as isize - 1, &mut sorted2, &|i| c2[i]);
    }
    for (i, &s) in sorted1.iter().enumerate() {
        fs.segment1[s].number = i as i32;
    }
    for (i, &s) in sorted2.iter().enumerate() {
        fs.segment2[s].number = i as i32;
    }

    let mut crossscore = vec![vec![0.0f64; count + 2]; count + 2];
    for i in 0..count {
        crossscore[fs.segment1[i].number as usize + 1][fs.segment2[i].number as usize + 1] =
            fs.segment1[i].score;
        fs.cut1[i + 1] = fs.segment1[sorted1[i]].center;
        fs.cut2[i + 1] = fs.segment2[sorted2[i]].center;
    }
    crossscore[0][0] = 10000000.0;
    fs.cut1[0] = 0;
    fs.cut2[0] = 0;
    crossscore[count + 1][count + 1] = 10000000.0;
    fs.cut1[count + 1] = len1 as i32;
    fs.cut2[count + 1] = len2 as i32;
    count += 2;

    block_align2(
        &mut fs.block,
        cn.penalty,
        &mut fs.cut1,
        &mut fs.cut2,
        &crossscore,
        &mut count,
    );

    // --- align each block between consecutive cuts ---
    let mut result1: Vec<Vec<u8>> = vec![Vec::new(); clus1];
    let mut result2: Vec<Vec<u8>> = vec![Vec::new(); clus2];
    let mut totalscore = 0.0f64;
    *fftlog = -1;
    for i in 0..count - 1 {
        *fftlog += 1;
        let headgp = if i == 0 { outgap } else { true };
        let tailgp = if i == count - 2 { outgap } else { true };

        let sgap1 = if fs.cut1[i] != 0 {
            getkyokaigap(group1, fs.cut1[i] as usize - 1)
        } else {
            vec![b'o'; clus1]
        };
        let sgap2 = if fs.cut2[i] != 0 {
            getkyokaigap(group2, fs.cut2[i] as usize - 1)
        } else {
            vec![b'o'; clus2]
        };
        let egap1 = if fs.cut1[i + 1] != len1 as i32 {
            getkyokaigap(group1, fs.cut1[i + 1] as usize)
        } else {
            vec![b'o'; clus1]
        };
        let egap2 = if fs.cut2[i + 1] != len2 as i32 {
            getkyokaigap(group2, fs.cut2[i + 1] as usize)
        } else {
            vec![b'o'; clus2]
        };

        let mut block1: Vec<Vec<u8>> = group1
            .iter()
            .map(|s| s[fs.cut1[i] as usize..fs.cut1[i + 1] as usize].to_vec())
            .collect();
        let mut block2: Vec<Vec<u8>> = group2
            .iter()
            .map(|s| s[fs.cut2[i] as usize..fs.cut2[i + 1] as usize].to_vec())
            .collect();

        if clus1 == 1 && clus2 == 1 {
            let mut s1 = std::mem::take(&mut block1[0]);
            let mut s2 = std::mem::take(&mut block2[0]);
            totalscore += ast.g_align11(
                cn,
                n_dynamicmtx,
                cn.penalty,
                cn.penalty_ex,
                &mut s1,
                &mut s2,
                headgp,
                tailgp,
                newgapstr,
            );
            block1[0] = s1;
            block2[0] = s2;
        } else {
            let (sc, _) = ast.a_align(
                cn,
                n_dynamicmtx,
                cn.penalty,
                cn.penalty_ex,
                &mut block1,
                &mut block2,
                eff1,
                eff2,
                headgp,
                tailgp,
                Some((&sgap1, &sgap2, &egap1, &egap2)),
                None,
                -1,
                -1,
                None,
                None,
                0.0,
                0.0,
                false,
                newgapstr,
            );
            totalscore += sc;
        }

        for (j, b) in block1.iter().enumerate() {
            result1[j].extend_from_slice(b);
        }
        for (j, b) in block2.iter().enumerate() {
            result2[j].extend_from_slice(b);
        }
    }

    for (j, r) in result1.into_iter().enumerate() {
        group1[j] = r;
    }
    for (j, r) in result2.into_iter().enumerate() {
        group2[j] = r;
    }

    totalscore
}

/// `seq_vec_3()` with `n20or4or2 == 1`: accumulate one member's contribution to
/// the indicator vector of alphabet index 0.
fn seq_vec_3(result: &mut [Cpx], incr: f64, seq: &[u8], cn: &Constants) {
    for (p, &c) in seq.iter().enumerate() {
        if cn.amino_n[c as usize] == 0 {
            result[p].r += incr;
        }
    }
}

/// `seq_vec_5()`: accumulate one member's contribution to the signal.
fn seq_vec_5(
    result: &mut [Cpx],
    score1: &[f64; 0x100],
    score2: &[f64; 0x100],
    incr: f64,
    seq: &[u8],
    cn: &Constants,
) {
    for (p, &c) in seq.iter().enumerate() {
        let n = cn.amino_n[c as usize];
        if n > 20 {
            continue;
        }
        result[p].r += incr * score1[n as usize];
        result[p].i += incr * score2[n as usize];
    }
}
