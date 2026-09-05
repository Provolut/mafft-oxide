//! FASTA input/output, transliterated from `core/io.c`.
//!
//! MAFFT's reader is character-oriented (`getc`/`ungetc`/`myfgets`) and its
//! quirks are load-bearing: names are truncated at `B-2` bytes, sequence
//! characters outside `[A-Za-z*.-]` are dropped, `*` becomes `-`, and residues
//! are upper-cased.  [`CFile`] reproduces the stdio behaviour the C relies on.

/// Name buffer size (`B` in mltaln.h).  Names are stored as `'=' + text`.
pub const B: usize = 256;

/// A byte buffer with C stdio `getc`/`ungetc`/`feof` semantics.
pub struct CFile {
    buf: Vec<u8>,
    pos: usize,
    eof: bool,
}

impl CFile {
    pub fn new(buf: Vec<u8>) -> CFile {
        CFile {
            buf,
            pos: 0,
            eof: false,
        }
    }

    pub fn rewind(&mut self) {
        self.pos = 0;
        self.eof = false;
    }

    /// `getc()`: returns `None` for EOF and latches the end-of-file flag.
    pub fn getc(&mut self) -> Option<u8> {
        if self.pos < self.buf.len() {
            let c = self.buf[self.pos];
            self.pos += 1;
            Some(c)
        } else {
            self.eof = true;
            None
        }
    }

    /// `ungetc()`: pushing back EOF is a no-op, as in C.
    pub fn ungetc(&mut self, c: Option<u8>) {
        if c.is_some() && self.pos > 0 {
            self.pos -= 1;
            self.eof = false;
        }
    }

    pub fn feof(&self) -> bool {
        self.eof
    }
}

/// `countKUorWA()`: number of lines starting with `>`.
fn count_ku_or_wa(fp: &mut CFile) -> usize {
    let mut value = 0usize;
    let mut b = b'\n';
    while let Some(c) = fp.getc() {
        if b == b'\n' && c == b'>' {
            value += 1;
        }
        b = c;
    }
    fp.rewind();
    value
}

/// `searchKUorWA()`: position the cursor on the first line-initial `>`.
fn search_ku_or_wa(fp: &mut CFile) {
    let mut b = b'\n';
    loop {
        let c = fp.getc();
        let is_stop = match c {
            Some(b'>') => b == b'\n',
            None => b == b'\n',
            _ => false,
        };
        if is_stop {
            fp.ungetc(c);
            return;
        }
        b = c.unwrap();
    }
}

/// `myfgets()`: read at most `l` bytes up to (and consuming) the newline.
fn myfgets(fp: &mut CFile, l: usize) -> Vec<u8> {
    let mut s = Vec::new();
    if fp.feof() {
        return s;
    }
    let mut c: Option<u8> = None;
    let mut i = 0usize;
    while i < l {
        c = fp.getc();
        match c {
            Some(b'\n') => break,
            Some(ch) => s.push(ch),
            None => break,
        }
        i += 1;
    }
    if c != Some(b'\n') {
        // Skip to the end of the line.
        loop {
            match fp.getc() {
                Some(b'\n') | None => break,
                _ => {}
            }
        }
    }
    s
}

/// `onlyAlpha_upper()` / `onlyAlpha_lower()` followed by `kake2hiku()`.
///
/// The C picks the case from the global `dorp`, which is still unset while
/// `getnumlen()` scans the file, so the first pass always upper-cases.
fn filter_seq(raw: &[u8], lower: bool) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len());
    for &ch in raw {
        if ch.is_ascii_alphabetic() || ch == b'-' || ch == b'*' || ch == b'.' {
            out.push(if ch == b'*' {
                b'-'
            } else if lower {
                ch.to_ascii_lowercase()
            } else {
                ch.to_ascii_uppercase()
            });
        }
    }
    out
}

/// How the reader cleans up a record's body. `--text` (`nblosum == -2`) keeps
/// every byte and only drops whitespace, where the ordinary modes keep letters
/// and case-fold.
#[derive(Clone, Copy, PartialEq)]
pub enum SeqFilter {
    Upper,
    Lower,
    /// `charfilter()`.
    Text,
}

/// `charfilter()`: drop LF, space and CR; keep everything else, whatever it is.
/// `=`, `<` and `>` are refused because the reader uses them as markers.
fn charfilter(raw: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(raw.len());
    for &c in raw {
        if c == b'=' || c == b'<' || c == b'>' {
            eprintln!();
            eprintln!("Characters '= < >' can be used only in the title lines in the --anysymbol or --text mode.");
            eprintln!();
            std::process::exit(1);
        }
        if c != 0x0a && c != 0x20 && c != 0x0d {
            out.push(c);
        }
    }
    out
}

/// `load1SeqWithoutName_realloc()`: read residues up to the next record.
fn load1_seq_without_name(fp: &mut CFile) -> Vec<u8> {
    load1_seq_without_name_case(fp, SeqFilter::Upper)
}

fn load1_seq_without_name_case(fp: &mut CFile, filter: SeqFilter) -> Vec<u8> {
    let mut raw = Vec::new();
    let mut b = b'\n';
    loop {
        let c = fp.getc();
        match c {
            None => {
                fp.ungetc(c);
                break;
            }
            Some(ch) => {
                if ch == b'>' && b == b'\n' {
                    fp.ungetc(c);
                    break;
                }
                raw.push(ch);
                b = ch;
            }
        }
    }
    match filter {
        SeqFilter::Text => charfilter(&raw),
        SeqFilter::Lower => filter_seq(&raw, true),
        SeqFilter::Upper => filter_seq(&raw, false),
    }
}

/// `countATGC()`: (nucleotide-ish characters, total alphabetic characters).
fn count_atgc(s: &[u8]) -> (i64, i64) {
    let mut n_atgc = 0i64;
    let mut n_char = 0i64;
    for &ch in s {
        let c = ch.to_ascii_lowercase();
        if c.is_ascii_alphabetic() {
            n_char += 1;
            if matches!(c, b'a' | b't' | b'g' | b'c' | b'u' | b'n') {
                n_atgc += 1;
            }
        }
    }
    (n_atgc, n_char)
}

/// `countATGCandN()`: (nucleotide-ish characters, N characters, total letters).
///
/// Separate from [`count_atgc`] rather than folded into it because the C keeps
/// two functions that differ only in whether they hand back the `n` tally, and
/// their callers are on different paths.
fn count_atgc_and_n(s: &[u8]) -> (i64, i64, i64) {
    let mut n_atgc = 0i64;
    let mut n_n = 0i64;
    let mut n_char = 0i64;
    for &ch in s {
        let c = ch.to_ascii_lowercase();
        if c.is_ascii_alphabetic() {
            n_char += 1;
            if matches!(c, b'a' | b't' | b'g' | b'c' | b'u' | b'n') {
                n_atgc += 1;
            }
            if c == b'n' {
                n_n += 1;
            }
        }
    }
    (n_atgc, n_n, n_char)
}

/// `countnogaplen()`: residues that are not `-`.
fn count_nogap_len(seq: &[u8]) -> usize {
    seq.iter().filter(|&&c| c != b'-').count()
}

/// What `getnumlen_nogap_countn()` measures about an input.
///
/// Attributes are the four values the C leaves in globals plus the two it
/// returns through pointers, which is why they travel together here: `countlen`
/// prints all of them in one line and nothing else uses them.
pub struct LengthCounts {
    /// `njob`: number of records.
    pub njob: usize,
    /// `nlenmax`: longest ungapped record.
    pub nlenmax: usize,
    /// The shortest ungapped record (the C's `*nlenminpt`).
    pub nlenmin: usize,
    /// `dorp`: `b'd'` for nucleotide, `b'p'` for protein, decided by the
    /// nucleotide-character fraction over the first 100000 sites.
    pub dorp: u8,
    /// The C's `*nfreq`: N characters per nucleotide-ish character.
    pub nfreq: f64,
}

/// `getnumlen_nogap_countn()`: measure an input without aligning it.
///
/// The sampling cap is the C's and is load-bearing rather than an optimisation:
/// only the first 100000 sites feed the nucleotide-fraction decision, so a file
/// whose head and tail differ in composition is typed from its head, and a port
/// that read all of it would disagree with the C on exactly those inputs.
///
/// Args are the raw file bytes and the `-D`/`-P` override, if one was given.
///
/// Returns:
///     The record count, the longest and shortest ungapped lengths, the
///     inferred sequence type, and the N frequency.
pub fn getnumlen_nogap_countn(bytes: Vec<u8>, dorp_in: Option<u8>) -> LengthCounts {
    let mut fp = CFile::new(bytes);
    let njob = count_ku_or_wa(&mut fp);
    search_ku_or_wa(&mut fp);
    let mut nlenmax = 0usize;
    let mut nlenmin = 99_999_999usize;
    let mut atgcnum = 0i64;
    let mut total = 0i64;
    let mut nnum = 0i64;
    // The C's `nN` is a function-scope local that only `countATGCandN` writes,
    // while `nnum += nN` sits *outside* the `total < 100000` guard. So once the
    // sampling cap is reached, every remaining record adds the last sampled
    // record's N count again. Transliterated rather than tidied: `nfreq` is
    // printed, and correcting the C here would make the two disagree on any
    // input longer than 100000 sites -- which is most of them.
    //
    // The zero here is a decision, not a transliteration, and it is the one
    // place this function knowingly differs from the C. The C leaves `nN`
    // uninitialised, and `countATGCandN` returns early for an empty sequence
    // body writing `*total` but not `*countN` -- so an input whose first record
    // has an empty body reaches the `nfreq` division with a stack value.
    // Valgrind names the allocation; upstream issue #8 has the four-line input
    // and the trace. There is no bit-exact answer to reproduce, so this picks
    // the defined one: zero, which agrees with the C on every input where the
    // C itself is defined. `verify_utils.sh` scores that input as
    // indeterminate rather than as a pass or a failure.
    let mut n_n = 0i64;
    for _ in 0..njob {
        let _name = myfgets(&mut fp, 5_000_000 - 1);
        let seq = load1_seq_without_name(&mut fp);
        let tmp = count_nogap_len(&seq);
        if tmp > nlenmax {
            nlenmax = tmp;
        }
        if tmp < nlenmin {
            nlenmin = tmp;
        }
        if total < 100_000 {
            let (n_atgc, sampled_n, n_site) = count_atgc_and_n(&seq);
            atgcnum += n_atgc;
            total += n_site;
            n_n = sampled_n;
        }
        nnum += n_n;
    }
    let atgcfreq = atgcnum as f64 / total as f64;
    let nfreq = nnum as f64 / atgcnum as f64;
    let dorp = dorp_in.unwrap_or(if atgcfreq > 0.75 { b'd' } else { b'p' });
    LengthCounts {
        njob,
        nlenmax,
        nlenmin,
        dorp,
        nfreq,
    }
}

/// One input record: the MAFFT-style name (with its leading `'='`) and residues.
pub struct Input {
    pub njob: usize,
    /// Longest input sequence; the C sizes its buffers from this.
    #[allow(dead_code)]
    pub nlenmax: usize,
    pub dorp: u8,
    pub name: Vec<Vec<u8>>,
    pub seq: Vec<Vec<u8>>,
}

/// `getnumlen()` + `readData_pointer()`.
///
/// `dorp_in` is the sequence type forced on the command line, if any.
pub fn read_data(bytes: Vec<u8>, dorp_in: Option<u8>, outnumber: bool, text: bool) -> Input {
    let mut fp = CFile::new(bytes);

    // --- getnumlen() ---
    let njob = count_ku_or_wa(&mut fp);
    search_ku_or_wa(&mut fp);
    let mut nlenmax = 0usize;
    let mut atgcnum = 0i64;
    let mut total = 0i64;
    for _ in 0..njob {
        let _ = myfgets(&mut fp, 5_000_000 - 1);
        // `getnumlen` runs before `constants()`, so `dorp` is still unset and
        // the ordinary path upper-cases; only `--text` differs here.
        let tmpseq = if text {
            load1_seq_without_name_case(&mut fp, SeqFilter::Text)
        } else {
            load1_seq_without_name(&mut fp)
        };
        if tmpseq.len() > nlenmax {
            nlenmax = tmpseq.len();
        }
        if total < 1_000_000 {
            let (a, t) = count_atgc(&tmpseq);
            atgcnum += a;
            total += t;
        }
    }
    let atgcfreq = atgcnum as f64 / total as f64;
    let dorp = match dorp_in {
        Some(d) => d,
        None => {
            if atgcfreq > 0.75 {
                b'd'
            } else {
                b'p'
            }
        }
    };

    // --- readData_pointer() ---
    fp.rewind();
    search_ku_or_wa(&mut fp);
    let mut name = Vec::with_capacity(njob);
    let mut seq = Vec::with_capacity(njob);
    for _ in 0..njob {
        fp.getc(); // consume '>'
        let mut nm = vec![b'='];
        nm.extend_from_slice(&myfgets(&mut fp, B - 2));
        name.push(nm);
        let filter = if text {
            SeqFilter::Text
        } else if dorp == b'd' {
            SeqFilter::Lower
        } else {
            SeqFilter::Upper
        };
        seq.push(load1_seq_without_name_case(&mut fp, filter));
    }

    if outnumber {
        for (i, nm) in name.iter_mut().enumerate() {
            let tail: Vec<u8> = match find_sub(nm, b"_numo_e_") {
                Some(p) => nm[p + 8..].to_vec(),
                None => nm[1..].to_vec(),
            };
            let mut buf = Vec::new();
            buf.push(b'=');
            buf.extend_from_slice(format!("_numo_s_{:08}_numo_e_", i + 1).as_bytes());
            buf.extend_from_slice(&tail);
            buf.truncate(B - 1);
            *nm = buf;
        }
    }

    Input {
        njob,
        nlenmax,
        dorp,
        name,
        seq,
    }
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    hay.windows(needle.len()).position(|w| w == needle)
}

/// `writeData_pointer()` with the default `LineLengthInFASTA < 0`.
pub fn write_data(out: &mut Vec<u8>, name: &[Vec<u8>], aseq: &[Vec<u8>]) {
    for i in 0..name.len() {
        out.push(b'>');
        out.extend_from_slice(&name[i][1..]);
        out.push(b'\n');
        out.extend_from_slice(&aseq[i]);
        out.push(b'\n');
    }
}

/// `gappick0()`: strip `-` from a sequence.
pub fn gappick0(seq: &[u8]) -> Vec<u8> {
    seq.iter().copied().filter(|&c| c != b'-').collect()
}

/// `WriteFloatHat2_pointer_halfmtx()`: the `hat2` distance-matrix file.
///
/// The line breaks are not cosmetic -- `dvtditr` parses this back, and the C
/// wraps after every 12th entry *or* at the end of a row, whichever comes
/// first.
pub fn write_hat2(njob: usize, name: &[Vec<u8>], mtx: &[Vec<f64>]) {
    use std::fmt::Write as _;
    let mut max = 0.0f64;
    for i in 0..njob - 1 {
        for j in 1..njob - i {
            if mtx[i][j] > max {
                max = mtx[i][j];
            }
        }
    }
    let mut s = String::new();
    let _ = writeln!(s, "{:5}", 1);
    let _ = writeln!(s, "{:5}", njob);
    let _ = writeln!(s, " {:6.3}", max * 2.5);
    for (i, nm) in name.iter().enumerate().take(njob) {
        let _ = writeln!(s, "{:4}. {}", i + 1, String::from_utf8_lossy(nm));
    }
    for i in 0..njob {
        for j in i + 1..njob {
            let _ = write!(s, "{:6.3}", mtx[i][j - i]);
            let ijsa = j - i;
            if ijsa % 12 == 0 || ijsa == njob - i - 1 {
                s.push('\n');
            }
        }
    }
    std::fs::write("hat2", s).expect("write hat2");
}

/// `readsubalignmentstable()`: `--merge`'s `_subalignmentstable`.
///
/// One group per line, sequences given as 1-origin numbers. A negative number
/// anywhere on the line sets that group's `preservegaps` flag (and its
/// magnitude is the sequence). Lines starting with `#`, and lines whose first
/// token reads as 0, are skipped.
///
/// Returns the groups (0-origin, no `-1` terminator) and their `preservegaps`
/// flags.
pub fn read_subalignments_table(nseq: usize) -> (Vec<Vec<i32>>, Vec<bool>) {
    let text = match std::fs::read("_subalignmentstable") {
        Ok(t) => t,
        Err(_) => {
            eprintln!("Cannot open _subalignmentstable");
            std::process::exit(1);
        }
    };
    let mut groups: Vec<Vec<i32>> = Vec::new();
    let mut preservegaps: Vec<bool> = Vec::new();
    let mut seen = vec![false; nseq];

    // `fgets` then `if( line[strlen(line)-1] != '\n' )`: every line the C
    // accepts ends in a newline, so a file whose last line does not is an
    // error -- except that `feof` breaks the loop first, so a missing final
    // newline just drops that line.
    for raw in text.split_inclusive(|&c| c == b'\n') {
        if raw.last() != Some(&b'\n') {
            break;
        }
        if raw.first() == Some(&b'#') {
            continue;
        }
        if atoi_bytes(raw) == 0 {
            continue;
        }
        let (mem, pg) = read_a_subalignment(raw);
        for &p in mem.iter() {
            if p < 0 || p as usize > nseq - 1 {
                eprintln!(
                    "Sequence {} does not exist in the input sequence file.",
                    p + 1
                );
                std::process::exit(1);
            }
            if seen[p as usize] {
                eprintln!("\nSequence {} appears in different groups.", p + 1);
                eprintln!("Hierarchical grouping is not supported.\n");
                std::process::exit(1);
            }
            seen[p as usize] = true;
        }
        groups.push(mem);
        preservegaps.push(pg);
    }
    (groups, preservegaps)
}

/// `readasubalignment()`: the numbers on one line, tabs counting as spaces.
fn read_a_subalignment(line: &[u8]) -> (Vec<i32>, bool) {
    let mut out = Vec::new();
    let mut preservegaps = false;
    let mut at_space = true;
    for (i, &c) in line.iter().enumerate() {
        if c == b' ' || c == b'\t' {
            at_space = true;
            continue;
        }
        if at_space {
            if c == b'\n' || c == b'#' {
                break;
            }
            at_space = false;
            let v = atoi_bytes(&line[i..]);
            if v == 0 {
                eprintln!("Format error? Sequences must be specified as 1, 2, 3...");
                std::process::exit(1);
            }
            if v < 0 {
                preservegaps = true;
            }
            out.push(v.abs() - 1);
        }
    }
    (out, preservegaps)
}

/// C's `atoi` on a byte slice: leading space, optional sign, digits, stop.
fn atoi_bytes(s: &[u8]) -> i32 {
    let mut i = 0usize;
    while i < s.len() && (s[i] == b' ' || s[i] == b'\t') {
        i += 1;
    }
    let neg = if i < s.len() && (s[i] == b'-' || s[i] == b'+') {
        let n = s[i] == b'-';
        i += 1;
        n
    } else {
        false
    };
    let start = i;
    let mut v: i64 = 0;
    while i < s.len() && s[i].is_ascii_digit() {
        v = v * 10 + (s[i] - b'0') as i64;
        if v > i32::MAX as i64 + 1 {
            v = i32::MAX as i64 + 1;
        }
        i += 1;
    }
    if i == start {
        return 0;
    }
    if neg {
        (-v) as i32
    } else {
        v as i32
    }
}

#[cfg(test)]
mod subaln_tests {
    use super::*;
    use std::io::Write as _;

    fn with_table<T>(tag: &str, contents: &str, f: impl FnOnce() -> T) -> T {
        let _g = crate::cwd_lock();
        let d =
            std::env::temp_dir().join(format!("mafft-rs-subaln-{}-{}", tag, std::process::id()));
        let _ = std::fs::create_dir_all(&d);
        let mut fh = std::fs::File::create(d.join("_subalignmentstable")).unwrap();
        fh.write_all(contents.as_bytes()).unwrap();
        drop(fh);
        let old = std::env::current_dir().unwrap();
        std::env::set_current_dir(&d).unwrap();
        let r = f();
        std::env::set_current_dir(old).unwrap();
        let _ = std::fs::remove_dir_all(&d);
        r
    }

    #[test]
    fn groups_are_read_one_per_line_and_made_zero_origin() {
        let (g, p) = with_table("plain", "1 2 3\n4 5\n", || read_subalignments_table(6));
        assert_eq!(g, vec![vec![0, 1, 2], vec![3, 4]]);
        assert_eq!(p, vec![false, false]);
    }

    /// A negative number anywhere on the line sets that group's
    /// `preservegaps`, and its magnitude is still the sequence.
    #[test]
    fn a_negative_number_means_preservegaps() {
        let (g, p) = with_table("neg", "1 2\n-3 4\n", || read_subalignments_table(5));
        assert_eq!(g, vec![vec![0, 1], vec![2, 3]]);
        assert_eq!(p, vec![false, true]);
    }

    #[test]
    fn comments_and_trailing_comments_are_skipped() {
        let (g, _) = with_table("cmt", "# a header\n1 2 # rest ignored\n3 4\n", || {
            read_subalignments_table(5)
        });
        assert_eq!(g, vec![vec![0, 1], vec![2, 3]]);
    }

    /// `fgets` + `feof`: a last line with no newline is dropped, not read.
    #[test]
    fn a_final_line_without_a_newline_is_dropped() {
        let (g, _) = with_table("nonl", "1 2\n3 4", || read_subalignments_table(5));
        assert_eq!(g, vec![vec![0, 1]]);
    }

    #[test]
    fn tabs_count_as_separators() {
        let (g, _) = with_table("tab", "1\t2\t3\n", || read_subalignments_table(4));
        assert_eq!(g, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn atoi_bytes_matches_c() {
        assert_eq!(atoi_bytes(b"  -12abc"), -12);
        assert_eq!(atoi_bytes(b"7 8"), 7);
        assert_eq!(atoi_bytes(b"x9"), 0);
        assert_eq!(atoi_bytes(b""), 0);
    }
}

#[cfg(test)]
mod charfilter_tests {
    use super::*;

    /// `--text` keeps every byte the ordinary reader would drop.
    #[test]
    fn the_text_filter_keeps_punctuation_and_digits() {
        let fasta = b">a\n!$%&@01\n>b\n{}[]~+/\n".to_vec();
        let r = read_data(fasta, Some(b'p'), false, true);
        assert_eq!(r.njob, 2);
        assert_eq!(r.seq[0], b"!$%&@01".to_vec());
        assert_eq!(r.seq[1], b"{}[]~+/".to_vec());
    }

    /// Only LF, space and CR go; case is left alone, unlike every other mode.
    #[test]
    fn the_text_filter_drops_only_whitespace_and_keeps_case() {
        let fasta = b">a\nAb Cd\tEf\n".to_vec();
        let r = read_data(fasta, Some(b'p'), false, true);
        // The tab is not in the C's drop list, so it survives.
        assert_eq!(r.seq[0], b"AbCd\tEf".to_vec());
    }

    /// The ordinary reader would have dropped all of that.
    #[test]
    fn the_ordinary_filter_still_drops_non_letters() {
        let fasta = b">a\nAb Cd!1\n".to_vec();
        let r = read_data(fasta, Some(b'p'), false, false);
        assert_eq!(r.seq[0], b"ABCD".to_vec());
    }
}

/// `input_new()`: one fixed-width number from a `hat2` matrix.
///
/// `D` is 6, and the field really is positional -- the C reads six characters
/// and calls `atof` on them, so a value that needs more room silently loses its
/// tail on both sides. One leading newline is skipped, which is what lets the
/// writer break the matrix into lines.
fn input_new(fp: &mut CFile, d: usize) -> f64 {
    let c = fp.getc();
    if c != Some(b'\n') {
        fp.ungetc(c);
    }
    let mut s = String::new();
    for _ in 0..d {
        match fp.getc() {
            Some(ch) => s.push(ch as char),
            None => break,
        }
    }
    // `atof` on a partly-numeric field takes the leading number and ignores the
    // rest, and returns 0.0 when there is none.
    let t = s.trim_start();
    let end = t
        .find(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == '+' || c == 'e' || c == 'E'))
        .unwrap_or(t.len());
    t[..end].parse::<f64>().unwrap_or(0.0)
}

/// Consume one line, as `fgets`/`myfgets` do.
fn skip_line(fp: &mut CFile) {
    while let Some(c) = fp.getc() {
        if c == b'\n' {
            return;
        }
    }
}

/// `readhat2_pointer()`: the pairwise distance matrix `dvtditr` builds its tree
/// from.
///
/// Three header lines, then one line per sequence (the names, which the C reads
/// and discards -- its name check is `#if 0`), then the upper triangle.
///
/// Args:
///     bytes: The `hat2` file.
///     nseq: How many sequences it must describe.
///
/// Returns:
///     `mtx[i][j]` for `i < j`; the rest is zero, as the C leaves it.
pub fn read_hat2(bytes: Vec<u8>, nseq: usize) -> Vec<Vec<f64>> {
    let mut fp = CFile::new(bytes);
    skip_line(&mut fp);
    let mut second = Vec::new();
    while let Some(c) = fp.getc() {
        if c == b'\n' {
            break;
        }
        second.push(c);
    }
    // `b[5] = 0` before `atoi`: the count is the first five characters.
    second.truncate(5);
    let nseq0: usize = String::from_utf8_lossy(&second).trim().parse().unwrap_or(0);
    if nseq != nseq0 {
        eprintln!("hat2 is wrong.");
        std::process::exit(1);
    }
    skip_line(&mut fp);
    for _ in 0..nseq {
        skip_line(&mut fp);
    }
    let mut mtx = vec![vec![0.0f64; nseq]; nseq];
    for i in 0..nseq.saturating_sub(1) {
        for j in i + 1..nseq {
            mtx[i][j] = input_new(&mut fp, 6);
        }
    }
    mtx
}

/// `readlocalhomtable2_half()`: the `hat3` constraint table.
///
/// One segment per line: `i j overlapaa opt start1 end1 start2 end2 infor`,
/// with `i < j` required. Several segments between the same pair chain onto one
/// another.
///
/// The `opt` stored is not the `opt` read: the C rescales it by `/ 5.8 * 600`.
///
/// Args:
///     bytes: The `hat3` file.
///     njob: Number of sequences.
///     table: The half table to fill, from `localhom::new_table`.
///
/// Returns:
///     Which sequences carry a structural (`korh == 'k'`) constraint.
pub fn read_localhomtable2_half(
    bytes: Vec<u8>,
    njob: usize,
    table: &mut crate::localhom::LocalHomTable,
) -> Vec<bool> {
    let mut kozoarivec = vec![false; njob];
    for line in bytes.split(|&c| c == b'\n') {
        if line.is_empty() {
            continue;
        }
        let text = String::from_utf8_lossy(line);
        let f: Vec<&str> = text.split_whitespace().collect();
        if f.len() < 8 {
            continue;
        }
        let i: usize = f[0].parse().unwrap_or(0);
        let j: usize = f[1].parse().unwrap_or(0);
        let overlapaa: i32 = f[2].parse().unwrap_or(0);
        let opt: f64 = f[3].parse().unwrap_or(0.0);
        let start1: i32 = f[4].parse().unwrap_or(0);
        let end1: i32 = f[5].parse().unwrap_or(0);
        let start2: i32 = f[6].parse().unwrap_or(0);
        let end2: i32 = f[7].parse().unwrap_or(0);
        let korh = f.get(8).and_then(|s| s.bytes().next()).unwrap_or(b'h');
        if korh == b'k' {
            kozoarivec[i] = true;
            kozoarivec[j] = true;
        }
        if j <= i || i >= njob || j >= njob {
            eprintln!("Check hat3.  The first sequence must be younger than the second one.");
            std::process::exit(1);
        }
        let e = crate::localhom::LocalHom {
            start1,
            end1,
            start2,
            end2,
            opt: opt / 5.8 * 600.0,
            overlapaa,
            korh,
            ..crate::localhom::LocalHom::default()
        };
        let chain = &mut table[i][j - i];
        // The C tracks `nokori` to decide between filling the head slot and
        // appending; the head is the one `new_table` made.
        if chain[0].opt == -1.0 && chain.len() == 1 {
            chain[0] = e;
        } else {
            chain.push(e);
        }
    }
    kozoarivec
}
