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

/// `load1SeqWithoutName_realloc()`: read residues up to the next record.
fn load1_seq_without_name(fp: &mut CFile) -> Vec<u8> {
    load1_seq_without_name_case(fp, false)
}

fn load1_seq_without_name_case(fp: &mut CFile, lower: bool) -> Vec<u8> {
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
    filter_seq(&raw, lower)
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
pub fn read_data(bytes: Vec<u8>, dorp_in: Option<u8>, outnumber: bool) -> Input {
    let mut fp = CFile::new(bytes);

    // --- getnumlen() ---
    let njob = count_ku_or_wa(&mut fp);
    search_ku_or_wa(&mut fp);
    let mut nlenmax = 0usize;
    let mut atgcnum = 0i64;
    let mut total = 0i64;
    for _ in 0..njob {
        let _ = myfgets(&mut fp, 5_000_000 - 1);
        let tmpseq = load1_seq_without_name(&mut fp);
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
        seq.push(load1_seq_without_name_case(&mut fp, dorp == b'd'));
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
