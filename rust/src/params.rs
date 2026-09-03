//! Command-line parsing, transliterated from `arguments()` in `core/disttbfast.c`.
//!
//! The parser accepts the same argument vector as the C `disttbfast` binary,
//! including its clustered single-letter form (`-Ot` == `-O -t`).  Options that
//! select a code path this port does not implement are accepted by the parser
//! and rejected later by [`Args::check_supported`], so that an unsupported run
//! fails loudly instead of quietly producing a different alignment.

/// `TMorJTT` values from `mltaln.h`.
pub const JTT: i32 = 0;
pub const TM: i32 = 1;

pub const NOTSPECIFIED: i32 = 100009;

/// `atof()` semantics: parse the longest numeric prefix, 0.0 on failure.
pub fn atof(s: &str) -> f64 {
    let b = s.as_bytes();
    let mut end = 0usize;
    let mut i = 0usize;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    let mut seen_digit = false;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        seen_digit = true;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            seen_digit = true;
        }
    }
    if seen_digit {
        end = i;
        // Optional exponent.
        if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
            let mut j = i + 1;
            if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
                j += 1;
            }
            let mut expdigit = false;
            while j < b.len() && b[j].is_ascii_digit() {
                j += 1;
                expdigit = true;
            }
            if expdigit {
                end = j;
            }
        }
    }
    if !seen_digit {
        return 0.0;
    }
    s[..end].parse::<f64>().unwrap_or(0.0)
}

/// `atoi()` semantics: parse the longest integer prefix, 0 on failure.
pub fn atoi(s: &str) -> i32 {
    let b = s.as_bytes();
    let mut i = 0usize;
    while i < b.len() && (b[i] as char).is_whitespace() {
        i += 1;
    }
    let start = i;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
    }
    if i == start {
        return 0;
    }
    s[start..i].parse::<i32>().unwrap_or(0)
}

/// Every global that `arguments()` sets.
///
/// A few fields are parsed but unused: they belong to C code paths this port
/// refuses to run, and are kept so the parser stays a faithful mirror of
/// `arguments()` rather than a subset that silently ignores arguments.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct Args {
    pub nthread: i32,
    pub nthreadpair: i32,
    pub nthreadtb: i32,
    pub outnumber: bool,
    /// Never set: the C's `case 'V': topin = 1` sits inside `#if 0`, and its
    /// branch in `disttbfast()` prints "--topin has been disabled" and exits.
    #[allow(dead_code)]
    pub topin: bool,
    pub treein: bool,
    pub treeout: i32,
    pub distout: u8,
    pub noalign: bool,
    pub nevermemsave: bool,
    pub inputfile: Option<String>,
    pub nadd: i32,
    pub addprofile: bool,
    pub constraint: bool,
    pub nblosum: i32,
    pub fmodel: i32,
    pub use_fft: bool,
    pub useexternalanchors: bool,
    pub oneiteration: bool,
    pub force_fft: bool,
    pub fftscore: bool,
    pub weight: i32,
    pub tbrweight: i32,
    pub tbutree: bool,
    pub disp: bool,
    pub outgap: bool,
    pub alg: u8,
    pub scmtd: i32,
    pub treemethod: u8,
    pub sueff_global: f64,
    pub scoremtx: i32,
    /// `'p'`, `'d'` or [`NOTSPECIFIED`] cast to a char slot.
    pub dorp: Option<u8>,
    pub ppenalty_dist: i32,
    pub ppenalty: i32,
    pub ppenalty_ex: i32,
    /// `-O` in the pairwise-stage parser: the generalised aligner's cost for
    /// opening a pair of unaligned runs.
    pub ppenalty_op: i32,
    pub penalty_shift_factor: f64,
    pub poffset: i32,
    pub kimura_r: i32,
    pub pam_n: i32,
    /// `JTT` (0) or `TM` (1); `-j` selects the first, `-m` the second.
    pub tm_or_jtt: i32,
    pub fft_win_size: i32,
    pub fft_threshold: i32,
    pub scoreout: bool,
    pub spscoreout: bool,
    pub tuplesize: i32,
    pub subalignment: bool,
    pub subalignmentoffset: i32,
    pub legacygapcost: bool,
    pub specificityconsideration: f64,
    pub nguidetree: i32,
    pub sparsepickup: i32,
    pub keeplength: bool,
    pub mapout: i32,
    pub smoothing: bool,
    pub nwildcard: bool,
    pub terminalmargin: i32,
    pub cut: f64,
}

impl Default for Args {
    fn default() -> Self {
        // Mirrors the initialisation block at the top of arguments().
        Args {
            nthread: 1,
            nthreadpair: 1,
            nthreadtb: 1,
            outnumber: false,
            topin: false,
            treein: false,
            treeout: 0,
            distout: 0,
            noalign: false,
            nevermemsave: false,
            inputfile: None,
            nadd: 0,
            addprofile: true,
            constraint: false,
            nblosum: 62,
            fmodel: 0,
            use_fft: false,
            useexternalanchors: false,
            oneiteration: false,
            force_fft: false,
            fftscore: true,
            weight: 3,
            tbrweight: 3,
            tbutree: true,
            disp: false,
            outgap: true,
            alg: b'A',
            scmtd: 5,
            treemethod: b'X',
            sueff_global: 0.1,
            scoremtx: 1,
            dorp: None,
            ppenalty_dist: NOTSPECIFIED,
            ppenalty: -1530,
            ppenalty_ex: NOTSPECIFIED,
            ppenalty_op: NOTSPECIFIED,
            penalty_shift_factor: 1000.0,
            poffset: -123,
            kimura_r: NOTSPECIFIED,
            pam_n: NOTSPECIFIED,
            tm_or_jtt: JTT,
            fft_win_size: NOTSPECIFIED,
            fft_threshold: NOTSPECIFIED,
            scoreout: false,
            spscoreout: false,
            tuplesize: 6,
            subalignment: false,
            subalignmentoffset: 0,
            legacygapcost: false,
            specificityconsideration: 0.0,
            nguidetree: 1,
            sparsepickup: 0,
            keeplength: false,
            mapout: 0,
            smoothing: false,
            nwildcard: false,
            terminalmargin: 0,
            cut: 0.0,
        }
    }
}

/// The unsupported-path report produced by [`Args::check_supported`].
pub struct Unsupported(pub String);

impl Args {
    /// Parses `argv` (without the program name) the way `arguments()` does.
    pub fn parse(argv: &[String]) -> Result<Args, String> {
        let mut a = Args::default();
        let mut i = 0usize;
        // C: while( --argc > 0 && (*++argv)[0] == '-' )
        'outer: while i < argv.len() && argv[i].starts_with('-') {
            let word: Vec<char> = argv[i].chars().collect();
            let mut k = 1usize; // skip the leading '-'
            while k < word.len() {
                let c = word[k];
                k += 1;
                // Helper closure state: options taking a value consume argv[i+1].
                macro_rules! val {
                    () => {{
                        i += 1;
                        if i >= argv.len() {
                            return Err(format!("option -{} needs an argument", c));
                        }
                        argv[i].clone()
                    }};
                }
                match c {
                    'i' => {
                        a.inputfile = Some(val!());
                        i += 1;
                        continue 'outer;
                    }
                    'I' => {
                        a.nadd = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'V' => {
                        a.ppenalty_dist = (atof(&val!()) * 1000.0 - 0.5) as i32;
                        i += 1;
                        continue 'outer;
                    }
                    'f' => {
                        a.ppenalty = (atof(&val!()) * 1000.0 - 0.5) as i32;
                        i += 1;
                        continue 'outer;
                    }
                    'Q' => {
                        a.penalty_shift_factor = atof(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'g' => {
                        a.ppenalty_ex = (atof(&val!()) * 1000.0 - 0.5) as i32;
                        i += 1;
                        continue 'outer;
                    }
                    'h' => {
                        a.poffset = (atof(&val!()) * 1000.0 - 0.5) as i32;
                        i += 1;
                        continue 'outer;
                    }
                    'k' => {
                        a.kimura_r = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'b' => {
                        a.nblosum = atoi(&val!());
                        a.scoremtx = 1;
                        i += 1;
                        continue 'outer;
                    }
                    'j' => {
                        a.pam_n = atoi(&val!());
                        a.scoremtx = 0;
                        a.tm_or_jtt = JTT;
                        i += 1;
                        continue 'outer;
                    }
                    'm' => {
                        a.pam_n = atoi(&val!());
                        a.scoremtx = 0;
                        a.tm_or_jtt = TM;
                        i += 1;
                        continue 'outer;
                    }
                    'C' => {
                        let v = val!();
                        a.nthreadpair = atoi(&v);
                        a.nthread = a.nthreadpair;
                        a.nthreadtb = match v.find('-') {
                            Some(p) => atoi(&v[p + 1..]),
                            None => a.nthread,
                        };
                        i += 1;
                        continue 'outer;
                    }
                    's' => {
                        a.specificityconsideration = atof(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'a' => a.fmodel = 1,
                    'K' => a.addprofile = false,
                    'y' => {
                        let v = val!();
                        a.distout = v.as_bytes().first().copied().unwrap_or(0);
                        if a.distout != b'c' && a.distout != b'h' {
                            return Err("Set -y c or -y h in v>=7.521.".into());
                        }
                        i += 1;
                        continue 'outer;
                    }
                    't' => a.treeout = 1,
                    '^' => a.treeout = 2,
                    'T' => a.noalign = true,
                    'r' => a.oneiteration = true,
                    'D' => a.dorp = Some(b'd'),
                    'P' => a.dorp = Some(b'p'),
                    'L' => a.legacygapcost = true,
                    'e' => a.fftscore = false,
                    'x' => {
                        a.terminalmargin = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'H' => {
                        a.subalignment = true;
                        a.subalignmentoffset = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'n' => a.outnumber = true,
                    'q' => {
                        a.sparsepickup = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'X' => {
                        a.treemethod = b'X';
                        a.sueff_global = atof(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'E' => {
                        a.nguidetree = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'A' => a.alg = b'A',
                    '&' => a.alg = b'a',
                    '@' => a.alg = b'd',
                    'N' => a.nevermemsave = true,
                    'M' => a.alg = b'M',
                    'S' => a.spscoreout = true,
                    'B' => {} // accepted and ignored, as in the C
                    'F' => a.use_fft = true,
                    'l' => {
                        // C falls through from 'l' into 'G'.
                        a.useexternalanchors = true;
                        a.use_fft = true;
                        a.force_fft = true;
                    }
                    'G' => {
                        a.use_fft = true;
                        a.force_fft = true;
                    }
                    'U' => a.treein = true,
                    'u' => {
                        a.weight = 0;
                        a.tbrweight = 0;
                    }
                    'v' => a.tbrweight = 3,
                    'd' => a.disp = true,
                    'O' => a.outgap = false,
                    'J' => a.tbutree = false,
                    'w' => {
                        a.fft_win_size = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'W' => {
                        a.tuplesize = atoi(&val!());
                        i += 1;
                        continue 'outer;
                    }
                    'Y' => a.keeplength = true,
                    'z' => a.mapout = 2,
                    'Z' => a.mapout = 1,
                    'p' => a.smoothing = true,
                    ':' => a.nwildcard = true,
                    other => return Err(format!("illegal option {}", other)),
                }
            }
            i += 1;
        }
        // C: a lone trailing non-option argument is the (unused) `cut` value.
        if i + 1 == argv.len() {
            a.cut = atof(&argv[i]);
            i += 1;
        }
        if i != argv.len() {
            return Err("options: Check source file !".into());
        }
        Ok(a)
    }

    /// Rejects every configuration whose C code path is outside this port.
    ///
    /// The port covers the progressive (`treebase`) path with the full
    /// group-to-group DP: `alg == 'A'`, no FFT anchoring, no constraints, no
    /// `--add`, no sub-alignments, no external anchors and no threads.
    pub fn check_supported(&self) -> Result<(), Unsupported> {
        let mut bad: Vec<&str> = Vec::new();
        if self.useexternalanchors {
            bad.push("-l (external anchors) is not ported");
        }
        if self.alg != b'A' {
            bad.push("only alg 'A' (-A) is ported; -M/-&/-@ are not");
        }
        if self.subalignment {
            bad.push("-H (sub-alignments) is not ported");
        }
        if self.oneiteration {
            bad.push("-r (one iteration) is not ported");
        }
        if self.keeplength && self.nadd == 0 {
            // The C reaches `restoreoriginalgaps( njob, bseq, NULL )` and
            // segfaults; there is no output to be identical to.
            bad.push("-Y (--keeplength) without -I <n> crashes the C too");
        }
        if self.nthread > 0 || self.nthreadpair > 0 || self.nthreadtb > 0 {
            bad.push("threads change the result; run with -C 0-0");
        }
        // For nucleotides `constants()` ignores -b and derives the matrix from
        // Kimura's model instead, so the BLOSUM check applies to protein only.
        if self.dorp != Some(b'd')
            && self.scoremtx == 1
            && !matches!(self.nblosum, 30 | 45 | 50 | 62 | 80 | 90 | 100 | -1)
        {
            // -1 is --aamatrix, a table read from `_aamtx`. -2 is --text,
            // which also switches the alphabet to 256 symbols.
            bad.push("only the built-in BLOSUM tables and --aamatrix are ported");
        }
        if self.scoremtx == 2 {
            bad.push("-b 2 (Miyata-Yasunaga) is not supported by the C either");
        }

        if bad.is_empty() {
            Ok(())
        } else {
            Err(Unsupported(bad.join("\n  ")))
        }
    }
}
