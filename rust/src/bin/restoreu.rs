//! `restoreu-rs` — a bit-exact Rust port of MAFFT's `restoreu`.
//!
//! The inverse of `replaceu`, and the last thing to touch the residues. It
//! takes the finished alignment, whose records are named `_os_<i>_oe_<name>`
//! and whose characters were folded to a known alphabet, and writes the
//! *original* characters back onto that alignment's gap pattern — restoring
//! case and every unusual symbol `replaceu` had mapped to `X` or `n`.
//!
//! So the alignment decides where the gaps are and the original decides what
//! the residues are. That split is why both files are read: `-a` the
//! alignment, `-i` the original.
#![forbid(unsafe_code)]

use std::io::Write;

use mafft_oxide::cargs::{self, Opt};
use mafft_oxide::io;

/// What the option loop leaves behind.
struct Args {
    inputfile: Option<String>,
    alignmentfile: Option<String>,
    /// `-Y`: re-insert the gaps `_deletelist` records, so the output keeps the
    /// input's length.
    keeplength: bool,
    /// `-z` / `-Z`: 2 and 1 respectively. Only 2 changes anything here, and
    /// only for records the delete list touched.
    mapout: i32,
}

/// Walk `arguments()`.
///
/// Returns:
///     The parsed options.
fn arguments(argv: &[String]) -> Args {
    let mut out = Args {
        inputfile: None,
        alignmentfile: None,
        keeplength: false,
        mapout: 0,
    };
    let opts = cargs::parse(argv, |c| match c {
        'i' | 'a' => Some(Opt::Value),
        'z' | 'Z' | 'p' | 'Y' => Some(Opt::Flag),
        _ => None,
    });
    for o in &opts {
        match o.c {
            'i' => out.inputfile.clone_from(&o.value),
            'a' => out.alignmentfile.clone_from(&o.value),
            'z' => out.mapout = 2,
            'Z' => out.mapout = 1,
            // `-p` exists only so `add2ndharfarg` can be passed through.
            'p' => {}
            'Y' => out.keeplength = true,
            _ => unreachable!("cargs only yields options this binary classified"),
        }
    }
    out
}

/// `gappick_samestring()`: drop every `-`, in place.
fn gappick_samestring(seq: &mut Vec<u8>) {
    seq.retain(|&c| c != b'-');
}

/// `strncmp( a, b, n ) == 0` over NUL-free byte strings.
///
/// The C compares ten bytes of the alignment's name against ten of the
/// original's, and stops early at a terminator on either side. Reproduced
/// rather than replaced with a slice equality because a name shorter than ten
/// bytes must still match.
fn strncmp_eq(a: &[u8], b: &[u8], n: usize) -> bool {
    for i in 0..n {
        let (x, y) = (a.get(i).copied(), b.get(i).copied());
        match (x, y) {
            (None, None) => return true,
            (Some(p), Some(q)) if p == q => {}
            // One string ended where the other did not, or they differ.
            _ => return false,
        }
    }
    true
}

/// `fillorichar()`: write the original residues onto the alignment's gaps.
///
/// Args:
///     oripos: For each aligned record, which original record it came from.
///     aseq: The aligned rows, edited in place.
///     seq: The ungapped originals.
fn fillorichar(oripos: &[usize], aseq: &mut [Vec<u8>], seq: &[Vec<u8>]) {
    for (i, row) in aseq.iter_mut().enumerate() {
        let src = &seq[oripos[i]];
        let mut k = 0usize;
        for ch in row.iter_mut() {
            if *ch != b'-' {
                // The C reads the original's NUL here and detects it by the
                // character it just wrote being 0. Running out is error (1).
                let Some(&c) = src.get(k) else {
                    eprintln!("ERROR!!");
                    eprintln!("alignment is inconsistent with the original sequences (1)");
                    std::process::exit(1);
                };
                *ch = c;
                k += 1;
            }
        }
        if k != src.len() {
            eprintln!("ERROR!!");
            eprintln!("alignment is inconsistent with the original sequences (2)");
            std::process::exit(1);
        }
    }
}

/// Find `needle` in `hay`.
fn find(hay: &[u8], needle: &[u8]) -> Option<usize> {
    hay.windows(needle.len()).position(|w| w == needle)
}

/// Read a file or die the way the C does.
fn read_or_die(path: &str) -> Vec<u8> {
    match std::fs::read(path) {
        Ok(b) => b,
        Err(_) => {
            eprintln!("Cannot open {path}");
            std::process::exit(1);
        }
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = arguments(&argv);

    let Some(inputfile) = args.inputfile.clone() else {
        // The C reaches `fopen(NULL)` only via stdin; with no -i it reads
        // stdin, and with no -a it refuses.
        let mut buf = Vec::new();
        use std::io::Read;
        let _ = std::io::stdin().read_to_end(&mut buf);
        restore(&args, buf);
        return;
    };
    let bytes = read_or_die(&inputfile);
    restore(&args, bytes);
}

/// The body, once the original's bytes are in hand.
fn restore(args: &Args, orig_bytes: Vec<u8>) {
    let Some(alignmentfile) = args.alignmentfile.clone() else {
        eprintln!("No alignment is given.");
        std::process::exit(1);
    };
    let align_bytes = read_or_die(&alignmentfile);

    let orig = io::read_data(orig_bytes, None, false, true);
    let njobs = orig.njob;
    let mut seq = orig.seq;
    let name = orig.name;

    // The alignment is read with the ordinary case-folding reader, as the C
    // does. Its residues are overwritten by `fillorichar` anyway; what is
    // taken from it is the gap pattern and the names.
    let align = io::read_data(align_bytes, None, false, false);
    let njoba = align.njob;
    let mut aseq = align.seq;
    let aname = align.name;

    // The C runs its record-count check *after* several loops that index the
    // original's arrays with the alignment's count, so a larger alignment
    // reads out of bounds before it ever reaches the check (this repo's issue
    // #10). Checked first here: same message, same exit status, without the
    // out-of-bounds read in between.
    if njoba != njobs {
        eprintln!("ERROR!!");
        eprintln!("In input file,");
        eprintln!("njob = {njobs}");
        eprintln!("but in alignment file,");
        eprintln!("njob = {njoba}");
        std::process::exit(1);
    }

    let mut insname = vec![false; njobs];
    for row in &mut seq {
        gappick_samestring(row);
    }

    if args.keeplength {
        let list = read_or_die("_deletelist");
        for line in String::from_utf8_lossy(&list).lines() {
            let mut it = line.split_whitespace();
            let (Some(i), Some(o), Some(l)) = (it.next(), it.next(), it.next()) else {
                continue;
            };
            let (i, mut o, l) = (
                cargs::myatoi(Some(i)) as usize,
                cargs::myatoi(Some(o)) as usize,
                cargs::myatoi(Some(l)),
            );
            for _ in 0..l {
                if let Some(c) = seq.get_mut(i).and_then(|r| r.get_mut(o)) {
                    *c = b'-';
                }
                o += 1;
            }
            if let Some(f) = insname.get_mut(i) {
                *f = true;
            }
        }
    }

    // Run again: a no-op without `-Y`, and what removes the just-inserted gaps
    // with it. Both passes are the C's.
    for row in &mut seq {
        gappick_samestring(row);
    }

    let mut oripos = vec![0usize; njoba];
    let mut oname: Vec<Vec<u8>> = Vec::with_capacity(njoba);
    for (i, an) in aname.iter().enumerate() {
        let Some(npt0) = find(an, b"_os_") else {
            eprintln!("Format error!");
            std::process::exit(1);
        };
        let o = cargs::myatoi(Some(&String::from_utf8_lossy(&an[npt0 + 4..]))) as usize;
        oripos[i] = o;
        let Some(oe) = find(an, b"_oe_") else {
            eprintln!("Format error!");
            std::process::exit(1);
        };
        // The C shifts the tail left over `_os_<n>_oe_`, in place, then
        // rebuilds the name from the prefix plus that tail -- which is the
        // same string. Written directly.
        let tail = an[oe + 4..].to_vec();
        let mut on = an[..npt0].to_vec();
        on.extend_from_slice(&tail);
        oname.push(on);

        // With `-z`, a record the delete list touched is matched from after
        // the first `|` instead of from the start of its name.
        let compare_from = if args.mapout == 2 && insname.get(o).copied().unwrap_or(false) {
            find(&tail, b"|").map_or(&tail[..], |p| &tail[p + 1..])
        } else {
            &tail[..]
        };
        let Some(original) = name.get(o) else {
            eprintln!("Format error!");
            std::process::exit(1);
        };
        if !strncmp_eq(compare_from, &original[1..], 10) {
            eprintln!("ERROR!!");
            eprintln!("In input file,");
            eprintln!("name[{}] = {}", o, String::from_utf8_lossy(original));
            eprintln!("but in alignment file,");
            eprintln!("oname[{}] = {}", i, String::from_utf8_lossy(&oname[i]));
            std::process::exit(1);
        }
    }

    fillorichar(&oripos, &mut aseq, &seq);

    let mut buf = Vec::new();
    io::write_data(&mut buf, &oname, &aseq);
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let _ = out.write_all(&buf);
    let _ = out.flush();
}
