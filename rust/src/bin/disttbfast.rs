//! `disttbfast-rs` — a bit-exact Rust port of MAFFT's `disttbfast` for the
//! progressive, non-FFT path (NW-NS-1 / NW-NS-2).
//!
//! It is a drop-in replacement for `libexec/mafft/disttbfast`: same argument
//! vector, FASTA on stdin, FASTA on stdout, and the same `order` file written
//! into the working directory.  Configurations whose C code path is not ported
//! are refused rather than approximated — see `Args::check_supported`.
#![forbid(unsafe_code)]

use std::io::{Read, Write};

use mafft_oxide::align::AlignState;
use mafft_oxide::constants::Constants;
use mafft_oxide::lalign::LalignState;
use mafft_oxide::params::Args;
use mafft_oxide::treebase::TreebaseOpts;
use mafft_oxide::{addfunctions, dist, falign, io, keeplength, score, tree, treebase};

/// `pickup()`: keep the longest sequence on each side of the last few tree
/// nodes and write the rest to `notused`.
///
/// The C's comment says it is not adapted to memsave, and it is not: with the
/// memsave `topol` each side holds a single member, so the "longest" search
/// runs over one candidate. Two more quirks come with it -- `m` counts loop
/// iterations while each iteration selects two sequences, and `m` is then used
/// as the record count for stdout, so the file can be shorter than the
/// selection.
fn pickup(
    n: usize,
    sparsepickup: i32,
    seqlen: &[usize],
    tr: &tree::Tree,
    name: &[Vec<u8>],
    seq: &[Vec<u8>],
) {
    let mut longestseq = vec![[-1i32, -1i32]; n.saturating_sub(1)];
    let mut select = vec![false; n];
    for i in 0..n - 1 {
        for side in 0..2 {
            // memsave: one member, then the -1 terminator, so the C's search
            // for the longest compares a single candidate against -1.
            let m = tr.topol[i][side];
            if m != -1 && seqlen[m as usize] as i64 > -1 {
                longestseq[i][side] = m;
            }
        }
    }
    let mut m = 1i32;
    for i in (0..n - 1).rev() {
        for side in 0..2 {
            let p = longestseq[i][side];
            if p >= 0 {
                select[p as usize] = true;
            }
        }
        m += 1;
        if m >= sparsepickup {
            break;
        }
    }

    let mut nameout = Vec::new();
    let mut seqout = Vec::new();
    let mut namenot = Vec::new();
    let mut seqnot = Vec::new();
    for i in 0..n {
        if select[i] {
            nameout.push(name[i].clone());
            seqout.push(seq[i].clone());
        } else {
            namenot.push(name[i].clone());
            seqnot.push(seq[i].clone());
        }
    }
    let mut out = Vec::new();
    let keep = (m as usize).min(nameout.len());
    io::write_data(&mut out, &nameout[..keep], &seqout[..keep]);
    std::io::stdout().write_all(&out).unwrap();
    let mut nu = Vec::new();
    let keepn = (n - m as usize).min(namenot.len());
    io::write_data(&mut nu, &namenot[..keepn], &seqnot[..keepn]);
    std::fs::write("notused", nu).expect("write notused");
}

/// `includemember()`: is every member of `mem` also in `cand`?
fn includemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.iter().all(|m| cand.contains(m))
}

/// `samemembern()`: `mem` has exactly `cand.len()` members and they are the
/// same set.
fn samemembern(mem: &[i32], cand: &[i32]) -> bool {
    mem.len() == cand.len() && mem.iter().all(|m| cand.contains(m))
}

/// `samemember()`: the same members, in any order. (The C has two definitions
/// of this; the one under `#if 0` compares position by position, the live one
/// compares as sets, and `topolorderz` does not hand them out in `addmem`'s
/// order.)
fn samemember(mem: &[i32], cand: &[i32]) -> bool {
    mem.len() == cand.len() && mem.iter().all(|m| cand.contains(m))
}

/// `varpairscore()`: how variable the pairwise scores of a random sample are.
///
/// `--treein auto` uses it to decide between a real tree and a chain: a low
/// coefficient of variation means the sequences are all about equally similar,
/// so the tree would not tell us much and a chain will do. Its `printf`s go to
/// *stdout*, ahead of the alignment, which is the C's behaviour and not a slip
/// here.
fn varpairscore(
    cn: &Constants,
    lst: &mut LalignState,
    nseq: usize,
    npick: usize,
    seq: &[Vec<u8>],
    seed: u32,
    out: &mut Vec<u8>,
) -> u8 {
    use std::io::Write as _;
    let nlenmax = seq.iter().map(|s| s.len()).max().unwrap_or(0);
    eprintln!("nseq = {}, nlenmax={}, seed={}", nseq, nlenmax, seed);

    let mut rng = mafft_oxide::cerand::CRand::new(seed);
    let mut slist: Vec<i32> = (0..nseq as i32).collect();
    mafft_oxide::cerand::stringshuffle(&mut slist, &mut rng);

    // The C indexes `slist` up to `npick` whatever `nseq` is; a larger npick
    // reads off the end of its array, so refuse rather than invent a result.
    if npick > nseq {
        eprintln!(
            "disttbfast-rs: npick={} exceeds the {} input sequences; the C reads past its array here",
            npick, nseq
        );
        std::process::exit(1);
    }
    let pickseq: Vec<Vec<u8>> = (0..npick)
        .map(|i| io::gappick0(&seq[slist[i] as usize]))
        .collect();
    // `pickseq = AllocateCharMtx( npick, nlenmax )` gives each row exactly
    // `nlenmax` bytes, but `gappick0` writes the residues *and* a NUL. A picked
    // sequence with no gaps and the maximal length therefore puts its
    // terminator one byte past its own row, in memory the allocator owns.
    //
    // That is not harmless. The very next `malloc` -- `G__align11_noalign`
    // allocating its own buffers on the first call -- writes a chunk header
    // over that byte, and from then on `strlen( pickseq[i] )` runs off the end
    // of the row into whatever follows. Measured on shuffled_domains.fa: the C
    // reports 280 before the loop and aligns 282 bytes inside it, so its score
    // for that pair is a function of the heap layout. Sometimes it aborts
    // instead, when the stray byte lands on a live chunk header.
    //
    // There is no defined answer to reproduce, so refuse.
    if let Some(i) = pickseq.iter().position(|p| p.len() >= nlenmax) {
        eprintln!(
            "disttbfast-rs: picked sequence {} is {} residues with nlenmax {}, so the C's\n  \
             gappick0 writes its terminator outside the pickseq row and then aligns\n  \
             whatever follows it in the heap (disttbfast.c:537,546). Refusing rather\n  \
             than guessing which bytes those were.",
            slist[i] + 1,
            pickseq[i].len(),
            nlenmax
        );
        std::process::exit(1);
    }

    let mut scoreav = 0.0f64;
    let mut scoreav2 = 0.0f64;
    let npair = npick * (npick - 1) / 2;
    for i in 1..npick {
        for j in 0..i {
            let score = lst.g_align11_noalign(
                cn,
                &cn.n_dis_consweight_multi,
                -1200,
                -60,
                &pickseq[i],
                &pickseq[j],
            );
            scoreav += score;
            scoreav2 += score * score;
            let _ = writeln!(out, "score = {}", score as i32);
        }
    }
    scoreav /= npair as f64;
    scoreav2 /= npair as f64;
    let scorevar = (scoreav2 - scoreav * scoreav) * npair as f64 / (npair - 1) as f64;
    let scorestd = scorevar.sqrt();
    let _ = writeln!(out, "av = {:.6}", scoreav);
    let _ = writeln!(out, "stddev = {:.6}", scorestd);
    let _ = writeln!(out, "cv = {:.6}", scorestd / scoreav);

    if scorestd / scoreav < 0.2 {
        b's'
    } else {
        b't'
    }
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = match Args::parse(&argv) {
        Ok(a) => a,
        Err(e) => {
            eprintln!("{}", e);
            std::process::exit(1);
        }
    };
    if let Err(u) = args.check_supported() {
        eprintln!("disttbfast-rs: unsupported configuration:\n  {}", u.0);
        std::process::exit(1);
    }

    let mut bytes = Vec::new();
    match &args.inputfile {
        Some(path) => {
            let mut f = std::fs::File::open(path).unwrap_or_else(|e| {
                eprintln!("Cannot open {}: {}", path, e);
                std::process::exit(1);
            });
            f.read_to_end(&mut bytes).expect("read input");
        }
        None => {
            std::io::stdin()
                .read_to_end(&mut bytes)
                .expect("read stdin");
        }
    }

    let input = io::read_data(bytes, args.dorp, args.outnumber, args.nblosum == -2);
    let njob = input.njob;
    let mut out = Vec::new();

    if njob < 2 {
        let ungapped: Vec<Vec<u8>> = input.seq.iter().map(|s| io::gappick0(s)).collect();
        io::write_data(&mut out, &input.name, &ungapped);
        std::io::stdout().write_all(&out).unwrap();
        eprintln!("Warning: Only {} sequence found.", njob);
        return;
    }

    // `constants( njob, seq )` sees the raw input, gaps and all; only
    // `--fmodel` reads it.
    let cn = Constants::new_seq(&args, input.dorp, &input.seq);

    // `-l` (`--anchors`): user anchors, read before `constants()` because the
    // C does, and indexed by the ungapped lengths of the input.
    let mut extanch: Option<Vec<mafft_oxide::anchors::ExtAnch>> = None;
    let mut anchindex: Vec<Vec<i32>> = Vec::new();
    if args.useexternalanchors {
        let nogaplen: Vec<usize> = input.seq.iter().map(|s| io::gappick0(s).len()).collect();
        eprintln!("reading anchors");
        let mut a = mafft_oxide::anchors::read_external_anchors(njob, &nogaplen);
        eprintln!("sorting anchors");
        anchindex = mafft_oxide::anchors::index_anchors(&mut a, njob);
        extanch = Some(a);
    }

    // `--merge` (-H): the groups that must stay monophyletic.
    let (subtable, preservegaps) = if args.subalignment {
        let (t, p) = io::read_subalignments_table(njob);
        eprintln!("nsubalignments = {}", t.len());
        (t, p)
    } else {
        (Vec::new(), Vec::new())
    };

    // seqcheck()
    for (i, s) in input.seq.iter().enumerate() {
        for (p, &c) in s.iter().enumerate() {
            if (c as usize) >= cn.charsize || cn.amino_n[c as usize] == -1 {
                eprintln!("Illegal character {}", c as char);
                eprintln!("=== Please check site {} in sequence {}.", p + 1, i + 1);
                std::process::exit(1);
            }
        }
    }

    if let Err(u) = args.check_blosum(cn.dorp) {
        eprintln!("disttbfast-rs: unsupported configuration:\n  {}", u.0);
        std::process::exit(1);
    }

    if cn.dorp == b'p' && args.tuplesize != 6 {
        eprintln!("tuplesize must be 6 for aa sequence");
        std::process::exit(1);
    }
    if cn.dorp == b'd' && args.tuplesize != 6 && args.tuplesize != 10 {
        eprintln!("tuplesize must be 6 or 10 for dna sequence");
        std::process::exit(1);
    }
    // `check_guidetreefile()`: the first line of `_guidetree` selects between a
    // real Newick tree and a family of generated ones.
    let mut treekind = b'l';
    let mut randomseed: u32 = 0;
    let mut npickx: usize = 200;
    // `compact` can clear `treein` after the fact, putting the run back on the
    // ordinary distance-matrix path.
    let mut treein_cleared = false;
    // `compacttree`: 1 = `compact`, 2 = `very compact` / `stepadd`,
    // 4 = `youngestlinkage`. All of them clear `treein`, so the run still goes
    // through "Making a distance matrix" -- but down its compact arm.
    let mut compacttree: i32 = 0;
    let mut maxdistmtxsize: f64 = 10.0 * 1000.0 * 1000.0 * 1000.0;
    if args.treein {
        let first = std::fs::read_to_string("_guidetree")
            .map(|t| {
                // `fgets` keeps the newline, and `strncmp` only looks at the
                // keyword, so a trailing newline is irrelevant.
                t.lines().next().unwrap_or("").to_string()
            })
            .unwrap_or_else(|_| {
                eprintln!("cannot open _guidetree");
                std::process::exit(1);
            });
        if let Some(rest) = first.strip_prefix("shuffle") {
            // `sscanf( string+7, "%d", seed )` -- absent leaves it at 0.
            randomseed = rest
                .split_whitespace()
                .next()
                .and_then(|t| t.parse().ok())
                .unwrap_or(0);
            eprintln!("shuffle, seed={}", randomseed);
            treekind = b's';
        } else if first.starts_with("pileup") {
            eprintln!("pileup.");
            treekind = b'p';
        } else if first.starts_with("nodepair") {
            // The C accepts the keyword and then falls off the end of its own
            // dispatch, so this is its behaviour, not a refusal of ours.
            eprintln!("Use nodepair.");
            eprintln!("Error. treein = {} or n", b'n');
            std::process::exit(1);
        } else if first.starts_with("auto") || first.starts_with("test") {
            let test = first.starts_with("test");
            let mut it = first[4..].split_whitespace();
            randomseed = it.next().and_then(|t| t.parse().ok()).unwrap_or(0);
            npickx = it.next().and_then(|t| t.parse().ok()).unwrap_or(200);
            if test {
                eprintln!("calc, seed={}, npick={}", randomseed, npickx);
            } else {
                eprintln!("auto, seed={}, npick={}", randomseed, npickx);
            }
            if npickx < 2 {
                eprintln!("Check npick");
                std::process::exit(1);
            }
            treekind = if test { b't' } else { b'a' };
        } else if first.starts_with("compact") {
            // `--initialramusage`: a size with a unit suffix. The C's first
            // test is `strchr(s,'k') || strchr(s,'k')` -- the same character
            // twice -- so an upper-case K is *not* a kilobyte here and falls
            // through to the error, while a lower-case k works. Kept.
            let sizestring = &first[7..];
            let has = |c: char| sizestring.contains(c);
            let tanni = if has('k') || has('k') {
                1e3
            } else if has('M') || has('m') {
                1e6
            } else if has('G') || has('g') {
                1e9
            } else if has('T') || has('t') {
                1e12
            } else {
                eprintln!("\nSpecify initial ram usage by '--initialramusage xGB'\n\n");
                std::process::exit(1);
            };
            let tmpd: f64 = {
                let t = sizestring.trim_start();
                let end = t
                    .find(|c: char| !(c.is_ascii_digit() || c == '.' || c == '-' || c == '+'))
                    .unwrap_or(t.len());
                t[..end].parse().unwrap_or(0.0)
            };
            maxdistmtxsize = tmpd * tanni;
            eprintln!("sizestring = {}", sizestring);
            eprintln!("Initial RAM usage = {:10.3}GB", maxdistmtxsize / 1e9);
            // `compacttree = 1; treein = 0;`, and then, once njob is known:
            // if the full matrix would fit in the stated budget the C drops
            // straight back to the ordinary builder.
            compacttree = 1;
            treekind = b'l';
            treein_cleared = true;
            if maxdistmtxsize > njob as f64 * (njob - 1) as f64 * 8.0 / 2.0 {
                eprintln!("Use conventional tree.");
                compacttree = 0;
            }
        } else if first.starts_with("very compact") {
            eprintln!("very compact.");
            compacttree = 2;
            treekind = b'l';
            treein_cleared = true;
        } else if first.starts_with("stepadd") {
            eprintln!("stepwise addition (disttbfast).");
            compacttree = 2;
            treekind = b'l';
            treein_cleared = true;
        } else if first.starts_with("youngestlinkage") {
            eprintln!("youngest linkage (disttbfast).");
            compacttree = 4;
            treekind = b'l';
            treein_cleared = true;
        } else {
            eprintln!("loadtree.");
        }
    }
    // `test` and `auto` both sample pairwise scores first. `test` then stops;
    // `auto` picks a chain or a real tree from the answer -- and either way
    // the sampled scores have already gone to stdout, ahead of everything.
    let mut treein = args.treein && !treein_cleared;
    let mut nguidetree_override: Option<usize> = None;
    if treekind == b't' || treekind == b'a' {
        let mut lst = LalignState::new();
        if treekind == b't' {
            varpairscore(
                &cn, &mut lst, njob, npickx, &input.seq, randomseed, &mut out,
            );
            std::io::stdout().write_all(&out).unwrap();
            std::process::exit(1);
        }
        if njob > 200000 {
            eprintln!("Chain?");
            treekind = b's';
            nguidetree_override = Some(1);
        } else if njob < 100
            || varpairscore(
                &cn, &mut lst, njob, npickx, &input.seq, randomseed, &mut out,
            ) == b't'
        {
            eprintln!("Tree!");
            treekind = b'l';
            treein = false;
            nguidetree_override = Some(2);
        } else {
            eprintln!("Chain!");
            treekind = b's';
            nguidetree_override = Some(1);
        }
    }

    // `shuffle` and `pileup` do not use the sequence weights at all: the C
    // clears `weight` and `tbrweight` at the call site.
    let tbrweight = if treekind == b's' || treekind == b'p' {
        0
    } else {
        args.tbrweight
    };

    // With a user tree there is no distance matrix at all; with a compact one
    // there is no *full* matrix, only nearest neighbours.
    let mut compact_mindist: Vec<f64> = Vec::new();
    let mut compact_mindistfrom: Vec<i32> = Vec::new();
    // `youngestlinkage` recomputes distances at every step, so it needs the
    // k-mer tables (first pass) or the MSA (second) kept around, not just the
    // nearest-neighbour summary the other keywords build from.
    let mut compact_kmer_selfscore: Vec<i32> = Vec::new();
    let mut compact_tuples: Option<mafft_oxide::dist::Tuples> = None;
    let mut compact_msa: Option<(Vec<Vec<u8>>, Vec<Vec<i32>>, Vec<i32>)> = None;
    let (mut mtx, nogaplen): (Vec<Vec<f64>>, Vec<usize>) = if treein {
        (Vec::new(), input.seq.iter().map(|s| s.len()).collect())
    } else if compacttree != 0 {
        eprintln!("\n\nMaking a distance matrix ..");
        let (md, mdf, ngl, _partmtx, kss) = mafft_oxide::compacttree::compact_distances(
            &cn,
            &input.seq,
            njob,
            args.tuplesize as usize,
            compacttree,
            maxdistmtxsize,
        );
        compact_mindist = md;
        compact_mindistfrom = mdf;
        compact_kmer_selfscore = kss;
        compact_tuples = Some(mafft_oxide::dist::build_tuples(
            &cn,
            &input.seq,
            njob,
            args.tuplesize as usize,
        ));
        (Vec::new(), ngl)
    } else {
        dist::distance_matrix(
            &cn,
            &input.seq,
            njob,
            args.tuplesize as usize,
            mafft_oxide::parallel::workers(args.nthreadpair),
        )
    };

    let mut nlen: Vec<i32> = input.seq.iter().map(|s| s.len() as i32).collect();
    // `-r` runs after the guide-tree loop, on whatever the last iteration
    // left, so the tree and the weights outlive it.
    let mut lasttree: Option<tree::Tree> = None;
    let mut lasteff: Vec<f64> = Vec::new();
    let mut st = AlignState::new();
    let mut fs = falign::FalignState::new();
    let mut bseq: Vec<Vec<u8>> = Vec::new();
    let mut originalgaps: Option<Vec<u8>> = None;
    let mut addbk: Option<Vec<Vec<u8>>> = None;
    let mut mapseq: Option<Vec<Vec<u8>>> = None;

    let nguidetree = match nguidetree_override {
        Some(n) => n,
        None if treein => 1,
        None => args.nguidetree.max(1) as usize,
    };
    for iguidetree in 0..nguidetree {
        let calcpairdists = iguidetree != nguidetree - 1;

        // `-y`: emit the k-tuple distance matrix. This happens before the tree
        // is built, and with `-T` (noalign) and no `-t` the run stops here --
        // before `order` is written.
        if args.distout != 0 && !args.treein {
            io::write_hat2(njob, &input.name, &mtx);
            if args.treeout == 0 && args.noalign {
                io::write_data(&mut out, &input.name, &input.seq);
                std::io::stdout().write_all(&out).unwrap();
                eprintln!();
                return;
            }
        }

        let mut subfull: Vec<[Vec<i32>; 2]> = Vec::new();
        let mut tr = if treekind == b's' || treekind == b'p' {
            tree::createchain(
                njob,
                if args.treeout != 0 {
                    Some(&input.name[..])
                } else {
                    None
                },
                args.outnumber,
                treekind == b's',
                randomseed,
            )
        } else if treein {
            match tree::load_tree(
                njob,
                if args.treeout != 0 {
                    Some(&input.name[..])
                } else {
                    None
                },
                args.outnumber,
            ) {
                Ok(t) => t,
                Err(e) => {
                    eprintln!("{}", e);
                    std::process::exit(1);
                }
            }
        } else if compacttree != 0 {
            eprintln!("Constructing a tree ... ");
            let mut return_tree_ct: Option<tree::Tree> = None;
            if compacttree == 4 {
                // `compacttree_memsaveselectable( ..., howcompact 2, memsave 1 )`.
                // The first pass scores k-mers, the second the alignment the
                // first produced -- the C picks on whether `seq` is non-NULL,
                // having freed `pointt` after the first pass and `bseq` after
                // the second.
                let src = match compact_msa.as_ref() {
                    Some((seq, skiptable, _)) => mafft_oxide::compacttree::CompactDist::Msa {
                        cn: &cn,
                        seq,
                        skiptable,
                    },
                    None => {
                        let t = compact_tuples.as_ref().unwrap();
                        mafft_oxide::compacttree::CompactDist::Kmer {
                            pointt: &t.pointt,
                            nlen: &t.nogaplen,
                            tsize: t.tsize,
                            lenfaca: t.lenfaca,
                            lenfacb: t.lenfacb,
                            lenfacc: t.lenfacc,
                            lenfacd: t.lenfacd,
                        }
                    }
                };
                let ss: &[i32] = match compact_msa.as_ref() {
                    Some((_, _, ss)) => ss,
                    None => &compact_kmer_selfscore,
                };
                return_tree_ct = Some(mafft_oxide::compacttree::compacttree_memsaveselectable(
                    njob,
                    &mut compact_mindistfrom,
                    &mut compact_mindist,
                    ss,
                    &src,
                    args.sueff_global,
                    if args.treeout != 0 {
                        Some(&input.name[..])
                    } else {
                        None
                    },
                    args.outnumber,
                    args.treeout != 0,
                ));
            }
            match return_tree_ct.take() {
                Some(t) => t,
                None => mafft_oxide::compacttree::compacttreegivendist(
                    njob,
                    &compact_mindist,
                    &compact_mindistfrom,
                    if args.treeout != 0 {
                        Some(&input.name[..])
                    } else {
                        None
                    },
                    args.outnumber,
                    args.treeout != 0,
                ),
            }
        } else if args.subalignment {
            eprintln!("Constructing a UPGMA tree ... ");
            let (t, f) = tree::upgma_treeout_constrained(
                njob,
                &mut mtx,
                args.sueff_global,
                &input.name,
                args.outnumber,
                &subtable,
            );
            subfull = f;
            t
        } else if args.treeout != 0 {
            tree::upgma_treeout(
                njob,
                &mut mtx,
                args.sueff_global,
                Some(&input.name),
                args.outnumber,
                args.treeout == 2,
            )
        } else {
            tree::upgma(njob, &mut mtx, args.sueff_global)
        };

        // `-q N`: keep a sparse subset and stop, before `order` is written.
        if args.sparsepickup != 0 && iguidetree == nguidetree - 1 {
            eprintln!("Sparsepickup! ");
            pickup(
                njob,
                args.sparsepickup,
                &nogaplen,
                &tr,
                &input.name,
                &input.seq,
            );
            eprintln!("done. ");
            return;
        }

        // The `order` file the mafft wrapper uses for --reorder / clustal out.
        {
            // For a chain the members are read straight off it, since
            // `topolorderz` would walk the comb one node at a time.
            let order: Vec<i32> = if treekind == b's' || treekind == b'p' {
                let mut o = vec![tr.topol[0][0]];
                for i in 1..njob {
                    o.push(tr.topol[i - 1][1]);
                }
                o
            } else {
                tree::topolorderz(&tr.topol, &tr.dep, njob - 2, 2)
            };
            let mut s = String::new();
            for v in order.iter().take(njob) {
                s.push_str(&format!("{}\n", v));
            }
            if let Err(e) = std::fs::write("order", s) {
                eprintln!("Cannot open 'order': {}", e);
                std::process::exit(1);
            }
        }

        if (args.treeout != 0 || args.distout != 0) && args.noalign {
            io::write_data(&mut out, &input.name, &input.seq);
            std::io::stdout().write_all(&out).unwrap();
            eprintln!();
            return;
        }

        let effarr: Vec<f64> = if tbrweight != 0 {
            tree::counteff_simple(njob, &mut tr)
        } else {
            vec![1.0; njob]
        };

        // `--add`: the first `njob - nadd` sequences are an existing alignment
        // and keep their columns; only the added ones are ungapped. Everything
        // else takes the plain path.
        let mut mergeoralign: Vec<u8> = vec![b'a'; njob.saturating_sub(1)];
        if args.nadd != 0 {
            let nadd = args.nadd as usize;
            let alignmentlength = input.seq[0].len();
            for i in 0..njob - nadd {
                if input.seq[i].len() != alignmentlength {
                    eprintln!("# The original {} sequences must be aligned", njob - nadd);
                    std::process::exit(1);
                }
            }
            let addmem: Vec<i32> = ((njob - nadd)..njob).map(|i| i as i32).collect();
            if args.addprofile {
                let alignmentlength = input.seq[njob - nadd].len();
                for i in (njob - nadd)..njob {
                    if input.seq[i].len() != alignmentlength {
                        eprintln!("# The {} additional sequences must be aligned", nadd);
                        std::process::exit(1);
                    }
                }
                let mut found = false;
                for i in 0..njob - 1 {
                    let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                    let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                    if samemember(&m0, &addmem) {
                        mergeoralign[i] = b'1';
                        found = true;
                    } else if samemember(&m1, &addmem) {
                        mergeoralign[i] = b'2';
                        found = true;
                    } else {
                        mergeoralign[i] = b'n';
                    }
                }
                if !found {
                    eprintln!("Not found the accurate branch.");
                    std::process::exit(1);
                }
            } else {
                // Plain `--add` starts every node at 'n' -- nodes inside the
                // existing alignment are already done -- and then two passes
                // mark the ones that touch an added sequence.
                for v in mergeoralign.iter_mut() {
                    *v = b'n';
                }
                for j in (njob - nadd)..njob {
                    let one = [j as i32];
                    for i in 0..njob - 1 {
                        let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                        let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                        if samemembern(&m0, &one) {
                            mergeoralign[i] = if mergeoralign[i] != b'n' { b'w' } else { b'1' };
                        } else if samemembern(&m1, &one) {
                            mergeoralign[i] = if mergeoralign[i] != b'n' { b'w' } else { b'2' };
                        }
                    }
                }
                for i in 0..njob - 1 {
                    let m0 = tree::topolorderz(&tr.topol, &tr.dep, i, 0);
                    let m1 = tree::topolorderz(&tr.topol, &tr.dep, i, 1);
                    let in0 = includemember(&m0, &addmem);
                    let in1 = includemember(&m1, &addmem);
                    if in0 && in1 {
                        mergeoralign[i] = b'w';
                    } else if in0 {
                        mergeoralign[i] = b'1';
                    } else if in1 {
                        mergeoralign[i] = b'2';
                    }
                }
            }
        }

        // `--keeplength`: remember the reference's all-gap columns before the
        // alignment can add any, and (for --mapout) keep an ungapped copy of
        // each added sequence to index the map by. Both read the input as it
        // was read, ahead of every commongappick below.
        if args.nadd != 0 && args.keeplength {
            let nadd = args.nadd as usize;
            originalgaps = Some(keeplength::record_original_gaps(njob - nadd, &input.seq));
            if args.mapout != 0 {
                addbk = Some(
                    input.seq[njob - nadd..]
                        .iter()
                        .map(|s| io::gappick0(s))
                        .collect(),
                );
            }
        }

        // `--merge`: every sub-alignment must already be aligned, and must
        // come out monophyletic in the tree; the nodes entirely inside one are
        // skipped, since it is kept as it is.
        let mut insubtable = vec![false; njob];
        if args.subalignment {
            for v in mergeoralign.iter_mut() {
                *v = b'a';
            }
            for (i, g) in subtable.iter().enumerate() {
                eprintln!("Checking subalignment {}:", i + 1);
                let alignmentlength = input.seq[g[0] as usize].len();
                for &m in g.iter() {
                    if m as usize >= njob {
                        eprintln!("No such sequence, {}.", m + 1);
                        std::process::exit(1);
                    }
                    if input.seq[m as usize].len() != alignmentlength {
                        eprintln!("# Subalignment {} must be aligned.", i + 1);
                        eprintln!(
                            "# {}. {:.10} -> {} letters (including gaps)",
                            g[0] + 1,
                            String::from_utf8_lossy(&input.name[g[0] as usize][1..]),
                            alignmentlength
                        );
                        eprintln!(
                            "# {}. {:.10} -> {} letters (including gaps)",
                            m + 1,
                            String::from_utf8_lossy(&input.name[m as usize][1..]),
                            input.seq[m as usize].len()
                        );
                        std::process::exit(1);
                    }
                    insubtable[m as usize] = true;
                }
                for j in 0..njob - 1 {
                    if includemember(&subfull[j][0], g) && includemember(&subfull[j][1], g) {
                        mergeoralign[j] = b'n';
                    }
                }
                let found = (0..njob - 1)
                    .any(|j| samemember(&subfull[j][0], g) || samemember(&subfull[j][1], g));
                if !found {
                    eprintln!(
                        "# Subalignment {} does not seem to form a monophyletic cluster",
                        i + 1
                    );
                    eprintln!(
                        "# in the guide tree ('GuideTree' in this directory) internally computed."
                    );
                    let _ = std::fs::copy("infile.tree", "GuideTree");
                    std::process::exit(1);
                }
                eprintln!(" -> OK");
            }
        }

        bseq = input.seq.iter().map(|s| io::gappick0(s)).collect();
        if args.subalignment {
            // Members of a sub-alignment keep their gaps; everyone else is
            // ungapped as usual.
            for (i, s) in input.seq.iter().enumerate() {
                if insubtable[i] {
                    bseq[i] = s.clone();
                }
            }
            for (i, g) in subtable.iter().enumerate() {
                if !preservegaps[i] {
                    let mut block: Vec<Vec<u8>> =
                        g.iter().map(|&m| bseq[m as usize].clone()).collect();
                    addfunctions::commongappick(&mut block);
                    for (j, s) in block.into_iter().enumerate() {
                        bseq[g[j] as usize] = s;
                    }
                }
            }
        }
        if args.nadd != 0 {
            let nadd = args.nadd as usize;
            // `--addprofile` keeps the added block's own columns: the C
            // collapses only its all-gap ones and copies it across verbatim,
            // where plain `--add` ungaps each added sequence (done above).
            if args.addprofile {
                let mut tail: Vec<Vec<u8>> = input.seq[njob - nadd..].to_vec();
                addfunctions::commongappick(&mut tail);
                // The C collapses these columns in `seq` itself, and the
                // compact --mapout indexes off `seq`, so it sees them gone.
                mapseq = Some(tail.clone());
                for (i, s) in tail.into_iter().enumerate() {
                    bseq[njob - nadd + i] = s;
                }
            }
            let mut head: Vec<Vec<u8>> = input.seq[..njob - nadd].to_vec();
            addfunctions::commongappick(&mut head);
            for (i, s) in head.into_iter().enumerate() {
                bseq[i] = s;
            }
        }

        let selfscore: Vec<i32> = if calcpairdists {
            input
                .seq
                .iter()
                .map(|s| score::naivepairscore11(&cn, s, s, cn.penalty_dist) as i32)
                .collect()
        } else {
            Vec::new()
        };

        eprintln!(
            "Progressive alignment {}/{}... ",
            iguidetree + 1,
            nguidetree
        );
        let opts = TreebaseOpts {
            extanch: extanch
                .as_ref()
                .map(|e| (e.as_slice(), anchindex.as_slice())),
            terminalmargin: args.terminalmargin,
            mergeoralign: if args.nadd != 0 || args.subalignment {
                Some(&mergeoralign)
            } else {
                None
            },
            nadd: args.nadd as usize,
            outgap: args.outgap,
            use_fft: args.use_fft,
            force_fft: args.force_fft,
            nevermemsave: args.nevermemsave,
            alg: args.alg,
            minimumweight: 0.0,
            constraint: None,
            fastathreshold: 0.0,
            use_dynamic_mtx: cn.specificityconsideration != 0.0,
        };
        treebase::treebase(
            &cn, &opts, njob, &mut bseq, &mut nlen, &tr, &effarr, &mut st, &mut fs,
        );

        lasteff = effarr;
        lasttree = Some(tr);

        if calcpairdists && compacttree != 0 {
            let (md, mdf) = mafft_oxide::compacttree::msa_compact_distances(
                &cn,
                &bseq,
                &selfscore,
                njob,
                compacttree,
                maxdistmtxsize,
            );
            compact_mindist = md;
            compact_mindistfrom = mdf;
            if compacttree == 4 {
                let skiptable: Vec<Vec<i32>> =
                    bseq.iter().map(|s| score::makeskiptable(s)).collect();
                compact_msa = Some((bseq.clone(), skiptable, selfscore.clone()));
                compact_tuples = None; // the C frees `pointt` after the first pass
            }
        } else if calcpairdists {
            let skiptable: Vec<Vec<i32>> = bseq.iter().map(|s| score::makeskiptable(s)).collect();
            // The O(n^2) half of an `--auto` build at these sizes: 7.3 million
            // pairs at 3831 sequences. Each cell reads only `bseq`, `skiptable`
            // and `selfscore`, none of which this loop writes, so the matrix is
            // the same for any worker count -- see `mafft_oxide::parallel`.
            mafft_oxide::parallel::fill_rows(
                &mut mtx,
                mafft_oxide::parallel::workers(args.nthreadpair),
                || (),
                |_, i, row| {
                    if i + 1 >= njob {
                        return;
                    }
                    let ssi = selfscore[i] as f64;
                    for j in i + 1..njob {
                        let ssj = selfscore[j] as f64;
                        let bunbo = if ssi < ssj { ssi } else { ssj };
                        let mut v = if bunbo == 0.0 {
                            2.0
                        } else {
                            let s = score::naivepairscorefast(
                                &cn,
                                &bseq[i],
                                &bseq[j],
                                &skiptable[i],
                                &skiptable[j],
                                cn.penalty_dist,
                            );
                            let mut t = (1.0 - s / bunbo) * 2.0;
                            if t > 10.0 {
                                t = 10.0;
                            }
                            t
                        };
                        if v < 0.0 {
                            eprintln!("WARNING: negative distance, iscoretmp = {}", v);
                            v = 0.0;
                        }
                        row[j - i] = v;
                    }
                },
            );
        }
    }

    if args.oneiteration {
        eprintln!("Iterative refinement (one vs others)");
        let opts = TreebaseOpts {
            // `-r` re-aligns without anchors, as the C's dooneiteration does
            // only when `useexternalanchors` is off; with anchors it uses them.
            extanch: extanch
                .as_ref()
                .map(|e| (e.as_slice(), anchindex.as_slice())),
            terminalmargin: args.terminalmargin,
            mergeoralign: None,
            nadd: 0,
            outgap: args.outgap,
            use_fft: args.use_fft,
            force_fft: args.force_fft,
            nevermemsave: args.nevermemsave,
            alg: args.alg,
            minimumweight: 0.0,
            constraint: None,
            fastathreshold: 0.0,
            use_dynamic_mtx: cn.specificityconsideration != 0.0,
        };
        treebase::dooneiteration(
            &cn,
            &opts,
            njob,
            &mut bseq,
            &mut nlen,
            lasttree.as_ref().unwrap(),
            &lasteff,
            &mut st,
            &mut fs,
        );
    }

    // `--keeplength`: cut the columns the alignment added to the reference
    // back out, record what that cost in `_deletelist` (and `_deletemap`),
    // then put the reference's own all-gap columns back.
    let mut names = input.name.clone();
    let mut ndeleted = 0usize;
    if args.keeplength {
        let nadd = args.nadd as usize;
        let (nd, deletelist) =
            keeplength::delete_new_insertions_whole(njob - nadd, nadd, &mut bseq);
        ndeleted = nd;
        let mut dlf = String::new();
        for (i, list) in deletelist.iter().enumerate() {
            for e in list.iter() {
                // 0-origin, as the C writes it.
                dlf.push_str(&format!("{} {} {}\n", njob - nadd + i, e.pos, e.len));
            }
        }
        std::fs::write("_deletelist", dlf).expect("Cannot open _deletelist");

        keeplength::restore_original_gaps(njob, &mut bseq, originalgaps.as_ref().unwrap());

        if args.mapout != 0 {
            let addbk = addbk.as_ref().unwrap();
            // The C indexes the residue-by-residue map by the *aligned*
            // sequences and the compact one by the input's own, so the two
            // read different arrays here.
            let map = if args.mapout == 1 {
                keeplength::reconstruct_delete_map(
                    addbk,
                    &deletelist,
                    &bseq[njob - nadd..],
                    &input.name[njob - nadd..],
                    args.outnumber,
                )
            } else {
                let mut tail: Vec<Vec<u8>> = names[njob - nadd..].to_vec();
                let realn: &[Vec<u8>] = match mapseq.as_ref() {
                    Some(v) => v,
                    None => &input.seq[njob - nadd..],
                };
                let m = keeplength::reconstruct_delete_map_compact(
                    addbk,
                    &deletelist,
                    realn,
                    &mut tail,
                    args.outnumber,
                );
                names[njob - nadd..].clone_from_slice(&tail);
                m
            };
            std::fs::write("_deletemap", map).expect("Cannot open _deletemap");
        }
    }

    io::write_data(&mut out, &names, &bseq);
    if ndeleted > 0 {
        eprintln!(
            "\nTo keep the alignment length, {} letters were DELETED.",
            ndeleted
        );
        if args.mapout != 0 {
            eprintln!("The deleted letters are shown in the (filename).map file.");
        } else {
            eprintln!("To know the positions of deleted letters, rerun the same command with the --mapout option.");
        }
    }
    std::io::stdout().write_all(&out).unwrap();
}
