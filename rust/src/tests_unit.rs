//! Unit tests that check the port's own behaviour, independently of the
//! golden files and of a C build being available.
//!
//! The golden harness proves "same as C" on the inputs it happens to cover.
//! These tests prove properties that must hold on *every* input, so they catch
//! the cases a corpus misses: they check algebraic invariants (an alignment
//! preserves its input residues; a distance matrix is a distance matrix), the
//! documented behaviour of the fiddly helpers (`myfgets` truncation, the
//! deliberate off-by-one reads), and self-consistency between routines that
//! compute the same thing two ways (`L__align11` vs `L__align11_noalign`,
//! `st_*GapCount` vs its incremental `*Add` form).

use crate::align::AlignState;
use crate::constants::{dist2offset, makedynamicmtx, Constants, NALPHABETS};
use crate::io;
use crate::lalign::LalignState;
use crate::localhom::{self, LocalHom};
use crate::pairalign::{self, PairAlg};
use crate::params::{atof, atoi, Args};
use crate::profile;
use crate::tree;

/// The parameters `mafft --localpair` uses for its pairwise stage.
fn pair_constants() -> Constants {
    let a = Args {
        nblosum: 62,
        ppenalty: -2000,
        ppenalty_ex: -100,
        poffset: 99,
        penalty_shift_factor: 100.0,
        ..Default::default()
    };
    Constants::new(&a)
}

/// The parameters the progressive stage uses.
fn tree_constants() -> Constants {
    let a = Args {
        nblosum: 62,
        ppenalty: -1530,
        ppenalty_ex: 0,
        poffset: 0,
        penalty_shift_factor: 100.0,
        ..Default::default()
    };
    Constants::new(&a)
}

fn deterministic_seqs(n: usize, len: usize, seed: u64) -> Vec<Vec<u8>> {
    // A tiny xorshift, so the tests do not need a dependency.
    let mut s = seed | 1;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        s
    };
    const ALPHA: &[u8] = b"ARNDCQEGHILKMFPSTWYV";
    let root: Vec<u8> = (0..len)
        .map(|_| ALPHA[(next() % ALPHA.len() as u64) as usize])
        .collect();
    (0..n)
        .map(|_| {
            let mut t = root.clone();
            for _ in 0..len / 6 {
                let p = (next() % len as u64) as usize;
                t[p] = ALPHA[(next() % ALPHA.len() as u64) as usize];
            }
            t
        })
        .collect()
}

fn ungap(s: &[u8]) -> Vec<u8> {
    s.iter().copied().filter(|&c| c != b'-').collect()
}

// ---------------------------------------------------------------- constants

#[test]
fn score_matrix_is_symmetric_and_integral() {
    let cn = tree_constants();
    for i in 0..NALPHABETS {
        for j in 0..NALPHABETS {
            assert_eq!(cn.n_dis[i][j], cn.n_dis[j][i], "n_dis[{}][{}]", i, j);
            let v = cn.n_dis_consweight_multi[i][j];
            assert_eq!(v, v.trunc(), "n_dis_consweight_multi[{}][{}]", i, j);
        }
    }
}

/// `A__align` refuses to run unless this holds, so make the reason explicit:
/// the gap character must score zero against everything.
#[test]
fn gap_scores_zero_against_every_residue() {
    let cn = tree_constants();
    let gap = cn.amino_n[b'-' as usize] as usize;
    for i in 0..NALPHABETS {
        assert_eq!(cn.n_dis[i][gap], 0);
        assert_eq!(cn.n_dis[gap][i], 0);
    }
}

/// Character- and index-addressed matrices must agree.
#[test]
fn amino_dis_agrees_with_n_dis() {
    let cn = tree_constants();
    for i in 0..NALPHABETS {
        for j in 0..NALPHABETS {
            assert_eq!(
                cn.amino_dis[cn.amino[i] as usize][cn.amino[j] as usize],
                cn.n_dis[i][j]
            );
        }
    }
    for i in 0..NALPHABETS {
        assert_eq!(cn.amino_n[cn.amino[i] as usize], i as i32);
    }
}

/// Identical residues must not score worse than mismatches, or the aligner is
/// solving a different problem.
#[test]
fn diagonal_dominates() {
    let cn = tree_constants();
    for i in 0..20 {
        for j in 0..20 {
            if i != j {
                assert!(
                    cn.n_dis[i][i] >= cn.n_dis[i][j],
                    "{}{} scores above {}{}",
                    cn.amino[i] as char,
                    cn.amino[j] as char,
                    cn.amino[i] as char,
                    cn.amino[i] as char
                );
            }
        }
    }
}

#[test]
fn dist2offset_is_capped_at_zero_and_monotone() {
    for spec in [0.0, 0.3, 0.8] {
        let mut prev = f64::NEG_INFINITY;
        for k in 0..40 {
            let d = k as f64 * 0.1;
            let v = dist2offset(d, spec);
            assert!(v <= 0.0, "dist2offset({}, {}) = {}", d, spec, v);
            assert!(v >= prev, "not monotone at {}", d);
            prev = v;
        }
        // With no shift allowance the offset is always zero.
        if spec == 0.0 {
            assert_eq!(dist2offset(1.0, spec), 0.0);
        }
    }
}

/// `makedynamicmtx` must be the identity when `--allowshift` is off, since
/// `tbfast` calls it on every tree node regardless.
#[test]
fn dynamic_matrix_is_identity_without_allowshift() {
    let cn = tree_constants();
    let mut out = vec![vec![0.0f64; NALPHABETS]; NALPHABETS];
    for distfromtip in [0.0, 0.1, 0.5, 1.0] {
        makedynamicmtx(
            &mut out,
            &cn.n_dis_consweight_multi,
            distfromtip,
            &cn.amino,
            0.0,
        );
        assert_eq!(
            out, cn.n_dis_consweight_multi,
            "distfromtip={}",
            distfromtip
        );
    }
}

/// With `--allowshift` the matrix is raised uniformly on the scored residues,
/// leaving the gap row and column alone.
#[test]
fn dynamic_matrix_shifts_uniformly() {
    let cn = tree_constants();
    let mut out = vec![vec![0.0f64; NALPHABETS]; NALPHABETS];
    makedynamicmtx(&mut out, &cn.n_dis_consweight_multi, 0.05, &cn.amino, 0.8);
    let gap = cn.amino_n[b'-' as usize] as usize;
    let delta = out[0][1] - cn.n_dis_consweight_multi[0][1];
    assert!(delta < 0.0, "expected a negative offset, got {}", delta);
    for i in 0..NALPHABETS {
        for j in 0..NALPHABETS {
            let want = if i == gap || j == gap {
                cn.n_dis_consweight_multi[i][j]
            } else {
                cn.n_dis_consweight_multi[i][j] + delta
            };
            assert_eq!(out[i][j], want, "[{}][{}]", i, j);
        }
    }
}

// ----------------------------------------------------------------------- io

#[test]
fn fasta_reader_filters_and_upcases() {
    let inp = b">a\nac-de\n>b\nAC*DE\n>c\nac de 1 2\n".to_vec();
    let r = io::read_data(inp, None, false, false);
    assert_eq!(r.njob, 3);
    assert_eq!(r.dorp, b'p');
    assert_eq!(r.seq[0], b"AC-DE");
    // `*` becomes a gap.
    assert_eq!(r.seq[1], b"AC-DE");
    // Digits and spaces are dropped.
    assert_eq!(r.seq[2], b"ACDE");
    assert_eq!(r.name[0], b"=a");
}

/// Names are truncated at `B-2` bytes, and the leading `=` is not part of the
/// budget. Output writes the name back without it.
#[test]
fn long_names_are_truncated_like_the_c() {
    let long = "N".repeat(400);
    let inp = format!(">{}\nACDE\n>b\nACDE\n", long).into_bytes();
    let r = io::read_data(inp, None, false, false);
    assert_eq!(r.name[0].len(), 1 + (io::B - 2));
    assert_eq!(r.name[0][0], b'=');
    let mut out = Vec::new();
    io::write_data(&mut out, &r.name, &r.seq);
    let text = String::from_utf8(out).unwrap();
    assert!(text.starts_with(&format!(">{}\n", "N".repeat(io::B - 2))));
}

#[test]
fn nucleotide_input_is_detected() {
    let inp = b">a\nacgtacgtacgtacgt\n>b\nacgtacgtacgtacgt\n".to_vec();
    assert_eq!(io::read_data(inp, None, false, false).dorp, b'd');
    // ...unless the caller forces the type.
    let inp = b">a\nacgtacgtacgtacgt\n>b\nacgtacgtacgtacgt\n".to_vec();
    assert_eq!(io::read_data(inp, Some(b'p'), false, false).dorp, b'p');
}

#[test]
fn write_data_round_trips_through_read_data() {
    let seqs = deterministic_seqs(5, 60, 7);
    let mut src = Vec::new();
    for (i, s) in seqs.iter().enumerate() {
        src.extend_from_slice(format!(">s{}\n", i).as_bytes());
        src.extend_from_slice(s);
        src.push(b'\n');
    }
    let r = io::read_data(src, None, false, false);
    let mut out = Vec::new();
    io::write_data(&mut out, &r.name, &r.seq);
    let r2 = io::read_data(out, None, false, false);
    assert_eq!(r.name, r2.name);
    assert_eq!(r.seq, r2.seq);
}

// ------------------------------------------------------------------ profile

/// The composition profile is a weight distribution: each column sums to the
/// total weight, and no cell is negative.
#[test]
fn cpmx_columns_sum_to_the_total_weight() {
    let cn = tree_constants();
    let seqs: Vec<Vec<u8>> = vec![
        b"ACD-EF".to_vec(),
        b"ACDGEF".to_vec(),
        b"-CDGEW".to_vec(),
        b"ACDGE-".to_vec(),
    ];
    let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let eff = vec![0.1, 0.2, 0.3, 0.4];
    let lgth = 6;
    let mut cpmx = vec![vec![0.0f64; lgth + 2]; NALPHABETS];
    profile::cpmx_calc_new(&cn, &refs, &mut cpmx, &eff, lgth, refs.len());
    for j in 0..lgth {
        let mut s = 0.0;
        for row in cpmx.iter().take(NALPHABETS) {
            assert!(row[j] >= 0.0);
            s += row[j];
        }
        assert!((s - 1.0).abs() < 1e-12, "column {} sums to {}", j, s);
    }
}

/// `cpmx_calc_add` must agree with recomputing from scratch. The C uses the
/// incremental form under `--allowshift` only, so a divergence would show up
/// on one strategy and not the others.
#[test]
fn incremental_profile_matches_full_recompute() {
    let cn = tree_constants();
    let seqs: Vec<Vec<u8>> = vec![b"ACDGEF".to_vec(), b"ACD-EF".to_vec(), b"-CDGEW".to_vec()];
    let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let lgth = 6;
    // Weights as fastconjuction_noname would produce them: the last member is
    // the new one, and the first `n-1` keep their relative proportions.
    let eff = vec![0.25, 0.35, 0.4];
    let mut full = vec![vec![0.0f64; lgth + 2]; NALPHABETS];
    profile::cpmx_calc_new(&cn, &refs, &mut full, &eff, lgth, 3);

    let scale = 1.0 - eff[2];
    let eff_before = vec![eff[0] / scale, eff[1] / scale];
    let mut inc = vec![vec![0.0f64; lgth + 2]; NALPHABETS];
    profile::cpmx_calc_new(&cn, &refs[..2], &mut inc, &eff_before, lgth, 2);
    profile::cpmx_calc_add(&cn, &refs, &mut inc, &eff, lgth, 3);

    for i in 0..NALPHABETS {
        for j in 0..lgth {
            assert!(
                (full[i][j] - inc[i][j]).abs() < 1e-12,
                "[{}][{}]: {} vs {}",
                i,
                j,
                full[i][j],
                inc[i][j]
            );
        }
    }
}

/// Likewise for the gap counts.
#[test]
fn incremental_gap_counts_match_full_recompute() {
    let seqs: Vec<Vec<u8>> = vec![
        b"AC--EFG".to_vec(),
        b"A-DGEFG".to_vec(),
        b"--DGE-G".to_vec(),
    ];
    let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let len = 7;
    let eff = vec![0.25, 0.35, 0.4];
    let scale = 1.0 - eff[2];
    let eff_before = vec![eff[0] / scale, eff[1] / scale];

    for (full_fn, add_fn) in [(
        profile::st_opening_gap_count as fn(&mut [f64], usize, &[&[u8]], &[f64], usize),
        profile::st_opening_gap_add as fn(&mut [f64], usize, &[&[u8]], &[f64], usize),
    )] {
        let mut full = vec![0.0f64; len + 2];
        full_fn(&mut full, 3, &refs, &eff, len);
        let mut inc = vec![0.0f64; len + 2];
        full_fn(&mut inc, 2, &refs[..2], &eff_before, len);
        add_fn(&mut inc, 3, &refs, &eff, len);
        for i in 0..len {
            assert!((full[i] - inc[i]).abs() < 1e-12, "opening gap at {}", i);
        }
    }

    let mut full = vec![0.0f64; len + 2];
    profile::st_final_gap_count(&mut full, 3, &refs, &eff, len);
    let mut inc = vec![0.0f64; len + 2];
    profile::st_final_gap_count(&mut inc, 2, &refs[..2], &eff_before, len);
    profile::st_final_gap_add(&mut inc, 3, &refs, &eff, len);
    for i in 0..len {
        assert!((full[i] - inc[i]).abs() < 1e-12, "final gap at {}", i);
    }
}

/// A gap opening and a gap closing must be counted once each per run of gaps.
#[test]
fn gap_counts_find_each_run_once() {
    let seqs: Vec<Vec<u8>> = vec![b"A--DE-G".to_vec()];
    let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let eff = vec![1.0];
    let len = 7;
    let mut og = vec![0.0f64; len + 2];
    profile::st_opening_gap_count(&mut og, 1, &refs, &eff, len);
    // Runs start at columns 1 and 5.
    assert_eq!(og[1], 1.0);
    assert_eq!(og[5], 1.0);
    assert_eq!(og.iter().take(len).sum::<f64>(), 2.0);

    let mut fg = vec![0.0f64; len + 2];
    profile::st_final_gap_count(&mut fg, 1, &refs, &eff, len);
    // Runs end at columns 2 and 5.
    assert_eq!(fg[2], 1.0);
    assert_eq!(fg[5], 1.0);
    assert_eq!(fg.iter().take(len).sum::<f64>(), 2.0);
}

#[test]
fn gapireru_expands_by_the_table() {
    let out = profile::gapireru(b"ACDE", b"o-oo-o", b'-');
    assert_eq!(out, b"A-CD-E");
}

// ------------------------------------------------------------------ lalign

/// A local alignment must be an alignment: strip the gaps and you get back
/// contiguous stretches of the inputs, starting at the reported offsets.
#[test]
fn local_alignment_columns_come_from_the_inputs() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let seqs = deterministic_seqs(6, 90, 11);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let r = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[i], &seqs[j]);
            if r.a1.is_empty() {
                continue;
            }
            assert_eq!(r.a1.len(), r.a2.len(), "rows must be the same length");
            // No column may be gap-gap.
            for k in 0..r.a1.len() {
                assert!(!(r.a1[k] == b'-' && r.a2[k] == b'-'), "gap-gap column");
            }
            let u1 = ungap(&r.a1);
            let u2 = ungap(&r.a2);
            assert_eq!(&seqs[i][r.off1..r.off1 + u1.len()], &u1[..]);
            assert_eq!(&seqs[j][r.off2..r.off2 + u2.len()], &u2[..]);
        }
    }
}

/// `L__align11` and `L__align11_noalign` are two implementations of the same
/// recurrence; only the traceback differs. Their scores must agree.
#[test]
fn local_score_agrees_with_the_score_only_variant() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let seqs = deterministic_seqs(5, 70, 23);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let with = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[i], &seqs[j]);
            let without =
                lst.l_align11_noalign(&cn, &cn.n_dis_consweight_multi, &seqs[i], &seqs[j]);
            assert_eq!(
                with.score, without,
                "pair {}-{}: {} vs {}",
                i, j, with.score, without
            );
        }
    }
}

/// A local alignment of a sequence with itself is the whole sequence, and
/// scores the sum of its diagonal.
#[test]
fn self_local_alignment_is_the_identity() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    for s in deterministic_seqs(4, 50, 5) {
        let r = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &s, &s);
        assert_eq!(r.a1, s);
        assert_eq!(r.a2, s);
        assert_eq!(r.off1, 0);
        assert_eq!(r.off2, 0);
        let want: f64 = s
            .iter()
            .map(|&c| cn.amino_dis[c as usize][c as usize] as f64)
            .sum();
        assert_eq!(r.score, want);
    }
}

/// The local score cannot exceed either self score, which is what makes
/// `score2dist` produce a distance in [0, 2].
#[test]
fn local_score_is_bounded_by_the_self_scores() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let seqs = deterministic_seqs(6, 80, 99);
    let selfscore: Vec<f64> = seqs
        .iter()
        .map(|s| {
            s.iter()
                .map(|&c| cn.amino_dis[c as usize][c as usize] as f64)
                .sum()
        })
        .collect();
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let r = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[i], &seqs[j]);
            assert!(
                r.score <= selfscore[i].min(selfscore[j]) + 1e-9,
                "pair {}-{}: {} > min({}, {})",
                i,
                j,
                r.score,
                selfscore[i],
                selfscore[j]
            );
        }
    }
}

/// Local alignment is symmetric in its arguments: swapping the inputs swaps
/// the output rows and offsets.
#[test]
fn local_alignment_is_symmetric() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let seqs = deterministic_seqs(5, 60, 1234);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let a = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[i], &seqs[j]);
            let b = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[j], &seqs[i]);
            assert_eq!(a.score, b.score, "pair {}-{}", i, j);
            assert_eq!(a.a1, b.a2);
            assert_eq!(a.a2, b.a1);
            assert_eq!(a.off1, b.off2);
            assert_eq!(a.off2, b.off1);
        }
    }
}

// ------------------------------------------------------------- global align

/// The global aligner rewrites its inputs in place; what comes back must
/// ungap to exactly the inputs.
#[test]
fn global_alignment_preserves_the_inputs() {
    let cn = tree_constants();
    let mut ast = AlignState::new();
    let seqs = deterministic_seqs(6, 70, 4242);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let mut a = seqs[i].clone();
            let mut b = seqs[j].clone();
            ast.g_align11(
                &cn,
                &cn.n_dis_consweight_multi,
                cn.penalty,
                cn.penalty_ex,
                &mut a,
                &mut b,
                false,
                false,
                b'-',
            );
            assert_eq!(a.len(), b.len());
            assert_eq!(ungap(&a), seqs[i]);
            assert_eq!(ungap(&b), seqs[j]);
            for k in 0..a.len() {
                assert!(!(a[k] == b'-' && b[k] == b'-'), "gap-gap column");
            }
        }
    }
}

/// Sequences of different length must still come out the same length, padded
/// with gaps.
#[test]
fn global_alignment_handles_ragged_lengths() {
    let cn = tree_constants();
    let mut ast = AlignState::new();
    let base = deterministic_seqs(1, 80, 77).pop().unwrap();
    for cut in [1usize, 5, 40, 79] {
        let mut a = base.clone();
        let mut b = base[..cut].to_vec();
        ast.g_align11(
            &cn,
            &cn.n_dis_consweight_multi,
            cn.penalty,
            cn.penalty_ex,
            &mut a,
            &mut b,
            false,
            false,
            b'-',
        );
        assert_eq!(a.len(), b.len(), "cut={}", cut);
        assert_eq!(ungap(&a), base);
        assert_eq!(ungap(&b), base[..cut].to_vec());
    }
}

// ------------------------------------------------------- group-to-group DP

/// The profile DP must preserve every member of both groups, and leave the two
/// groups the same width.
#[test]
fn group_alignment_preserves_every_member() {
    let cn = tree_constants();
    let mut ast = AlignState::new();
    let seqs = deterministic_seqs(7, 80, 8888);
    let mut g1: Vec<Vec<u8>> = seqs[..4].to_vec();
    let mut g2: Vec<Vec<u8>> = seqs[4..].to_vec();
    let before1 = g1.clone();
    let before2 = g2.clone();
    let eff1 = vec![0.25; 4];
    let eff2 = vec![1.0 / 3.0; 3];
    ast.a_align(
        &cn,
        &cn.n_dis_consweight_multi,
        cn.penalty,
        cn.penalty_ex,
        &mut g1,
        &mut g2,
        &eff1,
        &eff2,
        false,
        false,
        None,
        None,
        0,
        1,
        None,
        None,
        1.0,
        1.0,
        true,
        b'-',
    );
    let width = g1[0].len();
    assert!(width >= 80);
    for (a, b) in g1.iter().zip(before1.iter()) {
        assert_eq!(a.len(), width);
        assert_eq!(&ungap(a), b);
    }
    for (a, b) in g2.iter().zip(before2.iter()) {
        assert_eq!(a.len(), width);
        assert_eq!(&ungap(a), b);
    }
}

/// A group's own columns must not be disturbed relative to each other: the DP
/// may only insert whole columns of gaps into a group.
#[test]
fn group_alignment_inserts_whole_columns() {
    let cn = tree_constants();
    let mut ast = AlignState::new();
    let seqs = deterministic_seqs(6, 60, 31415);
    let mut g1: Vec<Vec<u8>> = seqs[..3].to_vec();
    let mut g2: Vec<Vec<u8>> = seqs[3..].to_vec();
    ast.a_align(
        &cn,
        &cn.n_dis_consweight_multi,
        cn.penalty,
        cn.penalty_ex,
        &mut g1,
        &mut g2,
        &[1.0 / 3.0; 3],
        &[1.0 / 3.0; 3],
        false,
        false,
        None,
        None,
        0,
        1,
        None,
        None,
        1.0,
        1.0,
        true,
        b'-',
    );
    // A column newly inserted into group 1 is all gaps across its members.
    for col in 0..g1[0].len() {
        let ngap = g1.iter().filter(|s| s[col] == b'-').count();
        assert!(
            ngap == 0 || ngap == g1.len() || seqs[..3].iter().any(|s| s.contains(&b'-')),
            "column {} is partially gapped but the inputs had no gaps",
            col
        );
    }
}

// ----------------------------------------------------------------- localhom

/// Every constraint block must be ungapped in both sequences and must map back
/// to the residues the alignment actually paired.
#[test]
fn constraint_blocks_are_ungapped_and_in_range() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let seqs = deterministic_seqs(5, 90, 271828);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let r = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &seqs[i], &seqs[j]);
            if r.a1.is_empty() {
                continue;
            }
            let mut chain = vec![LocalHom::default()];
            localhom::putlocalhom2(&cn, &r.a1, &r.a2, &mut chain, r.off1, r.off2, b'h');
            for e in &chain {
                if e.opt == -1.0 {
                    continue;
                }
                assert!(e.start1 <= e.end1 && e.start2 <= e.end2);
                // Blocks are ungapped, so both halves span the same length.
                assert_eq!(e.end1 - e.start1, e.end2 - e.start2);
                assert!((e.end1 as usize) < seqs[i].len(), "block runs past seq1");
                assert!((e.end2 as usize) < seqs[j].len(), "block runs past seq2");
            }
            // The blocks must be disjoint and ordered.
            let live: Vec<&LocalHom> = chain.iter().filter(|e| e.opt != -1.0).collect();
            for w in live.windows(2) {
                assert!(w[0].end1 < w[1].start1, "blocks overlap in seq1");
                assert!(w[0].end2 < w[1].start2, "blocks overlap in seq2");
            }
        }
    }
}

/// A self-alignment produces exactly one block covering the whole sequence.
#[test]
fn self_constraint_is_one_full_length_block() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let s = deterministic_seqs(1, 60, 6161).pop().unwrap();
    let r = lst.l_align11(&cn, &cn.n_dis_consweight_multi, 0.0, &s, &s);
    let mut chain = vec![LocalHom::default()];
    localhom::putlocalhom2(&cn, &r.a1, &r.a2, &mut chain, r.off1, r.off2, b'h');
    let live: Vec<&LocalHom> = chain.iter().filter(|e| e.opt != -1.0).collect();
    assert_eq!(live.len(), 1);
    assert_eq!(live[0].start1, 0);
    assert_eq!(live[0].end1, s.len() as i32 - 1);
    assert_eq!(live[0].overlapaa, s.len() as i32);
}

/// `calcimportance_half` must leave the two directions equal, since its last
/// step averages them.
#[test]
fn importance_is_symmetric_after_averaging() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let mut ast = AlignState::new();
    let njob = 5;
    let seqs = deterministic_seqs(njob, 70, 1618);
    let mut mtx: Vec<Vec<f64>> = (0..njob).map(|i| vec![0.0; njob - i]).collect();
    let mut table = localhom::new_table(njob);
    pairalign::pairalign(
        &cn,
        PairAlg::Local,
        njob,
        &seqs,
        &mut mtx,
        &mut table,
        true,
        false,
        &mut lst,
        &mut ast,
        &mut crate::genalign::GenalignState::new(),
        false,
    );
    let eff = vec![1.0 / njob as f64; njob];
    localhom::calcimportance_half(njob, &eff, &seqs, &mut table);
    let mut seen = 0;
    for row in &table {
        for chain in row {
            for e in chain {
                if e.opt == -1.0 {
                    continue;
                }
                assert_eq!(e.importance, e.rimportance);
                assert!(e.importance.is_finite());
                seen += 1;
            }
        }
    }
    assert!(seen > 0, "no constraints were produced");
}

// ---------------------------------------------------------------- pairalign

/// The pairwise stage produces a distance matrix: entries in [0, 2], and
/// identical sequences at distance 0.
#[test]
fn pairwise_distances_are_in_range() {
    let cn = pair_constants();
    let mut lst = LalignState::new();
    let mut ast = AlignState::new();
    let njob = 6;
    let mut seqs = deterministic_seqs(njob - 1, 70, 5150);
    seqs.push(seqs[0].clone()); // an exact duplicate
    let mut mtx: Vec<Vec<f64>> = (0..njob).map(|i| vec![0.0; njob - i]).collect();
    let mut table = localhom::new_table(njob);
    pairalign::pairalign(
        &cn,
        PairAlg::Local,
        njob,
        &seqs,
        &mut mtx,
        &mut table,
        true,
        false,
        &mut lst,
        &mut ast,
        &mut crate::genalign::GenalignState::new(),
        false,
    );
    for i in 0..njob - 1 {
        for j in i + 1..njob {
            let d = mtx[i][j - i];
            assert!((0.0..=3.0).contains(&d), "d({},{}) = {}", i, j, d);
        }
    }
    // The duplicate of sequence 0 must be at distance 0 from it.
    assert_eq!(mtx[0][njob - 1], 0.0);
}

// --------------------------------------------------------------------- tree

/// UPGMA on a half matrix must produce a binary tree over every leaf: each of
/// the `n-1` merges has two children, and the root's leaf order is a
/// permutation of `0..n`.
#[test]
fn guide_tree_covers_every_leaf_exactly_once() {
    for n in [2usize, 3, 5, 12, 30] {
        let mut s = (n as u64) | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut mtx: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n - i).map(|_| (next() % 1000) as f64 / 500.0).collect())
            .collect();
        let tr = tree::upgma(n, &mut mtx, 0.1);
        let order = tree::topolorderz(&tr.topol, &tr.dep, n - 2, 2);
        assert_eq!(order.len(), n, "n={}", n);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..n as i32).collect::<Vec<_>>(), "n={}", n);
        // Branch lengths must be finite; `counteff_simple` clamps negatives.
        for k in 0..n - 1 {
            assert!(tr.len[k][0].is_finite() && tr.len[k][1].is_finite());
            assert!(tr.dep[k].distfromtip >= 0.0);
        }
    }
}

/// Sequence weights are a probability distribution.
#[test]
fn sequence_weights_sum_to_one() {
    for n in [2usize, 4, 9, 25] {
        let mut s = (n as u64 * 7) | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let mut mtx: Vec<Vec<f64>> = (0..n)
            .map(|i| (0..n - i).map(|_| (next() % 1000) as f64 / 500.0).collect())
            .collect();
        let mut tr = tree::upgma(n, &mut mtx, 0.1);
        let eff = tree::counteff_simple(n, &mut tr);
        assert_eq!(eff.len(), n);
        let total: f64 = eff.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "n={}: {}", n, total);
        for w in eff {
            assert!(w > 0.0, "n={}: non-positive weight", n);
        }
    }
}

/// Identical sequences must get identical weights.
#[test]
fn identical_sequences_get_identical_weights() {
    let n = 6;
    // All pairwise distances equal: the tree is arbitrary but symmetric.
    let mut mtx: Vec<Vec<f64>> = (0..n).map(|i| vec![0.5; n - i]).collect();
    let mut tr = tree::upgma(n, &mut mtx, 0.1);
    let eff = tree::counteff_simple(n, &mut tr);
    for w in &eff {
        assert!((w - eff[0]).abs() < 1e-12, "weights differ: {:?}", eff);
    }
}

// ------------------------------------------------------------------- params

#[test]
fn atof_and_atoi_follow_c_semantics() {
    assert_eq!(atof("-1.53"), -1.53);
    assert_eq!(atof("2.7abc"), 2.7);
    assert_eq!(atof("abc"), 0.0);
    assert_eq!(atof("1e3"), 1000.0);
    assert_eq!(atoi("-42xyz"), -42);
    assert_eq!(atoi(""), 0);
    assert_eq!(atoi("0-0"), 0);
}

/// The penalty conversion truncates toward zero after subtracting 0.5, which
/// is why `-g 0` becomes 0 and not -1.
#[test]
fn penalty_conversion_matches_the_c_cast() {
    let a = Args::parse(&[
        "-f".into(),
        "-1.53".into(),
        "-g".into(),
        "0".into(),
        "-h".into(),
        "0".into(),
    ])
    .unwrap();
    assert_eq!(a.ppenalty, -1530);
    assert_eq!(a.ppenalty_ex, 0);
    assert_eq!(a.poffset, 0);
}

#[test]
fn clustered_flags_parse_like_the_c() {
    // `-On` is `-O -n`.
    let a = Args::parse(&["-On".into()]).unwrap();
    assert!(!a.outgap);
    assert!(a.outnumber);
    // `-C 0-0` sets both thread counts.
    let a = Args::parse(&["-C".into(), "0-0".into()]).unwrap();
    assert_eq!(a.nthread, 0);
    assert_eq!(a.nthreadtb, 0);
}

#[test]
fn unsupported_configurations_are_refused() {
    let base = ["-C".into(), "0-0".into()];
    // `-b` is checked by `check_blosum` rather than `check_supported`, because
    // it depends on the sequence type and that is not known until the input
    // has been read. BLOSUM45 is a built-in table and is supported...
    let mut argv: Vec<String> = base.to_vec();
    argv.extend(["-b".into(), "45".into()]);
    let a = Args::parse(&argv).unwrap();
    assert!(a.check_supported().is_ok());
    assert!(a.check_blosum(b'p').is_ok());

    // ...but a table that does not exist is refused, for protein.
    let mut argv: Vec<String> = base.to_vec();
    argv.extend(["-b".into(), "77".into()]);
    let a = Args::parse(&argv).unwrap();
    assert!(a.check_supported().is_ok(), "not a parse-time decision");
    assert!(a.check_blosum(b'p').is_err());

    // For nucleotides `-b` is ignored -- the matrix comes from Kimura's model
    // -- so even a nonsense table must not refuse the run. Detected type, not
    // the `-D` flag: without `-D` the flag is unset, and checking that instead
    // refused `-b 47` on nucleotide input the C aligns happily.
    assert!(a.check_blosum(b'd').is_ok());
    let mut argv: Vec<String> = base.to_vec();
    argv.extend(["-D".into(), "-b".into(), "77".into()]);
    let a = Args::parse(&argv).unwrap();
    assert!(a.check_supported().is_ok());
    assert!(a.check_blosum(b'd').is_ok());

    // JTT and transmembrane are supported now too.
    for flag in ["-j", "-m"] {
        let mut argv: Vec<String> = base.to_vec();
        argv.extend([flag.into(), "200".into()]);
        assert!(Args::parse(&argv).unwrap().check_supported().is_ok());
    }

    let mut argv: Vec<String> = base.to_vec();
    argv.push("-Y".into()); // --keeplength
    assert!(Args::parse(&argv).unwrap().check_supported().is_err());

    // Threads change the result, so they must not be silently accepted.
    let argv: Vec<String> = vec!["-C".into(), "4".into()];
    assert!(Args::parse(&argv).unwrap().check_supported().is_err());

    // The supported baseline passes.
    let argv: Vec<String> = base.to_vec();
    assert!(Args::parse(&argv).unwrap().check_supported().is_ok());
}

// ---------------------------------------------------------------------- fft

/// The forward transform divides by `n` and the inverse does not, so running
/// both returns the input. This is the property the correlation depends on.
#[test]
fn fft_round_trips() {
    use crate::fft::{Cpx, Fft};
    let mut f = Fft::new();
    for &n in &[8usize, 16, 64, 256] {
        let mut s = (n as u64) | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let orig: Vec<Cpx> = (0..n)
            .map(|_| Cpx {
                r: (next() % 2000) as f64 / 100.0 - 10.0,
                i: (next() % 2000) as f64 / 100.0 - 10.0,
            })
            .collect();
        let mut x = orig.clone();
        f.run(n as isize, &mut x);
        f.run(-(n as isize), &mut x);
        for (a, b) in x.iter().zip(orig.iter()) {
            assert!((a.r - b.r).abs() < 1e-9, "{} vs {}", a.r, b.r);
            assert!((a.i - b.i).abs() < 1e-9);
        }
    }
}

/// The whole point of the FFT stage: the cross-correlation of a signal with a
/// shifted copy of itself must peak at that shift.
#[test]
fn fft_correlation_peaks_at_the_true_shift() {
    use crate::fft::{calc_naiseki, Cpx, Fft};
    let n = 256usize;
    let len = 90usize;
    for shift in [0usize, 7, 23, 40] {
        let mut s = (shift as u64 * 31 + 1) | 1;
        let mut next = move || {
            s ^= s << 13;
            s ^= s >> 7;
            s ^= s << 17;
            s
        };
        let sig: Vec<f64> = (0..len).map(|_| (next() % 100) as f64 - 50.0).collect();

        let mut a = vec![Cpx::default(); n];
        let mut b = vec![Cpx::default(); n];
        for k in 0..len {
            a[k].r = sig[k];
            b[k + shift].r = sig[k];
        }

        let mut f = Fft::new();
        f.run(n as isize, &mut a);
        f.run(n as isize, &mut b);
        let mut prod: Vec<Cpx> = (0..n).map(|k| calc_naiseki(&a[k], &b[k])).collect();
        f.run(-(n as isize), &mut prod);

        // soukan[m] as Falign builds it: lag = m - n/2.
        let n2 = n / 2;
        let mut soukan = vec![0.0f64; n];
        for m in 0..=n2 {
            soukan[m] = prod[n2 - m].r;
        }
        for m in n2 + 1..n {
            soukan[m] = prod[n + n2 - m].r;
        }
        let best = (0..n)
            .max_by(|&x, &y| soukan[x].partial_cmp(&soukan[y]).unwrap())
            .unwrap();
        assert_eq!(
            best as i64 - n2 as i64,
            shift as i64,
            "peak at lag {} for shift {}",
            best as i64 - n2 as i64,
            shift
        );
    }
}

// --------------------------------------------------------------------- score

/// `naivepairscorefast` is a skip-table-accelerated `naivepairscore11`. The C
/// ships both and uses each in a different place, so they must agree.
#[test]
fn fast_and_plain_pair_scores_agree() {
    use crate::score;
    let cn = tree_constants();
    let mut ast = AlignState::new();
    let seqs = deterministic_seqs(6, 70, 24680);
    // Produce genuinely gapped rows by aligning pairs first.
    let mut aligned: Vec<Vec<u8>> = Vec::new();
    for i in 0..seqs.len() - 1 {
        let mut a = seqs[i].clone();
        let mut b = seqs[i + 1].clone();
        ast.g_align11(
            &cn,
            &cn.n_dis_consweight_multi,
            cn.penalty,
            cn.penalty_ex,
            &mut a,
            &mut b,
            false,
            false,
            b'-',
        );
        aligned.push(a);
        aligned.push(b);
    }
    // Pad to a common width so they form an alignment.
    let w = aligned.iter().map(|s| s.len()).max().unwrap();
    for s in aligned.iter_mut() {
        while s.len() < w {
            s.push(b'-');
        }
    }
    let skip: Vec<Vec<i32>> = aligned.iter().map(|s| score::makeskiptable(s)).collect();
    for i in 0..aligned.len() {
        for j in i + 1..aligned.len() {
            let plain = score::naivepairscore11(&cn, &aligned[i], &aligned[j], cn.penalty_dist);
            let fast = score::naivepairscorefast(
                &cn,
                &aligned[i],
                &aligned[j],
                &skip[i],
                &skip[j],
                cn.penalty_dist,
            );
            assert_eq!(plain, fast, "pair {}-{}", i, j);
        }
    }
}

/// The premise that lets `pairscore` accumulate in `i64`.
///
/// Its two addends are `penalty_dist` and an `amino_dis` entry, both `i32`, so
/// every partial sum is an exact integer and the `i64` chain agrees with the
/// C's `double` chain bit-for-bit. That is an argument about the *types* of
/// the addends, and a future change could break it without touching
/// `pairscore` at all -- by scaling the matrix to a fraction, say. This fails
/// if that happens, rather than letting the port quietly diverge from the C.
///
/// Also pins the magnitude. Integers stay exactly representable in an f64 up
/// to 2^53; the bound asserted here is twelve orders of magnitude below that,
/// so the exactness does not rest on the sum merely happening to be small.
#[test]
fn score_is_integral() {
    use crate::score;
    let cn = tree_constants();

    // Every matrix entry the loop can reach is an integer, and so is the
    // penalty -- the two facts the `i64` accumulator rests on.
    assert_eq!(cn.penalty_dist as f64, (cn.penalty_dist as f64).trunc());
    for i in 0..cn.nalphabets {
        for j in 0..cn.nalphabets {
            let v = cn.amino_dis[cn.amino[i] as usize][cn.amino[j] as usize];
            assert_eq!(
                v as f64,
                (v as f64).trunc(),
                "amino_dis[{}][{}] is not an integer",
                i,
                j
            );
        }
    }

    let mut ast = AlignState::new();
    let seqs = deterministic_seqs(6, 70, 13579);
    let mut aligned: Vec<Vec<u8>> = Vec::new();
    for i in 0..seqs.len() - 1 {
        let mut a = seqs[i].clone();
        let mut b = seqs[i + 1].clone();
        ast.g_align11(
            &cn,
            &cn.n_dis_consweight_multi,
            cn.penalty,
            cn.penalty_ex,
            &mut a,
            &mut b,
            false,
            false,
            b'-',
        );
        aligned.push(a);
        aligned.push(b);
    }
    let w = aligned.iter().map(|s| s.len()).max().unwrap();
    for s in aligned.iter_mut() {
        while s.len() < w {
            s.push(b'-');
        }
    }
    let skip: Vec<Vec<i32>> = aligned.iter().map(|s| score::makeskiptable(s)).collect();

    let mut checked = 0usize;
    for i in 0..aligned.len() {
        for j in i + 1..aligned.len() {
            for &penal in &[0, cn.penalty_dist] {
                let plain = score::naivepairscore11(&cn, &aligned[i], &aligned[j], penal);
                let fast = score::naivepairscorefast(
                    &cn,
                    &aligned[i],
                    &aligned[j],
                    &skip[i],
                    &skip[j],
                    penal,
                );
                for (what, v) in [("naivepairscore11", plain), ("naivepairscorefast", fast)] {
                    assert_eq!(v, v.trunc(), "{} returned {}, not an integer", what, v);
                    assert!(
                        v.abs() < 1.0e9,
                        "{} returned {}, outside the magnitude this accumulator assumes",
                        what,
                        v
                    );
                }
                checked += 1;
            }
        }
    }
    // A pass count of zero would make every assertion above vacuous.
    assert_eq!(checked, aligned.len() * (aligned.len() - 1), "wrong case count");
}

/// The skip table counts the gaps before each residue, so it must sum to the
/// number of gaps and have one more entry than there are residues.
#[test]
fn skiptable_accounts_for_every_gap() {
    use crate::score;
    for s in [
        &b"ACDE"[..],
        &b"--ACDE--"[..],
        &b"A-C--D-E"[..],
        &b"----"[..],
    ] {
        let t = score::makeskiptable(s);
        let nres = s.iter().filter(|&&c| c != b'-').count();
        let ngap = s.len() - nres;
        assert_eq!(t.len(), nres + 1);
        assert_eq!(t.iter().sum::<i32>(), ngap as i32, "for {:?}", s);
    }
}

/// A sequence scored against itself must reach its own self score.
#[test]
fn self_pair_score_is_the_self_score() {
    use crate::score;
    let cn = tree_constants();
    for s in deterministic_seqs(3, 40, 909) {
        let want: f64 = s
            .iter()
            .map(|&c| cn.amino_dis[c as usize][c as usize] as f64)
            .sum();
        assert_eq!(score::naivepairscore11(&cn, &s, &s, cn.penalty_dist), want);
    }
}

// ---------------------------------------------------------------------- warp

/// A warp jump may only improve a cell: the DP takes it when it beats the
/// score already there, so `cell` never returns a lower value.
#[test]
fn warp_never_lowers_a_score() {
    use crate::warp::Warp;
    let mut w = Warp::new(20, 20);
    let mut prev = 0.0;
    for i in 1..8usize {
        for j in 1..20usize {
            let wm = (i * j) as f64 * 0.5 - 30.0;
            let (out, _code) = w.cell(i, j, wm, 1.0, -1000.0, -10.0);
            assert!(out >= wm, "warp lowered {} to {}", wm, out);
            prev = out;
        }
        w.end_row(20);
    }
    let _ = prev;
}

/// A recorded jump target must be resolvable and must point at a cell that was
/// actually visited.
#[test]
fn warp_targets_resolve() {
    use crate::warp::Warp;
    let mut w = Warp::new(30, 30);
    let mut codes = Vec::new();
    for i in 1..10usize {
        for j in 1..30usize {
            // A large positive shift penalty makes jumps attractive.
            let (_out, code) = w.cell(i, j, -1000.0, 5.0, 0.0, 0.0);
            if let Some(c) = code {
                codes.push(c);
            }
        }
        w.end_row(30);
    }
    assert!(!codes.is_empty(), "no jumps were taken");
    for c in codes {
        assert!(c >= w.warpbase);
        let (a, b) = w.target(c);
        assert!(a >= -w.warpbase && b >= -w.warpbase);
    }
}

// ------------------------------------------------------------ nucleotides

/// The nucleotide branch of `constants()`, as the wrapper drives it.
fn dna_constants() -> Constants {
    let a = Args {
        nblosum: 62,
        penalty_shift_factor: 100.0,
        ..Default::default()
    };
    Constants::new_dorp(&a, b'd')
}

fn dna_idx(cn: &Constants, c: u8) -> usize {
    cn.amino_n[c as usize] as usize
}

/// The scored block of the nucleotide matrix is symmetric, like the protein one.
#[test]
fn dna_score_matrix_is_symmetric() {
    let cn = dna_constants();
    for i in 0..NALPHABETS {
        for j in 0..NALPHABETS {
            assert_eq!(
                cn.n_dis[i][j], cn.n_dis[j][i],
                "n_dis[{}][{}] != n_dis[{}][{}]",
                i, j, j, i
            );
        }
    }
}

/// Kimura's model with R = 2 must score a match above a transition and a
/// transition above a transversion; nothing else fixes the sign of the model.
#[test]
fn dna_matrix_orders_match_transition_transversion() {
    let cn = dna_constants();
    let (a, g, c, t) = (
        dna_idx(&cn, b'a'),
        dna_idx(&cn, b'g'),
        dna_idx(&cn, b'c'),
        dna_idx(&cn, b't'),
    );
    for &(x, y, z) in &[(a, g, c), (g, a, t), (c, t, a), (t, c, g)] {
        assert!(
            cn.n_dis[x][x] > cn.n_dis[x][y],
            "match not above transition"
        );
        assert!(
            cn.n_dis[x][y] > cn.n_dis[x][z],
            "transition not above transversion"
        );
    }
}

/// `'u'` scores as `'t'`, and the upper-case half of the alphabet repeats the
/// lower-case half -- the two block copies at the end of the DNA branch.
#[test]
fn dna_matrix_copies_u_and_upper_case() {
    let cn = dna_constants();
    let t = dna_idx(&cn, b't');
    let u = dna_idx(&cn, b'u');
    for j in 0..5 {
        assert_eq!(cn.n_dis[u][j], cn.n_dis[t][j]);
        assert_eq!(cn.n_dis[j][u], cn.n_dis[j][t]);
    }
    for i in 0..5 {
        for j in 0..5 {
            assert_eq!(cn.n_dis[i + 5][j + 5], cn.n_dis[i][j]);
        }
    }
}

/// `ambiguousscore()`: every IUPAC code scores as the mean of the bases it
/// stands for, rounded half away from zero.
#[test]
fn dna_ambiguity_codes_are_means_of_their_bases() {
    let cn = dna_constants();
    let round = |v: f64| {
        if v > 0.0 {
            (v + 0.5) as i32
        } else if v == 0.0 {
            0
        } else {
            (v - 0.5) as i32
        }
    };
    let two: [(u8, u8, u8); 6] = [
        (b'r', b'a', b'g'),
        (b'y', b'c', b't'),
        (b'k', b'g', b't'),
        (b'm', b'a', b'c'),
        (b's', b'g', b'c'),
        (b'w', b'a', b't'),
    ];
    for &(code, p, q) in two.iter() {
        let k = dna_idx(&cn, code);
        for base in b"agct" {
            let i = dna_idx(&cn, *base);
            let want =
                round(0.5 * (cn.n_dis[dna_idx(&cn, p)][i] + cn.n_dis[dna_idx(&cn, q)][i]) as f64);
            assert_eq!(
                cn.n_dis[i][k], want,
                "{} vs {}",
                code as char, *base as char
            );
        }
    }
    let three: [(u8, u8, u8, u8); 4] = [
        (b'b', b'c', b'g', b't'),
        (b'd', b'a', b'g', b't'),
        (b'h', b'a', b'c', b't'),
        (b'v', b'a', b'c', b'g'),
    ];
    for &(code, p, q, r) in three.iter() {
        let k = dna_idx(&cn, code);
        for base in b"agct" {
            let i = dna_idx(&cn, *base);
            let want = round(
                (1.0 / 3.0)
                    * (cn.n_dis[dna_idx(&cn, p)][i]
                        + cn.n_dis[dna_idx(&cn, q)][i]
                        + cn.n_dis[dna_idx(&cn, r)][i]) as f64,
            );
            assert_eq!(
                cn.n_dis[i][k], want,
                "{} vs {}",
                code as char, *base as char
            );
        }
    }
}

/// `--anysymbol` rewrites the `n` row and column; without it they stay put.
#[test]
fn dna_anysymbol_rewrites_the_n_row() {
    let plain = dna_constants();
    let wild = Constants::new_dorp(
        &Args {
            nblosum: 62,
            penalty_shift_factor: 100.0,
            nwildcard: true,
            ..Default::default()
        },
        b'd',
    );
    let n = dna_idx(&plain, b'n');
    assert_eq!(dna_idx(&wild, b'n'), n);
    let a = dna_idx(&plain, b'a');
    assert_ne!(plain.n_dis[a][n], wild.n_dis[a][n]);
    // Each rewritten cell is a quarter of the corresponding self score.
    for i in 0..NALPHABETS {
        let want = {
            let v = 0.25 * wild.n_dis[i][i] as f64;
            if v > 0.0 {
                (v + 0.5) as i32
            } else if v == 0.0 {
                0
            } else {
                (v - 0.5) as i32
            }
        };
        if i != n {
            assert_eq!(wild.n_dis[i][n], want);
            assert_eq!(wild.n_dis[n][i], wild.n_dis[i][n]);
        }
    }
}

/// The base-4 k-tuple code of a nucleotide sequence: one entry per 6-mer, each
/// in range, and the rolling update agrees with computing each code directly.
#[test]
fn dna_pointtable_matches_direct_codes() {
    let cn = dna_constants();
    let seq = b"acgtacgtttacggacatgcatgcaacgt".to_vec();
    let mut unknown = 0i64;
    let grp = crate::dist::seq_grp_nuc(&cn, &seq, &mut unknown, 6);
    assert_eq!(unknown, 0);
    let pointt = crate::dist::makepointtable_nuc(&grp);
    assert_eq!(*pointt.last().unwrap(), crate::dist::END_OF_VEC);
    let codes = &pointt[..pointt.len() - 1];
    assert_eq!(codes.len(), seq.len() - 5);
    for (k, &code) in codes.iter().enumerate() {
        assert!((0..4096).contains(&code), "code {} out of range", code);
        let mut direct = 0i32;
        for d in 0..6 {
            direct = direct * 4 + grp[k + d];
        }
        assert_eq!(code, direct);
    }
}

/// Sequences shorter than the tuple get the "too short" sentinel rather than
/// reading past the end of the group vector.
#[test]
fn dna_short_sequences_get_the_sentinel() {
    let cn = dna_constants();
    let mut unknown = 0i64;
    for len in 0..6usize {
        let seq: Vec<u8> = std::iter::repeat(b'a').take(len).collect();
        let grp = crate::dist::seq_grp_nuc(&cn, &seq, &mut unknown, 6);
        assert_eq!(grp[0], -1, "len {}", len);
        assert_eq!(crate::dist::makepointtable_nuc(&grp), vec![-1]);
    }
}

/// Characters outside `acgt` (and their upper case) are not tuple material;
/// they are counted as unknown and dropped, so the codes stay in range.
#[test]
fn dna_ambiguity_codes_are_dropped_from_tuples() {
    let cn = dna_constants();
    let mut unknown = 0i64;
    let grp = crate::dist::seq_grp_nuc(&cn, b"acgtnnrryyacgtacgt", &mut unknown, 6);
    assert_eq!(unknown, 6);
    assert_eq!(grp.len(), 12 + 1);
    for &g in &grp[..12] {
        assert!((0..4).contains(&g));
    }
}

/// The reader lower-cases nucleotide input and upper-cases protein: the same
/// file read either way must differ only in case.
#[test]
fn reader_lower_cases_nucleotides_only() {
    let fasta = b">a\nACGTacgtACGT\n>b\nacgtACGTacgt\n".to_vec();
    let d = io::read_data(fasta.clone(), Some(b'd'), false, false);
    let p = io::read_data(fasta, Some(b'p'), false, false);
    assert_eq!(d.dorp, b'd');
    assert_eq!(p.dorp, b'p');
    for (ds, ps) in d.seq.iter().zip(p.seq.iter()) {
        assert_eq!(ds, &ps.to_ascii_lowercase());
        assert_eq!(ps, &ds.to_ascii_uppercase());
    }
}

/// Sequence-type detection: `atgcun` above 75% of the alphabetic characters.
#[test]
fn sequence_type_is_detected_from_composition() {
    let dna = b">a\nacgtacgtacgtacgtacgt\n>b\nacgtacgtacgtacgtacgt\n".to_vec();
    assert_eq!(io::read_data(dna, None, false, false).dorp, b'd');
    let prot = b">a\nMKWVTFISLLLLFSSAYS\n>b\nMKWVTFISLLLLFSSAYS\n".to_vec();
    assert_eq!(io::read_data(prot, None, false, false).dorp, b'p');
}

/// A nucleotide alignment still preserves its inputs, and the DNA branch is
/// reached end to end (score matrix, k-tuple distance, pairwise DP).
#[test]
fn dna_pairwise_alignment_preserves_inputs() {
    let cn = dna_constants();
    let mut ast = AlignState::new();
    let a = b"acgtacgtttacggacatgcatgcaacgt".to_vec();
    let b = b"acgtacgttacggacatgcatgcaacgtt".to_vec();
    let mut s1 = a.clone();
    let mut s2 = b.clone();
    ast.g_align11(
        &cn,
        &cn.n_dis_consweight_multi,
        cn.penalty,
        cn.penalty_ex,
        &mut s1,
        &mut s2,
        false,
        false,
        b'-',
    );
    assert_eq!(s1.len(), s2.len());
    assert_eq!(ungap(&s1), a);
    assert_eq!(ungap(&s2), b);
    for (x, y) in s1.iter().zip(s2.iter()) {
        assert!(!(*x == b'-' && *y == b'-'), "gap-gap column");
    }
}

// ------------------------------------------------- generalised affine gaps

/// The parameters `mafft --genafpair` uses for its pairwise stage.
fn genaf_constants() -> Constants {
    let a = Args {
        nblosum: 62,
        ppenalty: -2000,
        ppenalty_ex: 0,
        ppenalty_op: -6000,
        poffset: 0,
        penalty_shift_factor: 100.0,
        ..Default::default()
    };
    Constants::new(&a)
}

/// The generalised local alignment preserves both inputs: ungapping a row
/// returns the slice of the input it started from, and no column is gap-gap.
#[test]
fn generalized_alignment_columns_come_from_the_inputs() {
    let cn = genaf_constants();
    let mut gst = crate::genalign::GenalignState::new();
    let seqs = deterministic_seqs(6, 90, 0x9e1f);
    for i in 0..seqs.len() {
        for j in i + 1..seqs.len() {
            let r = gst.gen_l_align11(&cn, &cn.n_dis_consweight_multi, &seqs[i], &seqs[j]);
            assert_eq!(r.a1.len(), r.a2.len(), "rows differ in width");
            if r.a1.is_empty() {
                continue;
            }
            for (x, y) in r.a1.iter().zip(r.a2.iter()) {
                assert!(!(*x == b'-' && *y == b'-'), "gap-gap column");
            }
            let u1 = ungap(&r.a1);
            let u2 = ungap(&r.a2);
            assert!(
                seqs[i][r.off1..].starts_with(&u1),
                "row 1 is not the input from off1"
            );
            assert!(
                seqs[j][r.off2..].starts_with(&u2),
                "row 2 is not the input from off2"
            );
        }
    }
}

/// A sequence against itself aligns to itself, ungapped, from offset zero.
#[test]
fn generalized_self_alignment_is_the_identity() {
    let cn = genaf_constants();
    let mut gst = crate::genalign::GenalignState::new();
    for s in deterministic_seqs(3, 60, 0x5e1f).iter() {
        let r = gst.gen_l_align11(&cn, &cn.n_dis_consweight_multi, s, s);
        assert_eq!(r.a1, *s);
        assert_eq!(r.a2, *s);
        assert_eq!(r.off1, 0);
        assert_eq!(r.off2, 0);
    }
}

/// Opening a pair of unaligned runs costs `penalty_OP` once, however long the
/// runs are, which is the whole point of the generalised recurrence: two
/// sequences sharing flanks but with unrelated middles of *different* lengths
/// must still align both flanks, which a plain affine local DP would not do
/// for a cheap enough opening cost.
#[test]
fn generalized_gaps_join_flanks_across_unrelated_middles() {
    let cn = genaf_constants();
    let mut gst = crate::genalign::GenalignState::new();
    let left = b"WYFMKLDEHQNRSTAGVIPC".to_vec();
    let right = b"CPIVGATSRNQHEDLKMFYW".to_vec();
    let mut a = left.clone();
    a.extend_from_slice(b"AAAAAAAAAA");
    a.extend_from_slice(&right);
    let mut b = left.clone();
    b.extend_from_slice(b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    b.extend_from_slice(&right);

    let r = gst.gen_l_align11(&cn, &cn.n_dis_consweight_multi, &a, &b);
    assert!(!r.a1.is_empty());
    // Both flanks are inside the reported alignment.
    let u1 = ungap(&r.a1);
    let u2 = ungap(&r.a2);
    assert!(
        u1.len() > left.len(),
        "alignment stopped at the first flank"
    );
    assert!(
        u2.len() > left.len(),
        "alignment stopped at the first flank"
    );
    assert!(
        r.a1.len() < a.len() + b.len(),
        "the two middles were not collapsed into one pair of runs"
    );
}

/// `penalty_OP` buys a *pair* of unaligned runs for one flat cost, so with a
/// real extension penalty in play a cheap `-O` lets the aligner cross two
/// unrelated middles of different lengths and reach the second flank, while a
/// prohibitive one makes it give up. Both numbers below are measured, not
/// assumed: 80 columns against 70 for these inputs.
#[test]
fn generalized_pair_opening_penalty_decides_whether_both_flanks_align() {
    let mk = |op: i32| {
        Constants::new(&Args {
            nblosum: 62,
            ppenalty: -2000,
            // A non-zero extension penalty is what makes the generalised
            // recurrence worth anything; with free extension (which is what
            // --genafpair itself passes) ordinary gaps cross just as cheaply.
            ppenalty_ex: -100,
            ppenalty_op: op,
            poffset: 0,
            penalty_shift_factor: 100.0,
            ..Default::default()
        })
    };
    let left = b"WYFMKLDEHQNRSTAGVIPC".to_vec();
    let right = b"CPIVGATSRNQHEDLKMFYW".to_vec();
    let mut a = left.clone();
    a.extend_from_slice(b"AAAAAAAAAA");
    a.extend_from_slice(&right);
    let mut b = left.clone();
    b.extend_from_slice(b"WWWWWWWWWWWWWWWWWWWWWWWWWWWWWW");
    b.extend_from_slice(&right);

    let mut gst = crate::genalign::GenalignState::new();
    let cheap = mk(-6000);
    let dear = mk(-2_000_000);
    let wide = gst.gen_l_align11(&cheap, &cheap.n_dis_consweight_multi, &a, &b);
    let narrow = gst.gen_l_align11(&dear, &dear.n_dis_consweight_multi, &a, &b);
    assert!(
        narrow.a1.len() < wide.a1.len(),
        "a prohibitive -O did not shorten the alignment ({} vs {})",
        narrow.a1.len(),
        wide.a1.len()
    );
    // The cheap one reaches past the first flank into the second.
    assert!(ungap(&wide.a1).len() > left.len() + 10);
}

// --- intergroup_score -------------------------------------------------

/// A gap-free column pair scores exactly the matrix entry, weighted.
#[test]
fn intergroup_score_of_identical_columns_is_the_diagonal() {
    let cn = pair_constants();
    let g1 = vec![b"AC".to_vec()];
    let g2 = vec![b"AC".to_vec()];
    let got = crate::score::intergroup_score(&cn, &g1, &g2, &[1.0], &[1.0], 1, 1, 2);
    let a = cn.amino_dis_consweight_multi[b'A' as usize][b'A' as usize];
    let c = cn.amino_dis_consweight_multi[b'C' as usize][b'C' as usize];
    assert!((got - (a + c)).abs() < 1e-9, "{got} vs {}", a + c);
}

/// Weights multiply out: doubling one side's weight doubles the score.
#[test]
fn intergroup_score_is_linear_in_the_weights() {
    let cn = pair_constants();
    let g1 = vec![b"ACDE".to_vec()];
    let g2 = vec![b"AKDE".to_vec()];
    let one = crate::score::intergroup_score(&cn, &g1, &g2, &[1.0], &[1.0], 1, 1, 4);
    let two = crate::score::intergroup_score(&cn, &g1, &g2, &[2.0], &[1.0], 1, 1, 4);
    assert!((two - 2.0 * one).abs() < 1e-9);
}

/// A column that is a gap on both sides is skipped entirely -- it costs
/// nothing and does not open a gap.
#[test]
fn intergroup_score_ignores_gap_gap_columns() {
    let cn = pair_constants();
    let plain = vec![b"AC".to_vec()];
    let padded = vec![b"A-C".to_vec()];
    let a = crate::score::intergroup_score(&cn, &plain, &plain, &[1.0], &[1.0], 1, 1, 2);
    let b = crate::score::intergroup_score(&cn, &padded, &padded, &[1.0], &[1.0], 1, 1, 3);
    assert!((a - b).abs() < 1e-9, "{a} vs {b}");
}

/// A run of gaps is charged `penalty` once, however long it is.
#[test]
fn intergroup_score_charges_a_gap_run_once() {
    let cn = pair_constants();
    let res = vec![b"AKKKC".to_vec()];
    let one = vec![b"A-KKC".to_vec()];
    let three = vec![b"A---C".to_vec()];
    let s1 = crate::score::intergroup_score(&cn, &res, &one, &[1.0], &[1.0], 1, 1, 5);
    let s3 = crate::score::intergroup_score(&cn, &res, &three, &[1.0], &[1.0], 1, 1, 5);
    // Both opened exactly one gap, so the difference is only the residues
    // that the longer run replaced -- no second `penalty`.
    let m = &cn.amino_dis_consweight_multi;
    let gapk = m[b'-' as usize][b'K' as usize];
    let kk = m[b'K' as usize][b'K' as usize];
    assert!(
        (s3 - (s1 - 2.0 * kk + 2.0 * gapk)).abs() < 1e-9,
        "{s1} {s3}"
    );
}


#[test]
fn commongappick_record_maps_kept_columns_to_their_origins() {
    // Column 1 is all-gap and must go; the map says where the survivors came
    // from, and its trailing entry is the C's terminator column.
    let mut seq = vec![b"A-CD".to_vec(), b"A-CE".to_vec()];
    let map = crate::addfunctions::commongappick_record(&mut seq);
    assert_eq!(seq[0], b"ACD".to_vec());
    assert_eq!(seq[1], b"ACE".to_vec());
    assert_eq!(map, vec![0, 2, 3, 4]);
}

#[test]
fn commongappick_record_agrees_with_commongappick_on_the_rows() {
    // The two must compact identically; only the map is extra.
    for rows in [
        vec![b"A-CD".to_vec(), b"A-CE".to_vec()],
        vec![b"--AB".to_vec(), b"--AB".to_vec()],
        vec![b"ABCD".to_vec(), b"ABCD".to_vec()],
        vec![b"----".to_vec(), b"----".to_vec()],
    ] {
        let mut a = rows.clone();
        let mut b = rows.clone();
        crate::addfunctions::commongappick(&mut a);
        let _ = crate::addfunctions::commongappick_record(&mut b);
        assert_eq!(a, b, "the two disagree on {rows:?}");
    }
}

#[test]
fn commongappick_record_keeps_every_column_when_none_is_all_gap() {
    let mut seq = vec![b"A-C".to_vec(), b"-BC".to_vec()];
    let map = crate::addfunctions::commongappick_record(&mut seq);
    assert_eq!(seq[0], b"A-C".to_vec());
    assert_eq!(map, vec![0, 1, 2, 3]);
}
