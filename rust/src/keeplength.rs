//! `--keeplength` (`-Y`), transliterated from `core/addfunctions.c`.
//!
//! With `--add`, `--keeplength` promises that the reference alignment comes
//! out exactly as wide as it went in. Any column the progressive alignment
//! inserted into it is therefore deleted again afterwards, taking whatever
//! residues of the added sequences sat in it -- so the option loses letters,
//! and `--mapout` exists to say which.
//!
//! Three passes, in the C's order: remember which columns of the reference
//! were all-gap before aligning (`recordoriginalgaps`), drop the columns the
//! alignment added (`deletenewinsertions_whole`, which also records what it
//! deleted), and put the remembered all-gap columns back
//! (`restoreoriginalgaps`).

use crate::io::B;

/// One deleted run: `pos` residues into the added sequence, `len` long.
/// The C's `GapPos`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GapPos {
    pub pos: i32,
    pub len: i32,
}

/// Read one byte the way C reads a NUL-terminated string: off the end is 0.
fn at(s: &[u8], i: usize) -> u8 {
    if i < s.len() {
        s[i]
    } else {
        0
    }
}

/// `recordoriginalgaps()`: `'-'` for each column that is a gap in all `n`
/// reference sequences, `'o'` for the rest. Sized from `s[0]`.
pub fn record_original_gaps(n: usize, s: &[Vec<u8>]) -> Vec<u8> {
    let len = s[0].len();
    (0..len)
        .map(|i| {
            if (0..n).all(|j| at(&s[j], i) == b'-') {
                b'-'
            } else {
                b'o'
            }
        })
        .collect()
}

/// `deletenewinsertions_whole()`: delete every column that is all-gap across
/// the `on` reference sequences, and report, per added sequence, the runs of
/// residues that went with them.
///
/// Returns the number of columns deleted. `seq` is `oseq` followed by `aseq`,
/// as the C's two pointers into one array are.
pub fn delete_new_insertions_whole(
    on: usize,
    an: usize,
    seq: &mut [Vec<u8>],
) -> (usize, Vec<Vec<GapPos>>) {
    let len = seq[0].len();

    // `eqseq`: '=' marks a column the alignment inserted into the reference.
    // Length len+1, the last entry the NUL the C's loops stop on.
    let mut eqseq = vec![0u8; len + 1];
    for (i, e) in eqseq.iter_mut().enumerate().take(len) {
        *e = if (0..on).all(|j| at(&seq[j], i) == b'-') {
            b'='
        } else {
            b'o'
        };
    }

    let mut deletelist: Vec<Vec<GapPos>> = Vec::with_capacity(an);
    for j in 0..an {
        let mut list: Vec<GapPos> = Vec::new();
        let mut q = 0usize;
        let mut i = 0usize;
        while i < len {
            let mut tmpc = at(&seq[on + j], i);
            let qstart = q;
            // Walk the whole run of inserted columns this sequence has
            // residues in. `i` can reach `len` here; `eqseq[len]` is the NUL,
            // so the loop stops there.
            while tmpc != b'-' && eqseq[i] == b'=' {
                i += 1;
                q += 1;
                tmpc = at(&seq[on + j], i);
            }
            if qstart != q {
                list.push(GapPos {
                    pos: qstart as i32,
                    len: (q - qstart) as i32,
                });
            }
            if tmpc != b'-' {
                q += 1;
            }
            i += 1;
        }
        deletelist.push(list);
    }

    // Squeeze out the inserted columns.
    let mut p = 0usize;
    for i in 0..len {
        if eqseq[i] != b'=' {
            if p != i {
                for j in 0..on + an {
                    seq[j][p] = seq[j][i];
                }
            }
            p += 1;
        }
    }
    for j in 0..on + an {
        seq[j].truncate(p);
    }

    (len - p, deletelist)
}

/// `restoreoriginalgaps()`: widen every sequence back out, putting a gap
/// wherever the reference had an all-gap column before the alignment.
pub fn restore_original_gaps(n: usize, seq: &mut [Vec<u8>], originalgaps: &[u8]) {
    let lenf = originalgaps.len();
    for s in seq.iter_mut().take(n) {
        let mut tmp = vec![0u8; lenf];
        let mut p = 0usize;
        for (j, t) in tmp.iter_mut().enumerate() {
            if originalgaps[j] == b'-' {
                *t = b'-';
            } else {
                *t = at(s, p);
                p += 1;
            }
        }
        *s = tmp;
    }
}

fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.len() > hay.len() {
        return None;
    }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// Where the printable part of a name starts: past the `'='` sentinel, past
/// the `--reorder` numbering if there is any, and past an `_oe_` marker.
///
/// Returned as an offset into `name` because the compact map needs to splice
/// its insertion summary in at exactly that point.
fn name_offset(name: &[u8], outnumber: bool) -> usize {
    let mut off = 1usize;
    if outnumber {
        // The C adds 8 to a `strstr` for the 7-byte "_numo_e", so it lands
        // past the trailing '_' as well.
        if let Some(p) = find_sub(&name[off..], b"_numo_e") {
            off += p + 8;
        }
    }
    if let Some(p) = find_sub(&name[off..], b"_oe_") {
        off += p + 4;
    }
    off
}

/// Paint the deleted runs onto a per-residue map of one added sequence.
fn gapped_of(len: usize, list: &[GapPos]) -> Vec<u8> {
    let mut gapped = vec![0u8; len + 1];
    for e in list.iter() {
        let mut p = e.pos as usize;
        for _ in 0..e.len {
            gapped[p] = b'-';
            p += 1;
        }
    }
    gapped
}

/// `reconstructdeletemap()`: one line per residue of each added sequence,
/// saying where it ended up in the reference alignment (or that it was cut).
pub fn reconstruct_delete_map(
    addbk: &[Vec<u8>],
    deletelist: &[Vec<GapPos>],
    realn: &[Vec<u8>],
    name: &[Vec<u8>],
    outnumber: bool,
) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    for i in 0..addbk.len() {
        let len = addbk[i].len();
        let gapped = gapped_of(len, &deletelist[i]);
        out.push(b'>');
        out.extend_from_slice(&name[i][name_offset(&name[i], outnumber)..]);
        out.push(b'\n');
        out.extend_from_slice(
            b"# letter, position in the original sequence, position in the reference alignment\n",
        );
        let mut p = 0usize;
        for j in 0..len {
            while at(&realn[i], p) == b'-' {
                p += 1;
            }
            if gapped[j] == b'-' {
                out.extend_from_slice(format!("{}, {}, -\n", addbk[i][j] as char, j + 1).as_bytes());
            } else {
                out.extend_from_slice(
                    format!("{}, {}, {}\n", addbk[i][j] as char, j + 1, p + 1).as_bytes(),
                );
                p += 1;
            }
        }
    }
    out
}

/// `reconstructdeletemap_compact()`: one line per deleted *run* rather than
/// per residue, and -- this is not cosmetic -- it also rewrites `name`, so the
/// added sequences carry their insertion summary into the alignment written
/// afterwards.
pub fn reconstruct_delete_map_compact(
    addbk: &[Vec<u8>],
    deletelist: &[Vec<GapPos>],
    realn: &[Vec<u8>],
    name: &mut [Vec<u8>],
    outnumber: bool,
) -> Vec<u8> {
    let mut out: Vec<u8> = Vec::new();
    out.extend_from_slice(b"# Insertion in added sequence > Position in reference\n");
    for i in 0..addbk.len() {
        let len = addbk[i].len();
        let gapped = gapped_of(len, &deletelist[i]);
        if deletelist[i].is_empty() {
            continue;
        }
        let off = name_offset(&name[i], outnumber);
        out.push(b'>');
        out.extend_from_slice(&name[i][off..]);
        out.push(b'\n');

        let mut status = -1i32;
        let mut insstr = String::new();
        let mut nins = 0i32;
        let mut p = 0usize;
        for j in 0..len {
            while at(&realn[i], p) == b'-' {
                p += 1;
            }
            if gapped[j] == b'-' {
                if status != 1 {
                    status = 1;
                    out.extend_from_slice(
                        format!("{}{} - ", j + 1, addbk[i][j] as char).as_bytes(),
                    );
                    if nins == 0 {
                        insstr = format!("{}{}-", j + 1, addbk[i][j] as char);
                    } else if nins == 1 {
                        insstr.push_str("etc,");
                    }
                    nins += 1;
                }
            } else {
                if status == 1 {
                    out.extend_from_slice(
                        format!("{}{} > {}v{}\n", j, addbk[i][j - 1] as char, p, p + 1).as_bytes(),
                    );
                    if nins == 1 {
                        insstr.push_str(&format!("{}{},", j, addbk[i][j - 1] as char));
                    }
                }
                status = 0;
                p += 1;
            }
        }
        if status == 1 {
            out.extend_from_slice(
                format!("{}{} > {}v{}\n", len, addbk[i][len - 1] as char, p, p + 1).as_bytes(),
            );
            if nins == 1 {
                insstr.push_str(&format!("{}{},", len, addbk[i][len - 1] as char));
            }
        }

        // `insstr[strlen(insstr)-1] = 0`: drop the trailing comma. There is
        // always at least one run here, so there is always one to drop.
        insstr.pop();
        let mut newname: Vec<u8> = name[i][..off].to_vec();
        newname.extend_from_slice(format!("{}ins:{}|", nins, insstr).as_bytes());
        newname.extend_from_slice(&name[i][off..]);
        // `newname[B] = 0`.
        newname.truncate(B);
        name[i] = newname;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(rows: &[&str]) -> Vec<Vec<u8>> {
        rows.iter().map(|r| r.as_bytes().to_vec()).collect()
    }
    fn s(rows: &[Vec<u8>]) -> Vec<String> {
        rows.iter()
            .map(|r| String::from_utf8(r.clone()).unwrap())
            .collect()
    }

    #[test]
    fn record_original_gaps_marks_all_gap_columns() {
        let a = v(&["A-C-", "A--Z"]);
        assert_eq!(record_original_gaps(2, &a), b"o-oo".to_vec());
    }

    #[test]
    fn record_original_gaps_sizes_from_the_first_row() {
        let a = v(&["A-C", "A--", "AB-"]);
        assert_eq!(record_original_gaps(3, &a).len(), 3);
    }

    /// The reference keeps exactly the width it came in with, which is the
    /// whole point of the option.
    #[test]
    fn the_reference_width_survives() {
        // Two reference rows 4 wide; the alignment has widened them to 6 by
        // inserting two columns that only the added sequence occupies.
        let mut seq = v(&["AB--CD", "AB--CE", "--XY--"]);
        let original = record_original_gaps(2, &v(&["ABCD", "ABCE"]));
        let (ndel, list) = delete_new_insertions_whole(2, 1, &mut seq);
        assert_eq!(ndel, 2);
        assert_eq!(s(&seq), vec!["ABCD".to_string(), "ABCE".into(), "----".into()]);
        // Both inserted residues went, as one run starting at residue 0.
        assert_eq!(list[0], vec![GapPos { pos: 0, len: 2 }]);
        restore_original_gaps(3, &mut seq, &original);
        assert_eq!(seq[0].len(), 4);
    }

    #[test]
    fn nothing_is_deleted_when_the_reference_gained_no_columns() {
        let mut seq = v(&["ABCD", "ABCE", "A-XD"]);
        let before = s(&seq);
        let (ndel, list) = delete_new_insertions_whole(2, 1, &mut seq);
        assert_eq!(ndel, 0);
        assert_eq!(s(&seq), before);
        assert!(list[0].is_empty());
    }

    #[test]
    fn a_run_split_by_a_gap_is_two_entries() {
        // Columns 1 and 3 are inserted; the added sequence has a residue in
        // each, separated by a column it is gapped in.
        let mut seq = v(&["A-B-C", "A-B-C", "-X-Y-"]);
        let (ndel, list) = delete_new_insertions_whole(2, 1, &mut seq);
        assert_eq!(ndel, 2);
        assert_eq!(
            list[0],
            vec![GapPos { pos: 0, len: 1 }, GapPos { pos: 1, len: 1 }]
        );
    }

    #[test]
    fn restore_original_gaps_is_the_inverse_of_dropping_them() {
        let full = v(&["A-C-D", "A--ZD"]);
        let map = record_original_gaps(2, &full);
        assert_eq!(map, b"o-ooo".to_vec());
        // Only column 1 is all-gap, so only column 1 comes out.
        let mut squeezed = v(&["AC-D", "A-ZD"]);
        restore_original_gaps(2, &mut squeezed, &map);
        assert_eq!(s(&squeezed), s(&full));
    }

    #[test]
    fn name_offset_skips_the_sentinel_and_the_markers() {
        assert_eq!(name_offset(b"=seq1", false), 1);
        assert_eq!(name_offset(b"=_numo_s_00000001_numo_e_seq1", true), 25);
        assert_eq!(name_offset(b"=x_oe_seq1", false), 6);
    }

    #[test]
    fn the_compact_map_renames_and_summarises() {
        let addbk = vec![b"ABCDE".to_vec()];
        let deletelist = vec![vec![GapPos { pos: 1, len: 2 }]];
        // The added sequence as the reference alignment holds it.
        let realn = vec![b"A--DE".to_vec()];
        let mut name = vec![b"=s1".to_vec()];
        let out = reconstruct_delete_map_compact(&addbk, &deletelist, &realn, &mut name, false);
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("# Insertion in added sequence > Position in reference\n"));
        assert!(text.contains(">s1\n"), "{text}");
        // One run: it opens at residue 2 (1-origin) and closes at 3.
        assert!(text.contains("2B - "), "{text}");
        assert_eq!(String::from_utf8(name[0].clone()).unwrap(), "=1ins:2B-3C|s1");
    }

    #[test]
    fn the_compact_map_skips_sequences_that_lost_nothing() {
        let addbk = vec![b"ABC".to_vec()];
        let deletelist = vec![vec![]];
        let realn = vec![b"ABC".to_vec()];
        let mut name = vec![b"=s1".to_vec()];
        let out = reconstruct_delete_map_compact(&addbk, &deletelist, &realn, &mut name, false);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "# Insertion in added sequence > Position in reference\n"
        );
        // ... and leaves its name alone.
        assert_eq!(name[0], b"=s1".to_vec());
    }

    #[test]
    fn the_full_map_has_one_line_per_residue() {
        let addbk = vec![b"ABCDE".to_vec()];
        let deletelist = vec![vec![GapPos { pos: 1, len: 2 }]];
        let realn = vec![b"A--DE".to_vec()];
        let name = vec![b"=s1".to_vec()];
        let out = reconstruct_delete_map(&addbk, &deletelist, &realn, &name, false);
        let text = String::from_utf8(out).unwrap();
        let body: Vec<&str> = text.lines().skip(2).collect();
        assert_eq!(body.len(), 5);
        // The third field is the *column* the residue occupies in the
        // reference alignment, gaps counted -- so D and E report 4 and 5,
        // not the 2 and 3 they would if it counted residues.
        assert_eq!(body[0], "A, 1, 1");
        assert_eq!(body[1], "B, 2, -");
        assert_eq!(body[2], "C, 3, -");
        assert_eq!(body[3], "D, 4, 4");
        assert_eq!(body[4], "E, 5, 5");
    }
}
