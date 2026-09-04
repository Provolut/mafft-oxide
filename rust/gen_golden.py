#!/usr/bin/env python3
"""Generate `tests/golden/matrix/` -- the bit-identity fixtures `cargo test` replays.

Every case runs the real C binary once here and the port once per `cargo test`.
The inputs are deliberately tiny: the whole matrix has to replay in seconds so
it can live in `cargo test` rather than in a nightly job.

    ./gen_golden.py <path-to-core-dir> [outdir]

Cases that need extra files in the working directory (`_externalanchors`,
`_guidetree`, `_aamtx`, a sub-alignment table) carry them as
`<key>.stage.<filename>` fixtures; the test copies them back before running.

Cases where the C exits non-zero are kept, with `expect=exit`: the port has to
refuse them too. That is where the refusal paths get their coverage, including
the three upstream crashes in issue #4.

Regenerate whenever a feature lands and commit the result -- a fixture that
changes is the port changing its answer.
"""
import os
import random
import shutil
import subprocess
import sys
import tempfile
import zlib

AA = "ACDEFGHIKLMNPQRSTVWY"
NUC = "acgt"


def mutate(s, rate, rng):
    out = list(s)
    for _ in range(int(len(s) * rate)):
        p = rng.randrange(len(out))
        out[p] = rng.choice(AA if out[p].isupper() else NUC)
    return "".join(out)


def corpus(rng):
    """Six inputs, each aimed at a different branch of the C."""
    c = {}
    base = "".join(rng.choice(AA) for _ in range(80))
    c["aa"] = [("s%d" % i, mutate(base, 0.15 * (i + 1), rng)) for i in range(6)]

    nbase = "".join(rng.choice(NUC) for _ in range(90))
    c["nuc"] = [("n%d" % i, mutate(nbase, 0.10 * (i + 1), rng)) for i in range(6)]

    # Ragged lengths: the terminal-gap handling and `getterminalmargin`.
    c["ragged"] = [("r%d" % i, mutate(base, 0.2, rng)[: 80 - 12 * i]) for i in range(5)]

    # Duplicates: every tie-break in the tree build and the tracebacks.
    d = mutate(base, 0.1, rng)
    c["dup"] = [("d%d" % i, d) for i in range(4)] + [("d4", mutate(base, 0.3, rng))]

    # Tandem repeats: the FFT anchor scan actually has to prune.
    unit = "".join(rng.choice(AA) for _ in range(30))
    c["tandem"] = [
        ("t%d" % i, "".join(mutate(unit, 0.1, rng) for _ in range(4 - (i % 2))))
        for i in range(5)
    ]

    # Already aligned, with gaps: --add, --merge and --keeplength read this.
    g = []
    for i in range(5):
        s = list(mutate(base, 0.15, rng))
        for p in range(10 + 3 * i, 16 + 3 * i):
            s[p] = "-"
        g.append(("g%d" % i, "".join(s)))
    c["gapped"] = g
    return c


def write_fasta(path, recs):
    with open(path, "w") as f:
        for n, s in recs:
            f.write(">%s\n%s\n" % (n, s))


D = "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
D1 = "-q 0 -E 1 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0"
TAIL = "-X 0.1 -x 100"
DF = f"{D} -F {TAIL}"

# ---------------------------------------------------------------------------
# Staging helpers. Each returns a dict of {filename: contents}; the generator
# writes them into the C's working directory and saves them as fixtures.


def stage_anchors(recs, rng, conflicting):
    lens = [len(s.replace("-", "")) for _, s in recs]
    n = len(recs)
    lines = ["# generated anchors"]
    width = 16
    for _ in range(6 if conflicting else 3):
        i = rng.randrange(n - 1)
        j = rng.randrange(i + 1, n)
        if min(lens[i], lens[j]) < width + 4:
            continue
        si = rng.randrange(1, lens[i] - width)
        sj = si if not conflicting else rng.randrange(1, lens[j] - width)
        lines.append(
            "%d %d %d %d %d %d %d"
            % (i + 1, j + 1, si, si + width - 1, sj, sj + width - 1, rng.randrange(1, 100))
        )
    return {"_externalanchors": "\n".join(lines) + "\n"}


HERE = os.path.dirname(os.path.abspath(__file__))


def newick(n, shape):
    """A guide tree in the C's `_guidetree` format.

    Delegated to `gen_treein.py` rather than reimplemented, so these fixtures
    and `verify.sh` feed the binaries the same trees.
    """
    out = subprocess.run(
        [sys.executable, os.path.join(HERE, "gen_treein.py"), str(n), shape],
        stdout=subprocess.PIPE,
        check=True,
    )
    return {"_guidetree": out.stdout.decode()}


def subalignments(aligned_path, spec):
    """`--merge`'s `_subalignmentstable`, via the generator `verify.sh` uses."""
    out = subprocess.run(
        [sys.executable, os.path.join(HERE, "gen_merge.py"), aligned_path] + spec.split(),
        stdout=subprocess.PIPE,
        check=True,
    )
    return {"_subalignmentstable": out.stdout.decode()}


def stage_keyword(kw):
    return {"_guidetree": kw + "\n"}


def stage_aamtx(rng):
    """A permuted-order user matrix, as `--aamatrix` reads it."""
    order = list(AA)
    rng.shuffle(order)
    out = [" ".join(order)]
    vals = {}
    for i, a in enumerate(order):
        for b in order[i:]:
            vals[(a, b)] = rng.randrange(-4, 12)
    for i, a in enumerate(order):
        row = [a] + ["%d" % vals[(a, b)] if (a, b) in vals else "%d" % vals[(b, a)] for b in order]
        out.append(" ".join(row))
    return {"_aamtx": "\n".join(out) + "\n"}


# ---------------------------------------------------------------------------
# The case matrix. Entries are (tag, argv, stage) or the long form
# (tag, argv, stage, applies_to, port_refuses):
#
#   stage         callable (recs, rng, input_path) -> {filename: contents}, or None
#   applies_to    predicate on the input's records; the case is skipped when false
#   port_refuses  the C succeeds but the port declines on purpose, so the
#                 fixture records `refuse` and the test requires exit != 0
ANY = lambda recs: True
ALIGNED = lambda recs: len({len(s) for _, s in recs}) == 1


def distt_cases(rng):
    """(tag, argv, stage) for `disttbfast`."""
    c = [
        ("nwns2", f"{D} {TAIL}", None),
        ("nwns1", f"{D1} {TAIL}", None),
        ("fftns2", DF, None),
        ("fftns1", f"{D1} -F {TAIL}", None),
        ("forcefft", f"{D} -G {TAIL}", None),
        ("nofreegap", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, None),
        ("unweighted", f"{D} {TAIL} -u", None),
        ("noscore", f"{D} {TAIL} -n", None),
        ("gapext", "-q 0 -E 2 -V -2.00 -s 0.0 -W 6 -O -C 0-0 -b 62 -g -0.05 -f -2.00 -Q 100.0 -h -0.1 -X 0.0 -x 100", None),
        ("dynmtx", f"{DF} -t", None),
        ("fftwin40", f"{D} -F -w 40 {TAIL}", None),
        ("anysymbol", f"{D} {TAIL} -:", None),
        ("blosum30", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 30 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, None),
        ("blosum80", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 80 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, None),
        ("jtt", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -j 200 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, None),
        ("fmodel", f"{D} {TAIL} -a", None),
        ("legacygap", f"{D} {TAIL} -L", None),
        ("legacygapfft", f"{DF} -L", None),
        ("distout", f"{DF} -y h", None),
        ("distonly", f"{DF} -y h -T", None),
        ("tuple6", f"{DF} -A 6", None),
        ("memsave", f"{D} {TAIL} -M -B", None),
        ("memsavefft", f"{DF} -M -B", None),
        ("aalign", f"{D} {TAIL} -&", None),
        ("aalignfft", f"{DF} -&", None),
        ("nomemsave", f"{D} {TAIL} -N", None),
        ("allowshift", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 62 -g 0 -f -1.53 -Q 100.0 -h 0 -X 0.1 -x 4", None),
        # -r: the C segfaults past the first node (issue #4); the port refuses.
        ("oneiteration", f"{DF} -r", None),
        # External anchors, consistent and deliberately conflicting.
        ("anchors", f"{D} {TAIL} -l", lambda r, g, p: stage_anchors(r, g, False)),
        ("anchorsx", f"{D} {TAIL} -l", lambda r, g, p: stage_anchors(r, g, True)),
        # A user guide tree, and the keywords that generate one instead.
        ("treeinlad", f"{DF} -U", lambda r, g, p: newick(len(r), "ladder")),
        ("treeinbal", f"{DF} -U", lambda r, g, p: newick(len(r), "balanced")),
        ("treeinbalt", f"{DF} -U -t", lambda r, g, p: newick(len(r), "balanced")),
        ("pileup", f"{DF} -U", lambda r, g, p: stage_keyword("pileup")),
        ("shuffle", f"{DF} -U", lambda r, g, p: stage_keyword("shuffle")),
        ("shuffle42", f"{DF} -U", lambda r, g, p: stage_keyword("shuffle 42")),
        ("nodepair", f"{DF} -U", lambda r, g, p: stage_keyword("nodepair")),
        ("autotree", f"{DF} -U", lambda r, g, p: stage_keyword("auto 0 4")),
        ("testtree", f"{DF} -U", lambda r, g, p: stage_keyword("test 1 4")),
        ("compactbig", f"{DF} -U", lambda r, g, p: stage_keyword("compact 10GB")),
        ("compactsm", f"{DF} -U", lambda r, g, p: stage_keyword("compact 3")),
        ("verycompact", f"{DF} -U", lambda r, g, p: stage_keyword("very compact")),
        ("stepadd", f"{DF} -U", lambda r, g, p: stage_keyword("stepadd")),
        ("youngest", f"{DF} -U", lambda r, g, p: stage_keyword("youngestlinkage")),
        # A user substitution matrix. `nblosum < 0` forces protein, so this runs
        # on the nucleotide input too and must still agree.
        ("aamatrix", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b -1 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, lambda r, g, p: stage_aamtx(g)),
        # --merge: sub-alignments whose columns are fixed.
        # --merge requires every sub-alignment to be aligned already, so these
        # only run on inputs whose rows are all the same length.
        ("merge2", f"{DF} -H 0", lambda r, g, p: subalignments(p, "2"), ALIGNED, False),
        ("merge3", f"{DF} -H 0", lambda r, g, p: subalignments(p, "3"), ALIGNED, False),
        ("mergeint", f"{DF} -H 0", lambda r, g, p: subalignments(p, "2 interleave"), ALIGNED, False),
        ("mergekeep", f"{DF} -H 0", lambda r, g, p: subalignments(p, "2 interleave preservegaps"), ALIGNED, False),
        # Configurations the port refuses up front: coverage for params.rs.
        # The port refuses threads rather than pretending -C 0-0: declared, not
        # discovered, so a future thread implementation has to update this too.
        ("threads", f"{DF} -C 2-2", None, ANY, True),
        ("badblosum", "-q 0 -E 2 -V -1.53 -s 0.0 -W 6 -O -C 0-0 -b 47 -g 0 -f -1.53 -Q 100.0 -h 0 " + TAIL, None),
    ]
    return c


TB_P = "-u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1"
TB_PG = "-u 0.0 -l 2.7 -C 0 -b 62 -g -0.10 -f -2.00 -Q 100.0 -h 0.1"
TB_PE = "-u 0.0 -l 2.7 -C 0 -b 62 -g 0.0 -f -2.00 -Q 100.0 -h 0.0 -O -6.00 -E -0.000 -N -Z"
TB_PEO = "-u 0.0 -l 2.7 -C 0 -b 62 -g -0.100 -f -2.00 -Q 100.0 -h 0.1 -O -6.00 -E -0.000 -N"
TB_T = "-W 0.00001 -V -1.53 -s 0.0 -C 0 -b 62 -f -1.53 -Q 100.0 -h 0 -F -l 2.7 -X 0.1"
TB_TE = "-W 0.00001 -V -1.53 -s 0.0 -O -C 0 -b 62 -f -1.53 -Q 100.0 -h 0.000 -F -l 2.7 -X 0.1"


def tbfast_cases(rng):
    return [
        ("linsi", f"_ {TB_P} -L _ -+ 2 {TB_T} -O", None),
        ("linsi1", f"_ {TB_P} -L _ -+ 0 {TB_T} -O", None),
        ("linsi16", f"_ {TB_P} -L _ -+ 16 {TB_T} -O", None),
        ("ginsi", f"_ {TB_PG} -A _ -+ 2 {TB_T}", None),
        ("ginsi1", f"_ {TB_PG} -A _ -+ 0 {TB_T}", None),
        ("einsi", f"_ {TB_PE} _ -+ 2 {TB_TE}", None),
        ("einsiold", f"_ {TB_PEO} _ -+ 2 {TB_TE}", None),
        ("linsidyn", f"_ {TB_P} -L _ -+ 2 {TB_T} -O -t", None),
        ("ginsiweight", f"_ {TB_PG} -A _ -+ 2 -W 0.001 -V -1.53 -s 0.0 -C 0 -b 62 -f -2.00 -Q 100.0 -h -0.1 -F -l 0.5 -X 0.5", None),
        ("linsitree", f"_ {TB_P} -L _ -+ 2 {TB_T} -O -U", lambda r, g, p: newick(len(r), "balanced")),
    ]


def run(binary, args, workdir, stdin_path):
    with open(stdin_path, "rb") as fh:
        p = subprocess.run(
            [binary] + args,
            stdin=fh,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            cwd=workdir,
        )
    return p.returncode, p.stdout


def main():
    core = os.path.abspath(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    out = (
        os.path.abspath(sys.argv[2])
        if len(sys.argv) > 2
        else os.path.join(here, "tests", "golden", "matrix")
    )
    if os.path.isdir(out):
        shutil.rmtree(out)
    os.makedirs(out)

    rng = random.Random(20260904)
    cases = corpus(rng)
    for name, recs in sorted(cases.items()):
        write_fasta(os.path.join(out, "%s.fa" % name), recs)

    # --add and --addprofile take a combined input: an existing alignment plus
    # the sequences to add. Build those by aligning a corpus entry with the C
    # and splicing the result, then treat each as an input in its own right.
    # `nadd` varies per input, so these carry their own case list.
    add_inputs = {}
    for inp in ("aa", "gapped", "ragged"):
        work = tempfile.mkdtemp()
        rc, aligned = run(
            os.path.join(core, "disttbfast"), (DF).split(), work, os.path.join(out, "%s.fa" % inp)
        )
        shutil.rmtree(work, ignore_errors=True)
        if rc != 0 or not aligned:
            continue
        apath = os.path.join(out, "%s.aligned.tmp" % inp)
        with open(apath, "wb") as fh:
            fh.write(aligned)
        n = aligned.count(b">")
        for nadd in (1, max(1, n // 2)):
            for mode in ("add", "profile"):
                name = "%s_%s%d" % (inp, mode, nadd)
                gen = subprocess.run(
                    [sys.executable, os.path.join(HERE, "gen_add.py"), apath, str(nadd), mode],
                    stdout=subprocess.PIPE,
                )
                if gen.returncode != 0 or not gen.stdout:
                    continue
                with open(os.path.join(out, "%s.fa" % name), "wb") as fh:
                    fh.write(gen.stdout)
                add_inputs[name] = (nadd, mode)
        os.remove(apath)

    manifest = []
    for inp in sorted(cases):
        recs = cases[inp]
        src = os.path.join(out, "%s.fa" % inp)
        for binname, caselist, outputs in (
            ("disttbfast", distt_cases(rng), ("order", "hat2", "infile.tree")),
            ("tbfast", tbfast_cases(rng), ("pre", "order", "hat2")),
        ):
            for case in caselist:
                tag, argstr, stage = case[0], case[1], case[2]
                applies = case[3] if len(case) > 3 else ANY
                port_refuses = case[4] if len(case) > 4 else False
                if not applies(recs):
                    continue
                key = "%s.%s" % (inp, tag)
                work = tempfile.mkdtemp()
                staged = []
                if stage is not None:
                    # One rng per case, seeded from the key, so regenerating is
                    # deterministic regardless of case order.
                    for fname, content in stage(recs, random.Random(zlib.crc32(key.encode())), src).items():
                        with open(os.path.join(work, fname), "w") as fh:
                            fh.write(content)
                        with open(os.path.join(out, "%s.stage.%s" % (key, fname)), "w") as fh:
                            fh.write(content)
                        staged.append(fname)
                rc, stdout = run(
                    os.path.join(core, binname), argstr.split(), work, src
                )
                if rc != 0:
                    manifest.append((key, binname, argstr, "", ",".join(staged), "exit"))
                    shutil.rmtree(work, ignore_errors=True)
                    continue
                if binname == "disttbfast":
                    with open(os.path.join(out, key + ".fa"), "wb") as fh:
                        fh.write(stdout)
                extras = []
                for e in outputs:
                    p = os.path.join(work, e)
                    if os.path.exists(p) and os.path.getsize(p):
                        with open(p, "rb") as fh:
                            data = fh.read()
                        with open(os.path.join(out, "%s.%s" % (key, e)), "wb") as g:
                            g.write(data)
                        extras.append(e)
                manifest.append(
                    (
                        key,
                        binname,
                        argstr,
                        ",".join(extras),
                        ",".join(staged),
                        "refuse" if port_refuses else "ok",
                    )
                )
                shutil.rmtree(work, ignore_errors=True)

    # `-K` is --add proper (the tail is ungapped); without it the tail is an
    # aligned profile. `-Y` (--keeplength) cuts the columns the alignment
    # inserted into the reference back out, and -Z/-z say what that cost.
    for name in sorted(add_inputs):
        nadd, mode = add_inputs[name]
        src = os.path.join(out, "%s.fa" % name)
        kflag = " -K" if mode == "add" else ""
        for tag, extra in (
            ("plain", ""),
            ("keeplength", " -Y"),
            ("mapout", " -Y -Z"),
            ("mapoutc", " -Y -z"),
            ("mapoutn", " -Y -Z -n"),
        ):
            key = "%s.%s" % (name, tag)
            argstr = "%s -I %d%s%s" % (DF, nadd, kflag, extra)
            work = tempfile.mkdtemp()
            rc, stdout = run(os.path.join(core, "disttbfast"), argstr.split(), work, src)
            if rc != 0:
                manifest.append((key, "disttbfast", argstr, "", "", "exit"))
                shutil.rmtree(work, ignore_errors=True)
                continue
            with open(os.path.join(out, key + ".fa"), "wb") as fh:
                fh.write(stdout)
            extras = []
            for e in ("order", "_deletelist", "_deletemap"):
                p = os.path.join(work, e)
                if os.path.exists(p) and os.path.getsize(p):
                    with open(p, "rb") as fh:
                        data = fh.read()
                    with open(os.path.join(out, "%s.%s" % (key, e)), "wb") as g:
                        g.write(data)
                    extras.append(e)
            manifest.append((key, "disttbfast", argstr, ",".join(extras), "", "ok"))
            shutil.rmtree(work, ignore_errors=True)

    with open(os.path.join(out, "MANIFEST"), "w") as fh:
        for row in manifest:
            fh.write("\t".join(row) + "\n")
    ok = sum(1 for r in manifest if r[5] == "ok")
    print("wrote %d cases (%d ok, %d refused) to %s" % (len(manifest), ok, len(manifest) - ok, out))


if __name__ == "__main__":
    main()
