# Optimisation attempts, 2026-09-04

The run GOAL.md describes: beat the port as it stood at `77dd5a6` by 10% on
**both** wall clock and peak RSS, or spend 90 minutes and write down what was
tried. This is that write-up. The target was not met; what follows is where
the time goes and which stones are now turned.

Output did not change at any point. Every kept change was measured with
`bench/profile.sh`, which refuses to record a row unless both binaries exit 0
and produce identical bytes.

> **This file is chronological, and a stale entry is not a contradiction.**
> Entries are appended and dated; an estimate or a ranking is superseded by a
> later entry that names and refutes it, never by editing the earlier one. So
> an item that looks wrong -- `## What is left` still ranks the scan fusion
> first at "perhaps 15%", which the next section measures and rejects -- is
> answered by reading forward, not by correcting it in place. That edit was
> drafted once and stopped only by reading on. It is worth saying why it would
> have been worse than the inconsistency it removed: it would have deleted the
> evidence that the estimate was tested and refuted, which is unrecoverable
> from the tree afterwards, and it would have left the document looking *more*
> coherent than the work it records. A log that reads cleanly because its
> refuted predictions were tidied away is no longer evidence of anything.

> **The 90 minutes here are an estimate, not a record.** GOAL.md asks for start
> and finish times in the write-up and they were never taken; the figure is the
> author's own tally of effort spent on optimisation attempts. Stated plainly
> because a budget nobody timed cannot be audited, and the honest fix after the
> fact is to say so rather than to reconstruct plausible timestamps. A later run
> against this clause should take the two readings at the time.

## Result

Not met. Both sides re-measured with the corrected harness (see *Measurement*),
`77dd5a6` rebuilt from a worktree so the two are measured identically:

| metric | baseline `77dd5a6` | now | target (<= 90%) |
|---|---|---|---|
| wall clock, 8 cases | 2372 ms | **2361 ms** (99.5%) | <= 2135 ms |
| mean peak RSS | 6084 KB | **6144 KB** (101.0%) | <= 5476 KB |

Per case, port against port:

```
                    instructions                       peak RSS        wall ms
sample.fftns2   1107938947 -> 1103611601  -0.4%      8864 ->  8940     112 -> 110
sample.nwns2    1515694127 -> 1515693464  -0.0%      8436 ->  8440     141 -> 138
sample.memsave  8950539338 -> 8283448864  -7.5%      5036 ->  5104    1096 -> 1078
sample.linsi    7606390975 -> 7606562578  +0.0%     13440 -> 13468    1001 -> 1014
aa.fftns2         16584195 ->   16545652  -0.2%      2956 ->  3348       5 ->   5
nuc.fftns2        19029295 ->   19006136  -0.1%      3364 ->  3316       6 ->   5
tandem.fftns2     19390787 ->   19222878  -0.9%      3268 ->  3288       6 ->   6
dup.nwns2         12795777 ->   12794559  -0.0%      3308 ->  3252       5 ->   5
```

One line of that is the whole story: **`sample.memsave` gave up 7.5% of its
instructions and 1.6% of its wall clock.** The work removed was real, and it
was nearly free to begin with. That case is not instruction-bound, and neither
is the total -- the port runs 1.33x the C's instructions at 1.015x its wall
clock.

Two of the eight cases carry the workload: `sample.memsave` (1078 ms) and
`sample.linsi` (1014 ms) are 88% of the total. The other six are 5-138 ms and
dominated by process startup, so nothing done to them can move the number. The
memory figure is `sample.linsi`'s: at 13.4 MB it is more than twice the mean,
and it did not move at all this session.

Against the frozen `bench/baseline-2026-09-04.tsv`, which was recorded with the
old harness, the same numbers read 99.7% and 105.9%. Its RSS column is not
comparable with anything measured since; see *Measurement*.

## What was tried

### Kept

1. **`prf1`/`prf2` hoisted out of the per-position loop** (`falign.rs`). The
   port allocated two `vec![0.0; na]` per position where the C re-zeroes one
   pair of buffers. Replaced by two buffers refilled with `fill(0.0)`.
   Allocations on `sample.memsave`: 625,467 -> 555,921. Instructions on that
   case: 8.95 G -> 8.28 G, a 7.5% cut.

2. **Lag-loop slice vectors hoisted** (`falign.rs`). The `Vec<&[u8]>` for each
   group was collected fresh per lag; now built once and `clear()`/`extend()`ed.
   Part of the same measurement as (1).

3. **Bounds check folded in the 26-element profile scan** (`falign.rs`).
   Slicing to `[..26]` up front lets the check fold. Worth **0%** wall clock,
   even though that scan is 21% of the case's instructions -- the checks were
   predicted perfectly. Kept only because it costs nothing.

4. **`FalignState`'s three `MAXSEG` arrays grow on demand.** They were
   allocated at 100000 entries x 24 B = 2.4 MB each regardless of input.
   First recorded as "no change", then as a 2.4 MB saving on `sample.linsi`.
   Both readings were artefacts of the measurement below. In a clean
   directory `sample.linsi` is unchanged by it to within 0.1% -- which is what
   should have been expected, since these arrays belong to `Falign`, the FFT
   path `disttbfast` uses, and the L-INS-i case never touches them. The
   allocation itself is still worth not making.

### Reverted or rejected

5. **`impmtx` walked as a rectangle** rather than per row. Measurably worse;
   reverted.

6. **`-C target-cpu=native`.** AVX-512 is emitted and the output stays
   byte-identical, and it is worth **0%** wall clock. That is not a surprise
   in hindsight: the DP inner loops carry an argmax with an index, which
   blocks vectorisation. `core/Lalign11.c` says exactly this in a comment
   above its own hand-written AVX-512 path -- and that path is compiled out in
   this build anyway, because `CFLAGS = -O3` never defines `__AVX512F__`. Both
   sides are running the scalar loop.

## Measurement, which was the real finding

Four conclusions in this session were wrong the first time, all from the same
metric: **peak RSS as `profile.sh` measured it was not measuring the program.**

The first symptom was spread. The same binary on the same input varied by 31%
between runs; `sample.linsi` clustered at 11.0, 13.4 and 14.4 MB while the C
sat inside 1%. On single samples that produced a change which only removes
allocations reading as a 7.6% memory *regression*, the same change reading as
"no difference", and a median compared against a single draw reading as 5.5%
worse. Sampling five times and taking the median fixed the spread but not the
cause.

The cause is the harness. Every measured run shared one working directory, and
`tbfast` and `disttbfast` write `hat2`, `hat3`, `order` and `pre` into the cwd
and read them back. So the second run of a case was not doing the work the
first one did, and each case inherited the previous case's files. Run in a
fresh empty directory, eleven repeats of `sample.linsi` land inside 1%:

```
port  14332 14356 14420 14420 14428 14432 14436 14448 14452 14468 14472
C      9160  9164  9176  9176  9176  9196  9216  9224  9236  9248  9340
```

`profile.sh` now gives every run -- callgrind, RSS poll and each timing
iteration -- its own directory. The wall-clock column was measured the same
polluted way, so it is re-measured too, on both sides.

`summarize.py` had a matching hole: a baseline whose rows all failed to
measure divided by zero and printed **TARGET MET**. It fired for real, on a
re-measurement pointed at a tree with no compiled C. It now refuses to compare
when either side has no rows marked `ok`, or when the two sides cover a
different number of cases.

The frozen `bench/baseline-2026-09-04.tsv` predates all of this. **It has been
left in place** rather than rewritten: re-measuring it raises the RSS target,
which is a change in my own favour, and that is not mine to make. The
re-measured file sits beside it and the comparison is reported against both.
Which one is the yardstick is a call for a human.

## One more thing the gates caught

The push was refused by `bench/coverage.sh`, which runs the suite a second
time. `blosum`, `io` and `tree` each had a test helper that swaps the
process-wide working directory under a `static LOCK` of its own -- four locks,
one directory, several test threads. About one run in eight, two helpers
interleaved and a test ran in the other's directory. Fixed with a single
`crate::cwd_lock()` taken by all four; no assertion changed. 20 consecutive
instrumented runs clean. Coverage 89.86% regions, 89.36% lines.

Worth recording here because the flake was invisible to the fast gate for as
long as only one run happened per push.

## Leads not taken

- **Three `commonIP` buffers against the C's one.** `align.rs`, `genalign.rs`
  and `lalign.rs` each keep their own traceback matrix with its own
  `common_alloc1`/`common_alloc2` growth. In the C, `commonIP` is a single TLS
  global shared by every aligner (`dp.h:15`, written from `Lalign11.c`,
  `Salignmm.c`, `Dalignmm.c`, `MSalign11.c`, `Galign11.c` and others). At
  `sample.linsi`'s sizes one such matrix is about 2.4 MB, and the L-INS-i run
  has a pairwise stage and a progressive stage that use different ones.
  The buffer is pure scratch -- every cell a traceback reads was written by
  the same call -- so sharing it should be output-neutral, and it is what the
  C does, so it is more faithful rather than less. **Not measured**: it needs
  the buffer threaded through the three states' call sites, which is more
  surgery than the remaining time allowed.

- **The argmax inside the DP row.** Both the port and the C's scalar tail
  update `maxwm`/`endali`/`endalj` inside the inner loop, and that is what
  stops either from vectorising. The C's AVX-512 path lifts it into a
  per-row reduction; doing the same in safe Rust means either a second pass
  over stored pre-clamp scores (extra stores, extra traffic) or a
  hand-blocked reduction. Worth trying, and untried here.

- **`sample.nwns2` at 1.60x the C's instructions** is the worst compute ratio
  left, but only 141 ms absolute; even halving it moves the total by 3%.

## What the ratios say

Instructions overstate the gap consistently. On `sample.memsave` the port
executes 1.23x the C's instructions and runs at 0.98x its wall clock; on
`sample.linsi`, 1.42x instructions for 1.07x wall clock. The port's extra work
is cache-friendly and pipelines behind the loads the C is also waiting on.
Optimising against `Ir` alone would therefore chase work that is already free.


# The L-INS-i path, 2026-09-04 (second session)

A different objective from the one above: not "beat the port by 10%", but
**reach the C on the workload a caller actually runs**. That caller is the
GEMME/ESCOTT evolutionary trace in provolut-ml, whose every sub-alignment
realignment is a `mafft --auto` that lands on L-INS-i. `sample.linsi` is the
only case in the database on that path.

Output did not change at any point.

## Result

The port is now **faster than the C** on `sample.linsi`, and the case's
instruction ratio fell by a fifth:

| `sample.linsi` | before | after |
|---|---|---|
| instructions vs C | 1.378x | **1.142x** |
| wall clock vs C | 1.018x | **0.967x** |
| peak RSS vs C | 1.633x | 1.663x |

On the caller's own corpus -- 10 sub-alignments drawn from the ProteinGym
P53_HUMAN alignment, replayed through `tbfast` alone with the argument vector
and staged input the `mafft` script hands it, 5 repetitions each -- the port
went from **0.940x the C to 0.986x**, every output byte-identical. Through the
whole `mafft --auto` pipeline, where `tbfast` is 78% of the call and the other
six binaries are still the C's, 0.950x to 0.988x.

The remaining gap on that corpus and the win on `sample.linsi` are the same
port: the two workloads differ in shape (51 sequences of 393 residues against
the 36-sequence `test/sample`), and the case that gains most is the one with
the most DP cells per call.

Whole-database, same host, before against after: wall clock 1280 ms -> 1223 ms
(95.5%), mean peak RSS 6505 KB -> 6527 KB (100.3%, inside the harness's noise).
The five `disttbfast` cases are unmoved to within noise, as expected -- nothing
here touches a function they call.

## What was done

All three changes remove bounds checks, and nothing else. Callgrind put 14% of
`l_align11`'s instructions in `core/src/slice/index.rs` inlined into it; after
the three, that line is gone from the profile. No arithmetic was reassociated,
no comparison reordered, no `unsafe` introduced.

1. **The DP row walked as zipped slices** rather than indexed. Slicing the seven
   buffers to `lastj` up front -- what the code already did, with a comment
   saying the checks would fold -- was not enough: seven indexed accesses per
   cell each kept a check. Zipping proves the ranges once. 0.940x -> 0.951x on
   the corpus.

2. **The `mi`/`mpi` prefix scan zipped** the same way, three buffers. The
   recurrence stays serial -- `miscan` feeds the next position -- so only the
   checks go. 0.951x -> 0.970x.

3. **`match_calc`'s row typed as `&[f64; 0x100]`** instead of a slice. The
   residue index is a `u8`, so with the length in the type the check on every
   cell of the row folds away. 0.970x -> 0.986x.

Measured one at a time, in that order, each on the 10-payload corpus.

## What is left, and why it was not taken

`if maxwm < w2v` is now the single hottest line at 10.8% of the function's
instructions. It is not a Rust-vs-C gap: `core/Lalign11.c` runs the same
in-loop argmax in its scalar path, and the branch is well predicted because
`maxwm` only ever increases. The C lifts it into a per-row reduction *only*
inside `LALIGN11_AVX512`, which `CFLAGS = -O3` never defines, so both sides run
the same construct. Matching it is the right outcome here; beating it means
vectorising the row, which is the anti-diagonal project the previous session
scoped and declined.

`ijp` is still `Vec<Vec<i32>>` against the C's single contiguous `commonIP`
(the "three `commonIP` buffers" lead above). Untried this session: rows are
already contiguous within themselves, so the gain is limited to row
transitions, and it is the change that would move `sample.linsi`'s peak RSS,
which is the metric already 1.66x the C's.

## Measurement note

`bench/profile.tsv` was re-recorded on a different host from the one that wrote
the committed copy -- `sample.linsi` takes 549 ms here against 1014 ms there --
so **the absolute `c_ms` / `r_ms` columns are not comparable with the previous
revision of that file; the ratio columns are**, and it is the ratios that carry
the claim. `bench/before-lalign11-bounds-checks.tsv` is this host's run of the
unchanged port, so before-and-after on one machine is comparable throughout.

## Gates

`cargo test --release` green (553-case golden matrix included). The wide
harnesses, run against the compiled `core/`: `verify.sh` pass=12018 fail=0,
`verify_tbfast.sh` pass=1690 fail=0, `verify_e2e.sh` pass=29 fail=0.
`bench/coverage.sh`: 89.88% regions, 89.37% lines.

---

# Optimisation attempts, 2026-09-05 — `dvtditr-rs` on the ESCOTT path

The port had just gained `dvtditr`, and it was the one binary measurably slower
than the C: 0.459 s against 0.354 s on a staged L-INS-i replay, which showed up
end-to-end as ESCOTT's `mafft --auto --preservecase` running 7.3% slower than
the all-C build. This is that pass.

Constraints, all three from the caller: bit-identical output, no memory blow-up,
and no `unsafe` (the port keeps `#![forbid(unsafe_code)]`).

## Result

| metric | before | after | C |
|---|---|---|---|
| instructions, L-INS-i replay | 12.71 G | **9.86 G** (-22.4%) | 8.35 G |
| wall clock, same replay (mean of 5) | 0.459 s | **0.374 s** (-18.5%) | 0.355 s |
| peak RSS | 18.9 MB | **18.9 MB** | 19.2 MB |
| 12 ESCOTT payloads, end to end | +7.3% vs C | **+3.9% vs C** | — |

Memory did not move, and the port still peaks *below* the C.

## What worked

**1. The DP inner loop's bounds checks — `-15.9%` (12.50 G -> 10.51 G).** By far
the largest. `part_a_align`'s innermost loop indexed ten `Vec` fields of `self`
per cell while writing through `&mut self`, so the optimiser had to treat every
length as a memory read that the loop's own stores might change, and kept all
ten checks. Binding each to a local slice (`&self.fgcp2[..lastj]`, …) before the
loop makes the lengths loop-invariant SSA values; the checks fold away.
`part_a_align` went 4.03 G -> 2.05 G, which is where nearly all of it came from.

**2. The sentinel in `match_calc`'s sparse walk — `-6.3%` (10.51 G -> 9.86 G).**
The C terminates each column's non-zero list with `-1` and re-tests for it on
every step of the innermost loop in the port's hottest function. Filling the
per-column buffers with `clear` + `push` instead makes the length carry that
information, and the walk becomes a plain zip over two slices. `clear` keeps the
capacity `ensure` allocated, so nothing reallocates.

**3. `scarr` on the stack, and its index masked — part of `-1.6%`.** A `vec!`
per call (~128 thousand calls) became a `[f64; 0x100]`, and `scarr[code & 0xff]`
lets the compiler drop a bounds check whose branch can never be taken.

**4. `fillimp`'s matrix zeroing via `fill` — the rest of that `-1.6%`.** The
element-by-element loop did not lower to a memset; `row[..lgth2].fill(0.0)`
does.

## What did not work

**Hoisting the row slices and the invariant product in `fillimp`: +0.07%.**
`seq1[i]`/`seq2[j]` re-indexed per residue and `e.importance * effij`
recomputed per matched cell both *look* like free wins. LLVM had already done
both, and spelling them out cost 6.5 M instructions. Reverted — a comment
claiming a win the measurement denies is worse than the original code.

## Where the remaining 1.5 G lives

```
match_calc        3.02 G (30.7%)
init_strict       2.06 G (20.9%)   C's fillimp + movereg + movereg_swap: 1.39 G
part_a_align      2.05 G (20.8%)
memset            1.21 G (12.3%)   C: 1.26 G -- parity, not a gap
intergroup_score  0.57 G ( 5.8%)
```

Two stones left unturned, both plausible and neither tried:

- **Flatten `doublework` / `intwork`.** They are `Vec<Vec<_>>`, so the hot loop
  chases a pointer and checks a bound once per *column* — 92 million column
  visits on this replay — before it reads any element. One flat buffer plus an
  offset table would remove that and improve locality; the per-column summation
  order, which is what must not move, is unaffected.
- **`init_strict` is 0.67 G over the C's equivalent.** Not yet localised beyond
  the function; `movereg` and the `impmtx[k1][k2]` double index are the
  candidates.

## Output did not change at any point

Every step was checked against the C before the next one started: the staged
L-INS-i replay byte-for-byte, then at the end `cargo test` 221, `verify_flh.sh`
98/98 plus its 60-call sequence, `verify_td.sh` 96/96, `verify_tdi.sh` 126/126,
`verify_dvtditr.sh` 32/32, `verify.sh`, `verify_tbfast.sh`, the 120-payload
ProteinGym corpus, and provolut-ml's GEMME/ESCOTT matrix test. `match_calc` is
shared with `tbfast` and `disttbfast`, so those two full differentials are not
optional here.

## Where ESCOTT's time actually is, and why `l_align11` resisted

Timing each binary of one `mafft --auto --preservecase` payload (51 sequences,
393 aa -- what the GEMME/ESCOTT trace hands MAFFT):

```
tbfast     820 ms
dvtditr    202 ms
version / f2cl / countlen   1 ms each
```

`tbfast` is 80% of the run, and **86% of `tbfast` is `l_align11`** -- so that
one function is roughly 69% of ESCOTT's whole aligner cost. Which makes the
`dvtditr` work above worth 3.9% end to end, and puts any further big saving in
`l_align11`.

It is already close to the wall. On the same input the C's scalar
`L__align11` retires 8.28 G instructions and the port's `l_align11` 8.97 G --
**8.4% apart**, with the port's loop already zipped (bounds checks gone),
branchless, and prefix-scanned.

**Tried: splitting the running maximum out of the elementwise body.**
`maxwm`/`endali`/`endalj` is a max *with an index*, which no vectoriser can
express, and it sat in the middle of an otherwise elementwise loop. Recording
each cell's score to a scratch row and reducing in a second pass removed it.
Result: instructions fell 4.6% (10.44 G -> 9.96 G on `tbfast`) and **wall clock
rose 1.3%** (0.782 s -> 0.791 s, best of 5, same input). The disassembly
explains the rest: still **zero packed arithmetic ops** -- LLVM did not
vectorise even with the reduction gone, so the change bought nothing but the
extra row's memory traffic. Reverted.

The lesson worth keeping: this loop is **memory/latency-bound, not
instruction-bound**. Instruction counts are the wrong metric here, which is why
the reverted change looked like a 4.6% win and was not one.

**So the remaining big saving in ESCOTT is vectorising `l_align11`, and only
explicitly.** That is what the C does -- `Lalign11.c`'s AVX-512 kernel behind
`#if defined(__AVX512F__) && defined(__AVX512VL__)`, which is the entire reason
`-march=native` runs 1.73x faster end to end (7.29 s vs 12.57 s over 12
payloads) while `-C target-cpu=native` buys the port 1.7%. The kernel's
operations are element-wise (`loadu_pd`, `cmp_pd_mask`, `mask_blend_pd`,
`add_pd`, `storeu_pd`) with no horizontal reduction, so a port of it would be
bit-exact by construction rather than by luck.

---

# The AVX-512 kernel for `l_align11`, 2026-09-05

The section above ended by saying the only big saving left was to vectorise
`l_align11` explicitly, as the C does. This is that, and it is the largest
single change the port has had.

## Result

One staged `tbfast` input (51 sequences, 393 aa — an ESCOTT realign payload),
best of five:

| build | time |
|---|---|
| `tbfast-rs`, scalar (`MAFFT_RS_NO_AVX512=1`) | 0.770 s |
| C `tbfast`, `-march=native` (its own AVX-512 kernel) | 0.342 s |
| **`tbfast-rs`, kernel** | **0.292 s** |

**2.63x over the port's scalar path, and 17% ahead of the C's own vectorised
build.** End to end over 12 ESCOTT payloads: C `-O3` 12.60 s, C `-march=native`
7.34 s, **port 6.86 s** — the port is now the fastest of the three, and unlike
the `-march=native` C it is a portable binary that decides at runtime.

Peak RSS is unchanged: alternating runs give 10.5–10.7 MB with the kernel on and
off alike, with no consistent direction. The kernel adds two 8-lane stack arrays
and nothing else.

## What it is

`src/lalign_avx512.rs`, a statement-for-statement transliteration of
`Lalign11.c`'s `LALIGN11_AVX512` block: eight cells at a time, `__m512d` for the
scores and `__m256i` for the traceback indices, with the per-lane running maximum
reduced once per row.

The C's own comment explains why neither compiler will do this automatically —
the loop mixes 64-bit scores with 32-bit indices, so there is no single
vectorisation factor. That matches what was measured here: LLVM emitted **zero**
packed arithmetic ops for this loop, and `-C target-cpu=native` bought 1.7% end
to end.

## Why it is bit-exact

Every lane runs the scalar computation on its own cell — the same comparisons
(`_CMP_GT_OQ` / `_CMP_LT_OQ`) on the same operands in the same order, the same
additions, no reassociation, and no horizontal sum anywhere. The only
cross-lane interaction is the running maximum, and the C's tie rule restores
scan order exactly: strict `>` keeps the earliest `j` within a lane, and the
smallest `j` wins across lanes.

That is an argument, not evidence, so `verify_lalign_kernel.sh` runs the same
alignments twice through one binary — once with `MAFFT_RS_NO_AVX512=1` — and
requires identical bytes. **This is the check the other harnesses structurally
cannot make**: they compare the port against the C, so they would pass if both
of the port's paths were wrong in the same way.

## The cost, stated plainly

The crate moved from `#![forbid(unsafe_code)]` to `#![deny(unsafe_code)]`, and
`lalign_avx512` is the one module that lifts it. `forbid` cannot be lifted
locally — that is the difference between the two lints — so this is the whole of
the licence, and it is visible in `lib.rs` rather than buried. Portable SIMD
(`std::simd`) would have avoided it, but it is still nightly-only with no
stabilisation date.


---

# The prefix scan, and two things dvtditr need not compute, 2026-09-06

With `L__align11`'s DP row vectorised, the profile changed shape enough that the
previous section's conclusions no longer held. Re-profiling found three items,
two of them algorithmic rather than instruction-level.

**Callgrind cannot see this work.** Valgrind masks AVX-512 out of `CPUID`, so
`lalign_avx512::available()` answers false under it and every profile shows the
scalar path — 85% in `l_align11`, which is the profile of a binary nobody runs.
Everything below is `perf record -F 4000` on the release build instead. If a
number here disagrees with `hotspots.sh`, `hotspots.sh` is profiling the
fallback.

## Where the time was

One staged `tbfast` payload (0.296 s) and one staged `dvtditr` payload (0.212 s),
an ESCOTT realign each:

| tbfast | | dvtditr | |
|---|---|---|---|
| `l_align11`, outside the kernel | 71% | `match_calc` | 27% |
| `lalign_avx512::row` | 13% | `init_strict` (via `fillimp`) | 24% |
| `a_align` | 7% | `part_a_align` | 19% |

## match_calc: 95% of a profile column is zeroes

`match_calc` built `scarr` by walking all 26 rows of `cpmx1`. A profile column is
how much of each residue type appears at one position, so it is nearly empty:
instrumented over the staged payload, **104125 of 2000024 entries are non-zero**,
1.35 of 26. A zero row contributes `scarr[l] += mtx[j][l] * 0.0` to every `l`.

Skipping it is exact, not approximately exact. `x + (±0.0)` is `x` for every `x`
except `-0.0` with a `+0.0` addend, and `scarr` cannot be `-0.0`: it starts at
`+0.0`, a sum is `-0.0` only when both addends are, and an exact cancellation
rounds to `+0.0` under round-to-nearest. So the skipped term could not have
changed a bit.

0.212 s → 0.191 s.

## fillimp: an answer that does not depend on the question

`map_range` answered "which column holds the n-th residue?" by counting from the
start of the sequence, and `fillimp` asked it twice per segment of every
constrained pair. The sequence does not change while it is being asked, so the
positions are computed once per member per call now.

The edges are the interesting part: a start past the last residue, a reversed
range, a negative start. None should arise from a real constraint, all are
reachable in the C, and each resolves as the final column for its own reason.
`tests::map_range_agrees` runs the old scan and the new lookup against random
gapped strings over the whole `-2..len + 2` range of endpoints — the check no
corpus can make, since the corpus only ever produces the ordinary case.

0.191 s → 0.167 s. Together: **dvtditr 1.27x.**

## The mi recurrence was 43% of tbfast

`l_align11` runs a serial recurrence down each row before the row itself:

```
miv[j] = m;  mpiv[j] = mp;
if( pw[j-1] > m ) { m = pw[j-1]; mp = j-1; }
m += fpenalty_ex;
```

Every step waits for the one before it — a compare plus an add of latency, about
eight cycles, for almost no work. That did not matter while the DP row cost more;
with the row vectorised it became the largest single item in the binary.
Measured directly by breaking the chain (`m = pw[idx] + e`, wrong answers, right
timing): **0.294 s → 0.167 s**. Replacing the triple `zip` with equal-length
slices changed nothing (0.296 → 0.294), which is what a latency-bound loop looks
like.

It is a max-plus scan, and those parallelise. Writing `d_s = p_s + (8-s)e` puts
every candidate on a common footing; a prefix maximum over `d` then answers eight
positions at once and only the carry between blocks stays serial.

**Two rewrites, one free and one conditional.** `max(a,b) + e == max(a+e, b+e)`
holds for every IEEE `a`, `b`, `e`, because rounding is monotone — so pushing the
`+ e` inside the maximum cannot move a value. Reaching a position with one
multiply, `m + t*e`, where the loop added `e` to itself `t` times, does not: they
agree only if nothing rounded. Whole numbers with sums under 2^53 give that, and
MAFFT's matrices and penalties are integers — but that is the input's property,
not the code's, so `l_align11` establishes it from the matrix and the penalties
and passes the answer down. When it fails, the serial loop runs.

The same condition settles the indices, which is the subtler half. Ties go to the
earliest position, and the two forms compare the same pair of candidates at
different points along the decay: if rounding could make two unequal values equal
later, they would disagree about which came first. Exact arithmetic keeps a
difference a difference.

One edge needed care. Position 0 must reproduce the incumbent untouched, and
`m + 0*e` does not: with a zero `e` the product is `+0.0`, which turns a `-0.0`
incumbent positive. Adding `-0.0` is the only identity that works for every
`f64`, and the unit test found this rather than the review.

**tbfast 0.296 s → 0.216 s.**

## Result

All three builds compiled and timed in one session on one machine, best of
seven per binary, outputs digested and identical:

| | C `-O3` | C `-march=native` | port before this | port |
|---|---|---|---|---|
| `tbfast` | 0.816 s | 0.363 s | 0.300 s | **0.216 s** |
| `dvtditr` | 0.187 s | 0.192 s | 0.210 s | **0.166 s** |
| 12 payloads, `mafft --auto` | 13.07 s | 7.75 s | 7.29 s | **5.63 s** |

**1.38x ahead of a C build that will not run on half the machines in the fleet**,
and 1.29x ahead of where the port stood before this round. The "port before
this" column is the AVX-512 kernel commit built from a worktree in the same
session, because comparing against a figure from an earlier day would have
measured the machine as much as the change: the C at `-march=native` alone
drifted 7.34 s to 7.75 s between the two.

Note the C at `-O3` beat the port's `dvtditr` until this round. That binary
never had a kernel to gain from -- `-march=native` does nothing for it either
(0.187 s against 0.192 s) -- so the two algorithmic wins above are the whole of
its lead.

## What is left

Measured, not guessed, and in the order they are worth doing:

1. **Fuse the scan into the row kernel** (tbfast). They walk `previousw` one
   after the other and the scan's output is consumed immediately by the row, so
   `miv`/`mpiv` need not be written to memory at all. Worth perhaps 15%; it is a
   real refactor, and the `trywarp` variant reads those arrays too.
2. **`part_a_align`** is 25% of `dvtditr` and has the same shape as
   `L__align11` — a serial `mi` with `+= fpenalty_ex`, then an element-wise row.
   The C has no kernel for it either. Hoisting the recurrence into a scan is the
   same transformation as above, with position-dependent `ogcp2[j]` coefficients
   in the candidate.
3. **`intergroup_score` (12%) and `alignable_reagion` (10%)** in `dvtditr` are
   unexamined; the two wins above were both algorithmic, so look for waste
   before reaching for intrinsics.
4. **Packing value and index into one `f64` key** in the scan would drop the
   `__m256i` half of the permutes, about 6 of 28 uops per block. It needs a
   second precondition (index under 2^21, values under 2^32) for a few percent.

---

# Fusing the scan into the row kernel does not pay, 2026-09-07

The previous section ranked this first, at "perhaps 15%". It was an estimate,
and measuring it refutes it: both variants are slower than the two separate
loops. Nothing shipped.

## Where ESCOTT's time actually is, measured through the caller

Earlier sections timed staged payloads. This one instruments the caller —
`pml.bioseq.gemme.pipeline`, whose `_run_aligner` is the single MAFFT
chokepoint — over a cold P53_HUMAN trace, 100 realign calls, Rust build, one
worker, no cache:

```
tbfast                       20.8 s   45%
dvtditr                      19.1 s   41%
mafft script + 5 binaries     3.5 s    8%
python trace + predict        2.6 s    6%
```

MAFFT is 96% of a cold GEMME/ESCOTT run and `tbfast` + `dvtditr` are 86% of
it. **`dvtditr` is 41% here, not the 35% the staged 12-payload figures imply**,
so the `part_a_align` item below was already worth more than this one's
ceiling.

## What was tried

One staged `tbfast` payload — 51 sequences of 393 aa, a real ESCOTT realign —
best of seven, each run in its own fresh directory, `pre` digested every time
and identical throughout.

| | wall | instructions | cycles | IPC |
|---|---|---|---|---|
| two loops (unfused) | **0.2089 s** | 3.610 G | 1.1247 G | 3.21 |
| fused, scan in step with the row | 0.2207 s (+5.7%) | 3.453 G (-4.4%) | 1.1876 G (+5.6%) | 2.91 |
| fused, scan one block ahead | 0.2136 s (+2.3%) | 3.620 G | 1.1511 G | 3.14 |

The block indices do line up exactly — with `jj = done + 1` the two loops read
the same eight values of `previousw`, the scan writes the `miv` slots the row
reads, and `done + 8 <= lastj - 1` is `jj + 8 <= lastj` — so the fusion is
sound, and its output is byte-identical on every check: kernel on/off 15/0,
`verify_tbfast` golden replay 1690/0, and the staged payload's digest.

It is simply not faster.

## Why

The first variant did what it was supposed to: 4.4% fewer instructions, from
`miv`/`mpiv` never being stored or reloaded and `previousw` being loaded once
instead of twice. Cycles went **up** anyway. Those arrays are about 3 KB at
this row length, so they were L1-resident and their traffic was already free in
a loop running at IPC 3.21 with spare load/store slots.

What the fusion costs is the reorder buffer. The scan's prefix tree is three
dependent permute-and-blend stages, and in the unfused arrangement that latency
was covered by the window holding many iterations of a very small loop at once.
A fused iteration is roughly twice the size, so half as many fit, and the chain
becomes visible — IPC 3.21 to 2.91.

The second variant confirms the diagnosis rather than fixing it. Issuing the
next block's scan before the current block's row gives the chain 28 uops of
independent work to hide behind, and it recovers most of the loss: IPC back to
3.14. But the staging costs the instructions the fusion saved (3.620 G against
3.610 G unfused), so the result is still slower than doing nothing.

This is the same lesson as the reverted running-maximum split two sections
above, and it is worth stating once more in the form that keeps being true
here: **this loop is latency- and throughput-bound, and instruction counts
predict the wrong sign.** Both fused variants would have looked like wins on
`Ir` alone.

## What is left

Replaces the ranking in the previous section; item 1 there is struck.

1. **`part_a_align`** in `dvtditr`, which is 41% of the aligner's cost rather
   than 35%. `perf record` on the staged ESCOTT payload puts the function at
   **20.4%** of that binary, not the 25% the previous ranking claimed, next to
   `match_calc` 25.2%, `intergroup_score` 13.0%, `init_strict` 11.4% and
   `alignable_reagion` 9.2%.

   Note it is **already fused**: the `mi` recurrence runs inline with the row,
   in one loop. So the work is to *split* it and vectorise both halves — which
   is the arrangement `l_align11` was already in when its kernel paid off, and
   the opposite of what the section above tried and rejected.

   Its `mi` recurrence is the same shape as `L__align11`'s with
   candidates `previousw[j - 1] + ogcp2[j] * gapfreq1[i - 1]` instead of
   `previousw[j - 1]` — and `partSalignmm.c` sets `USE_PENALTY_EX 0`, so there
   is no decay: it is a **pure prefix maximum**, exact for every `f64` with no
   integer precondition and no serial fallback. The row after it is elementwise
   in `j` once `miv`/`mpiv` exist, and carries no cross-lane reduction at all.
   Both halves are easier than the `l_align11` versions were.
2. **`intergroup_score` (12%) and `alignable_reagion` (10%)** in `dvtditr`,
   still unexamined.
3. **Packing value and index into one `f64` key** in the scan. Unchanged, and
   still low value.

---

# `part_a_align`, split and vectorised, 2026-09-07

The item the section above ranked first. It is the same two pieces as
`l_align11` — a scan and an element-wise row — with one difference that makes
it easier rather than harder, and one that makes it a split rather than a
fusion.

## What made it easier

`partSalignmm.c` opens with `#define USE_PENALTY_EX 0`, so the C never executes
the `mi += fpenalty_ex` that `Lalign11.c` does. The recurrence there is a plain
running maximum with no decay.

That matters because the scan's one conditional rewrite is reaching a position
with `m + t*e` where the loop added `e` to itself `t` times, which agree only
when nothing rounds. `l_align11` establishes that from its matrix and
penalties. Here `e` is zero — the `mafft` script never passes `-g` to
`dvtditr`, and the defaults are `DEFAULTGEP_B/N/J == 0` — so `t*e` is exactly
zero and the rewrite is exact for every `f64`, with no property of the
candidates required. The kernel is asked for only when `fpenalty_ex == 0.0`;
anything else runs the serial loop for the whole row.

The candidates differ from `l_align11`'s — `previousw[j-1] + ogcp2[j] *
gapfreq1[i-1]` rather than `previousw[j-1]` — but they enter the maximum and
nothing else, so writing them out as a row lets
`lalign_avx512::prefix_scan` serve this recurrence unchanged. The port does not
grow a second scan.

## What made it a split

`partA__align` runs its recurrence *inline* with the row, in one loop. So the
work is the opposite of the section above: separating them, not fusing them.
The separation is what makes the row element-wise in `j`, and an element-wise
row is the shape a kernel can take eight cells of at a time.

Separating costs something on its own — `miv`/`mpiv` become real memory and the
row is walked twice. Measured alone, before either kernel: **+6.2%**. That is
the price of admission, and it is why the previous section's result is not a
contradiction of this one: there, an already-separated pair was fused and lost;
here, an already-fused pair is separated and wins, because only one of the two
arrangements can be vectorised.

## The row kernel

`src/partalign_avx512.rs`, the port's second kernel and the second place it
contains `unsafe`. It is simpler than `lalign_avx512::row`: three comparisons
and four stores per block, and — because this row carries no running maximum —
**no cross-lane operation at all**. Every lane runs the scalar computation on
its own cell, so bit-exactness needs no tie-rule argument, only the absence of
reassociation.

The C has no kernel here either, so this is not catching up with it.

## Result

Staged ESCOTT `dvtditr` payload, best of seven, fresh directory per run, `pre`
digested every time:

| | dvtditr |
|---|---|
| before | 0.2081 s |
| split alone, no kernel | 0.2211 s |
| **split + both kernels** | **0.1911 s** (0.918x) |

`perf record` on the same payload, before and after:

```
                     before    after
match_calc            25.2%    27.4%
part_a_align          20.4%     8.7%   + 3.4% now in prefix_scan
intergroup_score      13.0%    14.2%
init_strict           11.4%    12.5%
alignable_reagion      9.2%    10.1%
```

End to end, the whole point: a cold GEMME/ESCOTT trace of P53_HUMAN, 100
realign calls, one worker, no cache, four repetitions each --

```
before   45.48  45.47  45.50  45.49 s
after    43.78  43.77  43.80  43.81 s
```

**3.7% off a cold run**, with the predicted matrices bit-identical (same digest
on all eight runs). `dvtditr` is 41% of the aligner and `part_a_align` was 20%
of it, so this is close to the whole of what the item was worth.

## Gates

`cargo test --release` 227 passed. Against the compiled `core/`: `verify_dvtditr`
32/32, `verify.sh` (disttbfast) 12018/0, `verify_tbfast.sh` 1690/0,
`verify_lalign_kernel.sh` 15/0. The staged payload's every output file
byte-identical to the previous build, and the end-to-end matrices digest-equal.

Of those, `pre-push` runs `cargo test`, the coverage floor and the two golden
replays; `verify_dvtditr` and `verify_lalign_kernel` were run by hand, the
first because it needs a compiled `core/` and the second because it is the only
one that can tell the port's two paths apart.

Note `verify_golden.sh` refuses to replay against a `core/` it did not record
against, so a rebuilt C needs the binary moved aside or the store re-recorded.
That guard fired here and is working as intended.

## One thing found on the way, not fixed here

The port applies `mi += fpenalty_ex` and `m[j] += fpenalty_ex` unconditionally
where the C has both inside `#if USE_PENALTY_EX`, which this file sets to 0 —
so the C never applies them. With `fpenalty_ex` zero the two agree (`x += 0.0`
differs from `x` only in the sign of a zero, which every comparison here treats
as equal), and zero is what every call through the `mafft` script passes. A
direct `dvtditr -g ...` would diverge.

Left alone deliberately: the differential corpus cannot reach the case, so the
fix could not be verified by the means this repository trusts, and an unverified
change to an untested path is worse than a recorded one.

## What is left

1. **`match_calc`, 27.4%** and the biggest single item in `dvtditr`. Already
   optimised twice; shared with `tbfast` and `disttbfast`, so a win there is
   worth more than its share here suggests.
2. **`intergroup_score` (14.2%) and `alignable_reagion` (10.1%)**, still
   unexamined. Look for algorithmic waste before intrinsics — every large win
   in this log has been algorithmic.
3. **`init_strict`, 12.5%**, still 0.67 G over the C's equivalent and still not
   localised beyond the function.


# Peak RSS, 2026-10-27 -- the instrument, not the metric

Where GOAL.md stands: wall clock is 1249 ms against the frozen baseline's 2367
(52.8%, met with room to spare, thanks to the AVX-512 kernels). Mean peak RSS
is 6725 KB against a 5220 KB target, and **not** met.

This section was first written to say that the four small cases are inherently
unmeasurable and that the metric therefore needed reweighting before anything
else. That was wrong. The cases are fine. Three instruments in a row were not,
and the record is kept in order because each failure is a different one and all
three are easy to repeat.

## Instrument one: polling VmHWM races the program

`profile.sh` polled `/proc/<pid>/status` for VmHWM while the process ran. The
four small cases finish in about a millisecond, so they exit before the first
read, and the poller reported whatever it caught during start-up -- or nothing.
`aa.fftns2` read as a 92-191% spread across repeats. `/bin/true` measured 0 KB,
which should have given it away at once.

This is a *different* fault from the one the 2026-09-04 section records. That
was the shared working directory, and it is fixed; these runs each had their
own directory and still spread 191%.

## Instrument two: fork copies the parent

The replacement read `ru_maxrss` from `wait4` in Python. `fork` copies the
parent's address space before `exec`, and `ru_maxrss` counts it, so the
interpreter put a floor of 8188 KB under every case. Five of the eight read
exactly 8188 with zero variance, which is how it showed.

## Instrument three: wait4 from a small parent

`bench/maxrss.c`, about thirty lines: fork, chdir, exec, `wait4`, print
`ru_maxrss`. The parent is a few hundred KB, so the floor is 1324 KB
(`/bin/true`), and the kernel's own high-water mark needs no sampling interval.
`profile.sh` builds it on demand; the bench path already needs `cc` for
`core/`.

Eleven runs per case, fresh directory each:

```
sample.fftns2    8888 (8800- 8944)      aa.fftns2      4208 (4084-4304)
sample.nwns2     8316 (8264- 8408)      nuc.fftns2     4212 (4140-4376)
sample.memsave   5308 (5184- 5372)      tandem.fftns2  4492 (4376-4592)
sample.linsi    14296 (14284-14396)     dup.nwns2      4048 (3960-4168)
```

Every case inside 5%, `sample.linsi` inside 0.8%. There was never anything
wrong with the small cases or with an unweighted mean; the harness simply could
not see a process that short.

## What the fixed instrument changes

The C was under-reported far worse than the port. Its small cases read 3260,
3472 and 3384 KB before; they are 22140, 22204 and 22212 KB. MAFFT allocates by
upper bound at start-up, and the poller never saw it happen.

```
mean peak RSS   port 6725 KB   C 16913 KB   (0.40x)
per-case ratio  mean 0.61x, worst 1.54x (sample.linsi)
```

So the port's memory advantage over the C is much larger than this database
has ever recorded, and the single case where it loses is `sample.linsi` at
14.3 MB against the C's 9.3 -- now measurable to under 1%, which is what makes
it worth working on.

**The baseline's RSS column came from instrument one.** Comparing 6725 KB
against its 5800 KB is not a comparison; `summarize.py` will keep saying 116%
of a number that was never measured. Correcting it relaxes the RSS target,
which `../CLAUDE.md` puts to a human, so it is left alone and flagged here and
in GOAL.md.

## A 1.6 MB red herring, turned twice

massif on `aa.fftns2` points straight at four eager allocations of 400,000 B
each:

- `cut1`, `cut2` -- `falign.rs:464-465`
- `result1`, `result2` -- `falign.rs:295-296` (`BlockAlignState::new`)

Four `vec![0; MAXSEG]` of `i32`, `MAXSEG = 100000`, allocated whatever the
input -- directly below the comment explaining that `segment`, `segment1` and
`segment2` were made grow-on-demand for this exact reason. 1.6 MB on a
510-byte input. An obvious win.

It is not a win. Rebuilt with `MAXSEG = 1000`, a 100x cut, and measured on both
instruments, because the first result rested on the broken one:

```
                 polled          wait4
aa.fa    100000  3316 KB         4176 KB
aa.fa      1000  3432 KB         4228 KB
sample   100000  8912 KB         8868 KB
sample     1000  8928 KB         8872 KB
```

Zero either way. `vec![0; n]` allocates through `alloc_zeroed`; nothing writes
those pages, so they never fault in and were never resident. **Allocated is not
resident.** Worth removing for allocation hygiene by anyone already in the
file; worth nothing for peak RSS. Recorded because massif will point the next
person at it just as invitingly.

## Where the RSS work actually is

`sample.linsi`: 14.3 MB, 1.54x the C, reproducible to 0.8%, and 21% of the
eight-case total on its own. Everything else is at or below the C, three cases
by a factor of five.

## The 5220 KB target is out of reach, with the floor measured

Measured on the fixed instrument, what a run costs before any alignment:

```
/bin/true                                 1324 KB   kernel + loader
version-rs (prints a string)              2092 KB   + Rust runtime
disttbfast-rs, no args (exits on usage)   2912 KB   + this binary mapped in
aa.fftns2 (6 seqs, a full alignment)      4228 KB   + 1316 KB of actual work
```

**2912 KB is spent before the aligner starts** -- 73% of the smallest case's
peak. No amount of work on the DP or the allocator moves it; only a smaller
binary would.

The arithmetic that follows:

- Today's mean is 6725 KB over eight cases, sum 53800. A 5220 KB mean needs the
  sum at or under 41760, a cut of 12040 KB.
- The four small cases cannot sum below 4 x 2912 = 11648 KB. Granting them
  their floor, the four heavy cases must sum to <= 30112 KB against 36832
  today, an 18% cut.
- The only case worse than the C is `sample.linsi`, 14348 against 9312.
  Bringing it all the way to C parity -- the most anyone could ask -- leaves
  the heavy sum at 31796, still **1684 KB over**.
- The remainder would have to come from `sample.fftns2`, `sample.nwns2` and
  `sample.memsave`, which already run at 0.35x, 1.13x and 0.24x of the C.

So condition (1) is not reachable by optimisation, and the reason is not that
the port is wasteful: it is that 5220 was derived as 90% of 5800, and 5800 came
from instrument one, which under-measured by three to seven times on half the
cases. The target is an artifact of a broken measurement, not a goal anyone set
on the evidence.

### The floor, measured rather than assumed

The paragraph above argued from "these cases already beat the C", which is not
a floor. So the floor was measured: massif's peak *live* heap per case, plus the
kernel-and-loader cost of any process at all, granting zero bytes for the
binary, the Rust runtime, malloc slack and the stack.

```
case              RSS    live heap   slack   generous floor
sample.fftns2    8864      7275       1589        8599
sample.nwns2     8296      6974       1322        8298
sample.memsave   5324      3390       1934        4714
sample.linsi    14348     12554       1794       13878
four small cases                                4 x 2912
                                    generous floor mean 5892 KB
                                    GOAL.md target      5220 KB
```

Malloc overhead is 46-65 KB per case, under 1%: the slack column is the binary
and the runtime, not fragmentation. These cases are **live-data bound**. There
is no waste to reclaim; the memory is the working set.

So the floor sits **672 KB above the target**, and that is with every case
granted a floor it cannot actually reach.

One structural redundancy is known and still present: `align.rs`, `genalign.rs`
and `lalign.rs` each keep their own `commonIP` traceback matrix where the C
shares one TLS global (`core/dp.h:15`), about 2.4 MB each at `sample.linsi`'s
sizes. Eliminating it *entirely* -- all of it, from every case -- gives:

```
  -2400 KB of live heap -> mean 5592 KB   still over
  -4800 KB of live heap -> mean 5292 KB   still over
```

Below that, the only remaining route is for the port to hold less live data
than a faithful transliteration holds. The C holds *more* on three of these
four cases -- 25668 KB against 8864 on `sample.fftns2` -- so this is not a
matter of the port being wasteful, and `../CLAUDE.md`'s "faithfulness beats
idiom" is the rule that would have to give.

### Timed attempt, 2026-10-27 19:46:46Z to 19:48:18Z (2 minutes)

Start and finish recorded this time, as GOAL.md asks. Two minutes, because the
attempt was to establish whether the one known structural win is real before
spending an hour on the refactor -- and the answer arrived from massif rather
than from editing.

**The two traceback matrices are genuinely live at the same time.** At
`sample.linsi`'s peak:

```
treebase::treebase (impmtx)   5,208,800 B   40%   <- largest single item
align::ensure_ijp             1,798,280 B   14%
lalign::l_align11 ijp         1,798,280 B   14%   <- same size, both live
lalign::l_align11 (f64)         522,240 B    4%
```

So consolidating them, as the C does with its one `commonIP` TLS global, saves
**1,798,280 B of live heap** -- measured as simultaneously live, not assumed.
`genalign.rs`'s pair does not appear: L-INS-i never takes the generalised
affine path, so only two of the three are ever in play here.

That is a real win and it is more faithful to the C, not less. It is also not
enough. Floor with it applied: 45381 KB, **mean 5673 KB against the 5220
target**. The conclusion above holds, now on a measured saving rather than the
2.4-4.8 MB that was guessed at.

**`impmtx`, the 5.2 MB that is 40% of the peak, is locked in.** `Salignmm.c:110`
allocates it with `AllocateFloatMtx( impalloclen, impalloclen )`, which looks
like a 2x saving sitting in plain sight -- until you read `mtxutl.c:174`, where
`AllocateFloatMtx` returns `double **` and allocates `sizeof( double * )`. The
name is a legacy misnomer; it allocates doubles. The port's `Vec<Vec<f64>>` is
exact, and narrowing it to `f32` would change every score it feeds. The largest
item in the heaviest case cannot be touched.

What a human needs to decide, and an agent must not:

1. Re-baseline against `bench/maxrss.c` and set the RSS number from numbers
   that exist. On today's database, 90% of the mean would be about 6050 KB.
2. Or drop the RSS half of the target and keep `sample.linsi` as the named
   memory item, since it is the one case that is both measurable and worse than
   the C.

Either way `summarize.py` will keep reporting "116% of 5800 KB" until someone
does, and that percentage is of a number that was never measured.


# Re-baselined, 2026-10-27

The RSS target was re-based on the fixed instrument, by the repository owner's
decision. What changed and what it means:

- The baseline is now `bench/baseline-rebased.tsv`: today's HEAD measured with
  `bench/maxrss.c` for peak RSS, with the wall-clock column carried over from
  2026-09-04 because the instrument fault was in RSS only.
- **The old target could not be reached.** 5220 KB was 90% of 5800 KB, and
  5800 came from the poller. The measured floor is 5892 KB, so the target sat
  672 KB *below* the floor.
- The new target is 6052 KB, 90% of a correctly measured 6725. It sits
  **above** the floor, so the RSS half is a real target for the first time:
  673 KB of mean to find against roughly 160 KB of headroom.

`77dd5a6`, the original baseline commit, could not be re-measured: the squash
merge of #5 collapsed it and the object is gone from the clone and from the
remote. Recorded because it is the reason the re-baseline had to use today's
HEAD rather than re-measuring the commit the old number came from, and because
it is a general lesson -- a squash merge destroys the ability to re-measure
anything the collapsed commits were the baseline for.

`bench/baseline-2026-09-04-remeasured.tsv` is deleted. It was the half-fix
(fresh working directory, still polling) and keeping a third baseline file
measured by the broken instrument is how today's confusion started.

## What is left on the RSS half

1. **`commonIP`, about 220 KB of the mean.** The two traceback matrices live
   simultaneously at `sample.linsi`'s peak, 1,798,280 B each, where the C
   shares one TLS global. Necessary, not sufficient. A `thread_local` pool
   with a take/put API is drafted and parked, not landed: it needs the
   hour-long differential harnesses behind it, and both aligners' restore
   points changed (`align.rs:1099`, `align.rs:1318`, `lalign.rs:515`).
2. **The other 453 KB is not identified.** `impmtx` -- 5.2 MB, 40% of
   `sample.linsi`'s heap -- is locked in: `Salignmm.c:110` allocates it with
   `AllocateFloatMtx`, which `mtxutl.c:174` shows returns `double **`. The
   name is a misnomer and the port's `f64` is exact.
3. The four small cases are 69% fixed cost (2912 KB before the aligner
   starts). Nothing in the DP moves them; only a smaller binary would.


# `commonIP` shared, 2026-10-27 20:33:48Z to 20:55:08Z (21 minutes)

Both clock readings taken, as the re-based GOAL.md now requires.

`AlignState` and `LalignState` each kept their own traceback matrix. At
`sample.linsi`'s peak both are live, 1,798,280 B apiece -- 28% of that case's
heap, for one matrix's worth of information. The C has had one for this all
along: `extern TLS int **commonIP` in `dp.h`, grown monotonically by whichever
aligner needs it biggest.

`src/commonip.rs` is a thread-local pool with `take(orlgth1, orlgth2)` and
`put(mtx)`. Five places changed, not the thirty a threaded parameter would
have needed: the two `ensure_ijp` bodies now borrow from the pool, and the
three restore points (`align.rs` twice, `lalign.rs` once) hand the matrix back
to it instead of to `self`. Both `common_alloc1`/`common_alloc2` pairs move
into the pool, which holds the running maximum the way the C's single
`commonAlloc1` does.

Why this cannot change a value: neither aligner reads a cell it has not written
in the same call -- which is what makes the C's sharing safe too -- and neither
re-zeroes the buffer on entry today, so handing over dirty contents is already
the behaviour. What changes is that the matrix is sized to the larger of the
two demands instead of allocated twice.

The buffer is lent out rather than borrowed in place, so no borrow is held
across a DP and a future reentrant path would cost memory rather than panic.
There is none today: `align.rs` never calls into `lalign.rs` or the reverse.

### Measured

```
sample.linsi peak RSS   14296 -> 12636 KB   (-1660 KB, median of 11, spread 0.8%)
mean peak RSS            6725 -> ~6517 KB   (target 6052)
```

Predicted 1756 KB from massif, got 1660: the difference is that not every page
of the second matrix was resident. Prediction high rather than low, which is
the right direction to be wrong in.

### Verified

```
verify.sh         (disttbfast, A__align)            pass=12018  fail=0
verify_tbfast.sh  (tbfast, L__align11 + A__align)   pass=1690   fail=0
verify_e2e.sh     (real mafft wrapper, L-INS-i in)  pass=28     fail=0  refused=1
cargo test        (199 unit + the 553-case matrix)  green
```

`verify_tbfast.sh` is the one that decides: `tbfast` is where both matrices are
live at once, so it exercises the sharing directly. The `refused` is by design,
a configuration outside the ported path.

`verify_e2e.sh` took three attempts to run at all -- twice for invocation
errors of mine (no arguments, then no `MAFFT=`), and once because `make` leaves
`core/mafft` without its execute bit so the harness's `-x` test rejected it.
Recorded because a harness that prints its usage message looks nothing like a
failure in a log, and the first two attempts produced no `pass=` line at all.

### Still open

~465 KB of mean. This was always necessary-but-not-sufficient; nothing about
the remaining gap is identified. `impmtx` -- 5.2 MB, 40% of `sample.linsi`'s
heap -- is locked in: `Salignmm.c:110` allocates it with `AllocateFloatMtx`,
and `mtxutl.c:174` shows that returns `double **`. The name is a legacy
misnomer and the port's `f64` is exact. The four small cases are 69% fixed cost
before the aligner starts, so only a smaller binary moves them.


# The small cases, 2026-10-27 21:04:14Z to 21:05:40Z

Timed attempts now total 23 + 2 = 25 minutes of the 90 the goal allows.

With `commonIP` landed the gap is 455 KB of mean, and the four small cases
looked like the place to find it: each sits about 1300 KB above the 2912 KB
floor, so 5200 KB of apparent headroom against 3640 KB needed.

It is not headroom. `aa.fftns2`, 4216 KB measured:

```
score matrix, 256x256 f64 (align.rs:1155)       522 KB
Sextet 6-mer tables, 6^6 x 4 B, two of them     373 KB
Constants::new_seq                              195 KB
accounted for                                  1090 KB   84% of the heap above floor
```

All three are fixed-size tables, and all three are the size the C makes them:

- **The score matrix is 0x100 x 0x100 because it is indexed by raw sequence
  bytes.** `Galign11.c` allocates the same, and that width is load-bearing in a
  way a Rust port makes sharper: the matrix is *written* only at `cn.amino[i]`
  positions but *read* at `seq[i]`, so a byte of 200 in the input reads
  `[200][...]`. In the C that is a calloc'd 0.0. In a 128-wide table it is a
  panic. Narrowing it would turn the C's silent zero on malformed input into a
  crash, which is a behaviour change however much memory it saves.
- **The sextet tables are 6^6 entries** because that is the 6-mer code space;
  `disttbfast.c` allocates them the same way. A six-sequence input does not
  make the code space smaller.

Four small cases at 522 KB of score matrix is 2088 KB of total, 261 KB of mean
-- more than half the remaining gap, sitting in a table whose dimensions are
not ours to choose.

## Where that leaves the RSS half

455 KB of mean, with no identified source. Each of the four largest items in
the database is now accounted for and locked:

| item | size | why it cannot move |
|---|---|---|
| `impmtx` | 5.2 MB | `AllocateFloatMtx` returns `double **`; the port's `f64` is exact |
| score matrix, x4 | 2.1 MB | 0x100 wide because raw bytes index it |
| sextet tables, x4 | 1.5 MB | 6^6 is the 6-mer code space |
| fixed process cost, x8 | 23.3 MB | 2912 KB before the aligner starts; only a smaller binary |

`commonIP` was the one piece of genuine duplication and it is spent. What is
left would need either a smaller binary, or a deviation from the C's own data
structures -- and `../CLAUDE.md` puts the second one out of reach without a
human deciding that memory now outranks faithfulness.


# Growing a matrix cost two of them, 2026-10-27 21:09:09Z to 21:36:37Z (27 min)

Started on `sample.nwns2`, the one case never broken down -- 8276 KB and 1.12x
the C, the second-worst ratio. massif found something better than a lead.

At its peak, **two** `commonip::take` allocations were live: 2,095,760 B and
1,798,280 B, 54% of that case's heap, for one matrix's worth of information.

My first reading was that the pool's fallback path had fired -- the one for "a
caller that asks while the matrix is already out", which `commonip.rs`'s own
comment says has no path today. Wrong. Reading the code instead of the comment:

```rust
p.mtx = vec![vec![0i32; ll2 + 10]; ll1 + 10];
```

An assignment evaluates its right-hand side before dropping the old value, so
during growth the old matrix and its replacement both exist. And growth happens
at exactly the moment the working set is largest, so the overlap lands *on* the
peak rather than beside it.

`p.mtx = Vec::new();` first. One line.

```
                      before    after    delta
sample.fftns2          8864     7596    -1268
sample.nwns2           8276     7220    -1056
sample.linsi          12560    12100     -460
sample.memsave         5332     5332        0
mean peak RSS          6507     6153   91.5% of baseline, target 6052
```

**2.8 MB, against the 1.66 MB that sharing `commonIP` was worth.** The
refactor's value turned out to be mostly that it put the allocation somewhere a
profiler would point at.

Not a regression introduced by that refactor, either: both `ensure_ijp`
implementations had the same `self.ijp = vec![..]` shape before they shared a
pool. The overlap has been there as long as the port has grown these matrices.

### Verified

```
verify.sh         pass=12018  fail=0
verify_tbfast.sh  pass=1690   fail=0
verify_e2e.sh     pass=28     fail=0  refused=1
cargo test        238 tests, including the 553-case matrix
```

### 101 KB short, and where to look next

The target is 6052 and the database says 6153. The gap is now small enough that
the next thing to try is obvious: **this pattern, everywhere else it occurs.**
Any `x = vec![...]` that replaces a live buffer pays for both copies at the
moment it is most expensive. It paid twice here. `grep` for assignments to a
field that already holds a `Vec` and check each one -- `align.rs`'s and
`lalign.rs`'s scratch buffers (`w1`, `w2`, `initverticalw`, `m`, `mp`) all grow
the same way, and so do `falign.rs`'s.

That is a mechanical sweep, not a redesign, and it is the first lead in this
log that is both unexplored and cheap.


# Target met, 2026-10-27 21:38:38Z to 22:06:32Z (28 minutes)

```
mean peak RSS   5870 KB   was 6725   87.3%   MET  (target <= 6052)
wall clock      1243 ms   was 2367   52.5%   MET
```

The sweep the previous section called "mechanical rather than a faithfulness
question", and it was. Thirteen large allocations had the shape
`x = vec![...]`, which evaluates its right-hand side before dropping the old
value -- so a growth holds two copies at the moment the working set is largest.

Almost all of the win was one of them:

```
sample.linsi    12100 -> 10236   -1864   impmtx, treebase.rs:230
sample.fftns2    7596 ->  7636      +40   (inside this case's spread)
sample.nwns2     7220 ->  7108     -112
```

`impmtx` is the 5.2 MB constraint matrix, regrown as the progressive alignment
lengthens. `sample.linsi`'s ratio against the C falls from 1.54x to 1.09x: it
was the only case meaningfully above the C on memory and now it is not.

For that one the fix is `drop(std::mem::take(&mut impmtx))` rather than
`impmtx = Vec::new()`. Both free the old matrix, but for a local the compiler
cannot see that the drop is the point and warns the value is never read; the
explicit form keeps a `-D warnings` build clean and the intent legible.

### Verified

```
verify.sh         pass=12018  fail=0
verify_tbfast.sh  pass=1690   fail=0
verify_e2e.sh     pass=28     fail=0  refused=1
verify_dvtditr.sh pass=32     fail=0
cargo test        238 tests, the 553-case matrix among them
```

`verify_dvtditr.sh` reported nothing on the first run, because it takes the
libexec *directory* and was handed the binary. That is the third harness
misinvocation today -- `verify_e2e.sh` took two -- and each printed a usage
line rather than a failure. **A harness that does not run looks like a harness
that passed, in a log.** Worth a guard: these scripts should exit non-zero in a
way the caller cannot mistake for silence, and the caller should require a
`pass=` line rather than its absence.

### What the whole run says

Five timed attempts, 2 + 21 + 2 + 27 + 28 = 80 minutes, every start and finish
recorded.

The two wins were both duplication, and neither was in the DP:

| change | worth |
|---|---|
| `commonIP` shared between the aligners, as the C shares it | 1.66 MB on `sample.linsi` |
| free a buffer before allocating its replacement, x13 sites | 2.8 MB, then 1.86 MB more |

Everything that looked like an algorithmic lead was locked: `impmtx`'s *size*
is the C's, the score matrices are `0x100` wide because raw bytes index them,
the sextet tables are 6^6 because that is the code space. The memory was never
in the arithmetic. It was in holding two of something while swapping one.

And the earlier conclusion in this log -- that the target was unreachable
because the measured floor was 5892 KB -- was right about the floor and wrong
about the ceiling: 5870 KB came in just under it, because the floor was
computed from *live heap at peak*, and the growth overlap was inflating that
very measurement. A floor derived from a number the bug was corrupting.


# Checking a consumer-side optimisation brief, 2026-10-28

A session working the consumer side (`provolut-ml`) sent a measured brief
ranking what to optimise, from 100 real realign payloads on a Zen 5. Three of
its four priorities rest on claims that do not reproduce here. Recorded with
what was measured, because the ranking changes.

Input: 51 unaligned sequences, 180-285 aa, synthesized to the brief's shape
(`countlen`: `51 x 285 - 180 p`, against the brief's `51 x 286 - 180 p`). This
machine has AVX-512 in hardware, so the kernel on/off split is meaningful.

## The brief's profiling recipe does not run the program

Its §5 gives the `tbfast` argv verbatim and says components read `infile` from
the current directory. Run exactly as given, `tbfast` prints
`Warning: Only 0 sequence found.`, **exits 0**, and finishes in 3 ms -- both
binaries "agree" at 3 ms and `pre` is 0 bytes.

The mechanism, corrected: the driver **redirects stdin** (`core/mafft:2525`
ends `< infile`; `core/mafft:2693` ends `< pre` for `dvtditr`). An earlier
version of this section said the driver appends `-i infile`; it does not, and
the brief's author was right about that. `-i <file>` is a second input path the
binary also accepts, and it is what the timings below used -- all three
variants wrote a byte-identical 26387-byte `pre`, so the measurements stand
even though the explanation of the driver was wrong.

A recipe that profiles nothing, exits 0, and reports plausible-looking
millisecond numbers is the same failure this log records three times already
under harness misinvocation. Anyone following §5 verbatim measures start-up.

## `dvtditr` already has a kernel (their priority 1)

The brief calls it "the biggest untouched win", "no AVX-512 kernel", "net gain
from the first instruction". `dvtditr-rs` carries **64 zmm references** --
`partalign_avx512.rs`, reached through `partalign.rs` for the restricted
iterative refinement. Measured, staged from one C `tbfast` run so every
variant sees identical inputs:

```
C dvtditr             294 ms
oxide, kernel on      245 ms   0.83x of C
oxide, NO_AVX512=1    291 ms   0.99x of C      kernel worth 16%
```

The port is already 17% faster than the C here, and the kernel is already
earning 16%. The brief's own table shows this faintly -- 94.8 ms on against
105.4 ms off -- and read a real 10% as absence.

Its author has since revised that to "the Rust code is ~7% behind C" scalar,
and deferred the kernel question ("You know your own code"). Neither half holds
on this machine: the kernel is there (64 zmm), and the scalar path is 0.99x of
C rather than 7% behind it.

## There is no 3.9x scalar deficit in `tbfast` (their priority 3)

```
                      this machine              brief's Zen 5
C                       581 ms                    121 ms
oxide, kernel on        154 ms  0.27x of C        119 ms  0.98x
oxide, NO_AVX512=1      525 ms  0.90x of C        473 ms  3.9x
```

The *Rust* columns nearly agree across the two machines (154/119, 525/473).
The **C** differs 4.8x. And `core/tbfast` here contains **0 zmm references**:
`Lalign11.c`'s kernel is behind `#if LALIGN11_AVX512`, which `CFLAGS = -O3`
never satisfies, so this repository's own C is scalar.

Against that scalar C the Rust scalar path is **0.90x -- marginally faster, at
parity**. The brief's 3.9x is a Rust-scalar-against-something-else comparison;
whatever its conda-forge build is doing, this repo's `gcc -O3` build does not
do it. So "every scalar-side win in tbfast is pure net speedup on top of a
kernel that already works" is still true as an aspiration, but the 3.9x hole it
was meant to fill is not there.

Why their C is 4.8x faster than ours is unresolved and cannot be settled from
here. It matters, because the whole ranking depends on which C is the baseline.

## `disttbfast` is the real target (their priority 2, confirmed)

```
            C zmm   Rust zmm
tbfast          0         68
dvtditr         0         64
disttbfast      0          0
dndpre          0          0
```

`disttbfast` has no kernel on either side, and the brief measures it at 99.3%
of a 17.8 s cold-build call. That is the one priority that survives checking,
and it is the largest single cost in the pipeline. Its inner loop is 6-mer
distances and two progressive cycles -- not what either existing kernel does,
so it needs its own.

## Bit-exactness

Unchanged and verified again on this shape: `tbfast` and `dvtditr` output
identical to the C with the kernel on **and** off. The `MAFFT_RS_NO_AVX512`
escape hatch proving both paths agree is what makes that checkable.

## The C baseline's provenance, and why the cross-machine ratio is still open

The consumer-side brief's author proposed that this repo's C baseline might be
unoptimised, which would make "Rust scalar at 0.90x of C" a property of the
baseline rather than of the port. Checked, because it is cheap and because it
underwrites every ratio in `profile.tsv`:

```
core/Makefile:23                  CFLAGS = -O3
readelf -p .comment core/tbfast   GCC: (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
objdump -d | grep -c zmm          0
objdump -d | grep -cE 'ymm|xmm'   19244      <- auto-vectorised at -O3
```

An optimised release build. The hypothesis is refuted and the baseline stands.

That does not reconcile the two machines, and the reason is worth keeping. The
brief reports, at n=100 on a 17121-byte payload on a 6-core Zen 5: C 121.2 ms,
Rust kernel-off 490.7 ms. This machine, on a 26387-byte payload on a 4-core
cloud Xeon: C 581 ms, Rust kernel-off 525 ms. The Rust arms agree to 1.07x; the
C arms differ by 4.8x.

Payload size and machine speed cannot explain that, because both would move the
two arms together. A 1.54x payload and a roughly 2x slower core predict the
single-threaded Rust number should be about 3x larger here; it is 1.07x larger.
The port is single-threaded on both sides -- no `std::thread` and no rayon
anywhere in `rust/src/` -- so that 1.07x is a clean per-core comparison, and it
says the two machines and the two payloads very nearly cancel. Whatever moves
the C arm by 4.8x therefore does not touch the Rust arm at all.

One candidate visible from here, unresolved: this build has threading compiled
in.

```
core/Makefile:9    ENABLE_MULTITHREAD = -Denablemultithread
core/Makefile:39   LIBS = -lm -lpthread
objdump -T core/tbfast | grep -c pthread_create   1
```

The C threads only when told to (`-thread n`) and the port never does, so if the
two invocations differ in that argument -- or in any other algorithm-selecting
flag -- the two C numbers are of different work and the ratio between them means
nothing. Settling it needs one shared payload and one pinned argv with
`-thread 1` forced, three arms each, not more sampling. Until then neither
`tbfast` scalar figure is actionable: 0.90x here against 4.05x behind there.

Recorded so the next person does not re-derive it: the disagreement is in the
argv, not in the optimiser.


# `disttbfast` was 1.6x the C on the `--auto` build path, 2026-10-29

Reported from the consumer side and reproduced here: `mafft --auto
--preservecase` was 1.65x the C's wall clock on a 3831x334 production payload
and 1.61x on a 790x509 payload generated by `gen/gencorpus.py`. Output was
byte-identical throughout; this was a speed regression only.

It had been invisible because of which cases `profile.tsv` holds. Its weighted
wall-clock total reads 0.68x of the C, but that figure is carried by
`sample.linsi` (`tbfast`, 0.247x) and three tiny `fftns2` cases that win on
cache footprint. The `disttbfast` rows at realistic size were *already* at
parity or slightly behind -- 1.208x, 1.384x and 1.243x on instructions -- and
the consumer's build step runs exactly that path. A headline that aggregates
over cases nobody runs will hide a regression in the one they do.

## Why this path and not the realign path

`core/mafft.tmpl` sends any input between 500 and 20000 sequences to
`ktuples, iterate=0, cycle=2`. The second cycle recomputes the whole O(n^2)
distance matrix over the *aligned* sequences (`core/disttbfast.c:2151`), so
`naivepairscorefast` runs once per pair: 311,655 calls at n=790, 7.33M at
n=3831. The realign path (L-INS-i, `tbfast` + `dvtditr`) is 3.4x *faster* than
the C and was never involved.

## Two hypotheses that measured out wrong

Recorded because both looked obviously right from reading the code, and the
cost of chasing them was a callgrind run each.

`amino_dis` is `Vec<Vec<i32>>` against the C's flat `int[0x100][0x100]`, so
every lookup is a double indirection into one of 256 scattered allocations.
That is real, and it is **0.7%** of the run: the line measures 202M of
28,217M. Flattening it would have bought nothing.

Memory behaviour was the second guess, from wall clock (1.61x) being further
behind than instructions (1.317x). cachegrind says the port's cache behaviour
is *better* than the C's: D1 miss rate 3.4% against 3.6%, LL misses 241,941
against 454,849.

## What it actually was

callgrind, per function, on the 790-sequence payload:

```
                                    C          port
naivepairscorefast family       5,885M      11,335M     1.93x
  of which Vec::push                 -       1,526M
  of which allocator plumbing        -         421M
alignment DP (A__align/a_align) 8,824M       5,414M     0.61x   <- port wins
total                          21,420M      28,217M     1.317x
```

The port returned two freshly allocated `Vec`s from
`commongappickpairfast` on every one of those calls. Fixed by writing into one
thread-local buffer pair grown monotonically -- the shape `commonip`'s matrix
pool already uses -- returning the column count, collapsing `p1`/`p2` into the
one cursor they always were, and hoisting the row length in both that function
and `pairscore` so the indexing drops its bounds checks.

```
                 before      after        C
instructions    28.22G      24.22G     21.42G     1.317x -> 1.131x
wall (790)       1326ms      1124ms      820ms     1.61x  -> 1.37x
hot family      11.34G       6.81G      5.89G     1.93x  -> 1.157x
```

Measured on the consumer's 3831x334 payload by their session, same two arms,
best of 3: 17.557 s -> 14.395 s against the C's 10.639 s, so 1.65x -> 1.35x
and 3.16 s saved. Their C arm reproduced to within 2 ms of its own earlier
independent run.

## Where the remaining 1.35x is

Instructions are now 1.131x but wall clock is 1.37x, and cachegrind rules out
memory. `--branch-sim=yes` locates it:

```
                        C          port
conditional branches  2.784G      4.748G     +1.96G
mispredict rate         2.5%        1.0%     port better
```

Nearly two billion extra conditional branches, predicted almost perfectly and
still costing branch-unit throughput -- at ~2 per cycle that is close to the
304 ms that remains. They are the bounds checks that are left, and callgrind
says where: `slice/index.rs` attributes 1,325M to `a_align`, 923M to the
distance-matrix loop in `disttbfast_rs::main`, 199M to `falign` and 123M to
`match_calc`. Closing it means iterator-based rewrites in those four loops,
each needing its own differential. It does **not** mean vectorising
`pairscore`, whose bounds checks are already gone.

## An f64 accumulator that did not need changing

`pairscore` accumulates in `f64` while both its addends are `i32`
(`penalty_dist`, and `amino_dis` via `n_dis[i][j] = n_distmp[i][j] as i32`),
so every partial sum is an exact integer and an `i64` accumulator is provably
bit-identical rather than merely close. It was built and measured: **1125 ms
against 1124 ms**, no change, because the addition was never the cost. Reverted
-- the C accumulates in `double`, and a deviation that buys nothing measurable
is what "faithfulness beats idiom" exists to prevent.

`score_is_integral` in `tests_unit.rs` was kept. The invariant is real, the C
depends on it implicitly, and it is the thing a future change would break
silently if anyone does want to vectorise this loop later.

## Verified

Byte-identical to the C on both payloads, and identical to the port's own
pre-patch output; `verify.sh` pass=12018 fail=0; `verify_tbfast.sh` pass=1690
fail=0; 243 tests; region coverage 91.30%.

## The FFT-NS-2 build path: 1.37x slower becomes 1.14x faster

The 1.35x left at the end of the previous entry was recorded as "bounds-check
branch throughput in the O(n²) distance matrix", with the work named as
iterator rewrites in `a_align`, the distance-matrix loop, `falign` and
`match_calc`. Three of those four would have been wasted effort. The entry also
said the warm path does not run the build step, which stopped being true: the
port became the default aligner on linux-64 for a consumer whose production
call is `--auto --preservecase` at 3831 sequences, so the build step is the
path that runs.

### Where the gap was, by subtraction rather than by profile

`-E` is the number of progressive cycles. `-E 1` does the ktuple matrix, the
guide tree and one progressive alignment; `-E 2` adds a second cycle and, in
between, the O(n²) distance matrix recomputed over the *aligned* sequences.
790 sequences, best of five, same binaries, same input:

```
arm     -E 1    -E 2    increment
C       357     823     466
port    357    1124     767
```

Everything `-E 1` touches is at parity to the millisecond. The whole of the gap
is in the increment, and the increment is dominated by the one function the
recompute calls per pair. That is a two-line measurement and it retired three
of the four candidate sites outright: `a_align` and `falign` run in the half
that already matches, and `match_calc` -- which costs 1.365G against the C's
0.361G, the worst *ratio* anywhere in the binary -- also runs there, so its
3.8x buys nothing. A ratio is not a priority.

### What was actually expensive, and the hypothesis that was wrong

A line-table build (`CARGO_PROFILE_RELEASE_DEBUG=line-tables-only` into a
separate target dir, so the shipped profile is untouched; output byte-identical
to the release binary) splits `naivepairscorefast`'s 6.81G:

```
commongappickpairfast  ~4.50G   skip1[k1] 555M, skip2[k2] 694M,
                                the two `!= '-'` branches 831M + 969M,
                                the scratch writes 415M + 138M
pairscore              ~1.55G   of which amino_dis 202M
```

So two thirds of it was the column *picking*, and `amino_dis` -- the line an
`i64` accumulator would have touched -- was 0.83% of the program. The previous
entry's revert of that accumulator was right, and right for this path and not
only for the refinement one.

The obvious reading of that table is that the scratch writes and the second
pass are the waste, so the two were fused: generate each kept column and score
it in the same walk, no buffers, columns read once. The addend chain is
untouched, so it is bit-identical by construction rather than by tolerance.

**It was worth 0.25%.** 24,224M instructions to 24,164M, and 1123 ms to 1081.
The second pass was walking a small cache-hot buffer with no carried
dependency, which is the cheapest kind of work there is; removing 1.7G of
instructions that retire several per cycle buys close to nothing.

### The real cost was the skip tables, consulted per column

What the fusion did not remove is the loop's serial dependency: `skip1[k1]`
loads, feeds a subtract, feeds a min, feeds the cursor, feeds the next load,
with two bounds checks in the middle. Dropping the tables entirely and walking
every column of the alignment -- same scoring, same order -- was measured as a
control:

```
best of 9, interleaved, 790 sequences, -E 2 -F, all byte-identical to the C
  C                    822 ms   1.000
  port, HEAD          1123 ms   1.366
  port, fused         1081 ms   1.315
  port, no tables      708 ms   0.861
```

1.59x, from deleting the accelerator. The per-column table traffic cost far
more than the columns it skipped.

That control cannot ship: it is O(alignment width) per pair where the table is
O(columns kept), and for two short rows in a wide alignment those diverge
without bound. That is the case `naivepairscorefast` was added upstream for.
Stated as a bound, not as a measurement -- a 300-row payload of 40- and
1200-residue sequences failed to show it, because at 300 rows the O(n·L²)
progressive step swamps the O(n²) matrix and all arms land within noise.

The shipped loop keeps the jump and pays for it only where it can fire. The
C computes `min(skip1[k1] - skipped1, skip2[k2] - skipped2)` at the top of
every column; that minimum is non-zero **exactly** when both rows are inside a
gap run, because a row sitting on a residue has consumed every gap recorded
before it and contributes 0. Both rows gapped is also exactly the condition for
dropping the column. So testing the two bytes first and reading the tables only
then makes the same jumps, with the loads and the dependency off the common
path:

```
  C                    822 ms   1.000
  port, hybrid         720 ms   0.876
```

1.56x against HEAD, and the port now beats the C on the build path. The hybrid
gives up 2% against the no-table control to keep the asymptotic bound, which is
the right trade for a function whose whole reason to exist is that bound.

### Verified

Byte-identical to the C, and to the port's own pre-patch output, on both
payloads and both `-E` arms. `verify.sh` pass=12018 fail=0;
`verify_tbfast.sh` pass=1690 fail=0; 244 tests; region coverage
91.33%, line 90.77%.

`fused_pair_score_matches_the_unfused_one_on_gap_heavy_rows` is the new test
and it generates its gaps rather than taking what an aligner produces, because
the paths the fusion introduced -- the skip-table jump, and a gap run in one
row spanning a column dropped for being gapped in both -- need both-gap columns
on purpose. It asserts bit equality against `naivepairscore11`, which drops
those columns by looking at them and consults no table, and it asserts both a
case count and that at least 200 of its row pairs contained a both-gap column,
so it cannot pass by having tested the old path only. Four mutations were run
against it: charging every gap column instead of once per run, failing to
terminate a run on a residue pair, and letting a dropped column terminate a run
are all caught; widening the jump test from "both rows gapped" to "row 2
gapped" is not, and should not be -- with only one row gapped the minimum is 0
and the `skip > 0` guard declines the jump, so that mutation is slower and not
wrong.

### Why `profile.tsv` could not have caught this

The database barely moves: `sample.fftns2` goes 923,130,040 instructions to
908,784,979, and the other rows likewise shift by about a percent. `test/sample`
is 36 sequences, so its cycle-2 recompute is 630 pairs where the production call
has 7.3 million, and the function that carried the whole gap is nearly absent
from every row. The database is still the right evidence for what it covers; it
just does not cover this, and the `-E 1`/`-E 2` subtraction above is what does.
`sample.memsave`'s `ms_ratio` reads 1.033 to 1.121 in the same regeneration
while its instruction count *falls*, which is wall-clock noise on a 848 ms case
and not a finding.

Those two counts are from a second run. The first one also reported 12018/0 and
1690/0, but `cargo` relinked `target/release/disttbfast-rs` two minutes and
forty-seven seconds into it, because the new test went into the library and the
binaries depend on the library. The mutation runs in between used
`cargo test --lib`, which does not build the binaries, and the timestamps say
the relink happened after the last mutation was reverted -- so the first run
almost certainly tested the right thing. "Almost certainly" is not what these
harnesses are for. They were re-run against binaries copied to a directory
nothing else writes to, with `sha256sum -c` after the run confirming neither
had changed while it was in progress.

### Confirmed at production size

Measured by the consumer's session on its own 3831 x 334 payload, best of 2,
fresh process, `--auto --preservecase`, output byte-identical at 4,542,015
bytes: C 10844 ms, the packaged `0.3.0 _2` 14682 ms (1.354x), and `cde06cc`
9054 ms (0.835x). 1.20x ahead of the C at production size against 1.14x at 790,
which is the direction the size ladder predicts -- the fix targets the O(n^2)
term, so its share grows with n.

No synthetic 3831 payload was built here to check it. Mutated copies of one
root align far more cleanly than real homologues, so the alignment width and
therefore the per-pair cost would both have been wrong, and the number would
have looked like evidence without being any.

## Threading the all-pairs matrices

`--thread n` now drives the three all-pairs matrices and nothing else. 790
sequences, `-E 2 -F`, best of 7, interleaved, every arm byte-identical to the
C's sequential run:

```
  C   -C 0-0     825 ms     port  -C 0-0     723 ms
  C   -C 8-0     569 ms     port  -C 2-2     594 ms
  C   -C 8-8     464 ms     port  -C 4-4     530 ms
                            port  -C 8-8     518 ms
                            port  -C 12-12   520 ms
```

1.59x against the C's unthreaded run, and ahead of the C's own pairwise-only
`-C 8-0`. Both stop scaling at about six workers: this box is a 6-core Ryzen
9600X with SMT, so `-C 8-8` already oversubscribes the physical cores and the
flat step from 8 to 12 is the hardware, not the partitioning.

### Why these loops and not the progressive step

Each cell of those matrices is a pure function of its own `(i, j)` and of data
the loop does not write. Nothing accumulates across cells and no cell reads one
another writes, so a thread count cannot change the result -- there is no
f64 ordering question here at all, which is the only reason this was reachable
without touching the one rule.

The C threads the same loops with a mutex-guarded `(i, j)` cursor every worker
pulls from (`disttbfast.c:1274-1296`), so its own assignment of pairs to
threads differs run to run; it reaches the same answer for the same reason.
Static round-robin by row is used here instead because it needs no lock, and
because a fixed assignment is reproducible under a debugger. Round-robin rather
than contiguous blocks: these are half matrices, so row `i` holds `njob - i`
cells and a block split would hand the first worker most of the work.

The progressive step is left sequential, and not out of caution. `alg` is a
process-wide global (`defs.c:80`) selecting between `A__align` and `MSalignmm`,
which are different recurrences, and `treebasethread` writes it with no lock at
`disttbfast.c:1929` when a node exceeds 30000 columns. Upstream's comment on
the next line -- "hoka no thread ga M ni shitakamo shirenainode", because
another thread may have set it to M -- documents that this is visible across
workers. So above 30000 columns **the C's threaded progressive step has no
fixed answer to be byte-identical to**, which is the `dvtditr` situation and not
a milder one. Part of the C's `-C 8-8` advantage is therefore an answer that is
not reproducible, and it is not a target.

### What the refusal used to say

`-C > 0` was refused outright with "threads change the result". That was never
measured. `core/disttbfast` built with `-Denablemultithread` gives
byte-identical stdout and side files at `-C 1-1`, `-C 8-8` and `-C 8-0` as at
`-C 0-0` on 1038 comparisons across four strategies and two payloads. The claim
was an assumption wearing an error message, and the cost of believing it was
that `--thread 1` -- semantically sequential -- exited 1 with no output.

### Verified

`verify_threads.sh` pass=2752 fail=0: 86 corpus inputs x 4 strategies x 8 `-C`
specs, each compared against the **C's** `-C 0-0` run rather than against the
port's own, because a port that was self-consistently wrong at every thread
count would pass the weaker comparison. The specs include 3 and 5, which do not
divide the row count, and 64, which exceeds it so some workers get no rows at
all; `-C 8-0` and `-C 0-8` catch a transposed parse of `-C a-b`, which every
symmetric spec would miss. Mutation-tested by dropping one worker's bucket:
exit 1, 24 of 64 failed.

`verify.sh` pass=12018 fail=0 and `verify_tbfast.sh` pass=1690 fail=0 unchanged,
against binaries frozen and `sha256sum -c`-checked across the run. 245 tests;
region coverage 91.37%, line 90.82%.

The golden matrix went 553 cases to 567. The seven `*.threads` fixtures were
declared refusals -- with a comment saying a future thread implementation would
have to update them -- and are now byte-identity cases against the C, which is
a stronger assertion than the refusal they replace; `threads3` (`-C 3-3`) and
`threadspair` (`-C 2-0`) were added. This is the permitted regeneration: a
feature landed and its fixtures came from the real C. All 546 non-thread
manifest entries are byte-identical across the regeneration, checked rather
than asserted.

## match_calc's row pointers, and two things that did not work

Once the matrices were threaded the progressive step became the dominant cost
at high thread counts, so `match_calc` was worth a second look. Two of its
hottest lines were addressing and nothing else:

```
282M   let codes = &intwork[j];
141M   let vals = &doublework[j];
```

423M instructions -- 2% of a build -- re-walking two outer `Vec`s once per
column, in a function called 578,166 times. Slicing all three row tables to one
length before the loop and zipping them removes it:

```
21,095M -> 20,534M instructions, -560M, -2.7%
wall clock at 790 sequences: 723 ms -> 722 ms
```

**The wall clock did not move**, and that is the honest headline. It is kept
anyway: the reduction is real and deterministic, the code is simpler than the
indexing it replaces, and a payload or a machine with less headroom than this
one will not necessarily hide it. But nobody should quote this change as a
speedup.

### The counted loop that was a wash

The inner loop over a column's non-zero entries was rewritten from
`codes.iter().zip(vals.iter())` to `for k in 0..m` over slices, on the theory
that building an adapter per column cost more than the one or two
multiplications it feeds -- the lists average 1.35 entries of 26. The theory was
wrong and the measurement says so precisely:

```
                   before      after
zip.rs              723M        144M
iter/range.rs          0        719M
match_calc total    2312M      2306M
```

579M left `zip.rs` and 719M appeared in `range.rs`, for a net 6M. `zip` over
two slices of equal length already compiles to the counted loop; the adapter was
never the cost. Reverted, and the row-pointer hoist re-measured on its own at
20,534,325,958 -- within 0.6M of the two changes together, so the whole 560M is
the hoist. That split is the only reason the attribution is known, and the
reason the two were not committed together.

### The AVX-512 kernel that was not written

Vectorising `naivepairscorefast` was the plan until the profile priced it. That
loop is now 3.80G against the C's 5.89G -- already 0.64x -- so the ceiling is
about 9% of instructions, and reaching it needs: an `unsafe` module, the f64
accumulator replaced by an integer one (lane-parallel f64 sums are
reassociation, which is the one thing that may not move, while an integer sum
is exact in any order -- `score_is_integral` already pins that both addends are
`i32`), a gather for `amino_dis`, and `vpcompressb` for the gap-run scan, whose
state is carried column to column. The row-pointer hoist above was the same
order of win with no unsafe and no ordering question. Priced and declined, not
attempted and abandoned.

### Noted, not acted on

On the `-INS-i` path `l_align11` is 71.2% of the run with a further 11.5% in
`zip.rs` and 0.8% in `enumerate.rs` -- 680M of 5.51G in adapters, a larger pool
than anything above. It is left alone for now for two reasons: the port is
already 2.32x faster than the C there, so it is upside rather than a gap; and
the wash recorded above is direct evidence that instructions leaving `zip.rs`
may simply reappear under another name. It needs the same split-and-measure
treatment, not a rewrite on the strength of a percentage.

### Verified

Byte-identical to the C on both `-E` arms and at `-C 0-0` and `-C 8-8`.
`verify.sh` pass=12018 fail=0; `verify_tbfast.sh` pass=1690 fail=0;
`verify_threads.sh` pass=2752 fail=0; `verify_dvtditr.sh` pass=32 fail=0 --
`match_calc` feeds `a_align`, which every strategy reaches, so the group and
refinement harnesses matter here and not just the build one. 245 tests; region
coverage 91.38%, line 90.83%. All four harnesses ran against binaries copied
to a directory nothing else writes to, with `sha256sum -c` after the last of
them confirming none had changed.

`verify_dvtditr.sh` first exited 2 with "no dvtditr in .../core/dvtditr": it
takes a libexec *directory* and was handed the binary's path. It refused rather
than reporting a vacuous pass, which is what that exit code is for.
