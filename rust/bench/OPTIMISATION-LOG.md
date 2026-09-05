# Optimisation attempts, 2026-09-04

The run GOAL.md describes: beat the port as it stood at `77dd5a6` by 10% on
**both** wall clock and peak RSS, or spend 90 minutes and write down what was
tried. This is that write-up. The target was not met; what follows is where
the time goes and which stones are now turned.

Output did not change at any point. Every kept change was measured with
`bench/profile.sh`, which refuses to record a row unless both binaries exit 0
and produce identical bytes.

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
