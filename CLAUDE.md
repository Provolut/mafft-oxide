# Instructions for Claude and any other agent working in this repository

## The one rule everything else serves

**Output must not change.** This repository exists to make MAFFT faster while
producing byte-identical alignments. A change that alters one residue of one
alignment is a bug, however fast it is, and however plausible the new answer
looks. `provolut/README.md` puts it as: verified byte-identical against the
unmodified upstream build, never *argued* to be equivalent.

## Two gates you may not lower

`.githooks/pre-push` refuses a push whose tests fail or whose region coverage
is under 90%. Both numbers are floors. **An agent may not lift, weaken, or
route around either one.**

The coverage floor tracks the suite's measured coverage a point or so below it,
rather than sitting far enough back to be unreachable. That is deliberate: it
catches code added ahead of its tests instead of only a collapse, so a module
and its tests belong in the same commit. `rust/README.md`'s harness table
carries what the suite currently measures; this file carries only the floor,
so the two cannot drift into disagreeing about the same number.

Concretely, all of the following are forbidden, whether or not they would make
a red run go green:

**On tests**

- Deleting a test, commenting one out, or marking it `#[ignore]`.
- Weakening an assertion so it passes: loosening an exact comparison to an
  approximate one, widening a tolerance, asserting a weaker property, or
  replacing a specific expected value with whatever the code currently
  produces.
- Editing a golden fixture in `rust/tests/golden/` to match new output. Those
  files are the C's answer. If the port's output no longer matches them, the
  port changed its answer and that is the finding — regenerate them with
  `gen_golden.py` **only** when the C itself changed or a new feature was added,
  never to accommodate a diff in the port.
- Removing cases from `gen_golden.py`, or narrowing the corpus, the argument
  vectors, or the `applies_to` predicates, so that fewer things are compared.
- Making `tests/matrix.rs` more permissive: relaxing `same()`, treating a
  mismatch as a warning, or reclassifying a failing case as `exit`/`refuse`
  when the C in fact succeeds.
- `git push --no-verify`, or `git config core.hooksPath` to something else, to
  get past the hook.

**On coverage**

- Lowering `FLOOR` in `.githooks/pre-push`.
- Changing `bench/coverage.sh` so it reports a higher number for the same
  code: dropping the two binaries from `llvm-cov`'s `--object` list in reverse,
  adding `--ignore-filename-regex` entries to hide uncovered files, counting a
  different metric, or measuring a subset of the crate.
- Deleting or excluding poorly-covered code to raise the percentage.
- Adding tests that execute code without asserting anything about it. Coverage
  is a proxy for "this is tested"; padding it with assertion-free exercise
  makes the proxy lie, which is worse than a low number.

If the gate fires, the answer is to fix the code or write a real test. If you
believe a gate is genuinely wrong — the floor is miscalibrated, a test asserts
something untrue, a fixture was generated from a broken build — **say so and
stop, and let a human decide.** Do not decide it yourself and proceed. Explain
what fired, why you think it is wrong, and what you would change.

A gate an agent can talk itself past is not a gate.

## Sources are UTF-8, and a hook keeps them that way

Upstream ships eight files whose Japanese comments are EUC-JP. That is not a
cosmetic problem: GNU grep in a UTF-8 locale treats such a file as binary and
reports **nothing** for it -- no match, no count, no warning -- so a "nothing
calls this function" or "this macro is unused" conclusion silently excludes
whichever files happen to be affected. Here that was `constants.c`, `io.c`,
`mltaln9.c`, `Lalign11.c`, `Lalignmm.c` and `MSalignmm.c`: all of I/O, the
memsave aligner, and the file `fillimp` lives in. It has already produced two
wrong statements in this repo's own documents.

They are converted, and `scripts/check-encoding.sh` (run by `pre-push`) refuses
a source that is not valid UTF-8. **After importing a new upstream release, run
`scripts/check-encoding.sh --fix`** -- the release will bring the EUC-JP files
back.

Converting is safe, and that was measured rather than assumed: every affected
byte is inside a comment, and `core/` built before and after the conversion
gives 34 of 34 binaries byte-identical.

## What to do when output diverges

Instrument a scratch copy of `core/`, dump the intermediates as f64 bit
patterns, dump the same from Rust behind an env var, and bisect. Every
mismatch in this port so far was found that way, never by reading the code and
reasoning about it.

## Faithfulness beats idiom

`rust/` transliterates C. Where the C's evaluation order, its persistent
`static` scratch buffers, or its edge-case behaviour look wrong, they are
usually load-bearing: `rust/README.md`'s "Notes for anyone extending this"
lists eleven cases where tidying the obvious way changes alignments. Do not
"clean up" C-shaped code without a differential run proving the output is
unchanged.

This applies to floating point in particular. **Never reassociate f64
arithmetic for speed.** Loop transposition that preserves each accumulator's
summation order is fine and has been done; changing the order in which values
are added is not, whatever the compiler flags allow.

## Where the ground truth lives

- `rust/tests/matrix.rs` — 553 cases whose expected bytes came from the real C.
  Runs in `cargo test`, needs no `core/` build. This is the fast gate.
- `rust/verify/verify.sh`, `verify_tbfast.sh`, `verify_dvtditr.sh`, `verify_e2e.sh`,
  `verify_gaplen.sh`, `verify_flh.sh`, `verify_td.sh`, `verify_tdi.sh` —
  differential harnesses against a compiled `core/`. Tens of thousands of
  comparisons, about an hour. These are the thorough gate; run them before
  claiming a change is safe, not just the fast one.
- `rust/bench/profile.tsv` — the profiling database, C and port side by side.
  Committed on purpose: it is the evidence for any performance claim, and a
  diff in it is reviewable.

Never report a harness as passing without running it. A harness that reports
`pass=0 fail=0` is vacuous, not green — this one has been fooled twice, and
both guards against it are in the scripts because of that.

## Performance claims

Measure, do not assert. One case is not the picture; `bench/summarize.py`
aggregates the database and prints the work-weighted ratio first, because an
unweighted mean of per-case ratios lets tiny inputs flatter the result by 60%.
Measure one change at a time — of the two `match_calc` fixes, the one that
looked obvious was worth 1% and the one that looked like a rounding detail was
worth 4%.

Watch for instrumentation that distorts what it measures. Putting line tables
on the release profile once read as a 9% memory regression; only the per-case
breakdown showed it was 750 KB of debug sections.
