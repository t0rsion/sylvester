# Gröbner engine benchmark: sylvester vs Singular, msolve, Macaulay2, Groebner.jl

Run 2026-08-19. sylvester 0.3.0, commit 517ca24.

## Protocol

- Field F_1073741827. Order: grevlex. Output: reduced Gröbner basis.
- Cores: `taskset -c 0-3,12-15` (the P-cores of the bench machine). Every
  child process runs inside a transient systemd scope with `MemoryMax=16G`
  and swap disabled.
- Deadline: 120 s per computation. A cell that exceeds it is DNF, and
  larger instances in the same family are then skipped for that engine.
- External tools (Singular, Macaulay2, Groebner.jl): median of 3 runs,
  timed with each tool's own internal timer around the compute call
  (Singular `rtimer`, M2 `elapsedTiming`, Julia `@elapsed`).
- msolve: timed in-process by `runner/bin/msolve-inproc`, which calls
  `core_msolve` directly instead of subtracting a process-startup baseline
  from the CLI's wall clock. It repeats the call until a 0.5 s time
  budget is met, with at least 3 repetitions and at most 5000, bounded by
  a 120 s single-repetition ceiling and a 130 s total deadline. Each
  repetition runs in its own forked child, so a leak in `core_msolve`'s
  `print_gb` path dies with the child instead of accumulating across
  repetitions. `results.json` records the repetition count and the
  `[min, max]` spread beside the median for every msolve cell: for
  example, katsura-9 ran 3 repetitions, median 0.229873543 s, spread
  [0.22986385, 0.248123346] s, and cyclic-4 ran 2768 repetitions, median
  0.000178806 s, spread [0.000174711, 0.000283503] s.
- Thread pinning: msolve `-t 1`, Julia `-t 1`, sylvester
  `ComputeOptions::threads(1)`, for every column except `f4/threads8`,
  which requests 8 threads.
- Engine versions, from `results.json`'s `_meta`: singular 4.4.1, msolve
  0.10.1, m2 1.26.05, julia 1.12.7, groebner.jl 0.10.3, sylvester
  0.3.0+517ca24. `_meta` carries no separate commit field; the commit is
  the build metadata suffix on the sylvester version string.

### Rational protocol

- Twelve instances at sizes smaller than the prime-field cells, named
  with a `-q` suffix: cyclic-4..6, katsura-4..7, noon-3..5, eco-8..9.
  Coefficient growth, not monomial count, drives cost over `Q`, so these
  are the sizes where the reference tools return in seconds. Input files
  add a `.sylq` format for sylvester, identical to `.syl` line for line
  except that a term's coefficient field also accepts
  `numerator/denominator`; Singular, msolve, and Groebner.jl each read
  their own rational input format (`.sing`, `.ms`, `.jl`).
- Reference tools: Singular 4.4.1, msolve 0.10.1, and Groebner.jl 0.10.3.
  Macaulay2 does not run a rational cell: `docs/rational-design.md` section
  1.4 names only these three as the rational design's references.
- msolve's rational timing is a plain CLI wall clock around
  `msolve -g 2`, not the in-process `msolve-inproc` runner the
  prime-field cells use above; there is no rational counterpart of that
  runner. The wall clock carries process startup, and it shows: the
  smallest rational cells cluster near 30 to 40 ms regardless of
  instance (cyclic-4-q 0.039 s, cyclic-5-q 0.042 s, katsura-4-q 0.036 s,
  noon-3-q 0.038 s, noon-4-q 0.027 s), well above sylvester's
  sub-millisecond time on the same cells, the signature of startup cost
  rather than compute cost.
- sylvester's rational run uses `ComputeOptions::threads(1)`. Every tool
  in this protocol runs on one thread: msolve `-t 1`, Julia `-t 1`, and
  sylvester as above. For the rational cells, `threads(1)` also caps
  `RationalOptions`' outer concurrency, so the multimodular driver runs
  one prime at a time instead of starting several prime runs together.
- Correctness is the gate here (`docs/rational-design.md` section 1.2); there
  is no speed gate over `Q`. Timing is recorded, not gated, because the
  per-prime cost this record shows is not yet the cross-prime
  learn-and-apply pipeline the design intends (section 3.9), which is
  later work.

## Correctness

### Prime field

For each instance, every engine that finished with a status of OK is
compared on basis size, on the sorted multiset of leading monomials (the
weaker check), and on the full canonical basis with coefficients
included (the check that matters; see `canon.py` and the "Full-output
comparison" section of `README.md` for the canonicalization). All 20
instances had every finisher's basis parsed successfully, so the
"finishers" column below is also the count compared on the full basis.

| instance | finishers | size/LM agreement | full basis agreement |
|---|---|---|---|
| cyclic-4 | 9 | AGREE | AGREE |
| cyclic-5 | 9 | AGREE | AGREE |
| cyclic-6 | 9 | AGREE | AGREE |
| cyclic-7 | 7 | AGREE | AGREE |
| cyclic-8 | 7 | AGREE | AGREE |
| katsura-4 | 9 | AGREE | AGREE |
| katsura-5 | 9 | AGREE | AGREE |
| katsura-6 | 9 | AGREE | AGREE |
| katsura-7 | 9 | AGREE | AGREE |
| katsura-8 | 8 | AGREE | AGREE |
| katsura-9 | 7 | AGREE | AGREE |
| katsura-10 | 5 | AGREE | AGREE |
| eco-8 | 9 | AGREE | AGREE |
| eco-9 | 7 | AGREE | AGREE |
| eco-10 | 7 | AGREE | AGREE |
| eco-11 | 5 | AGREE | AGREE |
| noon-3 | 9 | AGREE | AGREE |
| noon-4 | 9 | AGREE | AGREE |
| noon-5 | 9 | AGREE | AGREE |
| noon-6 | 8 | AGREE | AGREE |

A short script summed C(finishers with a parsed full basis, 2) over the
20 instances: 577 pairwise full-basis comparisons. No cell failed the
full comparison; every one of the 577 pairs agreed exactly, coefficients
included.

### Rational

The rational cells (`docs/rational-design.md` sections 1.2 and 11.4) compare
sylvester's multimodular engine against Singular, msolve, and Groebner.jl
over `Q`, on 12 instances across four families: cyclic-4..6,
katsura-4..7, noon-3..5, eco-8..9. Every finisher is compared on the full
canonical rational basis: exact fractions, monic, coefficients included,
not only leading monomials or size (`canon.py`'s `_q` functions build
this form once, the same way for every tool).

| instance | finishers | full basis agreement |
|---|---|---|
| cyclic-4-q | 4 | AGREE |
| cyclic-5-q | 4 | AGREE |
| cyclic-6-q | 4 | AGREE |
| katsura-4-q | 4 | AGREE |
| katsura-5-q | 4 | AGREE |
| katsura-6-q | 4 | AGREE |
| katsura-7-q | 4 | AGREE |
| eco-8-q | 4 | AGREE |
| eco-9-q | 4 | AGREE |
| noon-3-q | 4 | AGREE |
| noon-4-q | 4 | AGREE |
| noon-5-q | 4 | AGREE |

All 12 rational instances had 4 finishers with a parsed full output:
singular, msolve, groebner.jl, and sylvester (`results.json`'s `_meta`
records `rational_tools` the same way). A short script summed C(finishers
with a parsed full output, 2) over the 12 instances: 72 pairwise
full-basis comparisons. No cell failed the full comparison; every one of
the 72 pairs agreed exactly, exact fractions and coefficients included.

## The gate

Section 1.1 of `docs/f4-design.md` sets the release gate on four
cells, cyclic-7, katsura-9, eco-9, and noon-6, each run under 16 GB of
memory and a 120 s deadline, on the same thread count as the reference
engine (one thread, msolve's `-t 1`). It passes when the geometric mean
of the sylvester-to-msolve ratio over the four cells is at most 5x and
the worst single cell is at most 10x. The raw column is gated; the
certified column, reported below, is not. `docs/rational-design.md` section
1.2 adds no speed gate of its own: the design gates rational correctness only
(see above), so this section still measures the prime-field gate alone.

| cell | sylvester best raw | sylvester (s) | msolve (s) | ratio |
|---|---|---|---|---|
| cyclic-7 | f4/default | 0.074906 | 0.078739 | 0.95x |
| katsura-9 | f4/default | 0.250890 | 0.229874 | 1.09x |
| eco-9 | f4/default | 0.012994 | 0.013420 | 0.97x |
| noon-6 | f4/default | 0.011514 | 0.011741 | 0.98x |

Geometric mean: 1.00x (gate: at most 5.0x). Worst cell: 1.09x, at
katsura-9 (gate: at most 10.0x). Gate: PASS.

## Timings

Seconds, median of 3; for msolve, median of the in-process repetitions
described above. For the rational table, msolve's cell is a plain CLI
wall clock (see "Rational protocol").

### cyclic

| engine | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.002 | 0.079 | 1.228 |
| Groebner.jl | <0.001 | <0.001 | <0.001 | 0.023 | 0.590 |
| Singular | <0.001 | <0.001 | 0.006 | 0.769 | 22.0 |
| Macaulay2 | <0.001 | 0.001 | 0.015 | 2.083 | 77.6 |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.002 | 0.075 | 1.508 |
| sylvester f4 (8 threads) | <0.001 | <0.001 | 0.004 | 0.052 | 0.761 |
| sylvester classic | <0.001 | 0.002 | 0.228 | DNF | SKIP |

### katsura

| engine | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.002 | 0.007 | 0.040 | 0.230 | 1.641 |
| Groebner.jl | <0.001 | <0.001 | <0.001 | 0.002 | 0.027 | 0.149 | 0.398 |
| Singular | <0.001 | 0.001 | 0.007 | 0.067 | 0.554 | 5.110 | 38.0 |
| Macaulay2 | <0.001 | 0.002 | 0.020 | 0.201 | 2.110 | 26.8 | DNF |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.006 | 0.009 | 0.039 | 0.251 | 1.866 |
| sylvester f4 (8 threads) | <0.001 | <0.001 | 0.001 | 0.006 | 0.022 | 0.109 | 0.697 |
| sylvester classic | <0.001 | 0.006 | 0.032 | 0.348 | 4.204 | 54.1 | DNF |

### eco

| engine | 8 | 9 | 10 | 11 |
|---|---|---|---|---|
| msolve | 0.003 | 0.013 | 0.070 | 0.428 |
| Groebner.jl | 0.002 | 0.023 | 0.121 | 0.219 |
| Singular | 0.017 | 0.145 | 1.678 | 17.7 |
| Macaulay2 | 0.058 | 0.655 | 8.150 | DNF |
| sylvester f4 (1 thread) | 0.005 | 0.013 | 0.064 | 0.427 |
| sylvester f4 (8 threads) | 0.004 | 0.008 | 0.038 | 0.191 |
| sylvester classic | 4.802 | DNF | SKIP | SKIP |

### noon

| engine | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| msolve | <0.001 | 0.002 | 0.002 | 0.012 |
| Groebner.jl | <0.001 | <0.001 | 0.001 | 0.008 |
| Singular | <0.001 | <0.001 | 0.005 | 0.025 |
| Macaulay2 | <0.001 | <0.001 | 0.006 | 0.098 |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.002 | 0.012 |
| sylvester f4 (8 threads) | <0.001 | <0.001 | 0.002 | 0.015 |
| sylvester classic | <0.001 | 0.003 | 0.057 | 3.165 |

### rational

| engine | cyclic-4-q | cyclic-5-q | cyclic-6-q | katsura-4-q | katsura-5-q | katsura-6-q | katsura-7-q | eco-8-q | eco-9-q | noon-3-q | noon-4-q | noon-5-q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Singular | <0.001 | <0.001 | 0.026 | 0.001 | 0.005 | 0.032 | 0.295 | 0.030 | 0.433 | <0.001 | 0.001 | 0.005 |
| msolve | 0.039 | 0.042 | 0.031 | 0.036 | 0.016 | 0.041 | 0.048 | 0.034 | 0.060 | 0.038 | 0.027 | 0.041 |
| Groebner.jl | <0.001 | 0.001 | 0.004 | <0.001 | 0.002 | 0.005 | 0.030 | 0.008 | 0.039 | <0.001 | 0.001 | 0.006 |
| sylvester | <0.001 | 0.002 | 0.037 | 0.002 | 0.025 | 0.062 | 0.290 | 0.078 | 0.195 | 0.001 | 0.003 | 0.037 |

## Certified runs

Certification is defined only over `PrimeField`
(`Ideal<PrimeField>::groebner_basis_certified`); there is no certified
column over `Q`, so this section covers the prime-field cells alone. The
classic backend writes `sylv-gb-cert-v1` and F4 writes `sylv-gb-cert-v2`.
"split" is the certified total minus the same backend's raw seconds for
the same cell: certificate emission plus the first, bundled
verification, not a number the runner measures directly. `verify (s)` is
a second, standalone re-verification the runner times on its own, after
the certified call returns.

Certified katsura-9 does not finish inside 120 s under F4 (raw F4 itself
finishes katsura-9 in 0.250890 s; certifying it does not), and katsura-10
is then skipped for the F4 certified column. The classic
`sylv-gb-cert-v1` path does not finish katsura-8 and above: certified
katsura-8 is DNF, and katsura-9 and katsura-10 are then skipped for the
classic certified column.

| instance | backend | raw (s) | certified total (s) | split (s) | cert bytes | verify (s) | peak RSS raw (KB) | peak RSS certified (KB) |
|---|---|---|---|---|---|---|---|---|
| cyclic-4 | classic | <0.001 | <0.001 | 0.000110 | 3076 | 0.000045 | 3496 | 4076 |
| cyclic-4 | f4 | <0.001 | <0.001 | 0.000104 | 481 | 0.000043 | 4108 | 4028 |
| cyclic-5 | classic | 0.002 | 0.017 | 0.015005 | 166879 | 0.004171 | 3852 | 5312 |
| cyclic-5 | f4 | <0.001 | 0.002 | 0.001607 | 4529 | 0.001141 | 4000 | 4136 |
| cyclic-6 | classic | 0.228 | 3.219 | 2.990421 | 3968988 | 0.171528 | 4248 | 27804 |
| cyclic-6 | f4 | 0.002 | 0.025 | 0.023019 | 30041 | 0.019057 | 4184 | 4756 |
| cyclic-7 | classic | DNF | DNF | ? |  |  |  |  |
| cyclic-7 | f4 | 0.075 | 2.439 | 2.364454 | 784454 | 1.969116 | 9368 | 18288 |
| cyclic-8 | classic | SKIP | SKIP | ? |  |  |  |  |
| cyclic-8 | f4 | 1.508 | 36.4 | 34.855927 | 5549896 | 31.447833 | 62740 | 117196 |
| katsura-4 | classic | <0.001 | 0.003 | 0.002779 | 48198 | 0.001260 | 3692 | 4512 |
| katsura-4 | f4 | <0.001 | 0.001 | 0.001298 | 2725 | 0.001008 | 3916 | 4144 |
| katsura-5 | classic | 0.006 | 0.043 | 0.037786 | 376992 | 0.020118 | 3876 | 6340 |
| katsura-5 | f4 | <0.001 | 0.014 | 0.013872 | 10281 | 0.008358 | 4024 | 4344 |
| katsura-6 | classic | 0.032 | 1.025 | 0.992950 | 4289976 | 0.495490 | 4112 | 27104 |
| katsura-6 | f4 | 0.006 | 0.107 | 0.100732 | 48731 | 0.080612 | 4220 | 5488 |
| katsura-7 | classic | 0.348 | 24.5 | 24.144983 | 49045544 | 10.812658 | 4776 | 232460 |
| katsura-7 | f4 | 0.009 | 1.008 | 0.998863 | 235767 | 0.726632 | 4556 | 11056 |
| katsura-8 | classic | 4.204 | DNF | ? |  |  | 7348 |  |
| katsura-8 | f4 | 0.039 | 10.9 | 10.823677 | 1578004 | 7.604951 | 6752 | 45884 |
| katsura-9 | classic | 54.1 | SKIP | ? |  |  | 18472 |  |
| katsura-9 | f4 | 0.251 | DNF | ? |  |  | 15840 |  |
| katsura-10 | classic | DNF | SKIP | ? |  |  |  |  |
| katsura-10 | f4 | 1.866 | SKIP | ? |  |  | 47556 |  |
| eco-8 | classic | 4.802 | 32.5 | 27.690374 | 22000522 | 1.563540 | 4904 | 106676 |
| eco-8 | f4 | 0.005 | 0.318 | 0.312736 | 118909 | 0.226061 | 4512 | 7548 |
| eco-9 | classic | DNF | DNF | ? |  |  |  |  |
| eco-9 | f4 | 0.013 | 2.573 | 2.559821 | 535510 | 1.791814 | 5000 | 20508 |
| eco-10 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-10 | f4 | 0.064 | 25.3 | 25.236503 | 3466721 | 17.156467 | 7560 | 89164 |
| eco-11 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-11 | f4 | 0.427 | DNF | ? |  |  | 27472 |  |
| noon-3 | classic | <0.001 | <0.001 | 0.000642 | 12619 | 0.000298 | 3868 | 4104 |
| noon-3 | f4 | <0.001 | <0.001 | 0.000172 | 694 | 0.000101 | 3912 | 4084 |
| noon-4 | classic | 0.003 | 0.032 | 0.028904 | 329158 | 0.014712 | 3988 | 6224 |
| noon-4 | f4 | <0.001 | 0.008 | 0.007315 | 4337 | 0.002561 | 3996 | 4256 |
| noon-5 | classic | 0.057 | 3.204 | 3.146723 | 10855926 | 1.223008 | 4152 | 71624 |
| noon-5 | f4 | 0.002 | 0.042 | 0.040018 | 32875 | 0.026452 | 4524 | 5140 |
| noon-6 | classic | 3.165 | DNF | ? |  |  | 5728 |  |
| noon-6 | f4 | 0.012 | 0.758 | 0.746265 | 315964 | 0.492283 | 5400 | 11944 |

## Reading

At the gate's four cells, sylvester's `f4/default` config, one thread,
sits close to msolve, one thread: 0.95x on cyclic-7, 1.09x on katsura-9,
0.97x on eco-9, 0.98x on noon-6, a geometric mean of 1.00x against a 5x
budget and a worst cell of 1.09x against a 10x budget. The margin does not
hold at every size in the same families: at the largest instance of each,
sylvester trails msolve by more and Groebner.jl by more still. On
cyclic-8, sylvester is 1.508 s against msolve's 1.228 s (about 1.23x)
and Groebner.jl's 0.590 s (about 2.56x). On katsura-10, sylvester is
1.866 s against msolve's 1.641 s (about 1.14x) and Groebner.jl's 0.398 s
(about 4.69x); Macaulay2 is DNF there, and sylvester finishes. On
eco-11, sylvester is 0.427 s against msolve's 0.428 s (about even) and
Groebner.jl's 0.219 s (about 1.95x); Macaulay2 is DNF there too. Noon
tops out at noon-6, already a gate cell.

The `f4/threads8` column is not part of the gate: design 1.1 measures
one thread against msolve's one thread. Against `f4/default`, it
roughly halves the larger cells: cyclic-8 goes from 1.508 s to 0.761 s,
and katsura-9 from 0.251 s to 0.109 s.

Verifying a `sylv-gb-cert-v2` certificate costs far more than the raw
F4 computation it certifies, and the gap widens as the instance grows.
The verifier redoes every step with its own arithmetic, sharing no code
with the engine, so its cost tracks the size of the trace, not the
engine's. On the eco family, the standalone `verify (s)` column goes
from 0.226061 s on eco-8 (raw 0.005 s, about 45x) to 1.791814 s on eco-9
(raw 0.013 s, about 138x) to 17.156467 s on eco-10 (raw 0.064 s, about
268x). The limit is sharp: katsura-9 finishes raw F4 in 0.250890 s but
does not finish the certified call inside 120 s, and the classic
`sylv-gb-cert-v1` path stops finishing certified runs at katsura-8.

Since the prior record (`records/2026-08-19-f4/REPORT.md`), the engine code
moved from a prime-field-only representation to the sealed `D: Domain` parameter of
`docs/rational-design.md` section 2 and gained the `Budget` type behind
`ComputeOptions::budget`. Neither change touches the F4 kernel's
arithmetic or pair management, and the four gate cells show it: cyclic-7
moves from 0.075234 s to 0.074906 s, katsura-9 from
0.252204 s to 0.250890 s, eco-9 from 0.013555 s to 0.012994 s, and
noon-6 from 0.011349 s to 0.011514 s. The largest shift, on eco-9, is
about 4 percent; the rest are near 1 percent or less. That is run-to-run
noise, and the gate ratios above confirm it: 1.01x in the prior record,
1.00x here.

The rational cells are new in this record; the prior record covered the prime field
only. Groebner.jl is fastest on every rational cell. sylvester is the
slowest of the four reference tools on 4 of the 12 cells, cyclic-6-q
(0.037 s against Groebner.jl's 0.004 s, about 9x), katsura-5-q (0.025 s
against 0.002 s, about 13x), katsura-6-q (0.062 s against 0.005 s, about
12x), and eco-8-q (0.078 s against 0.008 s, about 10x), and
second-slowest, behind Singular, on 2 more: katsura-7-q (0.290 s against
Singular's 0.295 s) and eco-9-q (0.195 s against Singular's 0.433 s). On
the other 6, all the smallest cells in the set, msolve reads slowest,
but its number there is a CLI process's wall clock, and the "Rational
protocol" section above shows why that is not a fair read of its compute
time: those cells cluster near 30 to 40 ms regardless of instance size,
the mark of process startup, not work. Weighing the compute-bound cells
instead, sylvester's standing over `Q` is the weakest of the four tools
this record measures. Part of the reason is protocol:
`RationalOptions`' outer concurrency is capped at one thread
here, so the multimodular driver runs one prime at a time, and the
learn-and-apply cross-prime trace that would let later primes reuse the
first prime's structure is later work (`docs/rational-design.md` section
3.9). The design sets no speed gate over `Q` (section 1.2): the pipeline this
record measures is not the one the design intends.

## No new defects

This run found no new defect; `KNOWN_ISSUES.md` records the defects
found and fixed before it.
