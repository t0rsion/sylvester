# Gröbner engine benchmark: sylvester vs Singular, msolve, Macaulay2, Groebner.jl

Run 2026-08-19. sylvester 0.2.0, commit f413b0b.

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
  repetition runs in its own forked child. A leak in `core_msolve`'s
  `print_gb` path (about 650 KB a call, confirmed by a separate
  measurement) dies with the child instead of accumulating across
  repetitions. `results.json` records the repetition count and the
  `[min, max]` spread beside the median for every msolve cell: for
  example, katsura-9 ran 3 repetitions, median 0.232531 s, spread
  [0.230225106, 0.24535867] s, and cyclic-4 ran 872 repetitions, median
  0.00044973750000000003 s, spread [0.000176037, 0.001246908] s.
- Thread pinning: msolve `-t 1`, Julia `-t 1`, sylvester
  `ComputeOptions::threads(1)`, for every column except `f4/threads8`,
  which requests 8 threads.
- Engine versions, from `results.json`'s `_meta`: singular 4.4.1, msolve
  0.10.1, m2 1.26.05, julia 1.12.7, groebner.jl 0.10.3, sylvester
  0.2.0+f413b0b. `_meta` carries no separate commit field; the commit is
  the build metadata suffix on the sylvester version string.

## Correctness

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

## The gate

Section 1.1 of `docs/f4-design.md` sets the release gate on four
cells, cyclic-7, katsura-9, eco-9, and noon-6, each run under 16 GB of
memory and a 120 s deadline, on the same thread count as the reference
engine (one thread, msolve's `-t 1`). It passes when the geometric mean
of the sylvester-to-msolve ratio over the four cells is at most 5x and
the worst single cell is at most 10x. The raw column is gated; the
certified column, reported below, is not.

| cell | sylvester best raw | sylvester (s) | msolve (s) | ratio |
|---|---|---|---|---|
| cyclic-7 | f4/default | 0.075234 | 0.078256 | 0.96x |
| katsura-9 | f4/default | 0.252204 | 0.232531 | 1.08x |
| eco-9 | f4/default | 0.013555 | 0.013396 | 1.01x |
| noon-6 | f4/default | 0.011349 | 0.011732 | 0.97x |

Geometric mean: 1.01x (gate: at most 5.0x). Worst cell: 1.08x, at
katsura-9 (gate: at most 10.0x). Gate: PASS.

## Timings

Seconds, median of 3; for msolve, median of the in-process repetitions
described above.

### cyclic

| engine | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.002 | 0.078 | 1.215 |
| Groebner.jl | <0.001 | <0.001 | 0.001 | 0.023 | 0.595 |
| Singular | <0.001 | 0.001 | 0.010 | 0.812 | 21.9 |
| Macaulay2 | <0.001 | 0.001 | 0.015 | 2.087 | 78.1 |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.002 | 0.075 | 1.524 |
| sylvester f4 (8 threads) | <0.001 | <0.001 | 0.002 | 0.050 | 0.762 |
| sylvester classic | <0.001 | 0.002 | 0.234 | DNF | SKIP |

### katsura

| engine | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| msolve | <0.001 | 0.002 | 0.002 | 0.007 | 0.038 | 0.233 | 1.636 |
| Groebner.jl | <0.001 | <0.001 | <0.001 | 0.002 | 0.028 | 0.149 | 0.399 |
| Singular | <0.001 | 0.001 | 0.007 | 0.071 | 0.579 | 5.101 | 38.2 |
| Macaulay2 | <0.001 | 0.002 | 0.020 | 0.201 | 2.107 | 26.7 | DNF |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.001 | 0.011 | 0.042 | 0.252 | 1.862 |
| sylvester f4 (8 threads) | <0.001 | 0.002 | 0.001 | 0.009 | 0.021 | 0.110 | 0.711 |
| sylvester classic | <0.001 | 0.005 | 0.038 | 0.352 | 4.289 | 55.0 | DNF |

### eco

| engine | 8 | 9 | 10 | 11 |
|---|---|---|---|---|
| msolve | 0.003 | 0.013 | 0.069 | 0.429 |
| Groebner.jl | 0.002 | 0.022 | 0.120 | 0.215 |
| Singular | 0.015 | 0.145 | 1.673 | 17.7 |
| Macaulay2 | 0.058 | 0.654 | 8.053 | DNF |
| sylvester f4 (1 thread) | 0.003 | 0.014 | 0.066 | 0.436 |
| sylvester f4 (8 threads) | 0.007 | 0.011 | 0.038 | 0.201 |
| sylvester classic | 4.851 | DNF | SKIP | SKIP |

### noon

| engine | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.002 | 0.012 |
| Groebner.jl | <0.001 | <0.001 | 0.001 | 0.008 |
| Singular | <0.001 | <0.001 | 0.003 | 0.023 |
| Macaulay2 | <0.001 | <0.001 | 0.006 | 0.098 |
| sylvester f4 (1 thread) | <0.001 | <0.001 | 0.002 | 0.011 |
| sylvester f4 (8 threads) | <0.001 | 0.002 | 0.006 | 0.015 |
| sylvester classic | <0.001 | 0.002 | 0.063 | 3.201 |

## Certified runs

The classic backend writes `sylv-gb-cert-v1` and F4 writes
`sylv-gb-cert-v2`. "split" is the certified total minus the same
backend's raw seconds for the same cell: certificate emission plus the
first, bundled verification, not a number the runner measures directly.
`verify (s)` is a second, standalone re-verification the runner times on
its own, after the certified call returns.

Certified katsura-9 did not finish inside 120 s under F4 (raw F4 itself
finishes katsura-9 in 0.252204 s; certifying it does not). The classic
`sylv-gb-cert-v1` path does not finish katsura-8 and above: certified
katsura-8 is DNF, and katsura-9 and katsura-10 are then skipped for the
classic certified column.

| instance | backend | raw (s) | certified total (s) | split (s) | cert bytes | verify (s) | peak RSS raw (KB) | peak RSS certified (KB) |
|---|---|---|---|---|---|---|---|---|
| cyclic-4 | classic | <0.001 | <0.001 | 0.000129 | 3076 | 0.000049 | 3256 | 3412 |
| cyclic-4 | f4 | <0.001 | <0.001 | 0.000103 | 481 | 0.000038 | 3328 | 3392 |
| cyclic-5 | classic | 0.002 | 0.020 | 0.017711 | 166879 | 0.004165 | 3360 | 4648 |
| cyclic-5 | f4 | <0.001 | 0.002 | 0.001545 | 4529 | 0.001155 | 3364 | 3556 |
| cyclic-6 | classic | 0.234 | 3.272 | 3.038542 | 3968988 | 0.171516 | 3696 | 27244 |
| cyclic-6 | f4 | 0.002 | 0.026 | 0.023964 | 30041 | 0.019306 | 3620 | 4012 |
| cyclic-7 | classic | DNF | DNF | ? |  |  |  |  |
| cyclic-7 | f4 | 0.075 | 2.453 | 2.377364 | 784454 | 2.016108 | 8792 | 17708 |
| cyclic-8 | classic | SKIP | SKIP | ? |  |  |  |  |
| cyclic-8 | f4 | 1.524 | 36.3 | 34.794114 | 5549896 | 31.502089 | 62188 | 116580 |
| katsura-4 | classic | <0.001 | 0.008 | 0.007001 | 48198 | 0.001336 | 3188 | 3752 |
| katsura-4 | f4 | <0.001 | 0.007 | 0.006770 | 2725 | 0.002740 | 3388 | 3516 |
| katsura-5 | classic | 0.005 | 0.046 | 0.040982 | 376992 | 0.019876 | 3400 | 5712 |
| katsura-5 | f4 | <0.001 | 0.017 | 0.016720 | 10281 | 0.008537 | 3416 | 3736 |
| katsura-6 | classic | 0.038 | 1.037 | 0.999820 | 4289976 | 0.494796 | 3636 | 26680 |
| katsura-6 | f4 | 0.001 | 0.108 | 0.106423 | 48731 | 0.082451 | 3568 | 4916 |
| katsura-7 | classic | 0.352 | 24.7 | 24.392141 | 49045544 | 10.820062 | 4252 | 231644 |
| katsura-7 | f4 | 0.011 | 1.009 | 0.998060 | 235767 | 0.738677 | 3860 | 10496 |
| katsura-8 | classic | 4.289 | DNF | ? |  |  | 6876 |  |
| katsura-8 | f4 | 0.042 | 11.0 | 10.985526 | 1578004 | 7.723794 | 6080 | 45080 |
| katsura-9 | classic | 55.0 | SKIP | ? |  |  | 18056 |  |
| katsura-9 | f4 | 0.252 | DNF | ? |  |  | 15288 |  |
| katsura-10 | classic | DNF | SKIP | ? |  |  |  |  |
| katsura-10 | f4 | 1.862 | DNF | ? |  |  | 48392 |  |
| eco-8 | classic | 4.851 | 32.7 | 27.880073 | 22000522 | 1.562047 | 4344 | 107240 |
| eco-8 | f4 | 0.003 | 0.325 | 0.322393 | 118909 | 0.228960 | 3660 | 6984 |
| eco-9 | classic | DNF | DNF | ? |  |  |  |  |
| eco-9 | f4 | 0.014 | 2.606 | 2.592016 | 535510 | 1.807176 | 4600 | 20200 |
| eco-10 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-10 | f4 | 0.066 | 25.4 | 25.354439 | 3466721 | 17.292106 | 7084 | 88308 |
| eco-11 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-11 | f4 | 0.436 | DNF | ? |  |  | 27056 |  |
| noon-3 | classic | <0.001 | <0.001 | 0.000550 | 12619 | 0.000300 | 3184 | 3528 |
| noon-3 | f4 | <0.001 | <0.001 | 0.000467 | 694 | 0.000140 | 3372 | 3416 |
| noon-4 | classic | 0.002 | 0.033 | 0.030673 | 329158 | 0.014666 | 3400 | 5648 |
| noon-4 | f4 | <0.001 | 0.011 | 0.010248 | 4337 | 0.003307 | 3448 | 3600 |
| noon-5 | classic | 0.063 | 3.217 | 3.154439 | 10855926 | 1.215056 | 3668 | 71364 |
| noon-5 | f4 | 0.002 | 0.041 | 0.039269 | 32875 | 0.026579 | 3760 | 4540 |
| noon-6 | classic | 3.201 | DNF | ? |  |  | 4976 |  |
| noon-6 | f4 | 0.011 | 0.747 | 0.735250 | 315964 | 0.493821 | 4864 | 11308 |

## Reading

At the gate's four cells, sylvester's `f4/default` config, one thread,
sits close to msolve, one thread: 0.96x on cyclic-7, 1.08x on katsura-9,
1.01x on eco-9, 0.97x on noon-6, a geometric mean of 1.01x against a 5x
budget and a worst cell of 1.08x against a 10x budget. The gate compares
one thread to one thread, and the F4 engine's batched, degree-by-degree
matrix construction is what produces this speed. The margin does not
hold at every size in the same families: at the largest instance of each,
sylvester trails msolve by more and Groebner.jl by more still. On
cyclic-8, sylvester is 1.524 s against msolve's 1.215 s (about 1.25x)
and Groebner.jl's 0.595 s (about 2.56x). On katsura-10, sylvester is
1.862 s against msolve's 1.636 s (about 1.14x) and Groebner.jl's 0.399 s
(about 4.67x); Macaulay2 is DNF there, and sylvester finishes. On
eco-11, sylvester is 0.436 s against msolve's 0.429 s (about 1.02x) and
Groebner.jl's 0.215 s (about 2.03x); Macaulay2 is DNF there too. Noon
tops out at noon-6, already a gate cell.

The `f4/threads8` column is not part of the gate: design 1.1 measures
one thread against msolve's one thread. It is reported because it fits
the same `taskset -c 0-3,12-15` pinning. Against `f4/default`, it
roughly halves the larger cells: cyclic-8 goes from 1.524 s to 0.762 s,
and katsura-9 from 0.252 s to 0.110 s.

Verifying a `sylv-gb-cert-v2` certificate costs far more than the raw
F4 computation it certifies, and the gap widens as the instance grows.
The verifier redoes every step with its own arithmetic, sharing no code
with the engine, so its cost tracks the size of the trace, not the
engine's. On the eco family, the standalone `verify (s)` column goes
from 0.228960 s on eco-8 (raw 0.003 s, about 76x) to 1.807176 s on
eco-9 (raw 0.014 s, about 129x) to 17.292106 s on eco-10 (raw 0.066 s,
about 262x). The limit is sharp: katsura-9 finishes raw F4 in
0.252204 s but does not finish the certified call inside 120 s, and the
classic `sylv-gb-cert-v1` path stops finishing certified runs at
katsura-8.

Since the 2026-08-18 v0.1 record, the sylvester cells moved with the
new F4 backend, the default as of 0.2.0, which replaces the removed
matrix backend. On katsura-9, the v0.1 record's best backend was
matrix+parallel at 19.9 s; `f4/default` now takes 0.252 s, about 79x
faster. On noon-6, the v0.1 record's best backend was classic at
3.184 s; `f4/default` now takes 0.011 s, about 289x faster. Both
cyclic-7 and eco-9 were DNF for every v0.1 backend; both finish now, in
0.075 s and 0.014 s. The v0.1 record measured msolve by subtracting a
process-startup baseline from wall clock; this record measures msolve
in-process, forking a child per repetition (see Protocol). The two
methods are not the same measurement, so a msolve number in one record
is not comparable to a msolve number in the other, and neither is a
cross-record ratio to msolve.

The per-commit optimization measurements behind this v0.2 record are
not published as a record: the Macaulay2 timings taken during that
pass ran on the machine's E-cores, off the `taskset -c 0-3,12-15`
protocol this report uses, so they are not comparable to the numbers
above. That is a limit on what this report can show about where the
time went commit by commit, not on the record's own numbers.

## No new defects

This run found no new defect; `KNOWN_ISSUES.md` records the defects
found and fixed before it.
