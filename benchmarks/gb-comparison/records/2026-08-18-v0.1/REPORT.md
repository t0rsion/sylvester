# Gröbner engine benchmark: sylvester vs Singular, msolve, Macaulay2, Groebner.jl

Run 2026-08-18 on a 32-core box (Intel i9-13900KS), 188 GB RAM, no swap.

## Protocol

- Field F_p, p = 1073741827. Order: grevlex. Output: reduced Gröbner basis.
- Median of 3 runs per cell; a cell that exceeds 120 s is DNF, and larger
  instances in the same family are then skipped for that engine (timeouts are
  monotone in n across these families).
- Timing excludes parsing and startup wherever the tool exposes an internal
  timer: Singular `rtimer` around `std(I)`, M2 `elapsedTiming` around
  `groebnerBasis`, Julia `@elapsed` around `groebner()` in a warmed-up session,
  sylvester `Instant` around the compute call. msolve has no internal timer, so
  its figure is process wall clock minus a measured no-op baseline.
- Every child runs inside a transient systemd scope with `MemoryMax=16G` and
  swap disabled.

This run changes three things from the July protocol:

- Every process is pinned with `taskset -c 0-3,12-15`, the P-cores of the
  i9-13900KS.
- The msolve baseline is measured fresh for this run and printed by
  `report.py`: 0.00394 s. It replaces the 28.9 ms wrapped baseline of the
  July run.
- The sylvester runner configurations are `default` and `parallel` only.
  `montgomery` and `csr` no longer exist as features. The `parallel`
  configuration runs the matrix backend's rayon row elimination on the
  global pool over the 8 pinned logical CPUs.

The sylvester crate under test is v0.1.

## Correctness cross-check

For each instance, every engine that finished was compared on basis size and
on the full sorted multiset of leading monomials, exactly as in the July
report. `report.py` finds 20 instances, and every finisher agrees with every
other finisher on both.

Summing C(finishers, 2) over the 20 instances (the same count the July report
used) gives **307 pairwise comparisons, zero disagreements.** Where sylvester
finishes, it agrees exactly with all four external engines.

## Timings (seconds, median of 3)

### cyclic

| engine | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.002 | 0.086 | 1.245 |
| Groebner.jl | <0.001 | <0.001 | <0.001 | 0.024 | 0.577 |
| Singular | <0.001 | <0.001 | 0.006 | 0.771 | 22.0 |
| Macaulay2 | <0.001 | 0.001 | 0.015 | 2.089 | 78.1 |
| sylvester classic | <0.001 | 0.002 | 0.255 | DNF | skip |
| sylvester matrix | <0.001 | 0.004 | 0.272 | DNF | skip |
| sylvester matrix+parallel | <0.001 | 0.025 | 0.231 | DNF | skip |

### katsura

| engine | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| msolve | <0.001 | <0.001 | 0.001 | 0.007 | 0.053 | 0.247 | 1.810 |
| Groebner.jl | <0.001 | <0.001 | <0.001 | 0.002 | 0.027 | 0.144 | 0.401 |
| Singular | <0.001 | 0.001 | 0.007 | 0.067 | 0.573 | 5.119 | 38.0 |
| Macaulay2 | <0.001 | 0.002 | 0.020 | 0.203 | 2.099 | 26.9 | DNF |
| sylvester classic | <0.001 | 0.004 | 0.032 | 0.365 | 4.654 | 61.1 | DNF |
| sylvester matrix | <0.001 | 0.004 | 0.023 | 0.191 | 2.015 | 26.3 | DNF |
| sylvester matrix+parallel | 0.001 | 0.005 | 0.028 | 0.184 | 1.669 | 19.9 | DNF |

### eco

| engine | 8 | 9 | 10 | 11 |
|---|---|---|---|---|
| msolve | 0.029 | 0.029 | 0.082 | 0.473 |
| Groebner.jl | 0.002 | 0.022 | 0.115 | 0.214 |
| Singular | 0.016 | 0.146 | 1.685 | 17.8 |
| Macaulay2 | 0.058 | 0.654 | 8.215 | DNF |
| sylvester classic | 4.836 | DNF | skip | skip |
| sylvester matrix | DNF | skip | skip | skip |
| sylvester matrix+parallel | DNF | skip | skip | skip |

### noon

| engine | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| msolve | 0.002 | 0.001 | 0.004 | 0.030 |
| Groebner.jl | <0.001 | <0.001 | 0.001 | 0.008 |
| Singular | <0.001 | <0.001 | 0.005 | 0.023 |
| Macaulay2 | <0.001 | <0.001 | 0.006 | 0.098 |
| sylvester classic | <0.001 | 0.002 | 0.058 | 3.184 |
| sylvester matrix | <0.001 | 0.010 | 0.246 | DNF |
| sylvester matrix+parallel | <0.001 | 0.007 | 0.236 | 36.4 |

## Reading

sylvester still trails the external engines by orders of magnitude, and
neither backend finishes every instance the external engines finish.
Concretely:

- On katsura-9, the best sylvester configuration is matrix+parallel at
  19.9 s. msolve takes 0.247 s and Groebner.jl 0.144 s on the same instance,
  about 81x and 138x faster.
- The cyclic and eco families both DNF for every sylvester configuration
  from cyclic-7 and eco-9 onward; the smaller instances of each family
  finish.
- noon-6 finishes only in sylvester classic (3.184 s) and sylvester
  matrix+parallel (36.4 s). Plain matrix DNFs on noon-6.

### What changed since the July record

The July and this record used different core pinning (E-cores under
`taskset -c 16-31` for the July run, P-cores under `taskset -c 0-3,12-15`
here), so a cross-run ratio is approximate, not a clean speedup measurement.
With that caveat, the sylvester cells moved with the optimization pass
recorded in
`benchmarks/gb-comparison/records/2026-08-18-v0.1-optimizations.md`:

- katsura-8 matrix: 5.52 s in July, 2.015 s now.
- katsura-9: 88.1 s (matrix) and DNF (classic) in July; now 61.1 s classic,
  26.3 s matrix, 19.9 s matrix+parallel. Classic now finishes an instance it
  did not finish in July.
- noon-5 matrix: 1.61 s in July, 0.246 s now.
- eco-8 classic: 8.80 s in July, 4.836 s now.

## No new defects

This run found no new defect. The two resource-exhaustion defects the July
run found (unbounded memory in the removed `csr` path, stack exhaustion in
the classic backend on eco-9) are fixed and recorded in `KNOWN_ISSUES.md`.

## Artifacts

- `results.json`: every cell, with basis sizes and leading-monomial multisets.
- `results.csv`: the flattened form.
- `driver.py`: the harness; resumable, one cell at a time, memory-capped.
- `inputs/`: generated systems in each tool's syntax.
- `records/2026-07-25-archive/REPORT.md`: the July record, for historic
  numbers and the two resource-exhaustion defects it found.
- `records/2026-08-18-v0.1-optimizations.md`: the per-commit optimization
  measurements behind the sylvester numbers in this record.
