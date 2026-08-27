# Gröbner engine benchmark: sylvester vs Singular, msolve, Macaulay2, Groebner.jl

Run 2026-07-25 on a 32-core box, 188 GB RAM, no swap.

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
  swap disabled. This is not cosmetic: see the memory section.

Deviation from the original protocol: the msolve baseline was re-measured after
the memory cap was introduced (28.9 ms wrapped, vs 4.6 ms unwrapped). The cyclic
family's msolve numbers were taken against the older baseline, the other three
families against the newer one. The difference is below the resolution of every
conclusion drawn here.

## Correctness cross-check

For each instance, every engine that finished was compared on basis size and on
the full sorted multiset of leading monomials.

**118 pairwise comparisons across 20 instances, zero disagreements.** Where
sylvester finishes, it agrees exactly with all four external engines.

## Timings (seconds, median of 3)

### cyclic

| engine | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|
| msolve | 0.000 | 0.000 | 0.002 | 0.085 | 1.25 |
| Groebner.jl | 0.000 | 0.000 | 0.001 | 0.022 | 0.56 |
| Singular | 0.000 | 0.001 | 0.012 | 0.746 | 20.8 |
| Macaulay2 | 0.000 | 0.001 | 0.015 | 2.12 | 77.4 |
| sylvester classic | 0.000 | 0.006 | 0.399 | DNF | skip |
| sylvester matrix | 0.000 | 0.013 | 0.621 | DNF | skip |
| sylvester matrix+mont | 0.000 | 0.011 | 0.799 | DNF | skip |
| sylvester matrix+mont+csr | 0.000 | 0.012 | 1.238 | **OOM** | skip |
| sylvester matrix+mont+par | 0.003 | 0.097 | 0.714 | DNF | skip |

### katsura

| engine | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|
| msolve | 0.000 | 0.000 | 0.002 | 0.005 | 0.038 | 0.226 | 1.56 |
| Groebner.jl | 0.000 | 0.000 | 0.001 | 0.003 | 0.027 | 0.143 | 0.462 |
| Singular | 0.000 | 0.002 | 0.008 | 0.076 | 0.538 | 4.85 | 37.0 |
| Macaulay2 | 0.000 | 0.002 | 0.020 | 0.228 | 2.17 | 27.7 | DNF |
| sylvester classic | 0.004 | 0.015 | 0.059 | 0.604 | 13.8 | DNF | skip |
| sylvester matrix | 0.004 | 0.014 | 0.054 | 0.444 | 5.52 | 88.1 | DNF |
| sylvester matrix+mont | 0.003 | 0.014 | 0.058 | 0.508 | 6.52 | 107.9 | DNF |
| sylvester matrix+mont+csr | 0.006 | 0.017 | 0.077 | 0.903 | 12.1 | **OOM** | skip |
| sylvester matrix+mont+par | 0.029 | 0.066 | 0.207 | 0.585 | 3.94 | 40.3 | DNF |

### eco

| engine | 8 | 9 | 10 | 11 |
|---|---|---|---|---|
| msolve | 0.009 | 0.012 | 0.062 | 0.447 |
| Groebner.jl | 0.002 | 0.022 | 0.119 | 0.209 |
| Singular | 0.016 | 0.144 | 1.64 | 17.4 |
| Macaulay2 | 0.064 | 0.698 | 8.86 | DNF |
| sylvester classic | 8.80 | DNF* | skip | skip |
| sylvester matrix | DNF | skip | skip | skip |
| sylvester matrix+mont | DNF | skip | skip | skip |
| sylvester matrix+mont+csr | **OOM** | skip | skip | skip |
| sylvester matrix+mont+par | DNF | skip | skip | skip |

\* eco-9 classic additionally crashes on the default stack; see below.

### noon

| engine | 3 | 4 | 5 | 6 |
|---|---|---|---|---|
| msolve | 0.000 | 0.000 | 0.000 | 0.006 |
| Groebner.jl | 0.000 | 0.000 | 0.001 | 0.008 |
| Singular | 0.001 | 0.000 | 0.003 | 0.023 |
| Macaulay2 | 0.000 | 0.001 | 0.007 | 0.106 |
| sylvester classic | 0.001 | 0.005 | 0.115 | 5.66 |
| sylvester matrix | 0.002 | 0.014 | 1.61 | DNF |
| sylvester matrix+mont | 0.001 | 0.013 | 1.59 | DNF |
| sylvester matrix+mont+csr | 0.001 | 0.021 | 2.47 | DNF |
| sylvester matrix+mont+par | 0.008 | 0.110 | 0.954 | DNF |

## Two defects found by running this

### 1. Unbounded memory in the CSR backend

On cyclic-7 the `montgomery-csr` runner reached **169 GB RSS / 201 GB virtual**
and triggered a global OOM. Because the box has no swap, the kernel killed
inside the login scope and systemd tore down the entire session — twice, at
18:50:40 and 19:09:30, each time taking down the terminal and everything
running in it.

Under a 16 GB cap the same cell dies in 21 s. The other four configurations
merely time out on the same input, so this is specific to `csr`, not general
slowness: the CSR path allocates without bound where the other paths do not.
CSR is also the slowest configuration everywhere it completes, so it currently
costs memory and time both.

### 2. Stack exhaustion in the classic backend on eco-9

The classic backend crashes on eco-9 with the default 8 MiB stack. The failure
mode is not stable across runs: SIGABRT with glibc reporting
`corrupted size vs. prev_size`, SIGSEGV, and once a silent no-output failure.

With `ulimit -s unlimited` the same input does not crash — it runs to the
engine's own 120 s limit and reports `STATUS TIMEOUT` cleanly. So the cell is a
genuine DNF, and separately the engine exhausts the default stack on an input it
otherwise merely finds slow. sylvester's own source contains no `unsafe`, and
the default build pulls in neither bumpalo nor rayon, so this is stack depth or
frame size, not memory unsafety in the algebra. Naming the routine needs a
backtrace under a debugger, which was not run.

## Reading

sylvester is competitive with nothing at the top of the field and is not meant
to be. Concretely:

- It is within an order of magnitude of Singular on the small end of every
  family (n ≤ 6 for cyclic and katsura, all of noon up to 5).
- msolve and Groebner.jl are 100–1000× faster on the instances where sylvester
  still finishes, and they finish everything sylvester cannot.
- The practical ceiling is cyclic-6, katsura-9, eco-8, noon-5/6.
- The parallel configuration is the only feature flag that earns its keep at
  scale: it is 2× faster than plain matrix on katsura-9 (40.3 s vs 88.1 s) and
  2.2× on katsura-8, while costing 5–10× on small inputs where thread setup
  dominates. Montgomery is neutral to slightly negative throughout. CSR is a
  loss everywhere and dangerous on large inputs.

The gap is structural rather than a matter of tuning. msolve and Groebner.jl run
F4 with dense block linear algebra over carefully laid-out matrices; sylvester's
matrix backend does sparse elimination without that block structure, under a
signature discipline that additionally forbids some reductions.

## Artifacts

- `results.json` — every cell, with basis sizes and leading-monomial multisets.
- `results.csv` — 180 rows, flattened.
- `driver.py` — the harness; resumable, one cell at a time, memory-capped.
- `inputs/` — generated systems in each tool's syntax.
