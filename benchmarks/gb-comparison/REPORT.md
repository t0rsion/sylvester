# Gröbner engine benchmark: sylvester vs Singular, msolve, Macaulay2, Groebner.jl

Run 2026-09-18. The refreshed cells were measured from clean source commit
`68f75768b9acf255d0367f80fdad177fd3ed13d8` with sylvester version `0.4.0+68f7576`.
The release worktree may contain later interface and prose changes. Their
relationship to the measured source is recorded below.

## Scope and provenance

The record contains 228 cells. It remeasures 100 prime-field sylvester cells,
12 rational sylvester cells, and the four prime-field msolve gate cells. It
reuses 112 external cells from the frozen 2026-08-19 record. The retained cell
objects match the frozen record exactly. No frozen record was edited.

The 116 target cells contain 99 OK, 6 DNF, 7 SKIP, and 4 ERROR statuses. The
four ERROR statuses are writer-wrapped deadline results described below.

The measured source is detached commit `68f75768b9acf255d0367f80fdad177fd3ed13d8`. The source
record is `results.json.zst` from run `2026-08-19`, with sylvester
metadata `0.3.0+517ca24`. Its compressed SHA-256 is
`4ca4190eca309a58d039adbbfa2dca0987afdb2429bd395ac3d23312f141fb69` and its decoded JSON SHA-256 is
`8069269f2a67566916eb2b510dc068e0f535e396be16b0b275b9707ab43e5cb9`. The new compressed and decoded record hashes
are listed in `SHA256SUMS`. `_meta.provenance.source_engine_sha256` stores the
SHA-256 of every recorded engine, ring, polynomial, rational-check, quotient,
and verifier source file.

The final worktree audit found matching hashes for the compute, ring, polynomial,
normal-form, rational-check, and verifier paths. `src/lib.rs` and
`src/quotient/linear.rs` differ from the measured source because later changes
added result interfaces and moved a test module. Those changes are outside the
Gröbner benchmark hot path. The measured source commit remains the provenance
for every refreshed timing. The release runner also classifies writer timeouts
as DNF; the measured runner reported ERROR.

An unrelated CPU-bound process had affinity spanning the benchmark
CPUs during the full refresh. The full cell values remain preserved in this
record. The release gate uses the clean rerun in
`records/2026-09-18-clean-gate/`, completed after that process was stopped.

The current run metadata records Singular `4.4.1`,
msolve `0.10.1`, Macaulay2 `1.26.05`,
Julia `1.13.0`, no installed Groebner.jl package, and
sylvester `0.4.0+68f7576`. The retained external cells,
including Groebner.jl, retain the historical source versions: Julia
`1.12.7` and Groebner.jl
`0.10.3`. This distinction prevents the
current environment metadata from being read as the provenance of reused cells.

## Protocol

- Field: `F_1073741827`. Order: grevlex. Output: reduced Gröbner basis.
- The driver ran under `taskset -c 0-3,12-15`. Runner binaries were built on
  E-cores with `taskset -c 16-31`. The memory cap was 16 GB.
- Each computation had a 120 s deadline. DNF records the typed computation
  deadline. Larger members in the same family then receive SKIP.
- External cells are retained from the frozen record. Their protocol was a
  median of three internal-timer runs. msolve prime-field timings use the
  in-process runner and its repeated direct `core_msolve` calls.
- Sylvester raw cells use one thread except `f4/threads8`, which requests eight
  threads. Certified cells include computation, certificate writing, and the
  bundled verification. A later standalone verifier timing is recorded.
- The release gate is the `f4/default` result against msolve on the four cells
  named by section 14 of `docs/f4-design.md`. The rational cells have no speed
  gate.

## Correctness

The primary check compares complete canonical bases with coefficients included.
Basis size and leading monomials remain weaker diagnostics. Every compared
prime-field basis agrees.

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
| katsura-9 | 8 | AGREE | AGREE |
| katsura-10 | 5 | AGREE | AGREE |
| eco-8 | 9 | AGREE | AGREE |
| eco-9 | 7 | AGREE | AGREE |
| eco-10 | 7 | AGREE | AGREE |
| eco-11 | 5 | AGREE | AGREE |
| noon-3 | 9 | AGREE | AGREE |
| noon-4 | 9 | AGREE | AGREE |
| noon-5 | 9 | AGREE | AGREE |
| noon-6 | 8 | AGREE | AGREE |

Every compared rational basis also agrees. The rational comparison uses exact
fractions, monic normalization, coefficients, and the full canonical basis.

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

## Release gate

The gate uses `f4/default` and one thread. It compares sylvester against the
four msolve cells and applies the geometric-mean limit of 5.0x and the
single-cell limit of 10.0x. The values below come from the clean eight-cell
rerun in `records/2026-09-18-clean-gate/`; the full refresh values remain in
the raw configuration table.

| cell | sylvester raw | sylvester (s) | msolve (s) | ratio |
|---|---|---|---|---|
| cyclic-7 | f4/default | 0.076606215 | 0.078391977 | 0.98x |
| katsura-9 | f4/default | 0.252875262 | 0.229818113 | 1.10x |
| eco-9 | f4/default | 0.015791203 | 0.013231605 | 1.19x |
| noon-6 | f4/default | 0.011724404 | 0.011754756 | 1.00x |

geometric mean: 1.06x (gate: at most 5.0x)
worst cell: 1.19x (gate: at most 10.0x)
gate: PASS

## Timings

The generated tables below contain the retained external cells and the refreshed
sylvester cells. A sylvester value in the main table is the fastest one-thread
raw configuration for that instance.

### Main comparison

| instance | Singular | msolve | M2 | Groebner.jl | sylvester |
|---|---|---|---|---|---|
| cyclic-4 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (classic/default) |
| cyclic-5 | <0.001 | <0.001 | 0.001 | <0.001 | <0.001 (f4/default) |
| cyclic-6 | 0.006 | 0.002 | 0.015 | <0.001 | 0.002 (f4/default) |
| cyclic-7 | 0.769 | 0.078 | 2.083 | 0.023 | 0.075 (f4/default) |
| cyclic-8 | 22.0 | 1.228 | 77.6 | 0.590 | 1.507 (f4/default) |
| katsura-4 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (f4/default) |
| katsura-5 | 0.001 | <0.001 | 0.002 | <0.001 | <0.001 (f4/default) |
| katsura-6 | 0.007 | 0.002 | 0.020 | <0.001 | 0.001 (f4/default) |
| katsura-7 | 0.067 | 0.007 | 0.201 | 0.002 | 0.006 (f4/default) |
| katsura-8 | 0.554 | 0.040 | 2.110 | 0.027 | 0.038 (f4/default) |
| katsura-9 | 5.110 | 0.231 | 26.8 | 0.149 | 0.251 (f4/default) |
| katsura-10 | 38.0 | 1.641 | DNF | 0.398 | 1.895 (f4/default) |
| eco-8 | 0.017 | 0.003 | 0.058 | 0.002 | 0.003 (f4/default) |
| eco-9 | 0.145 | 0.013 | 0.655 | 0.023 | 0.012 (f4/default) |
| eco-10 | 1.678 | 0.070 | 8.150 | 0.121 | 0.064 (f4/default) |
| eco-11 | 17.7 | 0.428 | DNF | 0.219 | 0.434 (f4/default) |
| noon-3 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (f4/default) |
| noon-4 | <0.001 | 0.002 | <0.001 | <0.001 | <0.001 (f4/default) |
| noon-5 | 0.005 | 0.002 | 0.006 | 0.001 | 0.002 (f4/default) |
| noon-6 | 0.025 | 0.012 | 0.098 | 0.008 | 0.011 (f4/default) |

### Sylvester raw configurations

| instance | f4/default | f4/threads8 | classic/default |
|---|---|---|---|
| cyclic-4 | <0.001 | <0.001 | <0.001 |
| cyclic-5 | <0.001 | <0.001 | 0.003 |
| cyclic-6 | 0.002 | 0.002 | 0.304 |
| cyclic-7 | 0.075 | 0.050 | DNF |
| cyclic-8 | 1.507 | 0.758 | SKIP |
| katsura-4 | <0.001 | <0.001 | <0.001 |
| katsura-5 | <0.001 | <0.001 | 0.004 |
| katsura-6 | 0.001 | 0.001 | 0.041 |
| katsura-7 | 0.006 | 0.004 | 0.447 |
| katsura-8 | 0.038 | 0.020 | 5.370 |
| katsura-9 | 0.251 | 0.110 | 68.2 |
| katsura-10 | 1.895 | 0.685 | DNF |
| eco-8 | 0.003 | 0.003 | 5.310 |
| eco-9 | 0.012 | 0.008 | DNF |
| eco-10 | 0.064 | 0.037 | SKIP |
| eco-11 | 0.434 | 0.191 | SKIP |
| noon-3 | <0.001 | <0.001 | <0.001 |
| noon-4 | <0.001 | <0.001 | 0.002 |
| noon-5 | 0.002 | 0.002 | 0.060 |
| noon-6 | 0.011 | 0.010 | 3.224 |

### Rational comparison

| instance | Singular | msolve | Groebner.jl | sylvester |
|---|---|---|---|---|
| cyclic-4-q | <0.001 | 0.039 | <0.001 | <0.001 |
| cyclic-5-q | <0.001 | 0.042 | 0.001 | 0.002 |
| cyclic-6-q | 0.026 | 0.031 | 0.004 | 0.019 |
| katsura-4-q | 0.001 | 0.036 | <0.001 | 0.002 |
| katsura-5-q | 0.005 | 0.016 | 0.002 | 0.009 |
| katsura-6-q | 0.032 | 0.041 | 0.005 | 0.050 |
| katsura-7-q | 0.295 | 0.048 | 0.030 | 0.290 |
| eco-8-q | 0.030 | 0.034 | 0.008 | 0.040 |
| eco-9-q | 0.433 | 0.060 | 0.039 | 0.192 |
| noon-3-q | <0.001 | 0.038 | <0.001 | <0.001 |
| noon-4-q | 0.001 | 0.027 | 0.001 | 0.003 |
| noon-5-q | 0.005 | 0.041 | 0.006 | 0.018 |

## Certified runs

The certified table reports both certificate contracts. Classic writes
`sylv-gb-cert-v1`; F4 writes `sylv-gb-cert-v2`. `split` is certified total
minus the same backend's raw time. It includes certificate writing and bundled
verification. `verify` is a later standalone verification.

split = certified total - raw. verify is a later standalone check.

| instance | backend | raw (s) | certified total (s) | split (s) | cert bytes | verify (s) | peak RSS raw (KB) | peak RSS certified (KB) |
|---|---|---|---|---|---|---|---|---|
| cyclic-4 | classic | <0.001 | <0.001 | 0.000103 | 3076 | 0.000045 | 3480 | 3724 |
| cyclic-4 | f4 | <0.001 | <0.001 | 0.000083 | 481 | 0.000029 | 3684 | 3728 |
| cyclic-5 | classic | 0.003 | 0.017 | 0.014834 | 166879 | 0.003140 | 3720 | 5128 |
| cyclic-5 | f4 | <0.001 | 0.001 | 0.000967 | 4529 | 0.000529 | 3912 | 4068 |
| cyclic-6 | classic | 0.304 | 3.797 | 3.493243 | 3968988 | 0.109628 | 4268 | 27676 |
| cyclic-6 | f4 | 0.002 | 0.012 | 0.010315 | 30041 | 0.006298 | 4164 | 4672 |
| cyclic-7 | classic | DNF | DNF | ? |  |  |  |  |
| cyclic-7 | f4 | 0.075 | 1.222 | 1.146356 | 784454 | 0.735922 | 9172 | 18132 |
| cyclic-8 | classic | SKIP | SKIP | ? |  |  |  |  |
| cyclic-8 | f4 | 1.507 | 16.2 | 14.714220 | 5549896 | 11.082573 | 60304 | 120296 |
| katsura-4 | classic | <0.001 | 0.003 | 0.002640 | 48198 | 0.001006 | 3420 | 4160 |
| katsura-4 | f4 | <0.001 | <0.001 | 0.000675 | 2725 | 0.000337 | 3964 | 4084 |
| katsura-5 | classic | 0.004 | 0.041 | 0.036872 | 376992 | 0.014854 | 3848 | 6080 |
| katsura-5 | f4 | <0.001 | 0.005 | 0.005036 | 10281 | 0.002772 | 4092 | 4288 |
| katsura-6 | classic | 0.041 | 0.905 | 0.864169 | 4289976 | 0.339305 | 4052 | 27116 |
| katsura-6 | f4 | 0.001 | 0.053 | 0.051489 | 48731 | 0.027316 | 4120 | 5412 |
| katsura-7 | classic | 0.447 | 20.7 | 20.221176 | 49045544 | 6.603828 | 4624 | 233340 |
| katsura-7 | f4 | 0.006 | 0.524 | 0.517842 | 235767 | 0.263846 | 4636 | 11072 |
| katsura-8 | classic | 5.370 | ERROR | ? |  |  | 7180 |  |
| katsura-8 | f4 | 0.038 | 6.070 | 6.031715 | 1578004 | 2.898963 | 6272 | 45684 |
| katsura-9 | classic | 68.2 | DNF | ? |  |  | 17740 |  |
| katsura-9 | f4 | 0.251 | 66.5 | 66.218055 | 9459218 | 30.467932 | 15896 | 222224 |
| katsura-10 | classic | DNF | SKIP | ? |  |  |  |  |
| katsura-10 | f4 | 1.895 | ERROR | ? |  |  | 50800 |  |
| eco-8 | classic | 5.310 | 33.3 | 27.984966 | 22000522 | 1.041344 | 4712 | 107468 |
| eco-8 | f4 | 0.003 | 0.163 | 0.160476 | 118909 | 0.075703 | 4284 | 7508 |
| eco-9 | classic | DNF | DNF | ? |  |  |  |  |
| eco-9 | f4 | 0.012 | 1.418 | 1.405842 | 535510 | 0.630518 | 5116 | 20960 |
| eco-10 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-10 | f4 | 0.064 | 14.5 | 14.481155 | 3466721 | 6.376309 | 7320 | 88780 |
| eco-11 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-11 | f4 | 0.434 | ERROR | ? |  |  | 27476 |  |
| noon-3 | classic | <0.001 | <0.001 | 0.000542 | 12619 | 0.000267 | 3572 | 3904 |
| noon-3 | f4 | <0.001 | <0.001 | 0.000137 | 694 | 0.000075 | 3812 | 3932 |
| noon-4 | classic | 0.002 | 0.027 | 0.025394 | 329158 | 0.009988 | 3716 | 6008 |
| noon-4 | f4 | <0.001 | 0.002 | 0.001349 | 4337 | 0.000701 | 4024 | 4088 |
| noon-5 | classic | 0.060 | 2.603 | 2.543011 | 10855926 | 0.606831 | 4160 | 71304 |
| noon-5 | f4 | 0.002 | 0.023 | 0.020832 | 32875 | 0.009878 | 4236 | 4976 |
| noon-6 | classic | 3.224 | ERROR | ? |  |  | 5372 |  |
| noon-6 | f4 | 0.011 | 0.415 | 0.403881 | 315964 | 0.166808 | 5312 | 10744 |

## Certified deadline outcomes

Four certified cells reached the 120 s writer deadline. The runner emitted
`STATUS ERROR the certificate writer stopped: the computation passed its
deadline` for each cell:

| cell | raw status | cell wall (s) | frozen status |
|---|---|---:|---|
| katsura-8, classic certified | ERROR | 120.030 | DNF |
| katsura-10, F4 certified | ERROR | 120.031 | SKIP |
| eco-11, F4 certified | ERROR | 120.033 | DNF |
| noon-6, classic certified | ERROR | 120.050 | DNF |

The underlying error is `CertifyError::WriterExhausted(ComputeError::Timeout)`, a typed
resource deadline. The runner at the measured source maps engine timeouts and
verifier exhaustion to `STATUS TIMEOUT`, but its certified match does not yet
map the writer-wrapped timeout. The driver therefore preserves these four raw
`ERROR` records. They produced no basis and caused no correctness mismatch.
The frozen statuses show the same resource boundary in the previous runner's
classification. These cells are resource ceilings, not engine basis defects.
