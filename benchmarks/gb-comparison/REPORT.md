# v0.2 benchmark record

Run 2026-08-27. Sylvester 0.2.0, commit `caf634b`.

## Provenance

This run measures every sylvester cell and the four msolve gate cells.
Other external cells reuse the exact 2026-08-19 record at commit `f413b0b`.
The archived source is `records/2026-08-19-v0.2/results.json.zst`.

`results.json.zst` stores commands, counters, and complete bases. Its `_meta`
object records the reuse explicitly. `results.csv` stores the scalar cells.

## Protocol

- Field: $\mathbb F_{1073741827}$.
- Order: grevlex.
- Output: reduced Gröbner basis.
- Affinity: `taskset -c 0-3,12-15`.
- Memory: 16 GiB per process, no swap.
- Deadline: 120 s per computation.
- Gate threads: one for sylvester and msolve.
- Sylvester samples: median of three successful runs.
- External samples: each tool's recorded internal compute time.

Versions: Singular 4.4.1, msolve 0.10.1, Macaulay2 1.26.05, Julia 1.12.7,
and Groebner.jl 0.10.3.

## Correctness

The primary check compares each complete canonical basis, including
coefficients. Every complete basis agrees.

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

## Release gate

The gate compares `f4/default` with msolve on four cells. It requires a
geometric mean at most 5.0x and a worst cell at most 10.0x.

| cell | sylvester raw | sylvester (s) | msolve (s) | ratio |
|---|---|---|---|---|
| cyclic-7 | f4/default | 0.075122 | 0.076398 | 0.98x |
| katsura-9 | f4/default | 0.253319 | 0.230088 | 1.10x |
| eco-9 | f4/default | 0.013627 | 0.013244 | 1.03x |
| noon-6 | f4/default | 0.017326 | 0.011722 | 1.48x |

Geometric mean: 1.13x. Worst cell: 1.48x. Gate: PASS.

## Main results

Seconds. The sylvester column selects the faster one-thread raw backend.

| instance | Singular | msolve | M2 | Groebner.jl | sylvester |
|---|---|---|---|---|---|
| cyclic-4 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (classic/default) |
| cyclic-5 | 0.001 | <0.001 | 0.001 | <0.001 | <0.001 (f4/default) |
| cyclic-6 | 0.010 | 0.002 | 0.015 | 0.001 | 0.002 (f4/default) |
| cyclic-7 | 0.812 | 0.076 | 2.087 | 0.023 | 0.075 (f4/default) |
| cyclic-8 | 21.9 | 1.215 | 78.1 | 0.595 | 1.529 (f4/default) |
| katsura-4 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (f4/default) |
| katsura-5 | 0.001 | 0.002 | 0.002 | <0.001 | <0.001 (f4/default) |
| katsura-6 | 0.007 | 0.002 | 0.020 | <0.001 | 0.004 (f4/default) |
| katsura-7 | 0.071 | 0.007 | 0.201 | 0.002 | 0.017 (f4/default) |
| katsura-8 | 0.579 | 0.038 | 2.107 | 0.028 | 0.040 (f4/default) |
| katsura-9 | 5.101 | 0.230 | 26.7 | 0.149 | 0.253 (f4/default) |
| katsura-10 | 38.2 | 1.636 | DNF | 0.399 | 1.875 (f4/default) |
| eco-8 | 0.015 | 0.003 | 0.058 | 0.002 | 0.007 (f4/default) |
| eco-9 | 0.145 | 0.013 | 0.654 | 0.022 | 0.014 (f4/default) |
| eco-10 | 1.673 | 0.069 | 8.053 | 0.120 | 0.068 (f4/default) |
| eco-11 | 17.7 | 0.429 | DNF | 0.215 | 0.432 (f4/default) |
| noon-3 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (f4/default) |
| noon-4 | <0.001 | <0.001 | <0.001 | <0.001 | <0.001 (f4/default) |
| noon-5 | 0.003 | 0.002 | 0.006 | 0.001 | 0.002 (f4/default) |
| noon-6 | 0.023 | 0.012 | 0.098 | 0.008 | 0.017 (f4/default) |

## Sylvester raw results

| instance | f4/default | f4/threads8 | classic/default |
|---|---|---|---|
| cyclic-4 | <0.001 | <0.001 | <0.001 |
| cyclic-5 | <0.001 | <0.001 | 0.003 |
| cyclic-6 | 0.002 | 0.002 | 0.224 |
| cyclic-7 | 0.075 | 0.050 | DNF |
| cyclic-8 | 1.529 | 0.766 | SKIP |
| katsura-4 | <0.001 | <0.001 | <0.001 |
| katsura-5 | <0.001 | <0.001 | 0.004 |
| katsura-6 | 0.004 | 0.002 | 0.028 |
| katsura-7 | 0.017 | 0.007 | 0.303 |
| katsura-8 | 0.040 | 0.021 | 3.668 |
| katsura-9 | 0.253 | 0.110 | 44.4 |
| katsura-10 | 1.875 | 0.696 | DNF |
| eco-8 | 0.007 | 0.008 | 3.490 |
| eco-9 | 0.014 | 0.013 | DNF |
| eco-10 | 0.068 | 0.038 | SKIP |
| eco-11 | 0.432 | 0.194 | SKIP |
| noon-3 | <0.001 | <0.001 | <0.001 |
| noon-4 | <0.001 | <0.001 | 0.002 |
| noon-5 | 0.002 | 0.009 | 0.049 |
| noon-6 | 0.017 | 0.015 | 2.204 |

## Certified results

Classic writes `sylv-gb-cert-v1`. F4 writes `sylv-gb-cert-v2`. Certified
time covers computation, writing, and bundled verification. `verify` is a
later standalone replay. `split` is certified time minus raw time.

| instance | backend | raw (s) | certified total (s) | split (s) | cert bytes | verify (s) | peak RSS raw (KB) | peak RSS certified (KB) |
|---|---|---|---|---|---|---|---|---|
| cyclic-4 | classic | <0.001 | <0.001 | 0.000120 | 3076 | 0.000045 | 3148 | 3316 |
| cyclic-4 | f4 | <0.001 | <0.001 | 0.000106 | 481 | 0.000037 | 3384 | 3300 |
| cyclic-5 | classic | 0.003 | 0.019 | 0.015797 | 166879 | 0.003049 | 3336 | 4412 |
| cyclic-5 | f4 | <0.001 | 0.003 | 0.003036 | 4529 | 0.002705 | 3412 | 3632 |
| cyclic-6 | classic | 0.224 | 3.142 | 2.917750 | 3968988 | 0.107886 | 3740 | 27292 |
| cyclic-6 | f4 | 0.002 | 0.025 | 0.023440 | 30041 | 0.019028 | 3616 | 4112 |
| cyclic-7 | classic | DNF | DNF | ? |  |  |  |  |
| cyclic-7 | f4 | 0.075 | 2.440 | 2.364835 | 784454 | 1.991033 | 8756 | 17528 |
| cyclic-8 | classic | SKIP | SKIP | ? |  |  |  |  |
| cyclic-8 | f4 | 1.529 | 36.6 | 35.091894 | 5549896 | 31.742164 | 62156 | 116736 |
| katsura-4 | classic | <0.001 | 0.003 | 0.002363 | 48198 | 0.001000 | 3152 | 3744 |
| katsura-4 | f4 | <0.001 | 0.002 | 0.002160 | 2725 | 0.001249 | 3360 | 3536 |
| katsura-5 | classic | 0.004 | 0.041 | 0.036388 | 376992 | 0.014462 | 3352 | 5688 |
| katsura-5 | f4 | <0.001 | 0.017 | 0.016692 | 10281 | 0.008464 | 3424 | 3708 |
| katsura-6 | classic | 0.028 | 0.823 | 0.795111 | 4289976 | 0.334056 | 3552 | 25020 |
| katsura-6 | f4 | 0.004 | 0.103 | 0.098660 | 48731 | 0.080606 | 3548 | 4888 |
| katsura-7 | classic | 0.303 | 19.5 | 19.150033 | 49045544 | 6.510233 | 4296 | 231164 |
| katsura-7 | f4 | 0.017 | 0.960 | 0.943297 | 235767 | 0.724900 | 4056 | 10524 |
| katsura-8 | classic | 3.668 | DNF | ? |  |  | 6820 |  |
| katsura-8 | f4 | 0.040 | 10.2 | 10.199963 | 1578004 | 7.658739 | 6132 | 44864 |
| katsura-9 | classic | 44.4 | SKIP | ? |  |  | 18048 |  |
| katsura-9 | f4 | 0.253 | DNF | ? |  |  | 15424 |  |
| katsura-10 | classic | DNF | SKIP | ? |  |  |  |  |
| katsura-10 | f4 | 1.875 | SKIP | ? |  |  | 51764 |  |
| eco-8 | classic | 3.490 | 24.5 | 21.002242 | 22000522 | 1.015310 | 4264 | 106416 |
| eco-8 | f4 | 0.007 | 0.301 | 0.293971 | 118909 | 0.226576 | 3720 | 6996 |
| eco-9 | classic | DNF | DNF | ? |  |  |  |  |
| eco-9 | f4 | 0.014 | 2.404 | 2.390630 | 535510 | 1.784429 | 4556 | 20200 |
| eco-10 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-10 | f4 | 0.068 | 23.8 | 23.749439 | 3466721 | 17.143897 | 7052 | 88468 |
| eco-11 | classic | SKIP | SKIP | ? |  |  |  |  |
| eco-11 | f4 | 0.432 | DNF | ? |  |  | 27128 |  |
| noon-3 | classic | <0.001 | 0.001 | 0.001146 | 12619 | 0.000503 | 3236 | 3468 |
| noon-3 | f4 | <0.001 | <0.001 | 0.000298 | 694 | 0.000129 | 3456 | 3476 |
| noon-4 | classic | 0.002 | 0.027 | 0.025055 | 329158 | 0.009775 | 3420 | 5576 |
| noon-4 | f4 | <0.001 | 0.003 | 0.002524 | 4337 | 0.001498 | 3380 | 3648 |
| noon-5 | classic | 0.049 | 2.419 | 2.369741 | 10855926 | 0.599611 | 3648 | 69752 |
| noon-5 | f4 | 0.002 | 0.040 | 0.038349 | 32875 | 0.026313 | 3752 | 4536 |
| noon-6 | classic | 2.204 | DNF | ? |  |  | 5140 |  |
| noon-6 | f4 | 0.017 | 0.673 | 0.655898 | 315964 | 0.487429 | 4824 | 11040 |

## Limits

The benchmark is evidence, not a proof of engine correctness. F4 certification
does not finish `katsura-9` or `eco-11` within 120 s. Classic certification
does not finish `katsura-8`, `eco-9`, or `noon-6` within 120 s.

This run found no new defect. `KNOWN_ISSUES.md` records the repaired defects and
the remaining limits.
