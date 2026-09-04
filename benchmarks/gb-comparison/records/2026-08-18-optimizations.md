# v0.1 optimization record (2026-08-18)

Per-commit wall times of the two backends during the v0.1 optimization
pass. Every row is the runner `sylv-runner` on `inputs/<name>.syl`, p =
1073741827, one thread, `taskset -c 16-31` (E-cores of an i9-13900KS),
minimum of 3 in-process runs, release build. Output equality (`SIZE` and
every `LM` line) held against the pre-pass baseline at every commit.
Cores and profile differ from the comparison record, so the numbers here
are for relative change only.

## Classic backend and polynomial core

Commits, in order: hold the shifted monomial across `sub_scaled`
iterations (c1); widen the inline exponent array to eleven variables (c2);
keep field addition and small-modulus multiplication in 64 bits (c3);
build the reduction remainder descending (c4); check the monomial width
in debug builds only (c5); skip the coefficient multiply when the scalar
is one (c6); compare the F5 signature without building the product and
skip trivial inversions (c7); reject a divisor by total degree first (c8).

| input / backend | base | c1 | c2 | c3 | c4 | c5 | c6 | c7 | c8 | total |
|---|---|---|---|---|---|---|---|---|---|---|
| katsura-7 classic | 0.848 | 0.563 | 0.555 | 0.560 | 0.534 | 0.495 | 0.495 | 0.478 | 0.468 | 1.81x |
| katsura-7 matrix | 0.780 | 0.808 | 0.782 | 0.764 | 0.770 | 0.725 | 0.702 | 0.709 | 0.637 | 1.22x |
| katsura-8 classic | 23.584 | 16.936 | 6.833 | 6.811 | 6.609 | 6.112 | 6.100 | 5.952 | 5.750 | 4.10x |
| katsura-8 matrix | 10.762 | 11.167 | 10.058 | 9.598 | 10.002 | 9.209 | 9.088 | 9.248 | 8.258 | 1.30x |
| cyclic-6 classic | 0.567 | 0.399 | 0.397 | 0.393 | 0.387 | 0.352 | 0.350 | 0.335 | 0.324 | 1.75x |
| cyclic-6 matrix | 1.128 | 1.141 | 1.129 | 1.073 | 1.122 | 1.035 | 1.039 | 1.026 | 0.942 | 1.20x |
| noon-5 classic | 0.135 | 0.076 | 0.077 | 0.079 | 0.077 | 0.071 | 0.071 | 0.070 | 0.069 | 1.96x |
| noon-5 matrix | 2.757 | 2.928 | 2.904 | 2.832 | 2.823 | 2.744 | 2.750 | 2.627 | 2.653 | 1.04x |
| eco-8 classic | 11.774 | 6.917 | 6.849 | 6.846 | 6.714 | 6.324 | 6.300 | 6.234 | 6.079 | 1.94x |

Peak RSS on katsura-8: classic 9724 KB to 6848 KB, matrix 77404 KB to
47656 KB.

Inline width study (katsura-8 classic, after c1): `SmallVec<[u16; 8]>`
16.94 s, `[u16; 11]` 6.83 s, `[u16; 12]` 7.41 s, `[u16; 16]` 6.99 s.
Eleven is the largest width that keeps `size_of::<Monomial>()` at 40
bytes.

Measured and not landed: in-place `sub_scaled` into a reusable buffer,
0 to 2.9% on classic, inside build noise.

## Matrix backend

Measured on a branch from the pre-pass tree, before the classic commits
above were merged. Commits, in order: reduce rows in a dense accumulator
against a creation-ordered pivot list (m1); index the basis leads in
symbolic preprocessing (m2); build the remainder in order and compare
reducer signatures in place (m3); lower the parallel row threshold to 3
(m4, changes only the `parallel` build).

| input | base | m1 | m2 | m3 | final serial | final `parallel` |
|---|---|---|---|---|---|---|
| katsura-7 | 0.783 | 0.584 | 0.561 | 0.468 | 0.465 (1.68x) | 0.380 (2.06x) |
| katsura-8 | 10.842 | 7.368 | 7.124 | 6.271 | 6.117 (1.77x) | 4.661 (2.33x) |
| cyclic-6 | 1.125 | 0.814 | 0.774 | 0.716 | 0.715 (1.57x) | 0.493 (2.28x) |
| noon-5 | 2.775 | 0.635 | 0.607 | 0.565 | 0.564 (4.92x) | 0.450 (6.17x) |
| eco-8 | DNF | DNF | DNF | DNF | DNF | DNF |

Peak RSS on katsura-8: 77336 KB to 68176 KB serial, 70512 KB parallel.
Parallel threshold sweep over the four inputs: 3 beat 2, 4, 6, 8, and 16.

## Merged tree

After both streams merged and the helpers were shared (`cmp_shifted`,
`shift_monomial`, `var_mask`, `Signature::shifted_is_below`), plus the
classic degree-and-mask reducer filter:

| input | classic | matrix |
|---|---|---|
| katsura-7 | 0.459 | 0.321 |
| katsura-8 | 5.551 | 3.430 |
| cyclic-6 | 0.318 | 0.451 |
| noon-5 | 0.069 | 0.457 |
| eco-8 | 5.891 | DNF |

Against the pre-pass baseline: katsura-8 classic 23.584 to 5.551 (4.2x),
katsura-8 matrix 10.762 to 3.430 (3.1x), noon-5 matrix 2.757 to 0.457
(6.0x).

Remaining top symbols on katsura-8: classic `Polynomial::sub_scaled` 54%
and `SmallVec::extend` (monomial construction) 39%; matrix `reduce_row`
42% and post-elimination `f5_reduce` 39%. Both are the monomial
representation and the reduction structure that v0.2 replaces.
