# Gröbner basis comparison

Timings and full-output cross-checks for sylvester against Singular, msolve,
Macaulay2, and Groebner.jl, on cyclic-4..8, katsura-4..10, eco-8..11, and
noon-3..6 over F_1073741827 under grevlex. A smaller rational cell set,
against Singular, msolve, and Groebner.jl (no Macaulay2) over `Q`, checks
sylvester's multimodular engine the same way; see "Rational protocol"
below.

The families are:

- `cyclic-4` through `cyclic-8`;
- `katsura-4` through `katsura-10`;
- `eco-8` through `eco-11`;
- `noon-3` through `noon-6`.

`REPORT.md` is the current record. `results.csv` contains its scalar cells.
`results.json.zst` contains commands, counters, and complete bases. Frozen
records stay under `records/`.

The current record contains 228 cells. It remeasures 116 targets: every
sylvester cell and the four msolve gate cells. It reuses 112 external cells
from the frozen 2026-08-19 record. `results.json.zst` records this provenance
in `_meta`. The gate values in `REPORT.md` come from the clean eight-cell
rerun in `records/2026-09-18-clean-gate/`.

Decode and verify the exact raw record:

```text
zstd -d results.json.zst -o results.json
sha256sum -c SHA256SUMS
```

## Reproduce

```text
python3 gen.py
./build-runners.sh
taskset -c 0-3,12-15 python3 driver.py
python3 report.py > report-fragments.txt
zstd -19 --long=27 -T8 results.json -o results.json.zst
```

`driver.py` writes each cell immediately and skips existing cells. Stop and
restart it without losing completed work. Stored commands use paths relative to
this directory. They contain no machine-specific directory names.

- `GBBENCH_INSTANCES=cyclic-4,katsura-5` restricts `driver.py` to that
  comma-separated instance list.
- `GBBENCH_RESULTS=/path/to/file.json` points `driver.py` (and `report.py`,
  which reads the same variable) at a results file other than
  `results.json`. `report.py` writes the CSV next to that file. A smoke
  run therefore does not touch the real record.
- `GBBENCH_RESULTS_CSV=/path/to/file.csv` overrides the CSV path alone.

Use these variables for a smoke run:

```text
GBBENCH_INSTANCES=cyclic-4,katsura-5
GBBENCH_RESULTS=/tmp/sylvester-smoke.json
GBBENCH_RESULTS_CSV=/tmp/sylvester-smoke.csv
```

`GBBENCH_MEMMAX` changes the per-process memory cap. `GBBENCH_CORES` changes
the affinity string stored in the record. It does not set affinity.

## Correctness

The primary check compares the complete canonical basis. A basis-size or
leading-monomial check is insufficient: equal leading ideals do not imply equal
ideals.

`canon.py` normalizes each tool's output:

1. Parse every polynomial into exponent-vector terms.
2. Reduce coefficients modulo \(p\).
3. Divide each polynomial by its leading coefficient.
4. Sort terms by descending grevlex.
5. Sort the basis independently of engine output order.

The report also prints basis-size and leading-monomial agreement as a weaker
diagnostic.

Singular, Macaulay2, msolve, and Groebner.jl print polynomial expressions.
Sylvester prints numeric `POLY` blocks. Both paths use the same canonical
representation.

## Timing

Each raw cell excludes input parsing when the tool exposes an internal timer.
Successful cells run three times and report the median. msolve is the exception,
described below.

| Tool | Timed operation |
| --- | --- |
| Singular | `rtimer` around `std(I)` |
| Macaulay2 | `elapsedTiming` around `groebnerBasis` |
| Groebner.jl | `@elapsed` around `groebner` after warmup |
| sylvester | `Instant` around the compute call |
| msolve | repeated direct `core_msolve` calls |

Each computation has a 120 s limit. A timeout skips the larger members of that
family for the same tool and configuration.

### msolve

One untimed `msolve -t 1 -g 2` call supplies the complete basis. Its process
wall time is not a timing cell.

`runner/bin/msolve-inproc` supplies timing. It links against `libmsolve` and
calls `core_msolve`. It repeats the call until a 0.5 s time budget is met,
with at least 3 repetitions and at most 5000.

## Rational protocol

The rational cells (docs/rational-design.md sections 1.2 and 11.4) compare
sylvester's multimodular engine against Singular, msolve, and Groebner.jl
over `Q`, on the same four families at smaller sizes: cyclic-4..6,
katsura-4..7, noon-3..5, eco-8..9. Coefficient growth, not monomial count,
drives cost over `Q`, so these sizes are the ones where the reference
tools return in seconds; `gen.py` builds only these instances into
rational input files, distinct from the prime-field ones of the same
family and size. M2 is not part of the rational cells: `gen.py` writes no
`.m2` file for a rational instance, because Macaulay2 is not one of the
three reference implementations `docs/rational-design.md` section 1.4 names
for the rational design.

There is no gate on the rational cells. What is checked is that every
tool that ran agrees on the full canonical rational basis: exact
fractions, monic, sorted grevlex descending, sorted basis, the rational
counterpart of the prime-field full-output comparison above. Timing is
recorded for the main table, not gated (docs/rational-design.md 1.2), because
the per-prime cost the record shows is not the pipeline the design
intends: the learn-and-apply cross-prime trace of section 3.9 is later
work.

**Input files.** Four formats per rational instance, all under `inputs/`,
named `<instance>-q.<ext>`: `.ms` (msolve, characteristic line `0`),
`.sing` (Singular, `ring r = 0,(...),dp;`), `.jl` (Groebner.jl over `QQ`),
and `.sylq` (sylvester). `.sylq` is a separate extension from the
prime-field `.syl` format rather than a characteristic line added to
`.syl` itself, because `.syl` carries no ring at all (the domain comes
from `sylv-runner`'s own mode argument) and every existing `.syl` file
and its parsing in `runner/src/main.rs` stay untouched. `.sylq` is
otherwise identical to `.syl`, line for line and term for term, with the
coefficient field of a term also accepting `numerator/denominator`. The
family generators in `gen.py` write integers only. `sylv-runner` reads
the wider grammar (`docs/rational-design.md` section 8.3).

**Reference tool invocations, confirmed against the installed versions.**
msolve 0.10.1 emits a lifted rational basis with `-g 2` when the
characteristic line of the `.ms` file is `0`: content-cleared (every
denominator cleared, then the integer coefficients divided by their gcd),
not normalized to monic. Singular 4.4.1's `option(redSB)` alone
content-clears over `Q` too, unlike over `F_p`, where it is already
monic because every nonzero residue there has an inverse Singular
applies; reaching a monic basis over `Q` needs
`simplify(std(I), 1)` after it, which `gen.py`'s rational `.sing` script
runs. Groebner.jl 0.10.3's `groebner(sys, ordering=DegRevLex())` over
`QQ` returns a monic basis with no extra step. `canon.py`'s `monic_q`
divides every basis by its own leading coefficient, so all four tools'
output is compared under the same normalization.

**sylvester's rational run.** `sylv-runner`'s `rational` mode calls
`Ideal::<Rationals>::groebner_basis_with_report` with
`RationalOptions::new()` (the default stopping rule,
`RationalStop::Unchanged { extra: 2 }`) and prints the basis with exact
`numerator/denominator` coefficients, already monic and already sorted.
It also prints a
`LIFT` line: the fields of `sylvester::ModularLift` (primes consumed,
skipped, folded, and discarded, the confirming-prime count, the modulus
bit length, and `established`), which `driver.py` stores under the
cell's `"lift"` key.

Threads are 1 everywhere in this protocol, sylvester included: msolve
`-t 1`, Julia `-t 1`, and sylvester `ComputeOptions::threads(1)`. For the
prime-field cells that setting keeps the F4 kernel's row-reduction pool
at one worker; for the rational cells it does that and also caps
`RationalOptions`' outer concurrency, so the multimodular driver runs one
prime at a time rather than starting several prime runs at once. There is
no rational counterpart of the `f4/threads8` config: measuring
concurrent prime runs is future work, not this record.

**msolve's rational timing** is a CLI wall clock around `msolve -g 2`.
The prime-field cells use `msolve-inproc`, which calls `core_msolve`
directly. There is no rational counterpart of `msolve_inproc.c`, so the
rational msolve column includes process startup.

**Smoke run.** `GBBENCH_INSTANCES` restricts the rational families the
same way it restricts the prime-field ones, by the `-q` names:

```
GBBENCH_INSTANCES=cyclic-4-q,katsura-4-q GBBENCH_RESULTS=/tmp/smoke.json \
    taskset -c 0-3,12-15 python3 driver.py
```

`report.py`'s `"== RATIONAL: WHICH TOOLS PRODUCED A CELL =="` section
prints which tools finished at least one rational cell, from
`results.json`'s `"_meta"` key (`"rational_tools"`), which `driver.py`
writes after the rational cells run. It reads the cells and not the
version strings: a tool on `PATH` that ran nothing here is not listed.

The rational correctness table compares full canonical bases, and a
comparison takes two of them. A cell with fewer than two parsed full
outputs reads `NO COMPARISON`, and the agreement line at the end counts
only the cells that were compared.

## Rebaseline support

The installed `core_msolve` print path omits cleanup because the command-line
program exits after printing. The runner forks one child per repetition. The
timed region remains the `core_msolve` call. Child exit reclaims its
allocations.

`build-runners.sh` builds this runner only when `pkg-config` finds msolve.
Rebuild it after an msolve update. The `core_msolve` signature changes across
releases.

Five rows per instance, from one runner binary. The runner takes one mode:
`f4`, `classic`, `certified`, `f4-certified`, or `rational`
(see "Rational protocol"), and a thread count.
`f4/default` and `classic/default` are the raw one-thread configurations,
and the gate in docs/f4-design.md section 14 applies to `f4/default`.
`f4/threads8` is `f4` on 8 threads. `classic/certified` runs
`groebner_basis_certified` on the classic backend, which writes a
`sylv-gb-cert-v1` certificate, and `f4/certified` runs the same call on the
F4 backend, which writes a `sylv-gb-cert-v2` certificate.

The gate uses one thread:

- msolve uses `-t 1`;
- Julia uses `-t 1`;
- sylvester uses `ComputeOptions::threads(1)`;
- Singular and the plain Macaulay2 call use their serial paths.

`f4/threads8` is a diagnostic cell outside the gate.

Run the driver under `taskset -c 0-3,12-15`. Every child inherits that
affinity. The stored `cores` field must match the command.

## Sylvester modes

One binary, `runner/bin/sylv-runner`, provides five configurations:

| Result key | Runner mode | Threads | Certificate |
| --- | --- | ---: | --- |
| `sylvester-f4/default` | `f4` | 1 | none |
| `sylvester-f4/threads8` | `f4` | 8 | none |
| `sylvester-classic/default` | `classic` | 1 | none |
| `sylvester-classic/certified` | `certified` | 1 | v1 |
| `sylvester-f4/certified` | `f4-certified` | 1 | v2 |

A certified time covers computation, writing, and the bundled verification.
The runner also times one later call to `verify::verify` on the stored bytes.
The report labels this value `verify`. `split` is certified total minus raw
time. It bounds writing plus bundled verification. It does not measure either
part alone.

## Counters

Each result has a `counters` object. Sylvester F4 fills the fields of
`F4Counters` and `ComputeReport::threads_used`. Classic and external tools leave
engine counters null.

`peak_rss_kb` comes from `/proc/self/status` for sylvester. The msolve runner
uses child `getrusage` maxima. These values are diagnostics, not the engine's
tracked memory estimate.

## Containment

Each child runs in a transient user systemd scope:

```text
MemoryMax=16G
MemorySwapMax=0
```

This prevents one cell from exhausting the host. A cgroup kill records `OOM`.
A timeout records `DNF`.

## Release gate

The gate uses `f4/default` against msolve on `cyclic-7`, `katsura-9`, `eco-9`,
and `noon-6`. It requires:

\[
\operatorname{gmean}(t_{\mathrm{F4}}/t_{\mathrm{msolve}})\le5,
\qquad
\max(t_{\mathrm{F4}}/t_{\mathrm{msolve}})\le10.
\]

See section 14 of `docs/f4-design.md`. `report.py` computes the ratios from
the recorded cells.

## Files

- `gen.py`: deterministic input generator.
- `driver.py`: resumable runner.
- `canon.py`: full-basis parser and canonicalizer.
- `report.py`: CSV and Markdown table generator.
- `build-runners.sh`: sylvester and msolve runner build.
- `runner/`: runner sources and independent Cargo lock file.
- `inputs/`: generated inputs in each tool grammar.
- `records/`: frozen prior reports and raw data.
- `outputs/`: ignored transient tool output.
