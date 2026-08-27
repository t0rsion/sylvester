# Gröbner basis comparison

This harness compares sylvester with Singular, msolve, Macaulay2, and
Groebner.jl. It records timings and complete reduced bases over
\(\mathbb F_{1073741827}\) under grevlex.

The families are:

- `cyclic-4` through `cyclic-8`;
- `katsura-4` through `katsura-10`;
- `eco-8` through `eco-11`;
- `noon-3` through `noon-6`.

`REPORT.md` is the current record. `results.csv` contains its scalar cells.
`results.json.zst` contains commands, counters, and complete bases. Frozen
records stay under `records/`.

The current record remeasures every sylvester cell and the four msolve gate
cells. It reuses the remaining external cells from the frozen 2026-08-19 record.
`results.json.zst` records this provenance in `_meta`.

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

Required commands are `Singular`, `msolve`, `M2`, and `julia`. Groebner.jl must
exist in the Julia environment. A missing tool fails only its own cells.

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
calls `core_msolve`. It performs at least 3 and at most 5,000 repetitions until
the timing window is full.

The installed `core_msolve` print path omits cleanup because the command-line
program exits after printing. The runner forks one child per repetition. The
timed region remains the `core_msolve` call. Child exit reclaims its
allocations.

`build-runners.sh` builds this runner only when `pkg-config` finds msolve.
Rebuild it after an msolve update. The `core_msolve` signature changes across
releases.

## Threads

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

See section 14 of `docs/v0.2-design.md`. `report.py` computes the ratios from
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
