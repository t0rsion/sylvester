# Rational equality exploratory probe, 2026-09-18

This table is generated from the captured command output at measured source commit `68f75768b9acf255d0367f80fdad177fd3ed13d8`.
The checker uses `Polynomial<Rationals>`, `Budget`, and grevlex-v1.

| input | status | engine | exact check | process | raw basis elements | raw basis terms | final working bytes |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cyclic-5 | success | 1.618528ms | 15.667116835s | 15.734s | 46 | 670 | 7271232 |
| katsura-5 | success | 7.442002ms | 12.376897586s | 12.420s | 24 | 549 | 1320832 |
| noon-4 | success | 2.780568ms | 374.821413ms | 0.414s | 28 | 444 | 269040 |
| eco-8 | timeout | 32.871381ms | 300.002128755s | 300.071s |  |  |  |

The exact check elapsed time comes from the test output. The process time includes Cargo and test harness overhead.
The successful rows completed exact equality checks. The eco-8 row reached its typed timeout.
The four rows are single exploratory samples and are not a performance gate.
The lock serialized the probe with the release benchmark, but did not exclude unrelated host work.
The `final_working_bytes` value is a point estimate immediately before the temporary exact basis is dropped. It is not a running peak tracker or an RSS measurement.

The compressed JSON contains the portable commands, source and input hashes, limits, toolchain, CPU affinity, and sanitized captured output. `results.csv` carries the same parsed measurements.
