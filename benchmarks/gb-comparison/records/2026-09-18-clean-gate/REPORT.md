# Performance gate repeat, 2026-09-18

Run 2026-09-18 after an unrelated CPU-bound process was stopped.
This record supplies the release gate timing. The top-level record retains the full refresh and historical reference cells.

## Provenance

- Measured source commit: `68f75768b9acf255d0367f80fdad177fd3ed13d8`.
- Runner SHA-256: `1f514eb6a84f6e6d4f566c55ca269bbccc61f667776f783da5b1ab610208d763`.
- P-cores: `0-3,12-15`. Memory cap: `16G`.
- Started: `2026-09-18T00:45:36.849324+00:00`. Finished: `2026-09-18T00:45:42.292664+00:00`.
- Reason: Repeat performance gate after stopping an unrelated CPU-bound process.

## Gate

| cell | sylvester f4/default (s) | msolve (s) | ratio |
|---|---:|---:|---:|
| cyclic-7 | 0.076606215 | 0.078391977 | 0.98x |
| katsura-9 | 0.252875262 | 0.229818113 | 1.10x |
| eco-9 | 0.015791203 | 0.013231605 | 1.19x |
| noon-6 | 0.011724404 | 0.011754756 | 1.00x |

Geometric mean: 1.06x (limit: 5.0x).
Worst cell: 1.19x at eco-9 (limit: 10.0x).
Gate: PASS.

All eight cells returned OK. Their complete canonical bases agree exactly.
The four sylvester and four msolve cells were run under the shared P-core
lock with the requested affinity and memory cap.
