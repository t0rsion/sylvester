# sylv-gb-cert-v2 test vectors

The five certificates of `docs/certificate-v2.md` section 11, byte for
byte. `tiny.cert` is the listing of section 10. All five are over
`F_7[x, y]` with the variables in the order `x, y`.

| File | Input | Basis | Bytes |
|---|---|---|---|
| `tiny.cert` | `[x^2 + y, x]` | `[x, y]` | 93 |
| `square.cert` | `[x^2*y, x*y^2]` | `[x^2*y, x*y^2]` | 80 |
| `unit.cert` | `[3]` | `[1]` | 62 |
| `zero.cert` | `[0]` | `[]` | 49 |
| `empty.cert` | `[]` | `[]` | 47 |

The verifier accepts all five. The writer reproduces `square.cert`,
`unit.cert`, `zero.cert`, and `empty.cert` byte for byte, because section
9 fixes every node of those four traces. The writer's output for `tiny` is
not pinned to the listing: the listing is a hand-built derivation, and the
trace the writer records depends on how F4 batches the run. Section 11
states what must match instead.
