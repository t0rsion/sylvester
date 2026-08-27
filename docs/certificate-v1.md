# Certificate contract: sylv-gb-cert-v1

Status: frozen.

This contract belongs to the classic backend. Change it only by defining a new
schema. F4 uses `sylv-gb-cert-v2`, specified in `certificate-v2.md`.

## 1. Trust model

The engine and emitter are untrusted candidate producers. The verifier in
`src/verify` is the trust boundary. It accepts bytes only after it checks every
obligation below.

The verifier implements its own decoder, field arithmetic, monomials,
polynomials, grevlex order, and S-pair enumeration. It imports no engine,
emitter, ring, or public polynomial implementation.

Acceptance proves a result. It does not prove the engine correct. An engine may
still stop, exhaust a limit, or emit a rejected certificate.

## 2. Statement proved

Let \(F=(f_0,\ldots,f_{m-1})\) be `input`. Let
\(G=(g_0,\ldots,g_{k-1})\) be `basis`. Acceptance proves, over
\(\mathbb F_p\) under `grevlex-v1`:

1. \(g_j\in\langle F\rangle\) for every \(j\), by `origin`.
2. \(f_i\in\langle G\rangle\) for every \(i\), by `membership`.
3. Every non-coprime S-pair of \(G\) has a standard representation over
   \(G\), by `spairs`.
4. \(G\) is monic, minimal, interreduced, and canonically sorted.

Items 1 and 2 give \(\langle F\rangle=\langle G\rangle\). Item 3 and the
product criterion make \(G\) a Gröbner basis. Item 4 makes it the unique
reduced Gröbner basis of the input ideal.

## 3. Canonical JSON

The certificate is one UTF-8 JSON object.

- It contains no insignificant whitespace.
- Keys occur once, in the order listed below.
- Unknown, missing, duplicate, or reordered keys are invalid.
- Integers use unsigned decimal notation, with no leading zero except `0`.
- Floats, exponents, signs, escapes, and non-printable string bytes are invalid.
- Coefficients lie in \([1,p-1]\). A zero coefficient is invalid.
- Exponents lie in \([0,65535]\).

The emitter is deterministic. Equal inputs and options produce equal bytes.

## 4. Schema

The top-level keys occur in this order:

1. `schema`: the string `"sylv-gb-cert-v1"`.
2. `order`: the string `"grevlex-v1"`.
3. `modulus`: a prime \(p\), where \(2\le p\le 2^{31}-1\).
4. `nvars`: an integer \(n\), where \(0\le n\le256\).
5. `input`: the canonical input list \(F\). A zero polynomial remains `[]`.
6. `basis`: the claimed reduced basis \(G\).
7. `origin`: one cofactor list per basis element.
8. `membership`: one cofactor list per input polynomial.
9. `spairs`: the non-coprime S-pair witnesses.

A term is

```text
[coefficient,[e_0,...,e_(n-1)]]
```

A polynomial is an array of terms in strict descending grevlex order. The zero
polynomial is `[]`. Duplicate monomials are invalid.

Grevlex compares total degree first. At equal degree, scan exponents from the
last variable to the first. At the last difference, the smaller exponent is
the larger monomial.

`basis` is strictly descending by leading monomial. `[]` claims the zero ideal.
The one-element basis containing the constant polynomial 1 claims the unit
ideal.

## 5. Witnesses

### 5.1 Origin

`origin[j]` contains \(m\) polynomials \(c_{ji}\) and asserts

\[
g_j=\sum_{i=0}^{m-1}c_{ji}f_i.
\]

### 5.2 Membership

`membership[i]` contains \(k\) polynomials \(q_{ij}\) and asserts

\[
f_i=\sum_{j=0}^{k-1}q_{ij}g_j.
\]

### 5.3 S-pairs

One entry has the form `[i,j,h]`, where \(0\le i<j<k\) and `h` contains
\(k\) polynomials. It asserts

\[
S(g_i,g_j)=\sum_{l=0}^{k-1}h_lg_l.
\]

For every nonzero \(h_l\),

\[
\operatorname{lm}(h_lg_l)
\le \operatorname{lm}(S(g_i,g_j)).
\]

Entries use lexicographic `(i,j)` order and contain no duplicate pair. Every
non-coprime pair must occur. A coprime pair may be omitted. The verifier
enumerates all pairs itself.

## 6. Verification

The verifier performs these checks in order:

1. Decode the canonical JSON under the byte and count caps.
2. Check the schema, order, prime modulus, variable count, exponent width,
   coefficient range, and polynomial order.
3. Check that the basis is sorted, monic, minimal, and interreduced.
4. Check every origin identity.
5. Check every membership identity.
6. Enumerate all S-pairs and check coverage, identities, and leading bounds.

Every failure is a typed `VerifyError`. Exhaustion does not assert invalidity.
If a deadline passes, `DeadlineExceeded` outranks acceptance and invalidity. A
cap already reported remains `CapExceeded`.

## 7. Limits

The v1 verifier uses these fields of `verify::Limits`:

| Field | Meaning | Default |
| --- | --- | ---: |
| `max_bytes` | Certificate bytes | `64 << 20` |
| `max_polys` | Input, basis, and cofactor polynomials | 4,000,000 |
| `max_terms_per_poly` | Terms in one polynomial | `1 << 20` |
| `max_total_terms` | Terms in the certificate | 8,000,000 |
| `max_entries` | `origin`, `membership`, and `spairs` entries | 4,000,000 |
| `max_intermediate_bytes` | Live arithmetic buffers | `512 << 20` |
| `deadline` | Absolute stop instant | none |

One arithmetic term costs

```text
size_of::<Term>() + nvars * size_of::<Exp>()
```

The verifier charges buffer capacities before allocation. The byte estimate
excludes allocator metadata and is not a process RSS cap.

The decoder and arithmetic poll the deadline every 1,024 work steps. A caller
sets a deadline for untrusted bytes. Basis and pair checks are quadratic in the
basis count.

The remaining `Limits` fields belong to v2.

## 8. API

`Ideal::groebner_basis_certified(options)` follows `options.backend`:

- classic computes and writes v1;
- F4 computes and writes v2.

The method verifies its own bytes before returning. Its basis is reconstructed
from the accepted certificate. Emitter defects return `CertifyError::Emitter`.
Verifier rejection returns `CertifyError::Rejected`.

`verify::verify` dispatches by the first byte. `{` selects v1. The `SYLVGB`
magic selects v2. `verify::verify_with_limits` applies caller-supplied limits.

`Ideal::groebner_basis` returns the selected engine's result without independent
verification.

## 9. Emission

Classic tracks origin cofactors through computation and final interreduction.
The emitter derives membership and S-pair cofactors by recorded division over
the final basis.

The emitter charges input, basis, origin, generated cofactors, and output bytes
to the computation memory limit. It uses the same deadline as the engine and
the final verifier.

Emission never repairs a candidate. A wrong identity reaches the verifier and
is rejected.
