# Changelog

This file follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions follow [semantic versioning](https://semver.org/).

## 0.3.0 (2026-09-04)

A coefficient parameter, a multimodular rational engine, normal form and
membership, Hilbert series and dimension, a command line interface, and
Python bindings. Every item under "Changed" is a breaking change. There
is no compatibility alias.

### Added

- A sealed coefficient domain, `Domain`, with two implementations:
  `PrimeField` (the default type parameter of `PolynomialRing`,
  `Polynomial`, `Ideal`, and `GroebnerBasis`) and `Rationals`. `DomainOps`
  is the arithmetic each domain supplies; `PrimeOps` and `RationalOps`
  are its two implementations, both public with private fields.
- `PolynomialRing::rationals`, building a ring over `Q`. `Coefficient`
  and its `From` impls for every integer width, `BigInt`, and
  `BigRational` let one call build a polynomial in either domain.
- A `/` production in the text grammar: `1/2*x` parses in both domains,
  reduced modulo the prime over `F_p` and kept exact over `Q`.
  `ParseError::ZeroDenominator` and `RingError::{ZeroDenominator,
  CoefficientNotInvertible}` report a bad or unreadable fraction.
- `Felt`, public again, with `Felt::value` returning the residue as a
  `u64`.
- A multimodular rational engine: `Ideal::<Rationals>::groebner_basis`
  and `groebner_basis_with_report` clear denominators, run the chosen
  backend modulo a descending sequence of 31-bit primes, combine the
  runs by the Chinese remainder theorem, and lift the result by rational
  reconstruction. `RationalOptions` and `RationalStop::{Unchanged,
  ContainsInput}` control it. `GroebnerBasis::<Rationals>::lift` returns
  the `ModularLift` counters and `Established`, which names what the
  stopping rule observed: `Unchanged` establishes nothing about the
  ideal, `ContainsInput` establishes that the basis is the reduced
  Gröbner basis of an ideal containing the input, and neither
  establishes equality. `ComputeError::PrimesExhausted` reports a run
  that walked the whole prime sequence without stopping.
  `ComputeReport` gains `modular: Option<ModularLift>` and
  `modular_concurrency: Option<usize>`, filled on the rational path and
  `None` on the prime-field path.
- Normal form and ideal membership in both domains:
  `GroebnerBasis::normal_form` and `GroebnerBasis::contains`, under a new
  `Budget` type (`timeout`, `memory_limit`) and reporting
  `NormalFormError`. `ComputeOptions::budget` takes a `Budget` directly.
- A checked basis constructor, `GroebnerBasis::from_polynomials`, one
  method per domain, which accepts a caller-supplied list only after
  checking it is monic, sorted, interreduced, and a Gröbner basis of its
  own ideal, reporting `BasisError` for the first check that fails. This
  is a check, not a certificate.
- The Hilbert series of the quotient by the leading monomial ideal:
  `HilbertSeries`, `HilbertError`, `GroebnerBasis::hilbert_series`,
  `GroebnerBasis::krull_dimension` (`None` for the unit ideal), and
  `GroebnerBasis::is_homogeneous`. `Ideal::has_homogeneous_generators`
  tests the generators, which is sufficient and not necessary. The
  series, the dimension, and `contains` are about the ideal the basis
  generates; over `Q` what relates that ideal to the input is what
  `Established` says, and each docstring states it.
- The `sylvester-cli` crate, binary `sylv`: seven subcommands (`gb`,
  aliased `basis`; `certify`; `normal-form`, aliased `nf`; `member`;
  `hilbert`; `dim`; `verify`) over three input formats (`ms`, `syl`, and
  the crate's own `text`, plus `.sylq` for a rational `syl` instance),
  with `--certified` refused over `Q` and six exit codes for success,
  a rejected certificate, usage, a resource or structural limit, an
  internal defect, and an I/O failure. `--report` writes counters to
  standard error so the basis on standard output stays parseable.
  `hilbert` and `dim` over `Q` write an `# established:` comment above
  their output, because the value is of the ideal the lifted basis
  generates.
- The `sylvester-py` crate: Python bindings over pyo3 `=0.23.5`
  (abi3-py310), exposing the ring, the polynomial, the ideal, the basis,
  Hilbert series, certification, and `sylvester.verify`. Every exception
  subclasses `sylvester.SylvesterError` and a matching standard
  exception (`ValueError` for bad input, `RuntimeError` for an exhausted
  budget or a defect).
- `sylvester-cli` and `sylvester-py` as new workspace members. The CLI is
  publishable. The Python crate remains private and produces release
  artifacts.
- New dependencies of the library: `num-bigint`, `num-integer`,
  `num-rational`, `num-traits`. New dependencies of the new crates:
  `clap` (`sylvester-cli`) and `pyo3` (`sylvester-py`). The verifier
  gains none.
- CI checks Rust 1.88 and 1.92, the three supported operating systems,
  the slow oracle sweeps, the extracted source crate, and installed
  Python wheels and source distributions. `CITATION.cff` names the
  software citation.
- The CLI and Python binding default to the `ContainsInput` rational stop.
  Rust keeps `RationalStop::Unchanged` as its explicit low-cost default.

### Changed

- README files no longer name a release version. The changelog owns release
  numbers, while manifests and benchmark provenance retain machine versions.
- The F4 and rational design records use subject names instead of release
  numbers. Frozen benchmark directories use dates and subjects.
- The minimum Rust version is 1.88. The engine uses let chains, which
  stabilized in that release.
- `PolynomialRing`, `Polynomial`, `Ideal`, and `GroebnerBasis` take a
  sealed domain parameter `D: Domain`, defaulting to `PrimeField`, so a
  prime-field signature names no domain.
- `Polynomial::terms` yields `(&D::Coeff, &[u16])` instead of
  `(u64, &[u16])`. Over `PrimeField`, `Felt::value()` recovers the `u64`.
- `PolynomialRing::polynomial` takes `impl Into<Coefficient>` in place of
  `i64`, so every integer width still converts with no annotation.
- `ComputeOptions::timeout` and `memory_limit` now write into a `Budget`
  value the options hold; the setters are unchanged at the call site.
- `Ideal::groebner_basis_certified` is defined on `Ideal<PrimeField>`
  alone, and `CertifiedGroebnerBasis::basis` returns
  `&GroebnerBasis<PrimeField>`: certification over `Q` is a compile
  error, not a runtime one, because no isolated verifier checks the
  vote, the Chinese remainder step, or the rational reconstruction.
- The crate description reads "Gröbner bases over prime fields and the
  rationals: an F4 engine, a multimodular rational driver, a classic F5
  oracle, and independent certificate verifiers."

### Fixed

- The rational differential property now samples monomials of total
  degree at most 2. Its old per-variable bound admitted total degree 6,
  outside the scope the test states.
- The source crate excludes the benchmark records and carries its own
  `eco-9` stack-test fixture. Tests from the extracted crate no longer
  depend on the repository tree.
- With `GBBENCH_RESULTS` set, `report.py` writes its CSV next to that JSON
  file. A benchmark smoke run no longer overwrites the recorded CSV.

### Removed

Nothing beyond the signature changes above. `Backend`, `F4Counters`,
`verify::verify`, and both certificate contracts are untouched.

## 0.2.0 (2026-08-27)

New F4 engine, made the default backend, with a new binary certificate
format for it. The `sylv-gb-cert-v1` classic path stays.

### Added

- `Backend::F4`, now the default. It batches critical pairs by degree and
  reduces sparse matrices over interned monomials.
- `sylv-gb-cert-v2`, a binary trace certificate for F4. The frozen format is
  in `docs/certificate-v2.md`.
- An independent v2 verifier under `src/verify/v2`. It shares only public
  result and limit types with the v1 verifier.
- Certificate dispatch in `verify::verify` and `verify::verify_with_limits`.
  A leading `{` selects v1 JSON. The `SYLVGB` magic selects v2 binary.
- `ComputeOptions::threads`. Zero uses Rayon's global pool. One uses the
  calling thread. Certified runs use one thread.
- `Ideal::groebner_basis_with_report`, `ComputeReport`, and `F4Counters`.
- `ComputeError::ExponentLimit` for an F4 exponent above 65,535.
- `ComputeError::TableFull` for an exhausted F4 monomial table.
- `CertifyError::CapExceeded`, `CertificateCap`, and `TraceFault` for v2
  writer limits and defects.
- CI gates for the MSRV, extracted package, slow differential sweeps,
  cyclomatic complexity, dependency policy, rustdoc, Clippy, and rustfmt.

### Changed

- F4 replaces the signature-based matrix backend. Classic F5 remains an
  independent oracle and the v1 certificate backend.
- `Ideal::groebner_basis_certified` follows the selected backend. Classic
  writes v1. F4 writes v2.
- Rayon is an unconditional dependency. F4 uses parallel row reduction when
  a batch passes its internal work threshold.
- The memory budget charges F4 tables, matrices, worker accumulators, trace
  nodes, v2 input storage, and division traces before allocation.
- Deadline checks cover pair management, symbolic preprocessing, matrix
  construction, elimination, certificate writing, and verification.
- A passed verification deadline outranks acceptance and invalidity. A cap
  already reported remains a cap.
- Internal Rust functions now have cyclomatic complexity at most 10 under
  Lizard 1.24.0. Clippy uses a cognitive-complexity threshold of 10.
- Private F4 inspection methods compile only for tests. Unused private
  interfaces are removed.
- The one-thread F4 release gate passes against msolve at a 1.13x geometric
  mean and a 1.48x worst cell. All compared complete bases agree.

### Removed

- `Backend::Matrix` and `src/compute/matrix.rs`.
- The `parallel` Cargo feature. The crate has no feature flags.
- `CertifyError::NotCertifiable`, which named no reachable state.

## 0.1.0 (2026-08-21)

The first public alpha replaced the incorrect extraction preserved in
`main-local` history.

### Added

- `PolynomialRing` as the construction root for polynomials and ideals over
  `F_p`, with `2 <= p <= 2^31 - 1` and at most 256 variables.
- A sealed grevlex order over the ring's variable order.
- `Ideal::groebner_basis` with the repaired classic and matrix backends.
- `Ideal::groebner_basis_certified` for the classic backend.
- The canonical JSON contract `sylv-gb-cert-v1` and its independent verifier.
- Typed ring, parse, compute, certificate, and verification errors.
- Deadlines and memory budgets for computation, certificate writing, and
  verification.

### Fixed

- Non-regular critical pairs no longer record inflated syzygy signatures.
- Sig-redundant insertion stops duplicate leading-term and signature loops.
- Rewriter substitution replaces unsound pair deletion in the old matrix
  backend.
- Final interreduction runs under the computation budget.
- A critical pair past the exponent width returns `DegreeLimit`.
- `ComputeOptions::timeout(Duration::MAX)` no longer overflows an `Instant`.
- Classic runs `eco-9` on a 1 MiB test stack without stack exhaustion.

`KNOWN_ISSUES.md` gives the counterexamples, root causes, proof argument, and
evidence limits.
