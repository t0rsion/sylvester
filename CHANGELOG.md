# Changelog

This file follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Versions follow [semantic versioning](https://semver.org/).

## 0.2.0 (2026-08-27)

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
