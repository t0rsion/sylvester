# Rational coefficients, CLI, and Python design

## Release audit (2026-09-04)

The local release audit passes the default workspace suite, the ignored
oracle sweeps, Rust 1.88, the extracted source crate, and installed Python
wheels and source distributions. The complexity, dependency, workflow, and
privacy checks also pass. The manifests had named Rust 1.85 even though the
code needs 1.88.
The rational property had bounded each exponent by 2 instead of bounding
the monomial's total degree by 2. The source crate now carries the `eco-9`
fixture that its stack test reads.

## Status after implementation (2026-08-19)

The design is implemented in the chunks section 12 lays out. Each
item below names what a chunk shipped as designed and the deviations its
implementer recorded.

**A1, the sealed domain over the prime field.** Shipped as designed:
`PolynomialRing<D>`, `Polynomial<D>`, `Ideal<D>`, and `GroebnerBasis<D>`
all default to `D = PrimeField`, and `PolynomialRing` gained a content
`Hash` next to its content `PartialEq`. `RingError::ZeroDenominator` and
`RingError::CoefficientNotInvertible` landed in this chunk, ahead of
`Rationals` itself, because `Coefficient::normalized` needs both before
any domain reads a fraction. Two deviations: `DomainOps`, not only
`Domain`, is sealed, so a caller cannot write a third implementation of
either trait even indirectly; and the debug assertion in
`Polynomial::from_sorted_terms` is narrower than a full invariant check,
verifying the variable count and the strictly ascending order of the
terms and nothing about zero coefficients, which the domain's own
zero-is-`None` convention already rules out by construction.

**A2, the rational domain.** Shipped as designed: the parser is generic
over `D` and reads the `a/b` production of section 2.4 in both domains,
and `ParseError::ZeroDenominator` carries a `position` field like every
other `ParseError` variant. One deviation from the module list of
section 10: `Rationals` is declared in `src/ring.rs`, next to
`PrimeField`, rather than in `src/ring/rational.rs`, which holds
`RationalOps` and the arithmetic alone. Prime-field parsing also routes
every coefficient through `BigInt` first, the conversion the rational
domain needs anyway, rather than keeping a separate machine-integer path.

**A3, limits, budget, and the meters.** Shipped as designed: `ComputeLimits`
and `groebner_basis_with_limits` stay crate-private, as the design's own
body already says in section 3.8, even though the section 12 table lists
them under this chunk's surface; `Budget` is the public type, and
`ComputeOptions::budget` takes it directly. One deviation: cancellation
gains no representation in any public error type, not even indirectly. A
cancelled run reports `ComputeError::Timeout`, the same value an
exhausted deadline reports, since both say the same thing to a caller:
the run stopped before it produced a basis.

**B, normal form and the checked basis.** Shipped as designed:
`GroebnerBasis::normal_form` and `GroebnerBasis::contains` both take a
`Budget`; `GroebnerBasis::from_polynomials` is one method per domain,
reached by turbofish (`GroebnerBasis::<PrimeField>::from_polynomials`,
`GroebnerBasis::<Rationals>::from_polynomials`); `NormalFormError` and
`BasisError` carry the variants section 5 lists. `Monomial::checked_mul`
and `Polynomial::sub_scaled_checked` ship as a path separate from the
engines' own unchecked multiply, so the F4 and classic hot loops keep
their `expect`. One deviation: the Gröbner check inside
`from_polynomials` divides every pair's S-polynomial, with no skip for a
pair whose leading monomials are coprime. The engines' own Buchberger
product criterion, which does skip such pairs, was removed from the
engines during the initial repair because its interaction with signatures
had no soundness argument (`KNOWN_ISSUES.md`); reusing it here, on a
one-shot check rather than a hot loop, was judged not worth reopening
that question for.

**C, the multimodular engine.** Shipped as designed: `RationalOptions` and
`RationalStop` are public; `Ideal::<Rationals>::groebner_basis` and
`groebner_basis_with_report` exist; `ComputeError::PrimesExhausted` and
`ComputeReport::{modular, modular_concurrency}` ship as specified. Three
deviations: the CRT accumulator of section 3.4 is keyed by a `BTreeMap`,
not a hash map, so folding a run walks its monomials in a fixed order
with no separate sort before a candidate is built; a discarded run (one
folded into a category that later loses prevalence) keeps the sequence
position it was consumed at instead of being renumbered, so the counters
of section 4.2 read directly off the retained list; and the driver's
concurrency is `std::thread::scope`, not a rayon pool, because each prime
run already forces its own inner thread count to 1, leaving the driver to
join a bounded, already-spawned set of threads rather than schedule a
task graph. Decision 9 of section 13 is resolved: the installed msolve
0.10.1 does emit a lifted rational basis with `-g 2` over a `0`
characteristic line, so msolve stayed in the rational reference set. One
limit is narrower than section 3.8: the driver starts at most 256 prime
runs at once, whatever `ComputeOptions::threads` says, and
`ComputeReport::modular_concurrency` reports the largest number it
started rather than the number the caller asked for.

**D, Hilbert.** Shipped as designed: `HilbertSeries`, `HilbertError`,
`GroebnerBasis::{hilbert_series, krull_dimension, is_homogeneous}`, and
`Ideal::has_homogeneous_generators` all carry the signatures of section
6.1. Two deviations: there is no method returning the affine Hilbert
function `H(d)` itself, only the series coefficients that are its first
difference, since section 6.1 already states the two are not equal and a
second method inviting the reader to conflate them was judged to work
against that statement rather than for it. The Singular reference values
the tests compare against are committed as source constants, not as a
separate fixture file, because the set is small and fixed. The recursion
of section 6.2 runs on a heap allocated worklist, not on the calling
stack: the colon branch reaches one node per unit of the degree sum, and
the basis `[x^65535*y, z]` reaches 65,535 of them, which the calling
stack does not hold. The budget charges the worklist and every dense
coefficient vector before it is allocated.

**E, the CLI.** Shipped as designed: seven subcommands (`gb`/`basis`,
`certify`, `normal-form`/`nf`, `member`, `hilbert`, `dim`, `verify`),
`.sylq` as a `syl`-format extension for a rational instance, and the six
exit codes of section 8.4. `--certificate` implies `--certified`, and
`--certified` on a rational ring is a usage error (exit 2), reported
before any engine runs, rather than attempted and failed. `--report`
together with `--certified` is also a usage error, since the certified
path's single timed call has no separate counters to report. `--report`
alone writes counters to standard error, so a basis on standard output
stays parseable by the next tool in a pipeline. One deviation from the
output schema of section 8.2: over `Q`, `hilbert` and `dim` write an
`# established:` comment and one line of what it means above the value,
because the value is of the ideal the lifted basis generates and the
schema alone would read as a fact about the input ideal. A basis the
command reads with `--from-basis` carries no lift, so it carries no
comment either.

**F, Python.** Shipped as designed: every exception class carries two
bases, `sylvester.SylvesterError` and the matching standard exception
(`ValueError` for bad input, `RuntimeError` for an exhausted budget or a
defect), exactly as section 9.3 specifies. pyo3 is pinned to `=0.23.5` as
decided, built and tested against Python 3.14 through the `abi3-py310`
wheel, so the pin does not narrow which interpreter runs the tests. One
deviation: `num-bigint` and `num-rational` are direct dependencies of
`sylvester-py` itself, not only pyo3 features, because the binding
converts a `BigRational` coefficient by hand at the Python boundary
rather than relying on a pyo3 conversion alone. `__contains__` is not
implemented, as designed, so the sequence protocol's `in` iterates and
compares polynomials by value, the same fallback Python gives any
sequence that defines no `__contains__` of its own.

**G, the harness.** Shipped as designed: rational input files use the
`.sylq` extension and an instance name suffixed `-q` (`cyclic-4-q`, and
so on), which `GBBENCH_INSTANCES` matches the same way it matches a
prime-field instance name. msolve's content-cleared normalization
(denominators cleared, then the integer coefficients divided by their
gcd, not forced monic) is documented rather than assumed; `canon.py`'s
`monic_q` step is what makes every tool's output comparable regardless of
which normalization each one starts from.

Working spec for the implementation. It fixes the coefficient model,
the multimodular rational engine, the normal form and Hilbert series
surface, the `sylv` binary, the Python package, the test plan, and the
work split. Implementation agents follow this spec. A deviation needs a
recorded reason in the agent's final report.

The F4 engine of `docs/f4-design.md` does not change in its algebra or
its matrix kernel, and neither does what it computes over a prime field:
reduced Gröbner bases under grevlex, named `grevlex-v1`. What changes is
what can be given to it, what can be asked of the result, and its control
plane (section 2.2). Backward compatibility is not a goal.

Over a prime field this crate returns either an engine claim or a value
an isolated verifier accepted. Over the rationals it returns the result
of a heuristic modular computation, which no isolated verifier has
checked, and section 4 says why. Every rational result is weaker than
every prime-field certified result, and the types say so.

## 1. Goals, scope, and the gate

### 1.1 What the design adds

1. A coefficient parameter. `PolynomialRing<D>` carries a sealed domain
   `D`, either `PrimeField` or `Rationals`, and `Polynomial<D>`,
   `Ideal<D>`, and `GroebnerBasis<D>` follow it.
2. A multimodular rational engine: clear denominators, compute modulo
   many 31-bit primes with the existing F4 engine, combine by the Chinese
   remainder theorem, and lift by rational reconstruction.
3. Normal form and ideal membership against a computed basis, in both
   domains, and a checked public constructor for a basis the caller
   brings.
4. The Hilbert series of the quotient by the leading monomial ideal, and
   the Krull dimension read off it.
5. `sylv`, a command line binary in a new crate, over three input formats
   and the certificate files.
6. `sylvester`, a Python package in a new crate, built with pyo3 and
   maturin.

### 1.2 The gate

The design has a correctness gate and no speed gate.

- Every rational cell of the comparison harness agrees with the reference
  tools on the full canonical basis, coefficients included, not on basis
  size or leading monomials. Section 11.4 makes a working rational
  invocation for each tool a precondition of the gate rather than an
  assumption.
- The differential suite runs the rational path against an exact
  Buchberger oracle over `BigRational` on small systems, and against a
  fresh prime-field run modulo a prime the lift did not use on large
  ones.
- `cargo test` stays green, including `tests/known_defects.rs`, whose
  assertions are never edited.
- The prime-field path is unchanged in output: the committed golden
  existing certificates still verify byte for byte after the coefficient
  parameter lands.
- The Python smoke suite and the CLI golden tests pass on the release
  machine.

Rational timings are measured and recorded in
`benchmarks/gb-comparison/REPORT.md` under the same protocol as the
prime-field cells. They are reported, not gated. The ratios measure the
pipeline this release runs. They are not a test of the learn-and-apply
trace of section 3.9, which is later work and changes the per-prime cost.
A gate against a pipeline missing that optimization either rewards a slow
design or holds the release for out-of-scope work.

### 1.3 Non-goals

Not in this design: general monomial orders, elimination orders, FGLM, the
radical, saturation, solving, coefficient fields other than `F_p` and
`Q`, number fields, and modules. Certification over `Q` is a non-goal,
and section 4 says what that costs and why the near alternative is not
certification. The learn-and-apply modular trace is later work.

### 1.4 Reference implementations

The rational design borrows from three systems. Every statement below is
about one version, and the harness record carries the version, the
invocation, and the source location for each claim it relies on
(section 11.4). This document does not restate a mechanism it has not
read.

| system | what the design takes |
|---|---|
| Singular `modStd` (`modstd.lib`), 4.4.1 | the prevalence vote over leading monomial sets (section 3.3), Farey rational reconstruction, and the shape of the two exact tests over `Q` (section 3.7) |
| msolve 0.10.1 | the pipeline shape: clear denominators, run many primes, CRT, reconstruct. Its learn-and-apply trace is deferred (section 3.9) |
| Groebner.jl 0.10.3 | the shape of a confirmation step against a prime the lift did not use, in a deterministic form (section 3.7) |

For general inhomogeneous input, none of the three proves its rational
output in its default mode. Singular's `modStd` and Groebner.jl's
`certify=true` document a guarantee for homogeneous input under their own
hypotheses. This design implements neither homogeneous verification nor any
other proof over `Q`. Section 4 states the crate's own position.

## 2. The coefficient parameter

### 2.1 A sealed type parameter, not a runtime tag

The coefficient parameter is a type parameter, sealed to two domains
inside the crate.

```rust
// src/ring.rs
pub trait Domain: sealed::Sealed + Clone + Debug + Eq + Hash + Send + Sync + 'static {
    type Coeff: Clone + Debug + Eq + Hash + Send + Sync;
    type Ops: DomainOps<Coeff = Self::Coeff> + Clone + Debug + Eq + Hash + Send + Sync;
    /// What a basis of this domain records about its own origin.
    type BasisMeta: Clone + Debug + Eq + Send + Sync;
}

pub struct PrimeField;   // Coeff = Felt, Ops = PrimeOps, BasisMeta = ()
pub struct Rationals;    // Coeff = BigRational, Ops = RationalOps,
                         // BasisMeta = RationalMeta
```

`DomainOps` is the one trait the shared polynomial code calls. It is held
by value in the ring, and it replaces the `p: u64` argument that every
helper in `src/poly.rs` takes today. The name says "domain operations"
and not "arithmetic", because it also carries the conversion and the
formatting the generic parser and `Display` need:

```rust
pub trait DomainOps {
    type Coeff;
    fn one(&self) -> Self::Coeff;
    fn is_one(&self, a: &Self::Coeff) -> bool;
    fn add(&self, a: &Self::Coeff, b: &Self::Coeff) -> Option<Self::Coeff>;  // None is zero
    fn sub(&self, a: &Self::Coeff, b: &Self::Coeff) -> Option<Self::Coeff>;
    fn neg(&self, a: &Self::Coeff) -> Self::Coeff;
    fn mul(&self, a: &Self::Coeff, b: &Self::Coeff) -> Self::Coeff;
    fn inv(&self, a: &Self::Coeff) -> Self::Coeff;       // a is nonzero by the invariant
    fn convert(&self, value: &Coefficient) -> Result<Option<Self::Coeff>, RingError>;
    fn write(&self, a: &Self::Coeff, f: &mut fmt::Formatter<'_>) -> fmt::Result;
    fn is_negative(&self, a: &Self::Coeff) -> bool;
    fn heap_bytes(&self, a: &Self::Coeff) -> usize;
}
```

No stored coefficient is zero, so zero is `None` rather than a value:
`add` and `sub` return `None` when the result cancels, and `convert`
returns `Ok(None)` for an input that reduces to zero, which is what the
merge loops and constructors in `src/poly.rs` already branch on. `write`
and `is_negative` are what generic `Display` needs; over `F_p`
`is_negative` is always false, so the sign path is dead there and the
prime-field output is unchanged. `PrimeOps` and `RationalOps` are public types
with private fields, because a public trait's associated types cannot be
private.

`BasisMeta` is a crate-private field of `GroebnerBasis<D>`, read through
one accessor per domain rather than a generic one:
`GroebnerBasis<PrimeField>` has none, since `()` says nothing, and
`GroebnerBasis<Rationals>::lift` returns `Option<&ModularLift>` (section
4.2).

A runtime tag on the ring, with an enum of term storages, was the first
draft. The type parameter wins on three points.

1. It removes runtime errors that have no reason to exist.
   `groebner_basis_certified` is an inherent method on
   `Ideal<PrimeField>` and does not exist on `Ideal<Rationals>`, so "you
   cannot certify over Q" is a compile error, not a `CertifyError`
   variant. `PolynomialRing::modulus` stays `u64` on the prime-field ring
   instead of becoming `Option<u64>` for every caller.
2. The engines keep concrete types. `src/compute/classic.rs`,
   `src/cert/**`, and `src/compute/f4/**` take `Polynomial<PrimeField>`
   and go on using `Felt` and `p` directly. A runtime tag would put a
   domain match, or a debug assertion, on every one of those paths.
3. One generic kernel serves both domains for the shared work: parsing,
   `Display`, term merging, division, and interreduction. The alternative
   was two copies.

The cost lands in two places. Generic code gets a `D: Domain` bound in
every signature, and code that learns the domain at run time, which is
the CLI and the Python bindings, holds an enum over the two
instantiations and branches once per entry point (sections 8.1 and 9.2).
Section 13, item 1 records the runtime tag as the alternative the owner
can still take.

### 2.2 What changes in the existing code

The migration is mechanical and wide. The document lists it so no chunk
discovers it late.

| file | change |
|---|---|
| `src/ring.rs` | `PolynomialRing<D>`, `Domain`, `DomainOps`, `Coefficient`, `rationals`, `polynomial` over `Coefficient`, the parser over `D` |
| `src/ring/field.rs` | `Felt` becomes public and gains `value()`; `PrimeOps` holds `p` |
| `src/ring/rational.rs` | new: `RationalOps`, clear denominators, content division, the image of a cleared polynomial mod `p` |
| `src/poly.rs` | `Polynomial<D>`, `Term<D>`; `add`, `sub`, `scale_monomial`, `sub_scaled`, `make_monic`, `normal_form_refs`, `s_polynomial`, `push_term`, `check_operand` take `&D::Ops` in place of `p: u64`; `Monomial::checked_mul` and `sub_scaled_checked` (section 5); `Display` gains the sign path; `heap_bytes` counts coefficients |
| `src/ideal.rs` | `Ideal<D>`, `GroebnerBasis<D>` with `D::BasisMeta`; manual `PartialEq`; the compute methods split into one `impl` per domain; the public checked constructor of section 5 |
| `src/certificate.rs` | `CertifiedGroebnerBasis` holds `GroebnerBasis<PrimeField>`. It gains no parameter of its own, since rational certification does not exist |
| `src/compute/mod.rs` | signatures over the prime-field types; `ComputeLimits` (section 3.8); `verify::Limits` is imported as `VerifyLimits` and `cert::Budget` is renamed `WriterBudget`, because the public `Budget` of section 5 takes the short name |
| `src/compute/classic.rs`, `src/compute/interreduce.rs`, `src/cert/**` | prime-field signatures; bodies unchanged except the memory meters below |
| `src/compute/f4/mod.rs` | prime-field types at the boundary (`start`, `max_exponent`, `intern_generators`, `convert_out`), `Limits` folded into `ComputeLimits`, an internal thread count forced to 1 by the modular driver, and the cancellation check of section 3.8 in `Deadline::tick` |
| `src/compute/f4/kernel.rs` | `Deadline::fork` for the per-worker deadlines it builds, so a worker sees the cancellation flag. The reduction arithmetic is untouched |
| `src/compute/f4/{basis,pairs,symbolic,matrix,monomial,field,trace}.rs` | unchanged |
| `src/verify/**` | unchanged, and still importing nothing from the engine side |
| `tests/**`, `benches/**`, `benchmarks/gb-comparison/runner/**` | follow the signatures |

The F4 algebra, the monomial table, the matrix layout, and the reduction
arithmetic do not change. The control plane does, in three places: one
absolute deadline instead of a per-call one, a cancellation flag that
reaches the per-worker deadlines `kernel.rs` builds, and a thread count
the modular driver forces to 1.

`Monomial`, the sealed grevlex order, the exponent bound of 65,535, the
256-variable cap, and both certificate contracts do not change. A
monomial carries no coefficient, so nothing about the order or the
certificates is domain dependent.

### 2.3 Invariants

- Terms are sorted strictly ascending under grevlex, so the last term is
  the leading term. This remains unchanged in both domains.
- No stored coefficient is zero. A `Felt` is a reduced residue in
  `[1, p)`. A `BigRational` is in lowest terms with a positive
  denominator, which `num_rational::BigRational` maintains; the crate
  checks nonzero.
- Two rings are equal when they hold the same `D::Ops` value (the
  prime, or nothing) and the same variables in the same order.
  `PolynomialRing<D>` keeps content equality and gains a content `Hash`,
  which the Python bindings need (section 9.2). Two separate
  constructions of one ring stay equal and hash equally.
- `GroebnerBasis<D>` compares the ring and the ordered polynomial list.
  `D::BasisMeta` is not part of equality, so two computations of one
  basis with different prime counts stay equal. `PartialEq` is written by
  hand for that reason, and section 3.7 compares this algebraic value
  when it asks whether a lift changed.
- Memory accounting counts coefficients. `DomainOps::heap_bytes` returns 0
  for `Felt` and, for a `BigRational`, a documented estimate from the
  used bits of the numerator and the denominator rounded up to limbs.
  `num-bigint` does not expose its allocated capacity, so the number is a
  logical estimate and the docstring says so.

Three meters size polynomials today, and every one of them gains the
coefficient bytes: `src/compute/classic.rs` (origin cofactors),
`src/cert/mod.rs` (the writer budget), and `src/compute/interreduce.rs`
(`term_bytes`, which counts `size_of::<Term>()` alone and misses the
spilled exponents as well). The modular driver's ledger (section 3.8) and
the normal form temporaries are the two new meters.

### 2.4 Construction, text, and the term API

`PolynomialRing::rationals(variables)` mirrors `prime_field`: the same
name rules, the same 256-variable cap, the same `RingError` variants for
the name checks, and no modulus check. Both are inherent methods on
distinct instantiations, so `PolynomialRing::prime_field(32003, ["x"])`
infers `D = PrimeField` with no turbofish.

One input type serves both domains, and the ring converts it:

```rust
pub enum Coefficient {
    Small(i64),
    Integer(BigInt),
    Fraction { numerator: BigInt, denominator: BigInt },
    Rational(BigRational),
}
// From<i8, i16, i32, i64, i128, isize, u8, u16, u32, u64, u128, usize,
//      BigInt, BigRational>
```

Every ordinary integer width has an impl so that
`ring.polynomial([(1, [2, 0])])` compiles: an unsuffixed literal defaults
to `i32`, and inference does not pick `i64` out of a set of candidate
`Into<Coefficient>` impls. `From<u64>` and the wider integers use
`Small` when the value fits an `i64` and `Integer` otherwise.

`Fraction` is the raw pair, unnormalized, because `From` cannot fail and
a zero denominator must not panic. There is no `From<(i64, i64)>`; a
caller writes `Coefficient::Fraction { .. }` or builds a `BigRational`
itself.

`DomainOps::convert` normalizes before it does anything else, in both
domains: a zero denominator is `RingError::ZeroDenominator`, and
otherwise the fraction is reduced by the greatest common divisor and
given a positive denominator. Only then does the prime-field arm test the
denominator. Without that order, `3/3` over `F_3` would be rejected as
noninvertible while the equal value `1` is accepted. Named tests: `3/3`,
`0/3`, and `3/6` over `F_3`.

Over `F_p` a normalized integer reduces mod `p`, and a normalized
fraction becomes `a * b^{-1} mod p` unless `p` divides `b`, which is
`RingError::CoefficientNotInvertible { denominator: BigInt }`. Over `Q`
the normalized value is kept.

The text grammar gains one production and nothing else:

```text
term        := coefficient [ "*" powers ] | powers
coefficient := integer [ "/" integer ]
```

The denominator is a positive decimal integer with no sign and no
surrounding whitespace. A zero denominator is
`ParseError::ZeroDenominator`. `Display` writes `a/b` over `Q` when
`b > 1` and `a` otherwise, and separates terms with ` + ` or ` - ` by the
sign, printing the absolute value after a ` - `. A proptest asserts that
parse and `Display` round trip in both domains, over the coefficients
`1`, `-1`, `1/2`, `-1/2`, and a bare constant.

The term accessors keep their shape, with the coefficient type following
the domain:

```rust
impl<D: Domain> Polynomial<D> {
    pub fn terms(&self) -> impl ExactSizeIterator<Item = (&D::Coeff, &[u16])>;
    pub fn leading_term(&self) -> Option<(&D::Coeff, &[u16])>;
}
```

One slice iterator, no enum, no allocation, no per-term match.
`Felt::value() -> u64` recovers the residue a caller previously read directly.

## 3. The multimodular rational engine

### 3.1 The pipeline

`src/compute/modular/mod.rs` holds the driver. One run over the input
list `F` of rational polynomials:

1. **Normalize.** Drop the zero generators. An empty list after that is
   the zero ideal, whose reduced basis is empty; return it without
   running anything. Multiply each remaining generator by the least
   common multiple of its denominators, then divide by the greatest
   common divisor of the resulting integer coefficients. Scaling one
   generator by a nonzero rational leaves the ideal alone, so the step is
   exact and needs no record.
2. **Take the next prime** in the sequence of section 3.2 and reduce
   every generator modulo it. A prime that divides the leading
   coefficient of a cleared generator is skipped before any run starts,
   and counted.
3. **Run the backend** the options name over `F_p`, on one thread,
   getting `G_p`. Section 3.5 says what happens when the run fails.
4. **Vote** on the leading monomial set and pick the prevalent category
   (section 3.3).
5. **Combine** the runs of the prevalent category by CRT over the union
   of the monomials that occur (section 3.4).
6. **Lift** each combined residue to a rational number (section 3.6). A
   coefficient that does not lift sends the loop back to step 2.
7. **Stop** when the state machine of section 3.7 says so.

The result is monic, interreduced, and sorted descending by leading
monomial, the same shape a prime-field basis has.

### 3.2 Prime selection

Primes come from one deterministic descending sequence of odd primes,
starting at `2^31 - 1` and walking down the odd numbers. The prime 2 is
not in the sequence. Each prime of the range contributes about 31 bits to
the modulus and 2 contributes one, so including it would cost a full F4
run for one bit.

The sequence is a function of its position and nothing else, so two runs
over one input use the same primes in the same order at any thread count.

The sequence is finite. Exhausting it is
`ComputeError::PrimesExhausted`, a typed limit next to `DegreeLimit` and
`TableFull`. The variant exists because silently reusing a prime, or
silently returning an unlifted basis, is worse than a report that is
rarely seen.

The 31-bit bound is the one `PolynomialRing::prime_field` accepts, so
every modular image is a ring this crate can build and every prime run
uses the F4 field kernel as it stands. A rational ring has no such
bound of its own. Larger primes would cut the number of
runs and would need a wider kernel. That trade is outside this design.

Determinism has a price. The sequence is public, so an adversary who
picks the input knows which primes the confirmation of section 3.7 will
use. Section 13, item 2 asks the owner whether a seeded random sequence,
with the seed reported, is the better trade.

### 3.3 Unlucky primes and the prevalence vote

Start from what is true. For all but finitely many primes, the reduction
of the rational reduced basis modulo `p` is the reduced basis of the
modular ideal, so all but finitely many primes agree on the leading
monomial set. That is the whole justification for a vote.

What is not true, and what an earlier draft of this document claimed, is
that the input's leading coefficients identify the exceptions, or that a
strict containment between two observed leading ideals convicts one of
them. Both directions of failure occur. Take
`F = { x + y, x + (p + 1) y }` in `Q[x, y]` under grevlex. Both
generators are primitive with leading coefficient 1, so no leading
coefficient test skips `p`. Over `Q` the ideal is `(x, y)` and the
leading ideal is `(x, y)`. Modulo `p` the two generators coincide, the
ideal is `(x + y)`, and the leading ideal is `(x)`, strictly smaller. The
prime is unlucky, and it hides in a cofactor denominator that the input
does not show.

So the rule is a heuristic, named as one:

1. Group the consumed runs into categories by their leading monomial
   set. Every run is retained with its prime and its basis.
2. The prevalent category is the one with the most runs. A tie is broken
   by the lexicographically smallest sorted list of leading monomials.
   The tie break makes the choice deterministic and nothing more.
3. CRT accumulates only the runs of the prevalent category, in prime
   order.

Consuming one run is four steps, in this order, so two implementations
cannot disagree about a run that changes the winner:

1. Retain the result and add one to its category's count.
2. Recompute prevalence.
3. If the prevalent category changed, drop the accumulator, rebuild it
   from the retained runs of the new category in prime order, discard the
   held candidate, and reset the confirmation counter of section 3.7.
   The run just consumed is then already folded, or not, by the rebuild.
4. Otherwise fold the run if it is in the prevalent category, and
   otherwise leave it retained and unfolded.

A category can lose prevalence and win it back. Each transition is one
rebuild from what is retained, and the counter starts again at zero. No
further rule covers that case.

Retention costs the sum of every consumed modular output, over every
category, plus one accumulator. A speculative run discarded under section
3.8 is not retained, because it was never consumed. There is no
constant-factor bound relative to the accumulator: a run of `k` primes
retains `k` modular
bases whatever the vote does with them, and the term vectors, monomial
lists, and hash tables around them are counted too. The ledger of section
3.8 charges all of it, and section 13, item 4 offers a bounded pool with
recomputation as the alternative.

Monomials below the leading monomial take no part in the vote. A basis
coefficient divisible by `p` makes a monomial vanish in one run, which is
not evidence against that prime. Section 3.4 reads a missing monomial as
zero.

### 3.4 CRT

The accumulator holds one entry per basis element of the prevalent
category, keyed by the element's leading monomial, which the category
fixes. An entry maps a monomial to a residue modulo `M`, the product of
the folded primes, through a hash map, so a monomial first seen in the
fifth run needs no reindexing of the first four. The map is sorted into a
term list only when a candidate is built.

Folding a run extends `M` by one prime and updates each residue with the
standard incremental step. A monomial the run does not carry contributes
residue 0. A monomial the accumulator has not seen starts at residue 0
for the whole previous modulus, which is the correct value: the earlier
runs did say that coefficient was zero.

Runs are folded strictly in the order the prime sequence produced them,
never in completion order.

### 3.5 When a prime run fails

A prime run can stop with any `ComputeError`. The policy is fixed here,
because "unlucky prime, try another" is not a safe default for any of
them.

- `Timeout` stops the whole rational run and is returned as it is. The
  deadline is shared, so another prime would only spend what is gone.
- `MemoryLimitExceeded` under a per-run share is not global exhaustion.
  Section 3.8 says what the driver does before it reports it.
- `DegreeLimit`, `ExponentLimit`, and `TableFull` stop the whole rational
  run and are returned as they are. This is a policy, not an algebraic
  fact: these limits depend on the modular basis, on the pairs the run
  generates, and on its symbolic tables, so another prime can behave
  differently. The policy exists because the alternative is an unbounded
  retry loop that reports a structural limit as an unlucky prime, and
  because a run that reaches a structural limit says more about the
  ideal than about the prime, and nothing here measures how often the
  next prime would clear it.
  Section 13, item 5 offers a bounded retry as the alternative.

No error class is evidence about the prime. A prime is discarded only by
the vote of section 3.3, and the counts of skipped and discarded primes
are reported separately.

### 3.6 Rational reconstruction

For a residue `r` modulo `M`, find the unique reduced fraction `a / b`
with `b > 0`, `gcd(a, b) = 1`, `gcd(b, M) = 1`, `a ≡ r b (mod M)`,
`|a| <= A`, and `b <= B`, where `A = B = floor(sqrt((M - 1) / 2))`, so
`2 A B < M`. Uniqueness holds among reduced fractions inside those
bounds, and it needs `gcd(a, b) = 1` stated: without it, `0/1` and `0/2`
both satisfy the rest.

The search is the extended Euclidean algorithm on `(M, r)`, stopped at
the first remainder at or below `A`, taking the numerator from the
remainder and the denominator from the matching cofactor. The
implementation then verifies the result against every condition above,
including the reduction and the congruence, and returns "no lift" rather
than an unverified pair. `r = 0`, `b = 1`, and a negative numerator are
named cases in the tests.

`src/compute/modular/reconstruct.rs` owns this and nothing else. It is a
pure function over `num-bigint`, tested directly.

A lift that succeeds is the right answer only if `M` is large enough for
the true coefficients. Nothing inside the lift can know that. Section 3.7
is where it is addressed, and it does not settle it either.

### 3.7 Stopping: the state machine, and what it establishes

```rust
pub enum RationalStop {
    /// Stop when `extra` further primes of the prevalent category leave
    /// the lift unchanged. A heuristic. It establishes nothing about the
    /// ideal.
    Unchanged { extra: NonZeroUsize },        // the default, extra = 2
    /// Also run T1 and T2. The basis is then the reduced Gröbner basis
    /// of an ideal that contains the input ideal.
    ContainsInput { extra: NonZeroUsize },
}
```

`extra` is a `NonZeroUsize` because zero confirming primes is a stopping
rule that observes nothing, and a public field cannot be validated by a
builder.

The loop is one state machine over four pieces of state: the retained
runs, the accumulator, a held candidate, and a confirmation counter. A
fifth flag, `rejected`, records that the held candidate failed T1 or T2.
Every step consumes the next prime of the sequence, runs it, and applies
the four-step transition of section 3.3. The driver never chooses a
prime by the category it will land in, because it cannot know that before
the run.

1. Consume runs until the prevalent category holds two of them. A run in
   a losing category is retained and folded into nothing; it neither
   increments nor resets the counter, because the counter counts
   confirmations of the prevalent category.
2. Build the candidate: CRT over the category's runs, then lift. Clear
   `rejected`. A failed lift consumes the next prime and repeats this
   step.
3. Set the counter to zero. Consume further primes. For a run that lands
   in the prevalent category, fold it, lift again, and compare the new
   candidate with the held one as an algebraic value: same ring, same
   ordered polynomial list, ignoring `BasisMeta`. Equal increments the
   counter. A different candidate replaces the held one, clears
   `rejected`, and sets the counter back to zero. A failed lift holds no
   candidate at all: it drops the held one, clears `rejected`, sets the
   counter to zero, and leaves the machine waiting for a candidate, which
   the next successful lift becomes.
4. A change of prevalent category rebuilds as section 3.3 says and
   returns to step 2. Candidates from two categories are never compared:
   their leading monomial sets differ, so the comparison would mean
   nothing.
5. Stop when the counter reaches `extra`. Under
   `RationalStop::ContainsInput`, run T1 and T2 first, unless `rejected`
   is already set, in which case go straight back to step 3 without
   rerunning them. A candidate the tests reject sets `rejected` and
   returns to step 3, so the exact tests run again only after the lift
   itself changes. Without that flag, every further confirming prime
   would pay for the same rejected candidate again.
6. An exhausted budget, an exponent limit, or any other `ComputeError`
   raised inside T1 or T2 is not a rejection. It stops the whole rational
   run and is returned as itself. Only a nonzero remainder is evidence
   about the candidate.

So a stop takes at least `extra + 2` runs of one category. What a
confirming prime gives: folding prime `q` and getting the same lift means
the lifted basis agrees with the run over `q`, which the earlier lift did
not use. That is evidence of the same shape as the check Singular calls
`pTest_std` and the check Groebner.jl performs modulo a fresh prime, and
weaker than either, because those pick their test prime at random and
this sequence is fixed and public. No probability statement attaches to
it. The variant is named `Unchanged` after what it observes.

T1 and T2 are the exact tests over `Q`, run on the lifted basis `G` with
the division of section 5. They do not call the public `normal_form`,
which would start a fresh duration and a fresh memory scope. They call
the internal division entry point, which takes the run's absolute
deadline and the residual of the same ledger. Before they start, the
driver cancels and joins every outstanding prime run, so their working
memory is the whole residual and not a share of it. They run under a
fresh cancellation flag, per the generation rule of section 3.8.

- **T1.** Every generator in `F` reduces to zero modulo `G`. Then
  `<F>` is contained in `<G>`.
- **T2.** For every pair in `G`, the S-polynomial reduces to zero modulo
  `G`. Then `G` is a Gröbner basis of `<G>`, and being monic and
  interreduced, the reduced one.

Together: `G` is the reduced Gröbner basis of an ideal that contains
`<F>`. The containment of `<G>` in `<F>` is not tested and does not
follow. The gap is exactly the failure this repository already knows:
`G = {1}` passes T1 and T2 on every input. `tests/known_defects.rs`
exists because an output-only check accepts `{1}`, so the docstring of
`RationalStop::ContainsInput`, the README, and this section all say it.
Closing the gap needs cofactors over `Q`, which is section 4.

Cost. T2 enumerates every pair of `G` and divides over `Q`, with the
coefficient growth the multimodular pipeline exists to avoid. Whether it
dominates a run is a measurement the harness will make; no factor for it
appears in this document or in the README before that record exists.

### 3.8 Limits, threads, cancellation, determinism

**One absolute deadline.** `ComputeOptions::deadline()` builds
`Instant::now() + timeout` at every call, and `f4::Limits::of` does the
same, so calling `compute::groebner_basis` once per prime would restart
the timeout on every prime. The design promotes `f4::Limits` to
`compute::ComputeLimits` and adds an internal
`compute::groebner_basis_with_limits(ring, generators, options,
&ComputeLimits)`. The driver builds `ComputeLimits` once and passes the
same absolute deadline to every prime run. The public entry points build
it once at their top and call the same function, so there is one path.
The rename matters: `src/compute/mod.rs` already imports `verify::Limits`
and `cert::Budget`, which become `VerifyLimits` and `WriterBudget`, and
`f4::solve`, `f4::solve_recorded`, `groebner_basis`,
`groebner_basis_with_report`, and `groebner_basis_certified` all move to
the new type.

The deadline covers the driver's own work, not only the prime runs.
Normalization, the prime search, the modular images, the vote, the
accumulator rebuild, the CRT fold, and the reconstruction all check it in
their loops, once per element rather than once per phase. One big-integer
operation cannot be interrupted, so a check sits before and after each
one, and the granularity of the deadline is the cost of one such
operation. Without those checks, `ComputeOptions::timeout` would cover
the prime runs and not the whole rational computation.

**Cancellation.** `ComputeLimits` carries an optional cancellation flag,
an `&AtomicBool`. `Deadline::tick` in `src/compute/f4/mod.rs` and the
classic backend's deadline check read it, and `Deadline::fork` passes it
to the per-worker deadlines `src/compute/f4/kernel.rs` builds, so a
cancelled run stops inside a parallel batch too.

A flag is never reset. Each speculative wave gets a fresh
`Arc<AtomicBool>`, and so does each exact-test phase. Cancelling a wave
therefore cannot leak into the work that follows it, which resetting one
shared flag would risk: a worker that has not yet joined would clear a
flag the next phase depends on, or the next phase would start already
cancelled.

Cancellation stays below the public error boundary. `ComputeError` is a
public exhaustive enum and gains no variant: adding one would change
every downstream match, and a caller who never asked for cancellation
would have to handle it. The backends return a crate-private
`RunError::{Compute(ComputeError), Cancelled}`, the modular driver
consumes `Cancelled`, and the public entry points map the other arm
straight through. Nothing else in the crate sees the difference.

**Threads.** `ComputeOptions::threads(n)` runs up to `n` prime runs at
once. Each prime run is called with an internal options value whose
thread count is 1, so the outer concurrency is `n` and not `n^2`. This is
where the parallelism of the rational path lives. If the outer pool
cannot be built, the driver runs the primes on the calling thread and the
report says the concurrency it had.

**Speculation.** Runs are started in sequence order and consumed in
sequence order. When the state machine stops, the driver sets the
cancellation flag, joins every outstanding run, and discards its result
whole: it touches no vote, no accumulator, no lift, and no counter. A
cancelled run that crosses the deadline while it is winding down does not
turn the result into a `Timeout`, because the logical stop happened
first. That rule is stated so the outcome does not depend on how fast a
worker notices the flag.

**The ledger.** The memory limit charges, all by the meter of section
2.3: the cleared integer generators, the retained modular bases of every
category, the category tables, the CRT accumulator with its big-integer
bytes, the CRT and reconstruction scratch values, the modular ring and
generator images in flight, and the held candidate.

There is no separate reservation for an output whose size is not yet
known. The driver splits the whole residual into one share per concurrent
run. A run's working memory is charged to its share, and the share covers
both the engine structures and the output the run builds out of them,
which is what the engine's own accounting already does at `convert_out`.
When a run ends, the driver charges the measured size of its output to
the ledger and releases the share in the same step, so no window exists
in which either both or neither is counted.

A run that
reports `MemoryLimitExceeded` under that share has not shown that the
whole budget is gone. The driver cancels and joins the runs after it in
the sequence, discards their results, and retries the same prime alone
with the whole residual. Only a second failure is reported as
`ComputeError::MemoryLimitExceeded`. Speculation restarts after the
retry. The residual is that the ledger counts
tracked bytes and not process resident set size.

**Determinism, with its condition.** When the run succeeds, the basis and
every counter in `ModularLift` are identical at 1, 2, 8, and 16 threads,
because the prime sequence, the consumption order, the vote, the fold
order, and the lift are all deterministic and speculative results are
discarded. Whether a run succeeds is not thread-count independent: the
per-run memory share of `n` concurrent runs is smaller than of one, so a
thread count can be the difference between a completed run and
`MemoryLimitExceeded`. The report carries the concurrency, and the test
asserts equality across thread counts for runs that succeed.

### 3.9 What the design leaves out: learn and apply

msolve and Groebner.jl do not run full F4 per prime. A first run, or in
Groebner.jl's case a few of them, records the trace of the computation,
which pairs entered which batch and which rows reduced against which, and
later primes replay it: no symbolic preprocessing, no pair selection,
only the linear algebra.

This design does not implement it, for two reasons. The replay needs a failure
path for a pivot that was nonzero over the learning prime and vanishes
over a later one, which is a second notion of unluckiness with its own
detection and fallback. And the F4 engine's existing trace recorder
serves `docs/certificate-v2.md` section 9.4, whose contract is frozen and
is not the trace a replay needs; making one recorder serve both is a
design job of its own.

The consequence: the rational path runs a full F4 per prime, and the
prime runs are concurrent. How that compares
with msolve is what the record will say, and no ratio is claimed here.

## 4. Certification over Q

### 4.1 What a per-prime certificate would prove

Over `F_p`, `groebner_basis_certified` returns a value only after an
isolated verifier accepts bytes. `docs/certificate-v1.md` and
`docs/certificate-v2.md` section 2 state the facts that acceptance
establishes.

Over `Q` a run could write and verify one certificate per prime. A caller
holding `k` accepted certificates would know, for each prime, that `G_p`
is the reduced Gröbner basis of `F mod p`. That says nothing about `G`
over `Q`. The lift is outside every certificate: the vote, the CRT, and
the reconstruction are steps no verifier checked, and they are exactly
where a wrong answer comes from. A bundle of accepted certificates next
to an unverified lift is a certified-looking value with an uncertified
conclusion, which is worse than an honest unverified one.

### 4.2 The stopping rule

There is no certified path over `Q`, and the type system says so:
`groebner_basis_certified` is an inherent method on `Ideal<PrimeField>`,
so `Ideal<Rationals>` does not have it and no runtime error variant is
needed. `CertifiedGroebnerBasis` holds a `GroebnerBasis<PrimeField>` and
takes no parameter of its own. The CLI and the Python bindings, which
learn the domain at run time, report the refusal at their own boundary
(sections 8.4 and 9.3).

Every basis carries what its domain records about its own origin, through
`D::BasisMeta`:

```rust
pub enum RationalMeta {
    /// The multimodular driver produced it.
    Lifted(ModularLift),
    /// The caller supplied it and `from_polynomials` checked it
    /// (section 5). No modular run produced it.
    Checked,
}

impl GroebnerBasis<Rationals> { pub fn lift(&self) -> Option<&ModularLift>; }
```

The lift itself is counters and one proposition:

```rust
pub struct ModularLift {
    pub primes_consumed: usize,     // runs that finished and were consumed
    pub primes_skipped: usize,      // divided a leading coefficient, never run
    pub primes_folded: usize,       // in the final accumulator
    pub primes_discarded: usize,    // consumed but not in the final accumulator
    pub confirming_primes: usize,   // the counter of section 3.7 at the stop
    pub modulus_bits: u64,
    pub established: Established,
}
```

```rust
pub enum Established {
    /// The lift did not change over the confirming primes. Nothing about
    /// the ideal follows.
    Unchanged,
    /// Also T1 and T2 of section 3.7. The basis is the reduced Gröbner
    /// basis of an ideal that contains the input ideal. Equality is not
    /// established; the basis {1} passes both tests.
    ContainsInput,
}
```

Over `PrimeField`, `BasisMeta` is `()` and there is no accessor, because
`&()` reads as information and is not. Over `Rationals`, `lift` returns
`None` for a checked basis, which is the honest answer: the value exists
and no modular run produced it. Nothing fabricates counters for it.

The counters are exact, not indicative. `primes_consumed` counts the runs
the driver consumed in sequence order; a speculative run discarded under
section 3.8 is not one of them and appears in no counter.
`primes_folded + primes_discarded = primes_consumed` holds at the stop,
and a run folded into a category that later lost prevalence counts as
discarded, since the count is relative to the final accumulator.
`primes_skipped` counts primes rejected before any run and is not part of
that sum. `Established` names propositions; the counters sit beside it.

### 4.3 What a rational certificate would take

This is here because the answer shapes later work. Certification over `Q`
needs, for each `g` in `G`, an identity `g = sum_i h_i f_i` over `Q` that
a verifier can check by multiplying out.

Reconstructing the `h_i` multimodularly is not a matter of running the
existing cofactor tracking at many primes. Cofactors are not unique, and
the classic backend can take different reduction paths at different
primes, so per-prime cofactors need not be reductions of one rational
object. Coefficient-wise CRT over them has no guaranteed common rational
target: per-prime cofactors may align, and nothing in the engines makes
them. A rational cofactor lift therefore needs a derivation fixed across
primes, which is the
learn-and-apply trace of section 3.9 extended to carry a change matrix,
or a separate step that solves for a representation over `Q` once `G` is
known.

Size is the second obstacle, and it is a hypothesis here, not a
measurement: cofactor coefficients over `Q` are expected to be much
larger than basis coefficients, which would push both the certificate
bytes and the verifier's multiplication past what the caps in
`docs/certificate-v2.md` section 8 are set for. Before anyone designs a
v3 schema, that expectation needs a number from a real cofactor lift.

## 5. Normal form, membership, and a basis the caller brings

Division by a basis is one algorithm, in `src/normal_form.rs`, generic
over `D: Domain`, used by four callers: `GroebnerBasis::normal_form`, the
T1 and T2 tests of section 3.7, the checked constructor below, and the
CLI's `normal-form` subcommand. The verifiers keep their own division and
do not use it, which `tests/isolation.rs` continues to check.

```rust
pub struct Budget { timeout: Option<Duration>, memory_limit: Option<usize> }

impl Budget {
    pub fn new() -> Self;                              // no deadline, no limit
    pub fn timeout(self, timeout: Duration) -> Self;
    pub fn memory_limit(self, bytes: usize) -> Self;
}
impl ComputeOptions { pub fn budget(self, budget: Budget) -> Self; }

impl<D: Domain> GroebnerBasis<D> {
    pub fn normal_form(&self, f: &Polynomial<D>, budget: Budget)
        -> Result<Polynomial<D>, NormalFormError>;
    pub fn contains(&self, f: &Polynomial<D>, budget: Budget)
        -> Result<bool, NormalFormError>;
}
```

The checked constructor is one inherent method per domain, over one
private generic checker:

```rust
impl GroebnerBasis<PrimeField> {
    pub fn from_polynomials(ring: &PolynomialRing<PrimeField>,
        polynomials: Vec<Polynomial<PrimeField>>, budget: Budget)
        -> Result<Self, BasisError>;          // meta: ()
}
impl GroebnerBasis<Rationals> {
    pub fn from_polynomials(ring: &PolynomialRing<Rationals>,
        polynomials: Vec<Polynomial<Rationals>>, budget: Budget)
        -> Result<Self, BasisError>;          // meta: RationalMeta::Checked
}
```

The constructor is one per domain because `Domain` offers no way to build
a `D::BasisMeta`, and adding a `Default` bound to invent one would say
that every domain has a meaningful empty provenance. The checking work is
a private generic function both call.

`ComputeOptions::timeout` and `memory_limit` stay as setters that write
into the held `Budget`, so the common call is unchanged and the values
live in one place. `Budget` holds a duration; the deadline starts when
the call that receives it starts. `ComputeLimits` (section 3.8) is the
absolute form, and only the crate builds it. The budget charges the
working polynomial, the remainder, and any temporary the division holds,
by the same meter as section 2.3.

`contains` is `normal_form(...)?.is_zero()`, and it is a separate method
because that is the question callers ask.

The algorithm is full division: while the largest monomial of the working
polynomial is divisible by the leading monomial of some basis element,
subtract the matching multiple; a monomial no basis element divides moves
to the remainder. The basis is a Gröbner basis, so the remainder is the
unique normal form whatever order the divisors are picked in. The
implementation takes the first divisor in basis order, which makes it
deterministic as well as unique. Over `Q` the coefficients are exact and
can grow, which is what the budget is for.

`from_polynomials` is what the CLI and the Python bindings need to act on
a basis a user supplies. `GroebnerBasis::new` is crate-private because a
`GroebnerBasis` that is not one would make every method on it
meaningless. The constructor therefore checks, under the budget, that
the list holds no zero polynomial and is monic, sorted strictly
descending by leading monomial, duplicate-free, interreduced, and closed
under T2 of section 3.7.

```rust
pub enum BasisError {
    RingMismatch,
    ZeroPolynomial { index: usize },
    NotMonic { index: usize },
    NotSorted { index: usize },
    NotInterreduced { index: usize },
    NotGroebner { left: usize, right: usize },
    ExponentLimit { limit: u32 },
    Timeout,
    MemoryLimitExceeded,
}
```

`NotGroebner { left, right }` names the first index pair in ascending
lexicographic order, over the caller's own indices, whose S-polynomial
has a nonzero remainder. Naming the first pair makes the error a function
of the input and not of the order the checker happened to walk.

`ExponentLimit` is here for the same reason as in `NormalFormError`.
Building an S-polynomial multiplies the tail of each side by a quotient
monomial, and that product can pass 65,535 even when both leading
monomials and their least common multiple fit. The check uses
`Monomial::checked_mul` and `Polynomial::sub_scaled_checked` throughout,
so the case is an error and not a panic.

A checked basis is a check, not a certificate. Over `F_p` the certified
path is still the only one an isolated verifier stands behind, and over
`Q` the basis records `RationalMeta::Checked`, so `lift()` returns
`None`.

`NormalFormError` has four variants: `RingMismatch`,
`ExponentLimit { limit }`, `Timeout`, and `MemoryLimitExceeded`.

`ExponentLimit` is not theoretical, and it needs work in `src/poly.rs`.
Divide `x^65535 * y` by a basis holding `x^65535 + y^65535`. The
quotient monomial is `y`, and the tail multiple needs `y^65536`. Every
stored exponent fits a `u16`; the intermediate does not. Two functions
multiply monomials on this path today and both call `expect` on the
addition: `Monomial::mul`, correct for the engines because they bound a
pair's degree first, and `Polynomial::sub_scaled`, which multiplies the
multiplier into every term of the reducer. The design adds
`Monomial::checked_mul` and `Polynomial::sub_scaled_checked`, both
returning the exponent failure, and `src/normal_form.rs` uses only those.
The `expect` on the engine path stays where it is.

## 6. Hilbert series and dimension

### 6.1 What is computed, and what it equals

From the leading monomial ideal `L = <lm(g) : g in G>`, `src/hilbert.rs`
computes the Hilbert series of `k[x] / L` as a rational function

    N(t) / (1 - t)^n,     n = nvars,

with `N` a polynomial over `Z`. The claim is exact for `k[x] / L`. For
the ideal `I` the caller asked about, two statements hold and a third
does not.

- If `I` is homogeneous, the series of `k[x] / L` is the Hilbert series
  of `k[x] / I`.
- For any `I`, grevlex is a graded order, so the coefficient of `t^d` is
  the first difference of the affine Hilbert function of `I`, that is
  `H(d) - H(d - 1)` with `H(d) = dim_k R_{<=d} / (I ∩ R_{<=d})` and
  `H(-1) = 0`. The coefficients are not the affine Hilbert function
  itself. For the zero ideal in one variable they are 1 in every degree
  while `H(d) = d + 1`.
- The Krull dimension of `k[x] / I` equals the dimension read off this
  series. Passing to an initial ideal preserves dimension, and that is a
  theorem, so `krull_dimension` does not overclaim.

```rust
pub struct HilbertSeries { numerator: Vec<BigInt>, denominator_power: usize }

impl HilbertSeries {
    pub fn numerator(&self) -> &[BigInt];        // dense, low degree first
    pub fn denominator_power(&self) -> usize;    // n, before cancellation
    pub fn dimension(&self) -> Option<usize>;
    pub fn multiplicity(&self) -> Option<BigInt>;
    pub fn coefficient(&self, degree: u32) -> BigInt;  // dim of one graded piece
}
```

The basis is where a caller reaches it:

```rust
impl<D: Domain> GroebnerBasis<D> {
    pub fn hilbert_series(&self, budget: Budget) -> Result<HilbertSeries, HilbertError>;
    pub fn krull_dimension(&self, budget: Budget) -> Result<Option<usize>, HilbertError>;
    pub fn is_homogeneous(&self) -> bool;
}
impl<D: Domain> Ideal<D> { pub fn has_homogeneous_generators(&self) -> bool; }
```

`dimension` first tests whether `N` is the zero polynomial, which is the
unit ideal: the quotient ring is zero, its series is 0, and its dimension
is -1 by the usual convention. Both `dimension` and `multiplicity` return
`None` there, and the docstring names the convention. Cancelling `(1 - t)`
out of a zero numerator would not terminate, which is why the test comes
first. Otherwise `dimension` divides `N(t)` by `(1 - t)` while the
division is exact, `k` times, and returns `n - k`; `multiplicity` is the
value at 1 of what is left.

The two homogeneity methods are named for what they check.
`Ideal::has_homogeneous_generators` looks at the generators, which is
sufficient and not necessary: a homogeneous ideal can be given by
inhomogeneous generators. `GroebnerBasis::is_homogeneous` decides it,
because under a graded order an ideal is homogeneous exactly when its
reduced Gröbner basis is.

The series reads leading monomials only, so it is the same over `F_p` and
over `Q` for a basis of the same ideal, and the implementation is one
generic function with no coefficient arithmetic in it.

### 6.2 The algorithm

The numerator comes from the recursion on the minimal monomial generators
of `L` (Bigatti; Bayer and Stillman), not from counting standard
monomials. The recursion is the numerator form of the exact sequence

    0 -> (R / (L : x_j))(-1) --*x_j--> R / L -> R / (L + (x_j)) -> 0,

which gives

    N(L) = N(L + (x_j)) + t * N(L : x_j).

An earlier draft justified this by an ideal identity,
`L = (L : x_j) ∩ (L + (x_j))`, which is false: for `L = (x^2)` both sides
of the intersection are `(x)`. The exact sequence is the justification,
and the numerator identity itself is correct.

Base cases, tested in this order:

- A generator equal to 1, the unit ideal: `N = 0`.
- No generators: `N = 1`.
- One generator `m`: `N = 1 - t^{deg m}`.
- Every generator a pure power of a distinct variable: `N` is the
  product of `1 - t^{d_i}`, one factor per generator.

Otherwise pick the pivot variable `x_j` from a generator that is not a
pure power, taking the variable that occurs in the most generators among
those that appear in such a generator, and breaking a tie by the smallest
variable index. Correctness and termination do not depend on the tie
break; the traversal and the budget it spends do, so it is fixed here. That restriction is what makes the
recursion terminate: `x_j` is then not itself a minimal generator of `L`,
so `L + (x_j)` and `L : x_j` both have a strictly smaller sum of minimal
generator degrees than `L`, and that sum is the termination measure. A
pivot chosen without the restriction can return the parent problem: for
`L = (x, y z)` the variable `x` gives `L + (x) = L`. Generator sets are
interreduced to their minimal generators before each recursive call, and
canonical sets are memoized, so a subproblem the recursion reaches twice
is computed once.

The recursion branches, and no proved bound on its work is offered here.
So `hilbert_series` takes a `Budget` and returns `HilbertError::Timeout`
or `HilbertError::MemoryLimitExceeded`. The budget charges the live
generator sets, the memo table with its keys, every cached and live
numerator vector, and the temporaries of the polynomial additions,
counting each `BigInt` coefficient by the estimate of section 2.3. The
numerators are the large values here, so a budget that counted only
generator sets would report a limit that does not bound what the
recursion holds. Numerator coefficients are `BigInt` because a
fixed-width type would need an overflow story that no measurement here
supports.

## 7. Public API changes

Backward compatibility is not a goal before API stability. Everything below
is a break, disclosed in `CHANGELOG.md`, with no alias left behind.

| item | before | after |
|---|---|---|
| ring, polynomial, ideal, basis | `PolynomialRing` | `PolynomialRing<D>`, `Polynomial<D>`, `Ideal<D>`, `GroebnerBasis<D>` |
| `PolynomialRing::modulus` | `u64` | `u64`, on `PolynomialRing<PrimeField>` only |
| ring constructors | `prime_field` | plus `rationals` |
| `Polynomial::terms` | `(u64, &[u16])` | `(&D::Coeff, &[u16])`; `Felt::value()` gives the `u64` back |
| `PolynomialRing::polynomial` | `(i64, exps)` | `(impl Into<Coefficient>, exps)` |
| budget | fields of `ComputeOptions` | a `Budget` value held by `ComputeOptions` |
| `groebner_basis_certified` | on `Ideal` | on `Ideal<PrimeField>` only |
| `CertifiedGroebnerBasis::basis` | `&GroebnerBasis` | `&GroebnerBasis<PrimeField>` |

Added: `Domain`, `PrimeField`, `Rationals`, `DomainOps`, `PrimeOps`,
`RationalOps`, `Felt` (now public), `Coefficient`, `Budget`,
`RationalOptions`, `RationalStop`, `RationalMeta`, `ModularLift`,
`Established`, `HilbertSeries`, `HilbertError`, `NormalFormError`,
`BasisError`, `ComputeError::PrimesExhausted`,
`RingError::{CoefficientNotInvertible, ZeroDenominator}`,
`ParseError::ZeroDenominator`, `GroebnerBasis::{normal_form, contains,
from_polynomials, hilbert_series, krull_dimension, is_homogeneous}`,
`GroebnerBasis<Rationals>::lift`, and
`Ideal::has_homogeneous_generators`.

Removed: nothing beyond the signature changes above. `Backend`,
`F4Counters`, `verify::verify`, and both certificate contracts are
untouched.

The rational compute call takes its own options type, so no option on it
is inert:

```rust
pub struct RationalOptions { compute: ComputeOptions, stop: RationalStop }
impl RationalOptions {
    pub fn new() -> Self;                                  // stop: Unchanged { extra: 2 }
    pub fn compute(self, options: ComputeOptions) -> Self;
    pub fn stop(self, stop: RationalStop) -> Self;
}
impl Ideal<Rationals> {
    pub fn groebner_basis(&self, options: RationalOptions)
        -> Result<GroebnerBasis<Rationals>, ComputeError>;
    pub fn groebner_basis_with_report(&self, options: RationalOptions)
        -> Result<(GroebnerBasis<Rationals>, ComputeReport), ComputeError>;
}
```

`ComputeOptions` keeps the backend, the threads, and the budget, and
gains no rational field. `ComputeReport` gains
`modular: Option<ModularLift>` and `modular_concurrency: Option<usize>`,
and both fields are defined on both paths:

| field | prime-field path | rational path |
|---|---|---|
| `counters` | `Some(F4Counters)` under F4, `None` under classic | `None`: the counters describe one F4 run, a rational run has many, and summing them would name a run that never happened |
| `threads_used` | the engine pool the run had | 1, the thread count of each prime run |
| `modular` | `None` | `Some(ModularLift)` |
| `modular_concurrency` | `None` | the number of prime runs the driver ran at once |

`ComputeReport` stays `Copy`, so `ModularLift`, `Established`, and
`RationalMeta` derive `Copy` too; every field of all three is a counter
or a unit variant. Per-prime counters are not reported; the
harness reads what `ModularLift` carries.

## 8. The CLI

### 8.1 Crate and shape

`crates/sylvester-cli`, binary `sylv`, depends on the workspace library
and on `clap` with the derive feature. The path supports workspace
builds. The version makes the published package resolve on crates.io. `clap`
is a dependency of the binary alone.

The binary learns the domain at run time, so it holds one enum over the
two instantiations and branches once, right after it reads the ring:

```rust
enum AnyIdeal { Prime(Ideal<PrimeField>), Rational(Ideal<Rationals>) }
```

Every subcommand matches that enum once and then runs typed code. This is
the boundary cost of section 2.1, and it is confined to `main.rs` and the
format readers.

### 8.2 Subcommands

Each subcommand has its own grammar. `FILE` is a path, or `-` for
standard input.

```
sylv gb [FILE] [--in-format F] [--out-format F] [-o PATH]
        [--modulus P | --rationals] [--backend f4|classic] [--threads N]
        [--stop unchanged|contains-input] [--extra-primes N]
        [--timeout SECS] [--memory BYTES]
        [--certified] [--certificate PATH] [--report]
```

```
sylv normal-form --basis BASIS_FILE [--basis-format F]
        (--poly TEXT | --poly-file FILE [--in-format F])
        [--out-format F] [-o PATH] [--modulus P | --rationals]
        [--timeout SECS] [--memory BYTES]
```

```
sylv hilbert [FILE] [--in-format F] [-o PATH]
        [--modulus P | --rationals] [--from-basis]
        [--backend f4|classic] [--threads N]
        [--stop ...] [--extra-primes N] [--timeout SECS] [--memory BYTES]

sylv verify CERT [--max-bytes N] [--timeout SECS] [--out-format F] [--quiet]
```

`normal-form` takes the basis and the polynomial from separate sources,
each with its own format flag, and never from one positional argument.
Both must parse in the same ring, and a mismatch is a usage error. At
most one of the two sources reads standard input, since one stream cannot
carry two independently parsed inputs. `--poly-file` holds exactly one
polynomial; a second one is a usage error. The basis file goes through
`GroebnerBasis::from_polynomials` (section 5), so a file that is not a
reduced Gröbner basis is rejected with the check that failed, and not
silently used.

`hilbert` computes a basis from generators by default. With
`--from-basis` it reads a basis and checks it the same way, and then
`--backend`, `--threads`, `--stop`, and `--extra-primes` describe a
computation that does not happen, so passing any of them with
`--from-basis` is a usage error. Its output is not a polynomial system,
so it has no `--out-format` and writes its own schema: a `series:` line
with the numerator coefficients low degree first, a
`denominator_power:` line, a `dimension:` line, and a `multiplicity:`
line, with `dimension: none` for the unit ideal.

`--certificate` implies `--certified`; passing `--certified` without
`--certificate` runs the certified path and writes no file, which is
still a stronger claim about the printed basis.

Ring resolution is per command, not per file, because `normal-form` reads
two files and they must land in one ring. Three rules, in order:

1. Every source that carries a ring (`ms`, `text`) must name the same
   one. A disagreement is a usage error.
2. If any source carries a ring, that ring is the command's, a `syl`
   source inherits it, and `--modulus` or `--rationals` is a usage error.
3. If no source carries a ring, exactly one of `--modulus` and
   `--rationals` is required.

So a `syl` basis with a `text` polynomial is legal and takes its ring
from the `text` file, which the earlier draft made impossible to express.
Standard input with no `--in-format` defaults to `text`; an unknown
extension is a usage error rather than a guess. `--stop` and
`--extra-primes` describe the rational path only, so passing either with
a prime-field ring is a usage error, in the same way the Python bindings
reject a rational keyword on a prime-field ideal (section 9.2).

`--timeout` covers the whole command, not one library call. The binary
takes an instant when it starts and passes the remaining time to each
call it makes, so reading a basis, checking it, and reducing against it
share one deadline. A command whose remaining time reaches zero exits 3.

### 8.3 Formats

- `ms`, the msolve input format: a line of comma-separated variable
  names, a line with the characteristic (`0` for `Q`), then the
  polynomials, comma-separated. It is msolve's format, not Singular's:
  Singular reads the polynomial expressions only after the harness wraps
  them in a ring declaration, and section 11.4 owns that wrapper.
- `syl`, the benchmark format under `benchmarks/gb-comparison/inputs`: a
  line with the variable count, then one polynomial per line, terms
  separated by `;`, each term a coefficient and one exponent per
  variable, separated by `,`. It carries no ring. Variables are named
  `x1 .. xn`, which is the convention `gen.py` uses in every other
  format it writes; the sylvester runner under
  `benchmarks/gb-comparison/runner` names them `x0 .. x(n-1)` today and
  moves to `x1 .. xn` in chunk A1. Names change nothing about grevlex or
  the computed basis, only the text. A rational coefficient is written
  `a/b` in the coefficient field; the generator writes integers only, so
  reading stays a superset of what it writes.
- `text`, the crate's own syntax: a `# vars:` line, a `# modulus:` or
  `# coefficients: rationals` line, then one polynomial per line in the
  syntax `parse_polynomial` reads and `Display` writes. It is the default
  output format and it round trips through `sylv` itself.

### 8.4 Certificates, verification output, and exit codes

`sylv gb --certified --certificate out.cert` runs the certified path and
writes the accepted certificate bytes, both formats, with no encoding
step. Over `Q` there is no certified path, so the flag is a usage error
naming the domain (exit 2).

`sylv verify out.cert` runs `verify::verify_with_limits` on bytes it did
not produce and never loads an engine. `--max-bytes` applies before the
file is read: the reader takes at most `max + 1` bytes, from a file or
from standard input, and rejects anything longer without allocating past
the cap. What the command can print is limited by what a certificate
carries: `VerifiedGb` gives the modulus, the variable count, the input,
and the basis, and no variable names, because names are in neither
schema. So `verify` prints the basis under the synthetic names
`x1 .. xn`, and says on standard error that the names are synthetic. The
notice goes to standard error and not into the output, because a header
line would corrupt `ms` and `syl` output for the next tool that reads it.
`--quiet` prints nothing and reports through the exit code alone.

| code | meaning |
|---|---|
| 0 | success; for `normal-form` a remainder was printed, for `verify` the bytes were accepted |
| 1 | the verifier rejected untrusted certificate bytes |
| 2 | usage, input, or an operation the domain does not offer; a supplied basis that fails a shape or Gröbner check (`BasisError::{RingMismatch, ZeroPolynomial, NotMonic, NotSorted, NotInterreduced, NotGroebner}`), since that is bad input |
| 3 | a resource or structural limit: any `ComputeError`, `NormalFormError::{ExponentLimit, Timeout, MemoryLimitExceeded}`, `BasisError::{ExponentLimit, Timeout, MemoryLimitExceeded}`, `HilbertError`, a verifier cap or deadline, writer exhaustion |
| 4 | an internal defect: the crate wrote a certificate its own verifier rejected, an emitter fault, or an input mismatch on the certified path |
| 5 | an input or output failure: a file that cannot be read or written, or a broken pipe |

Code 1 is reserved for untrusted bytes. A freshly written certificate
that fails its own verifier is a defect in this crate, not an invalid
input, and it exits 4 so a script cannot confuse the two.

`--report` writes the counters to standard error, one per line, so the
basis on standard output stays parseable.

## 9. Python bindings

### 9.1 Crate, package, and layout

`crates/sylvester-py`, `publish = false`, `crate-type = ["cdylib"]`,
depending on `sylvester` and on pyo3 with features `extension-module`,
`abi3-py310`, `num-bigint`, and `num-rational`. The pyo3 requirement is
`"=0.23.5"`, an exact pin rather than the caret default, because the GIL
release call this document specifies is `Python::allow_threads`, which a
later pyo3 renames. Changing the pin is a decision with binding-wide
consequences (section 13, item 8), not a version bump.

`pyproject.toml` uses the maturin backend, with the build requirement
`maturin>=1.11,<2`, which is the sibling project's, the distribution name
`sylvester`, the crate's version,
`requires-python = ">=3.10"` to match `abi3-py310`, the license fields
and classifiers of the sibling `auslander-py`, and the mixed layout:

```toml
[tool.maturin]
python-source = "python"
module-name = "sylvester._sylvester"
```

`python/sylvester/__init__.py` re-exports from `sylvester._sylvester`
next to `__init__.pyi` and `py.typed`. The cost is one shim that lists
every export, and a test asserts the shim and the extension export the
same names.

Build and test, added to `AGENTS.md`:

```
cd crates/sylvester-py && maturin develop --release
python -m pytest tests
```

Release tags build wheels and a source distribution for the GitHub release.
The Rust crate has `publish = false`. The Python metadata carries the private
classifier, and the release process does not upload to PyPI.

### 9.2 The surface

The domain is known only at run time in Python, so each Rust pair collapses to
one Python class holding an enum, the same boundary cost as the CLI.

| Rust | Python |
|---|---|
| `PolynomialRing<D>` | `PolynomialRing.prime_field(p, vars)`, `PolynomialRing.rationals(vars)`, `ring.modulus` (`None` over `Q`), `ring.variables` |
| construction | `ring.parse(text)`, `ring.polynomial(terms)`, `ring.ideal(generators)` |
| `Polynomial<D>` | `Polynomial`: `terms()`, `degree()`, `__str__`, `__eq__`, `__hash__` |
| `Ideal<D>` | `Ideal`: `groebner_basis(**opts)`, `groebner_basis_certified(**opts)`, `has_homogeneous_generators()`, `generators` |
| `GroebnerBasis<D>` | `GroebnerBasis`: `__len__`, `__getitem__`, `__iter__`, `normal_form`, `contains`, `hilbert_series`, `krull_dimension`, `is_homogeneous`, `lift`, and the static `from_polynomials` |
| `HilbertSeries` | `HilbertSeries`: `numerator` (`list[int]`), `denominator_power`, `dimension()`, `multiplicity()`, `coefficient(d)` |
| `CertifiedGroebnerBasis` | `.basis`, `.certificate` (`bytes`) |
| `verify::verify` | `sylvester.verify(data) -> VerifiedGroebnerBasis` |

`Polynomial.terms()` returns `list[tuple[Coefficient, list[int]]]`,
largest monomial first, where `Coefficient` is `int` over `F_p` and
`fractions.Fraction` over `Q`. `ring.polynomial(terms)` accepts the same
shape and also `int`, `Fraction`, and `(numerator, denominator)` pairs
for the coefficient. Exponents cross as `list[int]`. numpy is not a
dependency: the arrays are short, ragged, and not numeric matrices.

`ComputeOptions` does not cross. The compute methods take keywords, all
defaulting to `None`, so an omitted argument is distinguishable from an
explicit one:
`backend=None`, `timeout=None` (seconds, float), `memory_limit=None`,
`threads=None`, `stop=None` (`"unchanged"` or `"contains_input"`), and
`extra_primes=None`. `stop=None` selects `ContainsInput`; the other omitted
values use the Rust defaults. A rational keyword passed to a prime-field
ideal is a `ValueError`, not a silent no-op. A wrong keyword value is a
`ValueError` naming the accepted set.

`GroebnerBasis` is sequence-like, so `f in G` would read as "is `f` one
of the listed polynomials". Ideal membership is `G.contains(f)`, and
`__contains__` is not implemented, because a wrong answer there would be
silent.

Semantics the classes need, none of which has an obvious default: every
wrapper is immutable; `__getitem__` accepts negative indices and slices,
and a slice returns a `list`; `__iter__` yields owned clones;
`__hash__` on `Polynomial` and `PolynomialRing` hashes content, so two
equal rings built separately hash equally; `lift` returns a dict of the
`ModularLift` fields with `established` as `"unchanged"` or
`"contains_input"`, returns `None` for a rational basis that came from
`from_polynomials`, and raises `ValueError` on a prime-field basis, which
is the runtime image of a method that does not exist in Rust.

`timeout` is validated before it becomes a `Duration`: a negative, an
infinite, or a NaN float is a `ValueError` naming the argument, not a
saturating conversion.

### 9.3 Errors and the GIL

Bad input is `ValueError`, an exhausted budget is a `BudgetExhausted`
subclass of `RuntimeError`, and a defect inside the crate is a distinct
`RuntimeError` subclass. The mapping follows `auslander-py`, and it
matches on wrapped errors rather than on the wrapper: `CertifyError`
carries `Engine(ComputeError)`, `WriterExhausted`, and
`VerifierExhausted(VerifyError)`, and each inner value decides the class.

| Rust | Python |
|---|---|
| `RingError`, `ParseError`, `NormalFormError::RingMismatch`, `BasisError::{RingMismatch, ZeroPolynomial, NotMonic, NotSorted, NotInterreduced, NotGroebner}` | `sylvester.RingError`, `sylvester.ParseError`, `sylvester.BasisError`, all `ValueError`, with `index` or `left` and `right` attached |
| `Timeout` and `MemoryLimitExceeded`, wherever they occur: `ComputeError`, `NormalFormError`, `BasisError`, `HilbertError`, `CertifyError::Engine`, `CertifyError::WriterExhausted`, and the deadline variants of `VerifyError` | `sylvester.Timeout`, `sylvester.MemoryLimitExceeded`, both `sylvester.BudgetExhausted(RuntimeError)` |
| `ComputeError::{DegreeLimit, ExponentLimit, TableFull, PrimesExhausted}`, `NormalFormError::ExponentLimit`, `BasisError::ExponentLimit`, `CertifyError::CapExceeded`, the cap variants of `VerifyError` | `sylvester.LimitExceeded(RuntimeError)`, with `limit` attached where the variant carries one |
| `CertifyError::{Rejected, Emitter, InputMismatch}` | `sylvester.InternalDefect(RuntimeError)`: the crate contradicted itself |
| the rejection variants of `VerifyError` from `sylvester.verify` | `sylvester.CertificateInvalid(ValueError)`: untrusted bytes are input |
| asking `Ideal(rationals).groebner_basis_certified` | `ValueError` naming the domain |

`sylvester.verify` returns `VerifiedGroebnerBasis`, a class distinct from
`GroebnerBasis`, with `input`, `basis`, `modulus`, and `nvars`. Its
polynomials use the synthetic names `x1 .. xn`, for the reason section
8.4 gives, and the class documents it.

Every compute call, every certified call, every `verify`, every normal
form, every basis check, and every Hilbert series releases the GIL around
the Rust work with `Python::allow_threads`. The rule around it: convert
and clone every input into an owned Rust value first, including copying
the bytes of `verify` out of the Python buffer, then release the GIL,
then convert results and errors into Python objects after it is back. A
`PyRef` or a borrowed buffer cannot cross into the released region. The
The existing invariants make the release sound: `PolynomialRing`, `Ideal`,
`Polynomial`, and `GroebnerBasis` are `Send + Sync`, the compute methods
take `&self` and return owned values, and no error borrows from the ring.
A static assertion in the binding crate fails the build if any of those
types stops being `Send + Sync`. A pytest case starts a computation in
one thread and passes a `threading.Event` back and forth with the main
thread, so it proves progress by synchronization and not by a wall clock
threshold.

## 10. Workspace, module tree, dependencies

The root package stays the library and becomes the workspace root. No
file under `src/` moves.

One exclusion is required, not optional. `benchmarks/gb-comparison/runner`
is a package (`sylv-runner`) with its own manifest and its own
`Cargo.lock`, and it depends on the library by path. Once the root is a
workspace, Cargo treats it as a member of that workspace unless the root
manifest says otherwise, and it loses its independent lock file, which
the benchmark protocol relies on. The root manifest therefore lists it
under `[workspace] exclude`.

```
Cargo.toml                     # [package] sylvester + [workspace] members and exclude
src/
  lib.rs  ring.rs  poly.rs  ideal.rs  certificate.rs
  ring/field.rs                # Felt public, PrimeOps
  ring/rational.rs             # new: RationalOps, clear denominators, mod-p image
  normal_form.rs               # new
  hilbert.rs                   # new
  compute/mod.rs  classic.rs  interreduce.rs  signature.rs  f4/*
  cert/*  verify/*             # unchanged apart from prime-field signatures
```

The modular driver is five files, one job each:

```
  compute/modular/mod.rs          # the driver of section 3.1
  compute/modular/primes.rs       # the sequence of section 3.2
  compute/modular/crt.rs          # the accumulator of section 3.4
  compute/modular/reconstruct.rs  # the lift of section 3.6
  compute/modular/check.rs        # T1 and T2 of section 3.7
crates/sylvester-cli/src/{main.rs, format/{ms,syl,text}.rs}
crates/sylvester-cli/tests/cli.rs
crates/sylvester-py/{src/lib.rs, python/sylvester/*, tests/*}
```

New dependencies of the library: `num-bigint`, `num-rational`,
`num-integer`, and `num-traits`. New dependencies of the new crates:
`clap` and `pyo3`. The verifier gains no dependency and imports none of
the new modules; `tests/isolation.rs` keeps checking that, and its
forbidden list grows by `crate::compute::modular`, `crate::normal_form`,
and `crate::hilbert`.

## 11. Test and benchmark plan

### 11.1 Rationals

- **Reconstruction, by contract.** Property tests assert the contract of
  section 3.6, not an intent: any returned pair is reduced, inside both
  bounds, and congruent; when the true fraction is inside the bounds, it
  is the one returned; when it is outside, a different in-bounds pair is
  a correct result, and only an out-of-contract pair fails. Named cases:
  `r = 0`, `b = 1`, negative numerators, and `M` on both sides of
  `2 A B`.
- **CRT.** Property test against a direct big-integer solve, including
  monomials that appear only in later runs and elements whose support
  changes between runs.
- **The vote.** `F = { x + y, x + (p + 1) y }` of section 3.3, with the
  prime `p` injected through a crate-private test hook, asserts that the
  category with leading ideal `(x)` loses once the later primes arrive,
  that the accumulator is rebuilt and the confirmation counter resets
  when prevalence changes, and that the skip and discard counters read
  what happened. A second case, built from a leading coefficient
  divisible by the injected prime, exercises the pre-run skip.
- **Exact oracle.** A Buchberger implementation over `BigRational` in
  `tests/rational.rs`, self-contained, run against the multimodular path
  on random small systems and on the small members of the fixed families,
  sized so at least one case makes the vote and the support union do real
  work.
- **Fresh-prime differential.** On the larger families, the lifted basis
  is reduced modulo a prime the run did not use and compared with a
  prime-field run over that prime. This repeats the production heuristic,
  so it is a regression test and not an oracle, and the plan says so.
- **Failure policy.** A prime run driven into `Timeout` and one driven
  into `ExponentLimit` both surface as themselves, not as a discarded
  prime. A run driven into the per-run memory share succeeds on the
  single retry of section 3.8.
- **Determinism.** For runs that succeed, byte-identical rational bases
  and identical `ModularLift` counters at 1, 2, 8, and 16 threads, and
  across two fresh processes.

### 11.2 Normal form, membership, Hilbert

- Normal form is idempotent, is zero exactly on ideal members, and agrees
  with the certified prime-field path where both apply.
- The `x^65535 * y` case of section 5 returns
  `NormalFormError::ExponentLimit` instead of a panic. The overflow is in
  the tail multiply alone: the quotient monomial comes from an exponent
  subtraction after divisibility holds, so it cannot overflow.
  `Monomial::checked_mul` is tested directly for its own overflow.
- `from_polynomials` rejects a non-monic list, an unsorted list, a
  non-interreduced list, and a list that fails T2, each with its own
  variant, and accepts a basis the certified path produced.
- Membership properties run against a basis an oracle or a certificate
  has confirmed, so a defect in the basis cannot pass as a property of
  membership. For random cofactors `h_i`, `sum h_i g_i` is a member; a
  random low-degree polynomial is compared with the oracle rather than
  assumed to be a non-member.
- Hilbert numerators on hand-checked monomial ideals, including
  `L = (x, y z)`, which loops under a pivot rule without the restriction
  of section 6.2, and against the count of standard monomials up to
  degree 8 for random monomial ideals in 3 and 4 variables.
- Dimension against Singular's `dim(std(I))` on the harness families,
  committed as fixtures, plus the unit ideal, the zero ideal, and a
  principal ideal.
- One homogeneous case where the series is the ideal's own Hilbert
  series, and one inhomogeneous case where it is not, asserted through
  the first-difference relation of section 6.1 rather than by equality.

### 11.3 CLI and Python

- Golden tests under `crates/sylvester-cli/tests/cli.rs`, so the binary
  is reachable through `CARGO_BIN_EXE_sylv` inside its own package. For
  each subcommand and each format: a committed input, a committed
  expected standard output, and the expected exit code. A `--certified`
  run writes a certificate that a second invocation of `sylv verify`
  accepts; a flipped byte makes that invocation exit 1. A rational
  `--certified` run exits 2. An unwritable output path exits 5.
- pytest: one module per surface area (ring and polynomial, ideal and
  basis, rationals, certificates, Hilbert), the error mapping of section
  9.3 asserted row by row including the wrapped variants,
  `sylvester.__version__` equal to the crate version, the stub and
  extension export lists equal, and the GIL case of section 9.3.

### 11.4 Benchmarks

`benchmarks/gb-comparison` gains rational cells for the small members of
the four families over `Q`, under the protocol of the prime-field cells:
full canonical output comparison, pinned cores, recorded repetitions,
recorded versions and full command lines.

One precondition comes first, and the gate of section 1.2 depends on it.
Each reference tool needs an invocation that emits a lifted rational
basis, and that has to be established rather than assumed. The harness
uses msolve today with `-g 2` over a prime field, and whether the
installed msolve 0.10.1 emits a rational basis is an open question
(section 13, item 9). Singular over `Q`, through both `std` and
`modStd`, and Groebner.jl over `Rational{BigInt}` are the references
known to produce one. If msolve cannot, it is dropped from the rational
cells and the record says so.

`sylv-runner` gains a rational mode. Every number in `README.md`,
`CHANGELOG.md`, or `KNOWN_ISSUES.md` comes from the record, copied,
never typed.

## 12. Work breakdown

An earlier draft ran four chunks in parallel behind one coefficient
chunk. The file dependencies do not allow it, the same way they did not
before. The coefficient migration alone is too large for one reviewable
commit, so it lands in three stages, each of which compiles and keeps
`cargo test` green.

| chunk | owns | depends on |
|---|---|---|
| A1. Sealed domain over the prime field | `src/ring.rs`, `src/ring/field.rs`, `src/poly.rs`, `src/ideal.rs`, `src/certificate.rs`, `src/lib.rs`, `src/compute/**`, `src/cert/**`, `tests/**`, `benches/**`, `benchmarks/gb-comparison/runner/{Cargo.toml, Cargo.lock, src/main.rs}`, root `Cargo.toml` and `Cargo.lock` | none |
| A2. The rational domain | `src/ring/rational.rs`, `src/ring.rs`, `src/poly.rs`, `src/ideal.rs`, `src/lib.rs`, `tests/roundtrip.rs`, `tests/isolation.rs` | A1 |
| A3. Limits, budget, and the meters | `src/compute/mod.rs`, `src/compute/f4/{mod,kernel}.rs`, `src/compute/{classic,interreduce}.rs`, `src/cert/mod.rs`, `src/lib.rs` | A2 |
| B. Normal form and the checked basis | `src/normal_form.rs`, `src/poly.rs` (checked multiply), `src/ideal.rs`, `src/lib.rs`, `tests/{normal_form,isolation}.rs` | A3 |
| C. Multimodular engine | `src/compute/modular/**`, `src/compute/mod.rs`, `src/ideal.rs`, `src/lib.rs`, `tests/{rational,isolation}.rs` | B |
| D. Hilbert | `src/hilbert.rs`, `src/ideal.rs`, `src/lib.rs`, `tests/{hilbert,isolation}.rs` | C |
| E. CLI | `crates/sylvester-cli/**`, root `Cargo.toml` and `Cargo.lock` | D |
| F. Python | `crates/sylvester-py/**`, root `Cargo.toml` and `Cargo.lock` | E |
| G. Harness | `benchmarks/gb-comparison/**` | C |
| H. Docs | `README.md`, `CHANGELOG.md`, `KNOWN_ISSUES.md`, `AGENTS.md`, `PLAN.md` status lines | all |

Three boundaries inside the A stages are set here, because each one is a
place where a naive split would not compile.

- A1 introduces `Coefficient` and the `num-*` dependencies, not A2.
  `DomainOps::convert` takes a `&Coefficient`, so the type has to exist
  with the trait. A1 also changes the public term accessor, so the tests,
  the benches, and the runner move with it, including the runner's own
  `Cargo.lock`, which the new dependencies change, and its variable
  naming, which goes from `x0 .. x(n-1)` to `x1 .. xn` in the same
  commit as the rest of its migration. A1 adds the workspace `exclude`
  entry for that package (section 10).
- A2 introduces `RationalMeta`, `ModularLift`, and `Established` as
  types, because `Rationals::BasisMeta` names `RationalMeta` the moment
  the domain exists. Chunk C is what first fills a `ModularLift` in.
  Until then the only rational basis a caller can build is a checked one
  from chunk B.
- A3 owns `src/compute/f4/kernel.rs`, because cancellation reaches the
  per-worker deadlines that file builds.

A1 is the sealed traits, the public operation types, and the generic
scaffold with `PrimeField` as the only domain, adapting every engine,
certificate, bench, and test signature. A2 adds `Rationals`, rational
arithmetic, the parser and `Display` paths, and the basis metadata types.
A3 extracts `Budget`, promotes `ComputeLimits`, adds cancellation, and
fixes the three memory meters. Nothing rational computes until C.

A1 through D are sequential, because each owns `src/ideal.rs` or
`src/lib.rs` in turn. E and F both touch the root manifest's workspace
members and the lock file, so F follows E rather than running beside it;
G shares no file with either and runs as soon as C lands. Every chunk
that adds a module also updates `tests/isolation.rs`. The existing golden
certificates must still verify byte for byte at the end of every stage.

Two rules hold unchanged. A chunk lands only with its tests. A
commit that claims a speed win carries a before and after number from the
harness record.

## 13. Implementation decisions

Each item carries the design's proposal. The owner may override any of
them. Items 1 to 5 change types or control flow, so they are answered
before chunk A1 starts.

1. **Coefficient parameter against a runtime tag.** Section 2.1 proposes
   the sealed type parameter `PolynomialRing<D>`, which makes "no
   certification over Q" a compile error and keeps the engines on
   concrete types, at the cost of an enum boundary in the CLI and the
   Python bindings. The alternative is one non-generic type with a domain
   tag and an enum of term storages, which keeps those two boundaries
   simple and moves the domain check to run time everywhere else. A
   related sub-decision: whether the public types get a default parameter
   (`PolynomialRing<D = PrimeField>`), which shortens prime-field
   signatures and hides the domain in the docs.
2. **Deterministic or seeded primes.** Section 3.2 proposes the fixed
   descending sequence, which makes the rational basis and its counters
   reproducible and gives an adversary the confirming primes in advance.
   The alternative is a seeded sequence with the seed in `ModularLift`,
   which supports a probability statement and gives up reproducibility
   across seeds.
3. **The default stopping rule.** Rust uses `RationalStop::Unchanged {
   extra: 2 }`. The CLI and Python binding default to `ContainsInput`,
   which adds T1 and T2 and still does not exclude `{1}`.
4. **How much a rational run may retain.** Section 3.3 keeps every
   completed run so the accumulator can be rebuilt when prevalence
   changes, and the retained bytes are the sum of every modular output.
   The alternative is a bounded pool with recomputation when the pool
   misses.
5. **Per-prime failure policy.** Section 3.5 stops the whole run on
   `DegreeLimit`, `ExponentLimit`, and `TableFull`. The alternative is a
   bounded retry with the retries counted.
6. **Hilbert resource limits.** Section 6.2 gives `hilbert_series` a
   `Budget` and an error type, because no proved bound on the recursion
   is offered. The alternative is an infallible method and a promise the
   design cannot back.
7. **The dimension of the zero quotient.** Proposed `Option<usize>` with
   `None` for the unit ideal. The alternative is a two-variant enum,
   which reads better at the call site and adds a public type.
8. **The pyo3 pin.** Proposed `=0.23.5`, matching `auslander-py`, with
   `Python::allow_threads`. A newer pin changes the GIL call and parts of
   the `Bound` API through the whole binding.
9. **The rational reference tool.** Section 11.4 makes a working
   rational invocation a precondition of the gate. If msolve 0.10.1
   cannot emit a lifted rational basis, confirm that Singular and
   Groebner.jl are enough and that msolve is dropped from the rational
   cells.
10. **New dependencies.** `num-bigint`, `num-rational`, `num-integer`,
    and `num-traits` in the library; `clap` in the CLI; `pyo3` in the
    bindings. The verifier gains none. Confirm.
11. **CLI input formats.** Proposed three: `ms`, `syl`, and `text`.
    `syl` carries no ring and exists for the harness inputs; dropping it
    removes a format and a flag interaction. The `syl` variable naming
    change (`x0` to `x1`) in the runner rides on this.
12. **The standalone verifier's output.** Section 8.4 has `sylv verify`
    print the basis under synthetic names, since certificates carry no
    variable names. The alternative is to report acceptance only.

## Review record

Codex (gpt-5.6-sol, xhigh) reviewed this design as an adversary over six
rounds and endorsed it at the end of the sixth. The record is the
disposition, not an endorsement of any earlier draft.

- Round 1 found the worst defect: the claim that a modular leading ideal
  always contains the rational one, used to convict primes as "proven
  unlucky". The counterexample `{ x + y, x + (p + 1) y }` is now in
  section 3.3. It also killed the false Hilbert identity
  `L = (L : x_j) ∩ (L + (x_j))`, the affine Hilbert function statement,
  the missing `gcd(a, b) = 1`, and the runtime domain tag, which the
  sealed type parameter replaced.
- Round 2 showed the Hilbert recursion could not terminate under the
  stated pivot rule (`L = (x, y z)` returns its parent), that
  `sub_scaled` multiplies unchecked so `Monomial::checked_mul` alone does
  not fix normal form, that `GroebnerBasis::new` is crate-private so the
  CLI could not act on a supplied basis, that speculative runs need a
  cancellation flag, and that `ComputeLimits` and `WriterBudget` collide
  with names `src/compute/mod.rs` already imports.
- Round 3 found that a checked rational basis has no honest
  `ModularLift` (now `RationalMeta::Checked`), that
  `ComputeError::Cancelled` cannot be internal-only on a public
  exhaustive enum (now a crate-private `RunError`), that `ContainsInput`
  could rerun T1 and T2 on a candidate it had already rejected, and that
  the A1 to A3 split did not compile at its boundaries.
- Round 4 corrected the order of normalization against the invertibility
  test (`3/3` over `F_3`), a generic constructor with no way to build
  `D::BasisMeta`, a cancellation flag that would have cancelled the work
  after it, an undefined memory reservation, and the workspace capture of
  the benchmark runner and its lock file.
- Round 5 closed the last three implementation defects: the deadline
  covers the driver's own loops and not only the prime runs, the CLI
  resolves one ring per command so a `syl` basis can meet a `text`
  polynomial, and three remaining claims were cut back to what is known.
- Round 6 was the endorsement and a style pass. It added the pivot tie
  break of section 6.2, limited the retention claim to consumed runs,
  fixed the scope of the 31-bit bound, and cut the remaining filler.
