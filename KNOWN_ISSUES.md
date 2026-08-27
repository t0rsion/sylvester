# Known issues

This file records correctness defects, their repairs, and the limits of the
evidence. The first commit in `main-local` preserves the incorrect extraction.
The current tree contains the repairs.

## Status

Both prime-field engines are empirically checked. Neither engine is proven
correct. An accepted certificate establishes the result through an independent
verifier. Neither verifier is machine-checked.

The F4 engine replaces the matrix engine in 0.2.0. The classic F5 engine remains
an independent oracle. The matrix defects below remain part of the record.

## Counterexamples

`tests/known_defects.rs` asserts both expected bases. Its checker is independent
of the engines. An output-only check would accept `{1}` and would miss the first
failure.

### Classic F5 returned a basis of the wrong ideal

Let the field be \(\mathbb F_2\), with variables \(x > y\) in grevlex:

\[
f_1=x^3+y^3,\qquad
f_2=x+x^2+x^2y,\qquad
f_3=y+x^2.
\]

The archived engine returned

\[
\{y^2,\ x+y\}.
\]

This set is a reduced Gröbner basis, but not for the input ideal. The normal
form of \(f_2\) is \(y\). The correct reduced basis is \(\{x,y\}\).

### Matrix F5 returned a non-Gröbner basis

Let the field be \(\mathbb F_3\), with variables \(x > y > z\) in grevlex:

\[
f_1=x^2+y^2,\qquad
f_2=xz+xy,\qquad
f_3=y+xy.
\]

The archived engine returned

\[
\{x^2+2yz,\ xy+y,\ y^2+yz,\ xz+2y\}.
\]

The S-polynomial of the first two elements has normal form \(y+yz^2\). The
correct reduced basis also contains \(yz^2+y\).

The archived internal check only proved that the output was reduced for the
ideal it generated. It did not prove equality with the input ideal.

## Root causes and repairs

### Inflated signatures from non-regular pairs

Both old engines constructed a pair when its component signatures were equal.
The common module monomial is not a trusted signature for the S-polynomial. If
the polynomial reduced to zero, the engine recorded an inflated syzygy
signature. That rule could then discard a necessary pair.

The repair rejects non-regular pairs at construction. Only regular pairs enter
the signature computation.

### Signature-blocked rows were dropped

The matrix engine dropped a surviving pivot row when an existing leading
monomial divided its leading monomial. Signature-safe elimination may forbid
that reduction. Such a row can therefore contain required information.

The repair completed every legal regular reduction. It dropped a result only
when an existing basis element divided both its signature and leading
monomial. The matrix engine was later removed.

### Duplicate insertions could prevent termination

The classic engine could repeatedly insert one signature and leading monomial.
Equal-signature reduction blocked each copy, while each copy created more
pairs.

The repair applies the same divisibility test before insertion. Fix one input
index. The accepted pairs

\[
(\operatorname{sig}(g),\operatorname{lm}(g))\in\mathbb N^{2n}
\]

form a Dickson-bad sequence: no earlier pair divides a later pair component by
component. Dickson's lemma makes this sequence finite. There are finitely many
input indices, so the basis and pair queue are finite. This is a paper proof,
not a machine proof.

### Rewriter pair deletion was order-dependent

The matrix engine processed pairs by least-common-multiple degree, not by
signature. Deleting a pair because a later element rewrote its signature could
remove a necessary S-polynomial.

The repair substituted the canonical rewriter's multiple. The replacement had
the same signature monomial. The matrix engine was later removed.

### The product criterion lacked a signature proof

The matrix engine used Buchberger's product criterion inside a signature
algorithm. No proof connected that deletion to the engine's signature coverage
argument. The repair removed the criterion. The classic engine never used it.

## F4 termination

The F4 engine accepts a nonunit candidate only if its leading monomial is
outside the current leading-monomial ideal. Each insertion therefore gives a
strict chain

\[
L_0\subsetneq L_1\subsetneq L_2\subsetneq\cdots
\]

of monomial ideals in \(\mathbb F_p[x_1,\ldots,x_n]\). The ascending chain
condition makes the number of insertions finite. Each insertion creates
finitely many pairs. A batch can restart each lane at most twice. Thus the F4
run terminates unless a typed resource or representation limit stops it.

This argument proves termination. It does not prove partial correctness.

## Resource defects

### Unbounded CSR memory

The old CSR matrix path allocated until the operating system killed the
process on `cyclic-7`. The plain sparse path did not show the same growth. The
repair removed the CSR path. Version 0.2.0 removes the whole matrix engine.
The crate has no feature flags.

The July comparison record preserves the measurements:
`benchmarks/gb-comparison/records/2026-07-25-archive/REPORT.md`.

### Classic stack exhaustion

The old classic engine exhausted the default stack on `eco-9`. The current
engine finishes or returns `ComputeError::Timeout` on a 1 MiB worker stack.
`tests/stack_depth.rs` asserts this property.

No isolated patch explains the repair. The API rebuild replaced the affected
control flow before a debugger identified one routine.

## Other fixed defects

- `Felt` is private and reduced at every construction site.
- Final interreduction runs under the computation deadline and memory budget.
- A classic pair whose least common multiple exceeds the exponent width returns
  `ComputeError::DegreeLimit`.
- An F4 monomial whose exponent exceeds the width returns
  `ComputeError::ExponentLimit`.
- An exhausted F4 monomial table returns `ComputeError::TableFull`.
- The crate contains no unsafe code.

## Evidence

The default differential suite compares both engines with a criterion-free
Buchberger oracle. It checks:

- 100 deterministic systems over each of \(\mathbb F_2\), \(\mathbb F_3\),
  and \(\mathbb F_5\);
- 25 larger systems over \(\mathbb F_7\);
- 25 systems over \(\mathbb F_{32003}\);
- both counterexamples in every generator order;
- the mutation-derived rewriter regression.

The ignored release suite adds all \(55^3=166{,}375\) ordered triples of
small monomials and binomials over \(\mathbb F_2\), plus 1,000 degree-reversal
systems over \(\mathbb F_3\). `tests/f4_engine.rs` also compares F4 with
classic on the `eco-8` family under the benchmark modulus.

These tests compare complete reduced bases, including coefficients. They
are evidence, not a proof.

The comparison harness checks complete bases against Singular, msolve, and
Groebner.jl. Its current record is `benchmarks/gb-comparison/REPORT.md`. Frozen
prior records stay under `benchmarks/gb-comparison/records/`.

## Certification boundary

Classic writes `sylv-gb-cert-v1`. F4 writes `sylv-gb-cert-v2`. Each verifier
implements its own arithmetic, monomials, polynomials, decoder, and checks.
Engine code does not enter either verifier tree.

Acceptance establishes input membership, basis membership, Buchberger's
criterion, and reduced canonical form. It therefore establishes the unique
reduced Gröbner basis of the input ideal under `grevlex-v1`.

The verifier is trusted code. Its tests include malformed data, cap boundaries,
deadline boundaries, truncation, and byte mutation. It is not machine-checked.

## Open limits

- Engine correctness is empirical and theory-backed, not proven.
- Verifier correctness is not machine-checked.
- Memory budgets charge explicit data structures. They do not cap allocator
  metadata, fragmentation, thread stacks, or process RSS.
- Verification caps bound work and memory estimates. A deadline bounds wall
  time with fixed polling intervals.
- F4 certification costs more than raw F4 computation because the verifier
  uses isolated arithmetic. The benchmark record measures that cost.
- Performance claims apply only to the recorded inputs, tools, limits, and
  machine.
- Input may be inhomogeneous. The differential suite includes inhomogeneous
  systems. This remains empirical coverage.

## Extraction note

The archived artifact is source-diff comparable with the old workspace. Lock
file updates make it not byte-identical.
