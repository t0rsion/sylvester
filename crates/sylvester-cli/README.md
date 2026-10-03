# sylvester-cli

`sylv`, the command line interface of the `sylvester` crate. It computes
reduced Gröbner bases under grevlex, over a prime field or over the
rational numbers. It also certifies a prime-field basis, checks certificate
bytes with the independent verifier, reduces a polynomial modulo a basis,
inspects finite quotient algebras, and reads the Hilbert series and the Krull
dimension of a basis.

Install the binary from crates.io:

```
cargo install sylvester-cli
```

Build it from the workspace root during development:

```
cargo build --release -p sylvester-cli
```

## Commands

```
sylv gb cyclic-7.ms
sylv gb system.syl --modulus 1073741827 --backend classic --threads 8
sylv gb system.text --rationals --stop contains-input --extra-primes 3
sylv certify system.ms --certificate system.cert
sylv verify system.cert
sylv normal-form --basis basis.text --poly "x^2*y - 1"
sylv member --basis basis.text --poly-file f.text
sylv hilbert system.ms
sylv dim system.ms --from-basis
sylv quotient system.text
sylv quotient system.text --matrix x --characteristic x --minimal x
```

`gb` takes `--certified` to run the certified path, and `--certificate
PATH` to write the accepted bytes. `certify` is `gb` with that path always
on. `nf` is the short name of `normal-form`, and `basis` of `gb`.
Over `Q`, the default stopping rule is `contains-input`. Pass `--stop
unchanged` only when the weaker claim is sufficient.

Over `Q`, `gb --check-equality` runs the bounded exact check after the
multimodular computation. It reports `# equality_check: passed` in text
output and records `equals_input` in JSON. The check is separate from the
multimodular `Established` claim. A quotient from `--from-basis` can use
the same option when its JSON record carries the original input.

`normal-form` and `member` take the basis and the polynomial from separate
sources. The basis goes through the checked constructor, so a list that is
not a reduced Gröbner basis is rejected with the check that failed. At most
one source reads standard input.

`normal-form --quotients` prints one labeled quotient for every basis element,
then the remainder. Parenthesized expressions, unary signs, and exact
fractions are accepted by `--poly` and by expression fields in input formats.

`hilbert` computes a basis from generators. With `--from-basis` it reads a
basis instead, and then `--backend`, `--threads`, `--stop`, and
`--extra-primes` describe a computation that does not happen, so passing
any of them is a usage error. `dim` prints the `dimension:` line of the
same output.

`quotient` computes a basis, or checks one with `--from-basis`, then builds
its finite quotient algebra. With no operation flag it prints the dimension
and standard monomials. `--standard-monomials` and `--dimension` select those
lines independently. `--matrix`, `--characteristic`, and `--minimal` apply a
polynomial to the quotient. A positive-dimensional quotient returns exit code
3 with a typed nonfinite error. The output uses labeled text lines.

`--timeout SECS` covers the whole command, not one library call. The binary
takes an instant when it starts and hands the remaining time to every call
it makes. `--memory BYTES` limits estimated live working data for the command,
including retained source and parsed input buffers. Bare decimal integers
remain valid; binary suffixes such as `512MiB` and decimal suffixes such as
`2MB` are also accepted. `--progress` writes named phases to standard error
without changing standard output. `--no-progress` suppresses them.

## Formats

`--in-format` and `--out-format` name one of four formats. A file with no
format option takes the format of its extension: `ms`, `syl`, `sylq`, `txt`,
`text`, or `json`. Standard input takes `text`. Any other extension is a
usage error rather than a guess.

- `ms`, the msolve input format: a line of comma-separated variable names,
  a line with the characteristic (`0` for `Q`), then comma-separated
  polynomials.
- `syl`, the benchmark format under `benchmarks/gb-comparison/inputs`: a
  line with the variable count, then one polynomial per line as terms
  separated by `;`, each term a coefficient and one exponent per variable,
  separated by `,`. It carries no ring, and a command that reads one names
  the variables `x1 .. xn` unless another source names them. A coefficient
  may be a fraction, `-1/3`. The `.sylq` files of the harness are the same
  grammar.
- `text`, the crate's own syntax: a `# vars:` line, a `# modulus:` or
  `# coefficients: rationals` line, then one polynomial per line in the
  syntax `PolynomialRing::parse_polynomial` reads and `Display` writes. It
  is the default output format, and it reads back through `sylv`.
- `json`, a `sylv-result-v1` computation record. It stores the ring, original
  input, basis, an untrusted provenance claim, and optional certificate bytes.
  `gb` writes this format with `--out-format json`; loading it does not trust
  the claim or certificate. `verify` rechecks an attached prime-field
  certificate and matches its full input and basis before accepting it.

`hilbert` and `dim` write their own schema: a `series:` line with the
numerator coefficients low degree first, a `denominator_power:` line, a
`dimension:` line, and a `multiplicity:` line. The unit ideal has
`dimension: none`.

## The ring of a command

The ring is resolved once per command, over every source it reads. Three
rules, in order:

1. Every source that names a part of the ring must name the same one. A
   disagreement is a usage error.
2. If a source names the coefficient domain, that domain is the command's,
   and `--modulus` and `--rationals` do not apply.
3. If no source names one, exactly one of `--modulus` and `--rationals` is
   required.

A source that names no variables takes the names of a source that does, so
a `syl` basis and a `text` polynomial work together and land in one ring.
When no source names any, the variables are `x1 .. xn`.

## Certificates

`sylv certify` runs the backend, writes a certificate for the run, and
returns the basis the verifier accepted, decoded from those bytes. The
format follows the backend: the classic backend writes `sylv-gb-cert-v1`
and F4 writes `sylv-gb-cert-v2`. There is no certified path over `Q`, so
`--certified` there is a usage error.

`sylv verify` reads bytes it did not produce and loads no engine.
`--max-bytes` applies before the file is read. A certificate carries no
variable names, so the basis prints under the synthetic names `x1 .. xn`,
and the notice about that goes to standard error, where it cannot corrupt
the output for the next tool. A JSON computation record carries variable names;
`verify` uses them only after independent certificate verification.

## Exit codes

| code | meaning |
|---|---|
| 0 | success |
| 1 | the verifier rejected untrusted certificate bytes |
| 2 | usage, input, or an operation the domain does not offer |
| 3 | a resource or structural limit: a deadline, a memory limit, a verifier cap |
| 4 | a defect in this crate: a certificate its own verifier rejected |
| 5 | a file that cannot be read or written, or a broken pipe |

A reader that closes the pipe, as `sylv gb system.ms | head` does, stops
the command with code 5 and no message on standard error. Every other fault
writes one line there, of the form `sylv: the modulus 12 is not prime`.
