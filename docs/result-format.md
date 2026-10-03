# Computation record format

`sylv-result-v1` stores a polynomial system and its computed basis as JSON.
It is a data format, not a certificate contract.

| Field | Meaning |
|---|---|
| `schema` | `sylv-result-v1` |
| `order` | `grevlex-v1` |
| `domain` | `{"kind":"prime_field","modulus":p}` or `{"kind":"rationals"}` |
| `variables` | Variable names in ring order |
| `input` | Original input polynomials as strings |
| `basis` | Claimed basis polynomials as strings |
| `claimed_provenance` | The producer's statement about the computation |
| `certificate` | Optional array of certificate bytes |

The provenance values are `unverified`, `supplied_basis`, `unchanged`,
`contains_input`, `equals_input`, and `certified`. Every value is an
untrusted claim on loading. Changing this field cannot establish a fact.

`unchanged`, `contains_input`, and `equals_input` apply to rational records.
`certified` requires a prime-field record with attached certificate bytes.
These shape checks do not validate the claims.

`unchanged` establishes nothing about the input ideal. `contains_input`
claims a reduced Gröbner basis of an ideal containing the input ideal.
It does not claim equality. `equals_input` claims that the producer checked
equality. None of these values constitutes independent certification.

To recover prime-field certification, verify the attached bytes and compare
the certificate's modulus, variable count, input, and basis with the record.
Variable names are local labels; their order fixes the coordinates.
`ResultEnvelope::verify_prime` performs the verification and comparison.
It applies both the caller's verifier caps and the operation budget.
Rational records have no certificate contract.

`ResultEnvelope::from_prime` attaches optional raw bytes with an `unverified`
claim. `ResultEnvelope::from_certified` verifies the certificate and its
association with the input and basis before recording `certified`.

`ResultEnvelope::from_checked_rational` records the input and basis retained
by `RationalEqualityCheck`. Its `equals_input` value is still an untrusted
claim after loading. `ResultEnvelope::check_rational` repeats the exact
check under a budget before returning the checked input, basis, and origins.

Loading checks the schema, monomial order, and ring metadata. Polynomial
strings remain data until a caller parses them under a budget. Checking a
loaded basis establishes properties of the ideal it generates. It does not
establish equality with the original input ideal.

Unknown fields, duplicate fields, and invalid JSON are rejected. Read and
write operations take a budget. The codec reads and writes bounded chunks and
polls cancellation and the deadline between chunks. The deadline covers
accounting, codec work, and header checks. Memory estimates cover JSON
storage; they do not cap process memory.
