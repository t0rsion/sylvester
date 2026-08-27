//! End-to-end tests for the `sylv-gb-cert-v2` verifier.
//!
//! The five fixtures are the byte strings of section 11 of
//! `docs/certificate-v2.md`, typed in as they stand there. The verifier
//! accepts all five and returns the modulus, the variable count, the input
//! list, and the basis of the table.
//!
//! The other certificates come from the small writer in this file. It
//! encodes the sections and computes the six section lengths, so a test can
//! change one field of a valid certificate and keep the bytes canonical.
//! The writer exists only to produce test bytes. It shares no code with
//! `src/cert`.

use std::time::{Duration, Instant};

use proptest::prelude::*;
use sylvester::verify::{
    BasisFault, BinaryFault, Cap, DivisionFault, DivisionSite, Limits, Location, NodeFault,
    PolyFault, PoolFault, VerifyError, WitnessFault, verify, verify_with_limits,
};

fn term_bytes(nvars: usize) -> usize {
    size_of::<sylvester::verify::Term>() + nvars * size_of::<sylvester::verify::Exp>()
}

/// The tiny fixture of section 11: `F = [x^2 + y, x]`, `G = [x, y]`.
#[rustfmt::skip]
const TINY: &[u8] = &[
    0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00, 0x0f, 0x73, 0x79, 0x6c,
    0x76, 0x2d, 0x67, 0x62, 0x2d, 0x63, 0x65, 0x72, 0x74, 0x2d, 0x76, 0x32,
    0x0a, 0x67, 0x72, 0x65, 0x76, 0x6c, 0x65, 0x78, 0x2d, 0x76, 0x31, 0x07,
    0x02, 0x0b, 0x09, 0x12, 0x03, 0x08, 0x01, 0x04, 0x00, 0x01, 0x01, 0x01,
    0x01, 0x00, 0x01, 0x01, 0x00, 0x02, 0x02, 0x02, 0x01, 0x03, 0x01, 0x01,
    0x01, 0x01, 0x02, 0x04, 0x00, 0x01, 0x00, 0x00, 0x02, 0x01, 0x01, 0x01,
    0x01, 0x02, 0x02, 0x01, 0x02, 0x00, 0x01, 0x01, 0x06, 0x02, 0x01, 0x03,
    0x02, 0x02, 0x00, 0x00, 0x01, 0x01, 0x00, 0x00, 0x00,
];

/// The square fixture: `F = [x^2*y, x*y^2]`, already the reduced basis.
#[rustfmt::skip]
const SQUARE: &[u8] = &[
    0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00, 0x0f, 0x73, 0x79, 0x6c,
    0x76, 0x2d, 0x67, 0x62, 0x2d, 0x63, 0x65, 0x72, 0x74, 0x2d, 0x76, 0x32,
    0x0a, 0x67, 0x72, 0x65, 0x76, 0x6c, 0x65, 0x78, 0x2d, 0x76, 0x31, 0x07,
    0x02, 0x0c, 0x07, 0x07, 0x03, 0x06, 0x02, 0x03, 0x00, 0x02, 0x00, 0x01,
    0x01, 0x02, 0x02, 0x00, 0x02, 0x01, 0x01, 0x02, 0x01, 0x01, 0x02, 0x01,
    0x01, 0x01, 0x02, 0x00, 0x01, 0x00, 0x00, 0x01, 0x01, 0x02, 0x00, 0x01,
    0x01, 0x00, 0x00, 0x01, 0x00, 0x01, 0x02, 0x00,
];

/// The unit fixture: `F = [3]`, `G = [1]`.
#[rustfmt::skip]
const UNIT: &[u8] = &[
    0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00, 0x0f, 0x73, 0x79, 0x6c,
    0x76, 0x2d, 0x67, 0x62, 0x2d, 0x63, 0x65, 0x72, 0x74, 0x2d, 0x76, 0x32,
    0x0a, 0x67, 0x72, 0x65, 0x76, 0x6c, 0x65, 0x78, 0x2d, 0x76, 0x31, 0x07,
    0x02, 0x02, 0x04, 0x08, 0x02, 0x03, 0x00, 0x01, 0x00, 0x01, 0x01, 0x03,
    0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x01, 0x00, 0x05, 0x01, 0x01, 0x01,
    0x00, 0x00,
];

/// The zero fixture: `F = [0]`, `G = []`.
#[rustfmt::skip]
const ZERO: &[u8] = &[
    0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00, 0x0f, 0x73, 0x79, 0x6c,
    0x76, 0x2d, 0x67, 0x62, 0x2d, 0x63, 0x65, 0x72, 0x74, 0x2d, 0x76, 0x32,
    0x0a, 0x67, 0x72, 0x65, 0x76, 0x6c, 0x65, 0x78, 0x2d, 0x76, 0x31, 0x07,
    0x02, 0x01, 0x02, 0x01, 0x01, 0x01, 0x00, 0x00, 0x01, 0x00, 0x00, 0x00,
    0x00,
];

/// The empty fixture: `F = []`, `G = []`.
#[rustfmt::skip]
const EMPTY: &[u8] = &[
    0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00, 0x0f, 0x73, 0x79, 0x6c,
    0x76, 0x2d, 0x67, 0x62, 0x2d, 0x63, 0x65, 0x72, 0x74, 0x2d, 0x76, 0x32,
    0x0a, 0x67, 0x72, 0x65, 0x76, 0x6c, 0x65, 0x78, 0x2d, 0x76, 0x31, 0x07,
    0x02, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
/// The five fixtures, by name.
fn fixtures() -> [(&'static str, &'static [u8]); 5] {
    [
        ("tiny", TINY),
        ("square", SQUARE),
        ("unit", UNIT),
        ("zero", ZERO),
        ("empty", EMPTY),
    ]
}

fn varint(value: u64, out: &mut Vec<u8>) {
    let mut left = value;
    loop {
        let byte = (left & 0x7f) as u8;
        left >>= 7;
        if left == 0 {
            out.push(byte);
            return;
        }
        out.push(byte | 0x80);
    }
}

fn string(text: &[u8], out: &mut Vec<u8>) {
    varint(text.len() as u64, out);
    out.extend_from_slice(text);
}

/// One trace node.
#[derive(Clone)]
enum Node {
    Input { index: u64, uses: u64 },
    Mul { src: u64, mono: u64, uses: u64 },
    Comb { steps: Vec<(u64, u64)>, uses: u64 },
    Scale { src: u64, scalar: u64, uses: u64 },
}

/// One pair witness.
#[derive(Clone)]
enum Witness {
    Coprime,
    Chain(u64),
    Reduce(Vec<(u64, u64)>),
}

/// A certificate, in the fields the byte format holds.
#[derive(Clone)]
struct Cert {
    magic: [u8; 8],
    schema: Vec<u8>,
    order: Vec<u8>,
    modulus: u64,
    nvars: u64,
    pool: Vec<Vec<(u64, u64)>>,
    input: Vec<Vec<(u64, u64)>>,
    trace: Vec<Node>,
    basis: Vec<u64>,
    membership: Vec<Vec<(u64, u64)>>,
    pairs: Vec<Witness>,
}

impl Cert {
    fn new(modulus: u64, nvars: u64) -> Self {
        Cert {
            magic: [0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00],
            schema: b"sylv-gb-cert-v2".to_vec(),
            order: b"grevlex-v1".to_vec(),
            modulus,
            nvars,
            pool: Vec::new(),
            input: Vec::new(),
            trace: Vec::new(),
            basis: Vec::new(),
            membership: Vec::new(),
            pairs: Vec::new(),
        }
    }

    fn pool_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        varint(self.pool.len() as u64, &mut out);
        for mono in &self.pool {
            varint(mono.len() as u64, &mut out);
            for &(variable, exponent) in mono {
                varint(variable, &mut out);
                varint(exponent, &mut out);
            }
        }
        out
    }

    fn input_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        varint(self.input.len() as u64, &mut out);
        for poly in &self.input {
            varint(poly.len() as u64, &mut out);
            for &(coeff, index) in poly {
                varint(coeff, &mut out);
                varint(index, &mut out);
            }
        }
        out
    }

    fn trace_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        varint(self.trace.len() as u64, &mut out);
        for node in &self.trace {
            match node {
                Node::Input { index, uses } => {
                    varint(0, &mut out);
                    varint(*uses, &mut out);
                    varint(*index, &mut out);
                }
                Node::Mul { src, mono, uses } => {
                    varint(1, &mut out);
                    varint(*uses, &mut out);
                    varint(*src, &mut out);
                    varint(*mono, &mut out);
                }
                Node::Comb { steps, uses } => {
                    varint(2, &mut out);
                    varint(*uses, &mut out);
                    varint(steps.len() as u64, &mut out);
                    let mut previous: Option<u64> = None;
                    for &(src, scalar) in steps {
                        match previous {
                            None => varint(src, &mut out),
                            Some(last) => varint(src - last - 1, &mut out),
                        }
                        varint(scalar, &mut out);
                        previous = Some(src);
                    }
                }
                Node::Scale { src, scalar, uses } => {
                    varint(3, &mut out);
                    varint(*uses, &mut out);
                    varint(*src, &mut out);
                    varint(*scalar, &mut out);
                }
            }
        }
        out
    }

    fn basis_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        varint(self.basis.len() as u64, &mut out);
        for &id in &self.basis {
            varint(id, &mut out);
        }
        out
    }

    fn steps(trace: &[(u64, u64)], out: &mut Vec<u8>) {
        varint(trace.len() as u64, out);
        for &(mono, element) in trace {
            varint(mono, out);
            varint(element, out);
        }
    }

    fn membership_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for trace in &self.membership {
            Cert::steps(trace, &mut out);
        }
        out
    }

    fn pairs_section(&self) -> Vec<u8> {
        let mut out = Vec::new();
        for witness in &self.pairs {
            match witness {
                Witness::Coprime => varint(0, &mut out),
                Witness::Chain(k) => {
                    varint(1, &mut out);
                    varint(*k, &mut out);
                }
                Witness::Reduce(trace) => {
                    varint(2, &mut out);
                    Cert::steps(trace, &mut out);
                }
            }
        }
        out
    }

    fn bytes(&self) -> Vec<u8> {
        let sections = [
            self.pool_section(),
            self.input_section(),
            self.trace_section(),
            self.basis_section(),
            self.membership_section(),
            self.pairs_section(),
        ];
        let mut out = Vec::new();
        out.extend_from_slice(&self.magic);
        string(&self.schema, &mut out);
        string(&self.order, &mut out);
        varint(self.modulus, &mut out);
        varint(self.nvars, &mut out);
        for section in &sections {
            varint(section.len() as u64, &mut out);
        }
        for section in &sections {
            out.extend_from_slice(section);
        }
        out
    }
}

/// The tiny fixture, in the writer's fields.
///
/// `F = [x^2 + y, x]` over `F_7[x, y]`, `G = [x, y]`. The pool is
/// `1, y, x, x^2`.
fn tiny() -> Cert {
    let mut cert = Cert::new(7, 2);
    cert.pool = vec![vec![], vec![(1, 1)], vec![(0, 1)], vec![(0, 2)]];
    cert.input = vec![vec![(1, 3), (1, 1)], vec![(1, 2)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 1 },
        Node::Input { index: 1, uses: 2 },
        Node::Mul {
            src: 1,
            mono: 2,
            uses: 1,
        },
        Node::Comb {
            steps: vec![(0, 1), (2, 6)],
            uses: 1,
        },
    ];
    cert.basis = vec![1, 3];
    cert.membership = vec![vec![(2, 0), (0, 1)], vec![(0, 0)]];
    cert.pairs = vec![Witness::Coprime];
    cert
}

/// `F = [x^2 - y, x*y - 1]` over `F_7[x, y]`.
///
/// The reduced basis is `G = [x^2 - y, x*y - 1, y^2 - x]`. The third
/// element comes from the S-polynomial of the first two:
/// `y * g_0 - x * g_1 = 6*y^2 + x`, scaled by 6. The trace holds one node
/// of every kind. The pair `(0, 2)` carries a `Chain` witness on `k = 1`,
/// and the other two pairs carry `Reduce` witnesses.
fn chained() -> Cert {
    let mut cert = Cert::new(7, 2);
    // 1, y, x, x*y, x^2
    cert.pool = vec![
        vec![],
        vec![(1, 1)],
        vec![(0, 1)],
        vec![(0, 1), (1, 1)],
        vec![(0, 2)],
    ];
    cert.input = vec![vec![(1, 4), (6, 1)], vec![(1, 3), (6, 0)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 2 },
        Node::Input { index: 1, uses: 2 },
        Node::Mul {
            src: 0,
            mono: 1,
            uses: 1,
        },
        Node::Mul {
            src: 1,
            mono: 2,
            uses: 1,
        },
        Node::Comb {
            steps: vec![(2, 1), (3, 6)],
            uses: 1,
        },
        Node::Scale {
            src: 4,
            scalar: 6,
            uses: 1,
        },
    ];
    cert.basis = vec![0, 1, 5];
    cert.membership = vec![vec![(0, 0)], vec![(0, 1)]];
    cert.pairs = vec![
        Witness::Reduce(vec![(0, 2)]),
        Witness::Chain(1),
        Witness::Reduce(vec![(0, 0)]),
    ];
    cert
}

/// The monomial ideal `[y^2*z^2, x*y^2, x*z^2]` over `F_7[x, y, z]`.
///
/// The three generators are their own reduced basis. Every S-polynomial of
/// the three is zero, so every pair takes an empty `Reduce` witness. The
/// leading monomials also make two `Chain` witnesses structurally valid in
/// both directions, which the cycle test uses.
fn monomials() -> Cert {
    let mut cert = Cert::new(7, 3);
    // 1, x*z^2, x*y^2, y^2*z^2
    cert.pool = vec![
        vec![],
        vec![(0, 1), (2, 2)],
        vec![(0, 1), (1, 2)],
        vec![(1, 2), (2, 2)],
    ];
    cert.input = vec![vec![(1, 3)], vec![(1, 2)], vec![(1, 1)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 1 },
        Node::Input { index: 1, uses: 1 },
        Node::Input { index: 2, uses: 1 },
    ];
    cert.basis = vec![0, 1, 2];
    cert.membership = vec![vec![(0, 0)], vec![(0, 1)], vec![(0, 2)]];
    cert.pairs = vec![
        Witness::Reduce(Vec::new()),
        Witness::Reduce(Vec::new()),
        Witness::Reduce(Vec::new()),
    ];
    cert
}

/// `F = [x]` with the claimed basis `[x^2, x]` over `F_7[x]`.
///
/// The list is monic, strictly descending, and duplicate-free. The one pair
/// has a zero S-polynomial, so an empty `Reduce` witness stands, and the
/// membership trace of `x` runs over the second element. Only minimality
/// rejects the certificate: `x` divides `x^2`, and the reduced basis is
/// `[x]`.
fn redundant() -> Cert {
    let mut cert = Cert::new(7, 1);
    cert.pool = vec![vec![], vec![(0, 1)]];
    cert.input = vec![vec![(1, 1)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 2 },
        Node::Mul {
            src: 0,
            mono: 1,
            uses: 1,
        },
    ];
    cert.basis = vec![1, 0];
    cert.membership = vec![vec![(0, 1)]];
    cert.pairs = vec![Witness::Reduce(Vec::new())];
    cert
}

fn accept(bytes: &[u8]) -> sylvester::verify::VerifiedGb {
    verify(bytes).expect("the certificate holds")
}

fn reject(cert: &Cert) -> VerifyError {
    verify(&cert.bytes()).expect_err("the verifier must reject the certificate")
}

/// The terms of a polynomial, as coefficients and exponent vectors.
fn terms(poly: &sylvester::verify::Poly) -> Vec<(u64, Vec<u32>)> {
    poly.terms()
        .iter()
        .map(|term| (term.coeff(), term.mono().exps().to_vec()))
        .collect()
}

fn with_byte(bytes: &[u8], at: usize, value: u8) -> Vec<u8> {
    let mut out = bytes.to_vec();
    out[at] = value;
    out
}

fn splice(bytes: &[u8], at: usize, take: usize, insert: &[u8]) -> Vec<u8> {
    let mut out = bytes[..at].to_vec();
    out.extend_from_slice(insert);
    out.extend_from_slice(&bytes[at + take..]);
    out
}

#[test]
fn the_five_contract_fixtures_are_accepted() {
    for (name, bytes) in fixtures() {
        let verified = verify(bytes).unwrap_or_else(|error| panic!("{name}: {error}"));
        assert_eq!(verified.modulus(), 7, "{name}");
        assert_eq!(verified.nvars(), 2, "{name}");
    }
}

#[test]
fn the_worked_example_decodes_to_the_stated_basis() {
    let verified = accept(TINY);
    assert_eq!(
        verified.input().iter().map(terms).collect::<Vec<_>>(),
        vec![
            vec![(1, vec![2, 0]), (1, vec![0, 1])],
            vec![(1, vec![1, 0])]
        ]
    );
    assert_eq!(
        verified.basis().iter().map(terms).collect::<Vec<_>>(),
        vec![vec![(1, vec![1, 0])], vec![(1, vec![0, 1])]]
    );
}

#[test]
fn the_square_fixture_decodes_to_its_own_input() {
    let verified = accept(SQUARE);
    assert_eq!(
        verified.basis().iter().map(terms).collect::<Vec<_>>(),
        vec![vec![(1, vec![2, 1])], vec![(1, vec![1, 2])]]
    );
    assert_eq!(
        verified.input().iter().map(terms).collect::<Vec<_>>(),
        verified.basis().iter().map(terms).collect::<Vec<_>>()
    );
}

#[test]
fn the_unit_fixture_decodes_to_the_unit_basis() {
    let verified = accept(UNIT);
    assert_eq!(
        verified.input().iter().map(terms).collect::<Vec<_>>(),
        vec![vec![(3, vec![0, 0])]]
    );
    assert_eq!(
        verified.basis().iter().map(terms).collect::<Vec<_>>(),
        vec![vec![(1, vec![0, 0])]]
    );
}

#[test]
fn the_zero_fixture_decodes_to_an_empty_basis() {
    let verified = accept(ZERO);
    assert_eq!(verified.input().len(), 1);
    assert!(verified.input()[0].is_zero());
    assert!(verified.basis().is_empty());
}

#[test]
fn the_empty_fixture_decodes_to_an_empty_input_and_basis() {
    let verified = accept(EMPTY);
    assert!(verified.input().is_empty());
    assert!(verified.basis().is_empty());
}

#[test]
fn the_v2_entry_point_accepts_the_same_bytes() {
    let verified = sylvester::verify::v2::verify(TINY).expect("the certificate holds");
    assert_eq!(verified.basis().len(), 2);
    let verified = sylvester::verify::v2::verify_with_limits(TINY, &Limits::default())
        .expect("the certificate holds");
    assert_eq!(verified.basis().len(), 2);
}

#[test]
fn the_test_writer_reproduces_the_tiny_fixture() {
    // The writer in this file encodes the same abstract certificate as the
    // listing of section 10. Every later test starts from a certificate it
    // built, so this equality pins the encoder the tests rely on.
    assert_eq!(tiny().bytes(), TINY);
}

#[test]
fn a_three_element_basis_with_a_chain_witness_is_accepted() {
    let verified = accept(&chained().bytes());
    assert_eq!(
        verified.basis().iter().map(terms).collect::<Vec<_>>(),
        vec![
            vec![(1, vec![2, 0]), (6, vec![0, 1])],
            vec![(1, vec![1, 1]), (6, vec![0, 0])],
            vec![(1, vec![0, 2]), (6, vec![1, 0])],
        ]
    );
}

#[test]
fn a_monomial_basis_over_three_variables_is_accepted() {
    let verified = accept(&monomials().bytes());
    assert_eq!(verified.nvars(), 3);
    assert_eq!(
        verified.basis().iter().map(terms).collect::<Vec<_>>(),
        vec![
            vec![(1, vec![0, 2, 2])],
            vec![(1, vec![1, 2, 0])],
            vec![(1, vec![1, 0, 2])],
        ]
    );
}

#[test]
fn a_chain_witness_may_depend_on_a_pair_that_carries_a_chain() {
    // Pair (0,1) takes Chain(2), whose dependencies are (0,2) and (1,2).
    // Pair (0,2) keeps its Reduce witness, so the graph is acyclic and the
    // reverse topological order validates (0,2) and (1,2) first.
    let mut cert = monomials();
    cert.pairs[0] = Witness::Chain(2);
    accept(&cert.bytes());
}

#[test]
fn a_wrong_first_byte_names_no_format() {
    assert_eq!(
        verify(&with_byte(TINY, 0, b'A')),
        Err(VerifyError::Format { first: Some(b'A') })
    );
    assert_eq!(verify(&[]), Err(VerifyError::Format { first: None }));
}

#[test]
fn a_wrong_magic_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 5, 0x00)),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Magic,
            offset: 0
        })
    );
}

#[test]
fn another_format_number_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 6, 0x03)),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Version {
                format: 3,
                revision: 0
            },
            offset: 6
        })
    );
    assert_eq!(
        verify(&with_byte(TINY, 7, 0x01)),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Version {
                format: 2,
                revision: 1
            },
            offset: 6
        })
    );
}

#[test]
fn a_wrong_schema_string_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 9, b't')),
        Err(VerifyError::Schema {
            found: "tylv-gb-cert-v2".to_string()
        })
    );
}

#[test]
fn a_wrong_order_string_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 25, b'l')),
        Err(VerifyError::Order {
            found: "lrevlex-v1".to_string()
        })
    );
}

#[test]
fn a_non_minimal_varint_is_a_rejection() {
    // The pool count 4 written in two bytes. The pool section grows by one
    // byte, so its declared length grows with it.
    let mut bytes = splice(TINY, 43, 1, &[0x84, 0x00]);
    bytes[37] = 0x0c;
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::NonMinimalVarint,
            offset: 43
        })
    );
}

#[test]
fn a_varint_longer_than_ten_bytes_is_a_rejection() {
    let long = [0xff; 11];
    let mut bytes = splice(TINY, 43, 1, &long);
    bytes[37] = 0x0b + 10;
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::VarintTooLong,
            offset: 43
        })
    );
}

#[test]
fn a_varint_above_the_range_of_u64_is_a_rejection() {
    let wide = [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x02];
    let mut bytes = splice(TINY, 43, 1, &wide);
    bytes[37] = 0x0b + 9;
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::VarintOverflow,
            offset: 43
        })
    );
}

#[test]
fn a_combination_delta_that_would_wrap_is_a_rejection() {
    // The second step of node 3 holds the delta 1. A delta of u64::MAX
    // would wrap the source id back to 0 without the checked addition.
    let wide = [0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x01];
    let mut bytes = splice(TINY, 79, 1, &wide);
    bytes[39] = 0x12 + 9;
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Overflow {
                what: "a combination source id"
            },
            offset: 89
        })
    );
}

#[test]
fn a_section_length_that_does_not_match_the_byte_count_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 37, 0x0c)),
        Err(VerifyError::MalformedBinary {
            reason: BinaryFault::SectionLengths,
            offset: 43
        })
    );
}

#[test]
fn a_pool_that_is_not_strictly_ascending_is_a_rejection() {
    // After y and x swap, the pool falls at index 2.
    let bytes = splice(TINY, 45, 6, &[0x01, 0x00, 0x01, 0x01, 0x01, 0x01]);
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::Pool {
            index: 2,
            fault: PoolFault::NotAscending
        })
    );
}

#[test]
fn a_duplicate_pool_entry_is_a_rejection() {
    // Entry 2 repeats entry 1, which strict ascent rejects.
    let bytes = splice(TINY, 48, 3, &[0x01, 0x01, 0x01]);
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::Pool {
            index: 2,
            fault: PoolFault::NotAscending
        })
    );
}

#[test]
fn a_zero_exponent_in_the_pool_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 47, 0x00)),
        Err(VerifyError::Pool {
            index: 1,
            fault: PoolFault::ExponentZero
        })
    );
}

#[test]
fn a_pool_variable_above_the_variable_count_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 46, 0x02)),
        Err(VerifyError::Pool {
            index: 1,
            fault: PoolFault::VariableOutOfRange { variable: 2 }
        })
    );
}

#[test]
fn an_unreferenced_pool_entry_is_a_rejection() {
    let mut cert = tiny();
    cert.pool.push(vec![(0, 3)]);
    assert_eq!(
        reject(&cert),
        VerifyError::Pool {
            index: 4,
            fault: PoolFault::Unreferenced
        }
    );
}

#[test]
fn a_pool_index_above_the_pool_count_is_a_rejection() {
    let mut cert = tiny();
    cert.input[1][0].1 = 9;
    assert_eq!(
        reject(&cert),
        VerifyError::IndexOutOfRange {
            what: "pool index",
            index: 9,
            bound: 4
        }
    );
}

#[test]
fn a_node_that_names_a_later_node_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 72, 0x02)),
        Err(VerifyError::Trace {
            node: 2,
            fault: NodeFault::SourceNotEarlier { src: 2 }
        })
    );
    assert_eq!(
        verify(&with_byte(TINY, 72, 0x03)),
        Err(VerifyError::Trace {
            node: 2,
            fault: NodeFault::SourceNotEarlier { src: 3 }
        })
    );
}

#[test]
fn a_wrong_use_count_is_a_rejection() {
    // Node 1 carries two references. A count of 3 leaves one behind.
    assert_eq!(
        verify(&with_byte(TINY, 68, 0x03)),
        Err(VerifyError::Trace {
            node: 1,
            fault: NodeFault::UseCountMismatch
        })
    );
    // A count of 1 frees the value before the basis section names it.
    assert_eq!(
        verify(&with_byte(TINY, 68, 0x01)),
        Err(VerifyError::Trace {
            node: 1,
            fault: NodeFault::UseCountMismatch
        })
    );
}

#[test]
fn a_zero_use_count_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 68, 0x00)),
        Err(VerifyError::Trace {
            node: 1,
            fault: NodeFault::UseCountZero
        })
    );
}

#[test]
fn an_unknown_node_kind_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 64, 0x04)),
        Err(VerifyError::Trace {
            node: 0,
            fault: NodeFault::UnknownKind { kind: 4 }
        })
    );
}

#[test]
fn a_combination_with_one_step_is_a_rejection() {
    let mut cert = tiny();
    cert.trace[3] = Node::Comb {
        steps: vec![(0, 1)],
        uses: 1,
    };
    assert_eq!(
        reject(&cert),
        VerifyError::Trace {
            node: 3,
            fault: NodeFault::CombTooShort { steps: 1 }
        }
    );
}

#[test]
fn a_scale_by_one_is_a_rejection() {
    let mut cert = chained();
    cert.trace[5] = Node::Scale {
        src: 4,
        scalar: 1,
        uses: 1,
    };
    assert_eq!(
        reject(&cert),
        VerifyError::Trace {
            node: 5,
            fault: NodeFault::ScalarOutOfRange { scalar: 1 }
        }
    );
}

#[test]
fn a_multiply_by_the_identity_monomial_is_a_rejection() {
    let mut cert = tiny();
    cert.trace[2] = Node::Mul {
        src: 1,
        mono: 0,
        uses: 1,
    };
    assert_eq!(
        reject(&cert),
        VerifyError::Trace {
            node: 2,
            fault: NodeFault::IdentityMultiplier
        }
    );
}

#[test]
fn an_input_node_after_another_kind_is_a_rejection() {
    let mut cert = tiny();
    cert.trace.push(Node::Input { index: 0, uses: 1 });
    assert_eq!(
        reject(&cert),
        VerifyError::Trace {
            node: 4,
            fault: NodeFault::InputOutOfOrder
        }
    );
}

#[test]
fn input_nodes_out_of_index_order_are_a_rejection() {
    let mut cert = tiny();
    cert.trace[0] = Node::Input { index: 1, uses: 1 };
    cert.trace[1] = Node::Input { index: 0, uses: 2 };
    assert_eq!(
        reject(&cert),
        VerifyError::Trace {
            node: 1,
            fault: NodeFault::InputOutOfOrder
        }
    );
}

#[test]
fn a_basis_node_id_above_the_node_count_is_a_rejection() {
    let mut cert = tiny();
    cert.basis[1] = 9;
    assert_eq!(
        reject(&cert),
        VerifyError::IndexOutOfRange {
            what: "basis node id",
            index: 9,
            bound: 4
        }
    );
}

#[test]
fn a_basis_element_that_is_not_monic_is_a_rejection() {
    // The second input polynomial becomes 2*x, and node 1 carries it into
    // the basis.
    assert_eq!(
        verify(&with_byte(TINY, 61, 0x02)),
        Err(VerifyError::Basis {
            index: 0,
            fault: BasisFault::NotMonic { lc: 2 }
        })
    );
}

#[test]
fn a_basis_that_is_not_sorted_is_a_rejection() {
    let bytes = splice(TINY, 82, 2, &[0x03, 0x01]);
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::Basis {
            index: 1,
            fault: BasisFault::NotDescending
        })
    );
}

#[test]
fn a_duplicate_basis_element_is_a_rejection() {
    // Two basis entries name node 0, which the use count of that node
    // covers. The elements are then equal, so the strict sort rejects them.
    let mut cert = monomials();
    cert.input.pop();
    cert.trace = vec![
        Node::Input { index: 0, uses: 2 },
        Node::Input { index: 1, uses: 1 },
    ];
    cert.basis = vec![0, 0, 1];
    cert.membership.pop();
    assert_eq!(
        reject(&cert),
        VerifyError::Basis {
            index: 1,
            fault: BasisFault::NotDescending
        }
    );
}

#[test]
fn a_basis_that_is_not_interreduced_is_a_rejection() {
    // G = [x^2 + y, y]. The tail of the first element is the leading
    // monomial of the second.
    let mut cert = Cert::new(7, 2);
    cert.pool = vec![vec![], vec![(1, 1)], vec![(0, 2)]];
    cert.input = vec![vec![(1, 2), (1, 1)], vec![(1, 1)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 1 },
        Node::Input { index: 1, uses: 1 },
    ];
    cert.basis = vec![0, 1];
    cert.membership = vec![vec![(0, 0)], vec![(0, 1)]];
    cert.pairs = vec![Witness::Coprime];
    assert_eq!(
        reject(&cert),
        VerifyError::Basis {
            index: 0,
            fault: BasisFault::TailReducible { term: 1, by: 1 }
        }
    );
}

#[test]
fn a_membership_trace_that_does_not_reach_zero_is_a_rejection() {
    let mut cert = tiny();
    cert.membership[0].pop();
    assert_eq!(
        reject(&cert),
        VerifyError::Division {
            at: DivisionSite::Membership { input: 0 },
            step: 1,
            fault: DivisionFault::NotZero
        }
    );
}

#[test]
fn a_membership_step_that_misses_the_lead_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 88, 0x00)),
        Err(VerifyError::Division {
            at: DivisionSite::Membership { input: 0 },
            step: 1,
            fault: DivisionFault::LeadMismatch
        })
    );
}

#[test]
fn a_coprime_witness_on_a_pair_that_shares_a_variable_is_a_rejection() {
    let mut cert = chained();
    cert.pairs[0] = Witness::Coprime;
    assert_eq!(
        reject(&cert),
        VerifyError::Witness {
            i: 0,
            j: 1,
            fault: WitnessFault::NotCoprime
        }
    );
}

#[test]
fn a_chain_whose_element_does_not_divide_the_lcm_is_a_rejection() {
    // The pair (0,1) has lcm x^2*y, and the third element leads with y^2.
    let mut cert = chained();
    cert.pairs[0] = Witness::Chain(2);
    assert_eq!(
        reject(&cert),
        VerifyError::Witness {
            i: 0,
            j: 1,
            fault: WitnessFault::ChainNotDividing { k: 2 }
        }
    );
}

#[test]
fn a_chain_that_names_its_own_pair_is_a_rejection() {
    let mut cert = chained();
    cert.pairs[1] = Witness::Chain(2);
    assert_eq!(
        reject(&cert),
        VerifyError::Witness {
            i: 0,
            j: 2,
            fault: WitnessFault::ChainNamesPair { k: 2 }
        }
    );
}

#[test]
fn a_chain_index_above_the_basis_count_is_a_rejection() {
    let mut cert = chained();
    cert.pairs[1] = Witness::Chain(9);
    assert_eq!(
        reject(&cert),
        VerifyError::IndexOutOfRange {
            what: "chain index",
            index: 9,
            bound: 3
        }
    );
}

#[test]
fn two_chain_witnesses_that_name_each_other_are_a_rejection() {
    // Pair (0,1) takes Chain(2), which depends on (0,2). Pair (0,2) takes
    // Chain(1), which depends on (0,1). Both pass the structural rules, so
    // only the cycle check stops them.
    let mut cert = monomials();
    cert.pairs[0] = Witness::Chain(2);
    cert.pairs[1] = Witness::Chain(1);
    assert_eq!(
        reject(&cert),
        VerifyError::Witness {
            i: 0,
            j: 1,
            fault: WitnessFault::Cycle
        }
    );
}

#[test]
fn a_reduce_witness_that_does_not_realize_the_s_polynomial_is_a_rejection() {
    let mut cert = chained();
    cert.pairs[0] = Witness::Reduce(vec![(0, 0)]);
    assert_eq!(
        reject(&cert),
        VerifyError::Division {
            at: DivisionSite::Pair { i: 0, j: 1 },
            step: 0,
            fault: DivisionFault::LeadMismatch
        }
    );
}

#[test]
fn a_reduce_witness_with_a_step_for_a_zero_s_polynomial_is_a_rejection() {
    let mut cert = monomials();
    cert.pairs[0] = Witness::Reduce(vec![(0, 0)]);
    assert_eq!(
        reject(&cert),
        VerifyError::Division {
            at: DivisionSite::Pair { i: 0, j: 1 },
            step: 0,
            fault: DivisionFault::ResidualZero
        }
    );
}

#[test]
fn a_missing_pair_is_a_rejection() {
    let mut cert = monomials();
    cert.pairs.pop();
    match reject(&cert) {
        VerifyError::MalformedBinary { reason, .. } => {
            assert_eq!(reason, BinaryFault::Truncated)
        }
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn an_extra_pair_is_a_rejection() {
    let mut cert = monomials();
    cert.pairs.push(Witness::Reduce(Vec::new()));
    match reject(&cert) {
        VerifyError::MalformedBinary { reason, .. } => {
            assert_eq!(reason, BinaryFault::SectionNotConsumed)
        }
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn a_missing_membership_trace_is_a_rejection() {
    let mut cert = tiny();
    cert.membership.pop();
    match reject(&cert) {
        VerifyError::MalformedBinary { reason, .. } => {
            assert_eq!(reason, BinaryFault::Truncated)
        }
        other => panic!("expected a malformed certificate, got {other}"),
    }
}

#[test]
fn a_composite_modulus_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 35, 0x09)),
        Err(VerifyError::Modulus {
            found: 9,
            composite: true
        })
    );
    // 1373653 is the smallest composite that passes the bases 2 and 3. The
    // verifier runs twelve bases, so it names it composite.
    let mut cert = tiny();
    cert.modulus = 1373653;
    assert_eq!(
        reject(&cert),
        VerifyError::Modulus {
            found: 1373653,
            composite: true
        }
    );
}

#[test]
fn a_modulus_outside_the_contract_range_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 35, 0x01)),
        Err(VerifyError::Modulus {
            found: 1,
            composite: false
        })
    );
    let mut cert = tiny();
    cert.modulus = 1 << 31;
    assert_eq!(
        reject(&cert),
        VerifyError::Modulus {
            found: 1 << 31,
            composite: false
        }
    );
}

#[test]
fn a_variable_count_above_the_contract_range_is_a_rejection() {
    // The check runs as soon as the header names the count, before any
    // exponent vector is allocated.
    let mut cert = tiny();
    cert.nvars = 257;
    assert_eq!(
        reject(&cert),
        VerifyError::Nvars {
            found: 257,
            max: 256
        }
    );
}

#[test]
fn a_coefficient_outside_the_field_is_a_rejection() {
    assert_eq!(
        verify(&with_byte(TINY, 56, 0x00)),
        Err(VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::CoefficientZero { term: 0 }
        })
    );
    assert_eq!(
        verify(&with_byte(TINY, 56, 0x07)),
        Err(VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::CoefficientOutOfRange { term: 0, coeff: 7 }
        })
    );
}

#[test]
fn input_terms_that_do_not_descend_are_a_rejection() {
    // The first input polynomial's terms rise, which the descent check rejects.
    let bytes = splice(TINY, 56, 4, &[0x01, 0x01, 0x01, 0x03]);
    assert_eq!(
        verify(&bytes),
        Err(VerifyError::Polynomial {
            at: Location::Input(0),
            fault: PolyFault::NotDescending { term: 1 }
        })
    );
}

#[test]
fn a_work_cap_below_the_replay_cost_reports_exhaustion() {
    let limits = Limits {
        max_work_units: 40,
        ..Limits::default()
    };
    let error = verify_with_limits(TINY, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::WorkUnits,
            limit: 40
        }
    );
}

#[test]
fn a_passed_deadline_reports_exhaustion() {
    let limits = Limits {
        deadline: Some(Instant::now() - Duration::from_secs(1)),
        ..Limits::default()
    };
    for (name, bytes) in fixtures() {
        let error = verify_with_limits(bytes, &limits).expect_err("the deadline holds");
        assert_eq!(error, VerifyError::DeadlineExceeded, "{name}");
        assert!(error.is_exhaustion(), "{name}");
    }
}

#[test]
fn every_count_cap_reports_exhaustion_before_the_check_runs() {
    let cases: [(Limits, Cap, usize); 6] = [
        (
            Limits {
                max_bytes: 16,
                ..Limits::default()
            },
            Cap::Bytes,
            16,
        ),
        (
            Limits {
                max_pool_monomials: 2,
                ..Limits::default()
            },
            Cap::PoolMonomials,
            2,
        ),
        (
            Limits {
                max_pool_entries: 1,
                ..Limits::default()
            },
            Cap::PoolEntries,
            1,
        ),
        (
            Limits {
                max_input_polys: 1,
                ..Limits::default()
            },
            Cap::InputPolynomials,
            1,
        ),
        (
            Limits {
                max_nodes: 3,
                ..Limits::default()
            },
            Cap::Nodes,
            3,
        ),
        (
            Limits {
                max_basis: 1,
                ..Limits::default()
            },
            Cap::BasisElements,
            1,
        ),
    ];
    for (limits, cap, limit) in cases {
        let error = verify_with_limits(TINY, &limits).expect_err("the cap holds");
        assert!(error.is_exhaustion());
        assert_eq!(error, VerifyError::CapExceeded { cap, limit });
    }
}

#[test]
fn the_division_step_caps_report_exhaustion() {
    let limits = Limits {
        max_division_steps: 1,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(TINY, &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::DivisionSteps,
            limit: 1
        })
    );
    let limits = Limits {
        max_total_division_steps: 2,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(TINY, &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::TotalDivisionSteps,
            limit: 2
        })
    );
}

#[test]
fn the_pair_cap_runs_before_anything_reads_the_pairs() {
    let limits = Limits {
        max_pairs: 2,
        ..Limits::default()
    };
    let cert = monomials();
    assert_eq!(
        verify_with_limits(&cert.bytes(), &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::Pairs,
            limit: 2
        })
    );
}

#[test]
fn the_live_byte_cap_reports_exhaustion() {
    let limits = Limits {
        max_live_bytes: term_bytes(2),
        ..Limits::default()
    };
    let error = verify_with_limits(&chained().bytes(), &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: term_bytes(2)
        }
    );
}

#[test]
fn every_rejection_apart_from_a_cap_reports_invalidity() {
    let corpus: Vec<Vec<u8>> = vec![
        with_byte(TINY, 0, b'A'),
        with_byte(TINY, 5, 0x00),
        with_byte(TINY, 35, 0x09),
        with_byte(TINY, 37, 0x0c),
        with_byte(TINY, 61, 0x02),
        with_byte(TINY, 68, 0x03),
        with_byte(TINY, 72, 0x02),
        with_byte(TINY, 88, 0x00),
    ];
    for bytes in corpus {
        let error = verify(&bytes).expect_err("the verifier must reject the change");
        assert!(!error.is_exhaustion(), "{error}");
    }
}

#[test]
fn no_byte_mutation_reaches_a_panic() {
    const REPLACEMENTS: [u8; 8] = [0x00, 0x01, 0x02, 0x03, 0x7f, 0x80, 0xfe, 0xff];
    let built = [
        chained().bytes(),
        monomials().bytes(),
        redundant().bytes(),
        shared_basis_node().bytes(),
    ];
    let mut mutations = 0;
    for bytes in fixtures()
        .iter()
        .map(|(_, bytes)| bytes.to_vec())
        .chain(built)
    {
        for position in 0..bytes.len() {
            for replacement in REPLACEMENTS {
                let mut mutated = bytes.clone();
                mutated[position] = replacement;
                let _ = verify(&mutated);
            }
            let mut dropped = bytes.clone();
            dropped.remove(position);
            let _ = verify(&dropped);
            let mut repeated = bytes.clone();
            repeated.insert(position, bytes[position]);
            let _ = verify(&repeated);
            mutations += REPLACEMENTS.len() + 2;
        }
    }
    assert!(mutations > 3_000, "the corpus ran {mutations} mutations");
}

#[test]
fn no_prefix_of_a_certificate_reaches_a_panic() {
    let built = [
        chained().bytes(),
        monomials().bytes(),
        redundant().bytes(),
        shared_basis_node().bytes(),
    ];
    for bytes in fixtures()
        .iter()
        .map(|(_, bytes)| bytes.to_vec())
        .chain(built)
    {
        for end in 0..bytes.len() {
            let error = verify(&bytes[..end]).expect_err("a prefix is not a certificate");
            assert!(!error.is_exhaustion(), "prefix of {end} bytes: {error}");
        }
    }
}

#[test]
fn the_v2_verifier_shares_no_code_with_the_engines() {
    let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("src/verify/v2");
    let forbidden = [
        "crate::ring",
        "crate::poly",
        "crate::compute",
        "crate::cert",
        "crate::ideal",
        "use super::super",
        "serde",
        "crate::verify::algebra",
        "crate::verify::json",
        "crate::verify::checks",
    ];
    let mut checked = 0;
    for entry in std::fs::read_dir(&root).expect("the v2 verifier directory is readable") {
        let path = entry.expect("the directory entry is readable").path();
        let text = std::fs::read_to_string(&path).expect("the file is readable");
        for needle in forbidden {
            assert!(!text.contains(needle), "{} holds {needle}", path.display());
        }
        checked += 1;
    }
    assert!(
        checked >= 5,
        "expected the v2 verifier module tree, found {checked} files"
    );
}

#[test]
fn a_pool_exponent_above_the_contract_range_is_a_rejection() {
    let mut cert = tiny();
    cert.pool[1] = vec![(1, 65536)];
    assert_eq!(
        reject(&cert),
        VerifyError::Pool {
            index: 1,
            fault: PoolFault::ExponentTooLarge { exponent: 65536 }
        }
    );
}

#[test]
fn a_product_above_the_exponent_bound_is_a_rejection() {
    // The bound holds for every monomial the verifier forms. The `Mul`
    // node squares x^65535.
    let mut cert = Cert::new(7, 2);
    cert.pool = vec![vec![], vec![(0, 65535)]];
    cert.input = vec![vec![(1, 1)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 1 },
        Node::Mul {
            src: 0,
            mono: 1,
            uses: 1,
        },
    ];
    cert.basis = vec![1];
    cert.membership = vec![vec![(0, 0)]];
    assert_eq!(reject(&cert), VerifyError::ExponentOverflow { max: 65535 });
}

#[test]
fn the_cap_on_one_arithmetic_result_reports_exhaustion() {
    let limits = Limits {
        max_intermediate_bytes: term_bytes(2),
        ..Limits::default()
    };
    let error = verify_with_limits(TINY, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit: term_bytes(2)
        }
    );
}

#[test]
fn the_combination_step_caps_report_exhaustion() {
    let limits = Limits {
        max_comb_steps: 1,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(TINY, &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::CombSteps,
            limit: 1
        })
    );
    let limits = Limits {
        max_trace_steps: 1,
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(TINY, &limits),
        Err(VerifyError::CapExceeded {
            cap: Cap::TraceSteps,
            limit: 1
        })
    );
}

#[test]
fn a_lead_that_divides_an_earlier_lead_is_a_rejection() {
    // Section 5.7: no leading monomial divides another leading monomial.
    // The divisor comes later here, because a monomial order never puts a
    // divisor above its multiple.
    assert_eq!(
        reject(&redundant()),
        VerifyError::Basis {
            index: 0,
            fault: BasisFault::LeadDivisible { by: 1 }
        }
    );
}

/// Raises the declared use count of one trace node.
fn bump_uses(node: &mut Node) {
    let uses = match node {
        Node::Input { uses, .. }
        | Node::Mul { uses, .. }
        | Node::Comb { uses, .. }
        | Node::Scale { uses, .. } => uses,
    };
    *uses += 1;
}

/// The five permutations of three elements that move at least one of them.
const SHUFFLES: [[usize; 3]; 5] = [[0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]];

proptest! {
    /// The reduced basis is unique, so only one order of a three-element
    /// basis can pass. Every other order breaks the strict descent.
    #[test]
    fn a_reordered_basis_is_a_rejection(choice in 0usize..SHUFFLES.len()) {
        let base = chained();
        let mut cert = base.clone();
        cert.basis = SHUFFLES[choice].iter().map(|&at| base.basis[at]).collect();
        let error = verify(&cert.bytes()).expect_err("only one order is sorted");
        prop_assert!(matches!(error, VerifyError::Basis { .. }), "{error}");
    }

    /// A multiple of an element is redundant wherever it is inserted. Its
    /// leading monomial carries the leading monomial of the element it
    /// multiplies, so the basis is no longer minimal.
    #[test]
    fn a_redundant_multiple_in_the_basis_is_a_rejection(
        source in 0usize..3,
        mono in 1u64..5,
        at in 0usize..4,
    ) {
        let mut cert = chained();
        let node = cert.basis[source];
        bump_uses(&mut cert.trace[node as usize]);
        let id = cert.trace.len() as u64;
        cert.trace.push(Node::Mul { src: node, mono, uses: 1 });
        cert.basis.insert(at, id);
        let error = verify(&cert.bytes()).expect_err("a multiple of an element is redundant");
        prop_assert!(matches!(error, VerifyError::Basis { .. }), "{error}");
    }
}

/// A certificate whose replay cost is almost all monomial comparisons.
///
/// `F = [f]` over `F_7` with 64 variables. `f` holds 8 terms, and every
/// monomial of `f` carries all 64 variables. `G = [f]`. Each of the 8
/// rounds scales `f` by 2 and adds `5 * f` to the result, which cancels to
/// zero and compares 8 pairs of monomials of support 64. The zero values
/// meet `f` in one final combination, which returns `f` and compares
/// nothing, because the zero values hold no term.
fn wide_comparisons() -> Cert {
    const VARS: u64 = 64;
    const TERMS: u64 = 8;
    const ROUNDS: u64 = 8;
    let mut cert = Cert::new(7, VARS);
    cert.pool = vec![Vec::new()];
    for term in 0..TERMS {
        cert.pool.push(
            (0..VARS)
                .map(|variable| (variable, if variable == 0 { term + 1 } else { 1 }))
                .collect(),
        );
    }
    cert.input = vec![(1..=TERMS).rev().map(|index| (1, index)).collect()];
    cert.trace = vec![Node::Input {
        index: 0,
        uses: 2 * ROUNDS + 1,
    }];
    for round in 0..ROUNDS {
        cert.trace.push(Node::Scale {
            src: 0,
            scalar: 2,
            uses: 1,
        });
        cert.trace.push(Node::Comb {
            steps: vec![(0, 5), (1 + 2 * round, 1)],
            uses: 1,
        });
    }
    let mut steps = vec![(0, 1)];
    steps.extend((1..=ROUNDS).map(|round| (2 * round, 1)));
    cert.trace.push(Node::Comb { steps, uses: 1 });
    cert.basis = vec![2 * ROUNDS + 1];
    cert.membership = vec![vec![(0, 0)]];
    cert
}

/// `G = [x^2 + y, x + 1, x + 1]`, where two basis entries name one node.
///
/// The format allows a shared basis node. The list then repeats an element,
/// which the strict descent of section 5.7 rejects. Before that check runs,
/// the basis section holds three elements of two terms each, so the live
/// total covers six terms.
fn shared_basis_node() -> Cert {
    let mut cert = Cert::new(7, 2);
    cert.pool = vec![vec![], vec![(1, 1)], vec![(0, 1)], vec![(0, 2)]];
    cert.input = vec![vec![(1, 3), (1, 1)], vec![(1, 2), (1, 0)]];
    cert.trace = vec![
        Node::Input { index: 0, uses: 1 },
        Node::Input { index: 1, uses: 2 },
    ];
    cert.basis = vec![0, 1, 1];
    cert.membership = vec![vec![(0, 0)], vec![(0, 1)]];
    cert.pairs = vec![
        Witness::Reduce(Vec::new()),
        Witness::Reduce(Vec::new()),
        Witness::Reduce(Vec::new()),
    ];
    cert
}

#[test]
fn the_comparison_charge_bounds_a_replay_over_wide_monomials() {
    // The eight rounds run 64 monomial comparisons, each over two supports
    // of 64 variables, so section 8.1 charges 8192 units for them alone.
    // The rest of the replay costs under 4000 units, so a cap of 8000 holds
    // only when the comparisons carry their support charge.
    let bytes = wide_comparisons().bytes();
    let limits = Limits {
        max_work_units: 8_000,
        ..Limits::default()
    };
    let error = verify_with_limits(&bytes, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::WorkUnits,
            limit: 8_000
        }
    );
    let verified = accept(&bytes);
    assert_eq!(verified.basis().len(), 1);
    assert_eq!(verified.basis()[0].terms().len(), 8);
}

#[test]
fn a_shared_basis_node_costs_one_live_copy_per_basis_entry() {
    // Three basis entries name two nodes of two terms each, so G holds six
    // terms. One copy per node would hold four and pass a cap of five.
    let bytes = shared_basis_node().bytes();
    let limit = 5 * term_bytes(2);
    let limits = Limits {
        max_live_bytes: limit,
        ..Limits::default()
    };
    let error = verify_with_limits(&bytes, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit
        }
    );
    // The format allows the shared node. The repeated element is what the
    // basis shape rejects.
    assert_eq!(
        reject(&shared_basis_node()),
        VerifyError::Basis {
            index: 2,
            fault: BasisFault::NotDescending
        }
    );
}

/// A certificate whose pool copies cost more than its bytes.
///
/// `F` holds eight copies of one polynomial `f` over `F_7` with 64
/// variables. `f` holds 8 terms, and every monomial of `f` carries all 64
/// variables. `G = [f]`, and each of the eight membership traces is one
/// step. The pool holds the 8 monomials once, so every copy of `f` reads
/// them again.
fn repeated_wide_terms() -> Cert {
    const VARS: u64 = 64;
    const TERMS: u64 = 8;
    const COPIES: usize = 8;
    let mut cert = Cert::new(7, VARS);
    cert.pool = vec![Vec::new()];
    for term in 0..TERMS {
        cert.pool.push(
            (0..VARS)
                .map(|variable| (variable, if variable == 0 { term + 1 } else { 1 }))
                .collect(),
        );
    }
    let poly: Vec<(u64, u64)> = (1..=TERMS).rev().map(|index| (1, index)).collect();
    cert.input = vec![poly; COPIES];
    cert.trace = vec![Node::Input { index: 0, uses: 1 }];
    cert.basis = vec![0];
    cert.membership = vec![vec![(0, 0)]; COPIES];
    cert
}

#[test]
fn the_pool_copy_charge_covers_the_support_it_copies() {
    // The eight copies of f read 64 pool monomials of support 64, so the
    // copies alone charge 4096 units. The rest of the check costs under
    // 18000 units, so a cap of 19000 holds only when a copy charges the
    // support it copies.
    let bytes = repeated_wide_terms().bytes();
    let limits = Limits {
        max_work_units: 19_000,
        ..Limits::default()
    };
    let error = verify_with_limits(&bytes, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::WorkUnits,
            limit: 19_000
        }
    );
    let verified = accept(&bytes);
    assert_eq!(verified.input().len(), 8);
    assert_eq!(verified.basis().len(), 1);
}

/// `F = [x0]` over `F_7` in 64 variables, `G = [x0]`. The single membership
/// trace carries 200 steps, each a distinct 64-variable multiplier from the
/// pool. The steps decode, and the trace holds one cloned monomial each,
/// before the replay runs.
fn long_membership_trace() -> Cert {
    const VARS: u64 = 64;
    const STEPS: u64 = 200;
    let mut cert = Cert::new(7, VARS);
    cert.pool = vec![vec![(0, 1)]];
    for step in 0..STEPS {
        let mut mono = vec![(0, step + 2)];
        for variable in 1..VARS {
            mono.push((variable, 1));
        }
        cert.pool.push(mono);
    }
    cert.input = vec![vec![(1, 0)]];
    cert.trace = vec![Node::Input { index: 0, uses: 1 }];
    cert.basis = vec![0];
    cert.membership = vec![(1..=STEPS).map(|index| (index, 0)).collect()];
    cert
}

#[test]
fn input_storage_counts_against_the_live_budget() {
    // repeated_wide_terms holds eight copies of an eight-term poly over 64
    // variables, so the input list is 64 resident terms. The trace and the
    // basis hold at most sixteen terms at once, so a cap of forty terms
    // admits the check only when the input storage goes uncharged. The
    // check charges it, so the cap reports exhaustion.
    let bytes = repeated_wide_terms().bytes();
    let limit = 40 * term_bytes(64);
    let limits = Limits {
        max_live_bytes: limit,
        ..Limits::default()
    };
    let error = verify_with_limits(&bytes, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit
        }
    );
    // The same certificate verifies under the default caps.
    assert_eq!(accept(&bytes).basis().len(), 1);
}

#[test]
fn a_membership_trace_counts_its_steps_against_the_live_budget() {
    // The 200 stored step monomials are 64 variables wide. A cap of 100
    // terms trips on them before the replay runs, so a division trace
    // cannot hold gigabytes under the step caps alone.
    let bytes = long_membership_trace().bytes();
    let limit = 100 * term_bytes(64);
    let limits = Limits {
        max_live_bytes: limit,
        ..Limits::default()
    };
    let error = verify_with_limits(&bytes, &limits).expect_err("the cap holds");
    assert!(error.is_exhaustion());
    assert_eq!(
        error,
        VerifyError::CapExceeded {
            cap: Cap::IntermediateBytes,
            limit
        }
    );
}

#[test]
fn a_passed_deadline_outranks_a_v2_acceptance() {
    // The empty certificate is valid. A passed deadline outranks the
    // acceptance, so the verifier reports the deadline, not the basis.
    let limits = Limits {
        deadline: Some(Instant::now() - Duration::from_secs(1)),
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(EMPTY, &limits),
        Err(VerifyError::DeadlineExceeded)
    );
}

#[test]
fn a_passed_deadline_outranks_a_v2_invalidity() {
    // A wrong schema string is an invalidity the header finds before any
    // polling loop runs. A passed deadline still outranks it, so the
    // verifier reports the deadline, not the schema fault.
    let mut cert = tiny();
    cert.schema = b"sylv-gb-cert-v9".to_vec();
    let limits = Limits {
        deadline: Some(Instant::now() - Duration::from_secs(1)),
        ..Limits::default()
    };
    assert_eq!(
        verify_with_limits(&cert.bytes(), &limits),
        Err(VerifyError::DeadlineExceeded)
    );
}
