//! Typed rejections. Every value names what failed.

use std::fmt;

/// A syntax fault in the certificate bytes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Syntax {
    /// The byte at this position is not the one the grammar requires.
    UnexpectedByte {
        /// The byte the decoder read.
        found: u8,
        /// The item the grammar requires here.
        expected: &'static str,
    },
    /// The bytes end inside the object.
    Truncated,
    /// Bytes follow the closing brace.
    TrailingBytes,
    /// Insignificant whitespace. The canonical encoding has none.
    Whitespace,
    /// A key appears twice.
    DuplicateKey(String),
    /// A key is not in the schema.
    UnknownKey(String),
    /// A key the schema requires is absent.
    MissingKey(&'static str),
    /// A key appears before the key the schema puts first.
    KeyOutOfOrder {
        /// The key the decoder read.
        found: String,
        /// The key the schema puts at this position.
        expected: &'static str,
    },
    /// An integer carries a leading zero.
    LeadingZero,
    /// A number carries a fraction part.
    Float,
    /// A number carries an exponent part.
    Exponent,
    /// A number carries a sign. Every value in the schema is positive or zero.
    Sign,
    /// An integer does not fit 64 bits.
    IntegerOverflow {
        /// The name of the value.
        what: &'static str,
    },
    /// An integer is outside the range the schema allows.
    IntegerOutOfRange {
        /// The name of the value.
        what: &'static str,
        /// The value the decoder read.
        value: u64,
        /// The largest value the schema allows.
        max: u64,
    },
    /// A string holds a byte the canonical encoding forbids.
    BadString,
    /// A term's exponent array holds more entries than the contract's
    /// variable width allows.
    ExponentWidth {
        /// The largest width the contract allows.
        max: u64,
    },
}

impl fmt::Display for Syntax {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Syntax::UnexpectedByte { found, expected } => {
                write!(f, "found byte {found:#04x}, expected {expected}")
            }
            Syntax::Truncated => write!(f, "the bytes end inside the object"),
            Syntax::TrailingBytes => write!(f, "bytes follow the object"),
            Syntax::Whitespace => write!(f, "whitespace is not canonical"),
            Syntax::DuplicateKey(key) => write!(f, "duplicate key \"{key}\""),
            Syntax::UnknownKey(key) => write!(f, "unknown key \"{key}\""),
            Syntax::MissingKey(key) => write!(f, "missing key \"{key}\""),
            Syntax::KeyOutOfOrder { found, expected } => {
                write!(f, "found key \"{found}\", expected key \"{expected}\"")
            }
            Syntax::LeadingZero => write!(f, "leading zero in an integer"),
            Syntax::Float => write!(f, "a fraction part is not an integer"),
            Syntax::Exponent => write!(f, "an exponent part is not an integer"),
            Syntax::Sign => write!(f, "a sign is not canonical here"),
            Syntax::IntegerOverflow { what } => write!(f, "{what} does not fit 64 bits"),
            Syntax::IntegerOutOfRange { what, value, max } => {
                write!(f, "{what} is {value}, the maximum is {max}")
            }
            Syntax::BadString => write!(f, "a string holds a forbidden byte"),
            Syntax::ExponentWidth { max } => {
                write!(f, "a term holds more than {max} exponents")
            }
        }
    }
}

/// The place a polynomial holds in the certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Location {
    /// An input polynomial, by index.
    Input(usize),
    /// A basis element, by index.
    Basis(usize),
    /// An origin cofactor.
    Origin {
        /// The index of the basis element.
        basis: usize,
        /// The index of the input polynomial.
        input: usize,
    },
    /// A membership cofactor.
    Membership {
        /// The index of the input polynomial.
        input: usize,
        /// The index of the basis element.
        basis: usize,
    },
    /// An S-pair cofactor.
    SpairCofactor {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
        /// The index of the basis element the cofactor multiplies.
        basis: usize,
    },
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Location::Input(index) => write!(f, "input polynomial {index}"),
            Location::Basis(index) => write!(f, "basis polynomial {index}"),
            Location::Origin { basis, input } => {
                write!(f, "origin cofactor of basis {basis} for input {input}")
            }
            Location::Membership { input, basis } => {
                write!(f, "membership cofactor of input {input} for basis {basis}")
            }
            Location::SpairCofactor { i, j, basis } => {
                write!(f, "S-pair ({i},{j}) cofactor for basis {basis}")
            }
        }
    }
}

/// A polynomial that the encoding rules reject.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PolyFault {
    /// The coefficient is zero. The encoding drops zero terms.
    CoefficientZero {
        /// The index of the term.
        term: usize,
    },
    /// The coefficient is not below the modulus.
    CoefficientOutOfRange {
        /// The index of the term.
        term: usize,
        /// The coefficient the decoder read.
        coeff: u64,
    },
    /// The exponent vector does not hold one entry per variable.
    ExponentCount {
        /// The index of the term.
        term: usize,
        /// The number of exponents the term holds.
        found: usize,
        /// The number of variables the certificate declares.
        expected: usize,
    },
    /// The term is not below the term before it.
    NotDescending {
        /// The index of the term.
        term: usize,
    },
    /// The term repeats the monomial of the term before it.
    DuplicateMonomial {
        /// The index of the term.
        term: usize,
    },
}

impl fmt::Display for PolyFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PolyFault::CoefficientZero { term } => {
                write!(f, "term {term} has coefficient 0")
            }
            PolyFault::CoefficientOutOfRange { term, coeff } => {
                write!(
                    f,
                    "term {term} has coefficient {coeff}, which the modulus excludes"
                )
            }
            PolyFault::ExponentCount {
                term,
                found,
                expected,
            } => write!(
                f,
                "term {term} holds {found} exponents, the certificate declares {expected} variables"
            ),
            PolyFault::NotDescending { term } => {
                write!(f, "term {term} is not below the term before it")
            }
            PolyFault::DuplicateMonomial { term } => {
                write!(f, "term {term} repeats the monomial of the term before it")
            }
        }
    }
}

/// A basis that is not the reduced basis shape the contract requires.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisFault {
    /// The element is the zero polynomial.
    Zero,
    /// The leading coefficient is not 1.
    NotMonic {
        /// The leading coefficient.
        lc: u64,
    },
    /// The leading monomial is not below the leading monomial before it.
    NotDescending,
    /// A non-leading term is divisible by the leading monomial of another
    /// element, so the basis is not interreduced.
    TailReducible {
        /// The index of the term.
        term: usize,
        /// The index of the element that reduces the term.
        by: usize,
    },
    /// The leading monomial is divisible by the leading monomial of another
    /// element, so the element is redundant.
    LeadDivisible {
        /// The index of the element that divides the leading monomial.
        by: usize,
    },
}

impl fmt::Display for BasisFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BasisFault::Zero => write!(f, "the element is zero"),
            BasisFault::NotMonic { lc } => {
                write!(f, "the leading coefficient is {lc}, not 1")
            }
            BasisFault::NotDescending => {
                write!(f, "the leading monomial is not below the one before it")
            }
            BasisFault::TailReducible { term, by } => write!(
                f,
                "term {term} is divisible by the leading monomial of element {by}"
            ),
            BasisFault::LeadDivisible { by } => write!(
                f,
                "the leading monomial is divisible by the leading monomial of element {by}"
            ),
        }
    }
}

/// A resource cap the verifier enforces before it allocates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cap {
    /// The size of the certificate.
    Bytes,
    /// The number of polynomials in the certificate.
    Polynomials,
    /// The number of terms in one polynomial.
    TermsPerPolynomial,
    /// The number of terms in the certificate.
    TotalTerms,
    /// The number of entries in the `origin`, `membership`, and `spairs`
    /// arrays together. One entry holds one list of cofactors.
    Entries,
    /// The number of bytes the verifier arithmetic holds at once. One
    /// prospective term costs `size_of::<Term>() + nvars * size_of::<Exp>()`.
    IntermediateBytes,
}

impl fmt::Display for Cap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Cap::Bytes => write!(f, "certificate bytes"),
            Cap::Polynomials => write!(f, "polynomial count"),
            Cap::TermsPerPolynomial => write!(f, "terms in one polynomial"),
            Cap::TotalTerms => write!(f, "total terms"),
            Cap::Entries => write!(f, "array entries"),
            Cap::IntermediateBytes => write!(f, "bytes of verifier arithmetic"),
        }
    }
}

/// Why the verifier rejects a certificate.
///
/// [`VerifyError::CapExceeded`] and [`VerifyError::DeadlineExceeded`] report
/// exhaustion. They say nothing about the certificate. Every other value
/// reports an invalid certificate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum VerifyError {
    /// The bytes are not a canonical encoding of the schema.
    Malformed {
        /// The syntax fault.
        reason: Syntax,
        /// The byte offset of the fault.
        offset: usize,
    },
    /// The `schema` string is not the one this verifier implements.
    Schema {
        /// The schema string in the certificate.
        found: String,
    },
    /// The `order` string is not the one this verifier implements.
    Order {
        /// The order string in the certificate.
        found: String,
    },
    /// The modulus is out of range or composite.
    Modulus {
        /// The modulus in the certificate.
        found: u64,
        /// True if the modulus is in range and composite.
        composite: bool,
    },
    /// The variable count is out of range.
    Nvars {
        /// The variable count in the certificate.
        found: u64,
        /// The largest variable count the contract allows.
        max: u64,
    },
    /// A polynomial encoding is not canonical.
    Polynomial {
        /// The place the polynomial holds.
        at: Location,
        /// The encoding rule the polynomial breaks.
        fault: PolyFault,
    },
    /// The basis is not a reduced basis.
    Basis {
        /// The index of the basis element.
        index: usize,
        /// The rule the basis breaks.
        fault: BasisFault,
    },
    /// The origin identity of a basis element does not hold.
    OriginIdentity {
        /// The index of the basis element.
        basis: usize,
    },
    /// The membership identity of an input polynomial does not hold.
    MembershipIdentity {
        /// The index of the input polynomial.
        input: usize,
    },
    /// A pair with non-coprime leading monomials has no S-pair entry.
    SpairMissing {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
    },
    /// Two S-pair entries name the same pair.
    SpairDuplicate {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
    },
    /// An S-pair entry does not follow the entry before it.
    SpairUnsorted {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
    },
    /// An S-pair entry names indices that are not in increasing order.
    SpairIndexOrder {
        /// The first index the entry names.
        i: usize,
        /// The second index the entry names.
        j: usize,
    },
    /// The S-pair identity does not hold.
    SpairIdentity {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
    },
    /// A summand of an S-pair identity is above the S-polynomial.
    SpairBound {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
        /// The index of the summand.
        summand: usize,
    },
    /// An index does not address an element.
    IndexOutOfRange {
        /// The name of the index.
        what: &'static str,
        /// The value of the index.
        index: u64,
        /// The number of elements the index must address.
        bound: usize,
    },
    /// An array holds the wrong number of entries.
    CountMismatch {
        /// The name of the array.
        what: &'static str,
        /// The index of the entry, when one entry holds the wrong count.
        index: Option<usize>,
        /// The number of entries the array holds.
        found: usize,
        /// The number of entries the contract requires.
        expected: usize,
    },
    /// A resource cap stopped the work.
    CapExceeded {
        /// The cap that stopped the work.
        cap: Cap,
        /// The value of the cap.
        limit: usize,
    },
    /// The deadline passed.
    DeadlineExceeded,
}

impl fmt::Display for VerifyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VerifyError::Malformed { reason, offset } => {
                write!(f, "malformed certificate at byte {offset}: {reason}")
            }
            VerifyError::Schema { found } => {
                write!(f, "schema is \"{found}\", expected \"sylv-gb-cert-v1\"")
            }
            VerifyError::Order { found } => {
                write!(f, "order is \"{found}\", expected \"grevlex-v1\"")
            }
            VerifyError::Modulus { found, composite } => {
                if *composite {
                    write!(f, "modulus {found} is composite")
                } else {
                    write!(f, "modulus {found} is out of range")
                }
            }
            VerifyError::Nvars { found, max } => {
                write!(f, "nvars is {found}, the maximum is {max}")
            }
            VerifyError::Polynomial { at, fault } => write!(f, "{at}: {fault}"),
            VerifyError::Basis { index, fault } => {
                write!(f, "basis element {index}: {fault}")
            }
            VerifyError::OriginIdentity { basis } => write!(
                f,
                "the origin identity of basis element {basis} does not hold"
            ),
            VerifyError::MembershipIdentity { input } => write!(
                f,
                "the membership identity of input polynomial {input} does not hold"
            ),
            VerifyError::SpairMissing { i, j } => write!(
                f,
                "pair ({i},{j}) has no entry and its leading monomials are not coprime"
            ),
            VerifyError::SpairDuplicate { i, j } => {
                write!(f, "pair ({i},{j}) has more than one entry")
            }
            VerifyError::SpairUnsorted { i, j } => {
                write!(f, "the entry for pair ({i},{j}) is out of order")
            }
            VerifyError::SpairIndexOrder { i, j } => {
                write!(f, "the entry names ({i},{j}); index {i} must be below {j}")
            }
            VerifyError::SpairIdentity { i, j } => {
                write!(f, "the identity of pair ({i},{j}) does not hold")
            }
            VerifyError::SpairBound { i, j, summand } => write!(
                f,
                "summand {summand} of pair ({i},{j}) is above the S-polynomial"
            ),
            VerifyError::IndexOutOfRange { what, index, bound } => {
                write!(f, "{what} is {index}, the count is {bound}")
            }
            VerifyError::CountMismatch {
                what,
                index,
                found,
                expected,
            } => match index {
                Some(index) => write!(
                    f,
                    "{what} {index} holds {found} entries, expected {expected}"
                ),
                None => write!(f, "{what} holds {found} entries, expected {expected}"),
            },
            VerifyError::CapExceeded { cap, limit } => {
                write!(
                    f,
                    "the certificate needs more than the cap on {cap}, which is {limit}"
                )
            }
            VerifyError::DeadlineExceeded => write!(f, "the deadline passed"),
        }
    }
}

impl std::error::Error for VerifyError {}

impl VerifyError {
    /// Report whether the value comes from a resource cap or a deadline.
    ///
    /// An exhausted verifier says nothing about the certificate.
    pub fn is_exhaustion(&self) -> bool {
        matches!(
            self,
            VerifyError::CapExceeded { .. } | VerifyError::DeadlineExceeded
        )
    }
}
