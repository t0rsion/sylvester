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
    /// The canonical encoding has no insignificant whitespace.
    Whitespace,
    /// A key appears twice.
    DuplicateKey(String),
    /// A key is not in the schema.
    UnknownKey(String),
    /// A key the schema requires is absent.
    MissingKey(&'static str),
    /// A key is not the one the schema puts at this position.
    KeyOutOfOrder {
        /// The key the decoder read.
        found: String,
        /// The key the schema puts at this position.
        expected: &'static str,
    },
    /// The canonical encoding has no leading zero except `0`.
    LeadingZero,
    /// The canonical encoding has no fraction part.
    Float,
    /// The canonical encoding has no exponent part.
    Exponent,
    /// Every value in the schema is positive or zero.
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
    /// The encoding drops zero terms.
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
    /// A reduced basis has no zero element.
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

/// A syntax fault in the v2 binary encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BinaryFault {
    /// The first eight bytes are not the v2 magic.
    Magic,
    /// The magic names another format number or another revision.
    Version {
        /// The format number in the magic.
        format: u8,
        /// The revision in the magic.
        revision: u8,
    },
    /// The bytes end inside an item, or an item crosses a section boundary.
    Truncated,
    /// The header size plus the six section lengths is not the byte count.
    SectionLengths,
    /// A section leaves bytes unread inside its own range.
    SectionNotConsumed,
    /// The canonical form has no trailing zero byte.
    NonMinimalVarint,
    /// A varint is longer than ten bytes.
    VarintTooLong,
    /// A varint does not fit 64 bits.
    VarintOverflow,
    /// A sum or a product of decoded integers does not fit 64 bits.
    Overflow {
        /// The name of the value.
        what: &'static str,
    },
}

impl fmt::Display for BinaryFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BinaryFault::Magic => write!(f, "the bytes do not carry the v2 magic"),
            BinaryFault::Version { format, revision } => write!(
                f,
                "the magic names format {format} revision {revision}, this verifier reads format 2 revision 0"
            ),
            BinaryFault::Truncated => write!(f, "the bytes end inside an item"),
            BinaryFault::SectionLengths => {
                write!(f, "the section lengths do not add up to the byte count")
            }
            BinaryFault::SectionNotConsumed => {
                write!(f, "the section leaves bytes unread")
            }
            BinaryFault::NonMinimalVarint => write!(f, "a varint carries a trailing zero byte"),
            BinaryFault::VarintTooLong => write!(f, "a varint is longer than ten bytes"),
            BinaryFault::VarintOverflow => write!(f, "a varint does not fit 64 bits"),
            BinaryFault::Overflow { what } => write!(f, "{what} does not fit 64 bits"),
        }
    }
}

/// A monomial pool entry that breaks a rule of the pool encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PoolFault {
    /// The support holds more entries than the certificate has variables.
    SupportTooLarge {
        /// The support size the decoder read.
        found: u64,
        /// The number of variables the certificate declares.
        max: usize,
    },
    /// The variable index is not below the variable count.
    VariableOutOfRange {
        /// The variable index.
        variable: u64,
    },
    /// The variable index is not above the one before it.
    VariableNotIncreasing,
    /// The sparse form drops zero exponents.
    ExponentZero,
    /// The exponent is above the contract range.
    ExponentTooLarge {
        /// The exponent the decoder read.
        exponent: u64,
    },
    /// The monomial is not above the monomial before it.
    NotAscending,
    /// No term, node, or step references the entry.
    Unreferenced,
}

impl fmt::Display for PoolFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PoolFault::SupportTooLarge { found, max } => write!(
                f,
                "the support holds {found} entries, the certificate declares {max} variables"
            ),
            PoolFault::VariableOutOfRange { variable } => {
                write!(f, "variable {variable} is not below the variable count")
            }
            PoolFault::VariableNotIncreasing => {
                write!(f, "the variable index is not above the one before it")
            }
            PoolFault::ExponentZero => write!(f, "the exponent is zero"),
            PoolFault::ExponentTooLarge { exponent } => {
                write!(f, "the exponent is {exponent}, the maximum is 65535")
            }
            PoolFault::NotAscending => {
                write!(f, "the monomial is not above the one before it")
            }
            PoolFault::Unreferenced => write!(f, "nothing references the entry"),
        }
    }
}

/// A trace node that breaks a rule of the trace encoding.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NodeFault {
    /// The node kind is not one the contract defines.
    UnknownKind {
        /// The kind the decoder read.
        kind: u64,
    },
    /// Every node has at least one reference.
    UseCountZero,
    /// The declared use count is not the number of references the verifier
    /// counted.
    UseCountMismatch,
    /// The source id is not below the node's own id.
    SourceNotEarlier {
        /// The source id.
        src: u64,
    },
    /// The `Input` nodes are not first, or their indices do not increase.
    InputOutOfOrder,
    /// The `Mul` monomial is the identity monomial.
    IdentityMultiplier,
    /// The `Scale` scalar is 1, or it is not below the modulus.
    ScalarOutOfRange {
        /// The scalar the decoder read.
        scalar: u64,
    },
    /// The `Comb` holds fewer than two steps.
    CombTooShort {
        /// The step count the decoder read.
        steps: u64,
    },
}

impl fmt::Display for NodeFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NodeFault::UnknownKind { kind } => write!(f, "node kind {kind} is not defined"),
            NodeFault::UseCountZero => write!(f, "the use count is zero"),
            NodeFault::UseCountMismatch => {
                write!(f, "the use count is not the number of references")
            }
            NodeFault::SourceNotEarlier { src } => {
                write!(f, "source {src} is not below the node's own id")
            }
            NodeFault::InputOutOfOrder => {
                write!(f, "the input nodes are not first in increasing input index")
            }
            NodeFault::IdentityMultiplier => {
                write!(f, "the multiplier is the identity monomial")
            }
            NodeFault::ScalarOutOfRange { scalar } => {
                write!(
                    f,
                    "the scalar is {scalar}, the range is 2 to the modulus minus 1"
                )
            }
            NodeFault::CombTooShort { steps } => {
                write!(f, "the combination holds {steps} steps, the minimum is 2")
            }
        }
    }
}

/// The place a division trace holds in the certificate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DivisionSite {
    /// The membership trace of an input polynomial.
    Membership {
        /// The index of the input polynomial.
        input: usize,
    },
    /// The `Reduce` witness of a pair.
    Pair {
        /// The first index of the pair.
        i: usize,
        /// The second index of the pair.
        j: usize,
    },
}

impl fmt::Display for DivisionSite {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DivisionSite::Membership { input } => {
                write!(f, "the membership trace of input polynomial {input}")
            }
            DivisionSite::Pair { i, j } => write!(f, "the reduce witness of pair ({i},{j})"),
        }
    }
}

/// A division trace that breaks a step rule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DivisionFault {
    /// The step runs on a zero residual.
    ResidualZero,
    /// The step's multiple does not lead with the residual's leading
    /// monomial.
    LeadMismatch,
    /// The residual after the last step is not zero.
    NotZero,
}

impl fmt::Display for DivisionFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DivisionFault::ResidualZero => write!(f, "the residual is already zero"),
            DivisionFault::LeadMismatch => write!(
                f,
                "the multiple does not lead with the leading monomial of the residual"
            ),
            DivisionFault::NotZero => write!(f, "the residual is not zero"),
        }
    }
}

/// A pair witness that breaks a rule.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WitnessFault {
    /// The witness kind is not one the contract defines.
    UnknownKind {
        /// The kind the decoder read.
        kind: u64,
    },
    /// A `Coprime` witness names a pair whose leading monomials share a
    /// variable.
    NotCoprime,
    /// A `Chain` witness names one of the two elements of its own pair.
    ChainNamesPair {
        /// The basis index the witness names.
        k: usize,
    },
    /// The leading monomial of the named element does not divide the least
    /// common multiple of the pair.
    ChainNotDividing {
        /// The basis index the witness names.
        k: usize,
    },
    /// The witness lies on a cycle of the dependency graph.
    Cycle,
}

impl fmt::Display for WitnessFault {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            WitnessFault::UnknownKind { kind } => write!(f, "witness kind {kind} is not defined"),
            WitnessFault::NotCoprime => {
                write!(f, "the leading monomials share a variable")
            }
            WitnessFault::ChainNamesPair { k } => {
                write!(f, "the chain names element {k} of its own pair")
            }
            WitnessFault::ChainNotDividing { k } => write!(
                f,
                "the leading monomial of element {k} does not divide the least common multiple"
            ),
            WitnessFault::Cycle => write!(f, "the witness lies on a dependency cycle"),
        }
    }
}

/// A resource cap the verifier enforces before it allocates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Cap {
    /// Certificate size in bytes.
    Bytes,
    /// Polynomials in the input, the basis, and every cofactor.
    Polynomials,
    /// Terms in one polynomial.
    TermsPerPolynomial,
    /// Terms in the certificate.
    TotalTerms,
    /// Entries in the `origin`, `membership`, and `spairs` arrays together.
    /// One entry holds one list of cofactors.
    Entries,
    /// Bytes of verifier arithmetic held at once. One prospective term
    /// costs `size_of::<Term>() + nvars * size_of::<Exp>()`.
    IntermediateBytes,
    /// Monomials in the v2 pool.
    PoolMonomials,
    /// Variable/exponent entries in the v2 pool.
    PoolEntries,
    /// Polynomials in the v2 input section.
    InputPolynomials,
    /// Nodes in the v2 trace.
    Nodes,
    /// Steps in one v2 combination node.
    CombSteps,
    /// Combination steps in the v2 trace.
    TraceSteps,
    /// Elements in the v2 basis.
    BasisElements,
    /// Pairs of v2 basis elements.
    Pairs,
    /// Steps in one v2 division trace.
    DivisionSteps,
    /// Division steps across the v2 certificate.
    TotalDivisionSteps,
    /// Work units charged by the v2 verifier.
    WorkUnits,
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
            Cap::PoolMonomials => write!(f, "monomials in the pool"),
            Cap::PoolEntries => write!(f, "entries in the monomial pool"),
            Cap::InputPolynomials => write!(f, "input polynomials"),
            Cap::Nodes => write!(f, "trace nodes"),
            Cap::CombSteps => write!(f, "steps in one combination"),
            Cap::TraceSteps => write!(f, "combination steps in the trace"),
            Cap::BasisElements => write!(f, "basis elements"),
            Cap::Pairs => write!(f, "basis pairs"),
            Cap::DivisionSteps => write!(f, "steps in one division trace"),
            Cap::TotalDivisionSteps => write!(f, "division steps in the certificate"),
            Cap::WorkUnits => write!(f, "verifier work units"),
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
    /// The first byte selects no supported certificate format.
    Format {
        /// The first byte, or `None` for empty input.
        first: Option<u8>,
    },
    /// The bytes are not a canonical v2 binary encoding.
    MalformedBinary {
        /// The binary syntax fault.
        reason: BinaryFault,
        /// The byte offset of the fault.
        offset: usize,
    },
    /// A v2 monomial pool entry breaks an encoding rule.
    Pool {
        /// The pool entry index.
        index: usize,
        /// The rule the entry breaks.
        fault: PoolFault,
    },
    /// A v2 trace node breaks an encoding rule.
    Trace {
        /// The node id.
        node: usize,
        /// The rule the node breaks.
        fault: NodeFault,
    },
    /// A division trace does not reduce its polynomial to zero.
    Division {
        /// The trace's place in the certificate.
        at: DivisionSite,
        /// The failing step, or the step count for a nonzero residual.
        step: usize,
        /// The rule the trace breaks.
        fault: DivisionFault,
    },
    /// A pair witness does not justify its pair.
    Witness {
        /// The first basis index.
        i: usize,
        /// The second basis index.
        j: usize,
        /// The rule the witness breaks.
        fault: WitnessFault,
    },
    /// A monomial formed by the verifier exceeds the exponent bound.
    ExponentOverflow {
        /// The largest permitted exponent.
        max: u64,
    },
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
    /// The variable count is outside the contract range.
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
    /// Entries contain no duplicate pair.
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
    /// An index does not address an element of the collection it names.
    IndexOutOfRange {
        /// The name of the index.
        what: &'static str,
        /// The value of the index.
        index: u64,
        /// The number of elements the index must address.
        bound: usize,
    },
    /// An array holds a different number of entries than the contract
    /// requires.
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
            VerifyError::Format { first } => match first {
                Some(byte) => write!(
                    f,
                    "the bytes start with {byte:#04x}, which selects no certificate format"
                ),
                None => f.write_str("the bytes are empty"),
            },
            VerifyError::MalformedBinary { reason, offset } => {
                write!(f, "malformed certificate at byte {offset}: {reason}")
            }
            VerifyError::Pool { index, fault } => write!(f, "pool entry {index}: {fault}"),
            VerifyError::Trace { node, fault } => write!(f, "trace node {node}: {fault}"),
            VerifyError::Division { at, step, fault } => {
                write!(f, "{at}, step {step}: {fault}")
            }
            VerifyError::Witness { i, j, fault } => {
                write!(f, "the witness of pair ({i},{j}): {fault}")
            }
            VerifyError::ExponentOverflow { max } => {
                write!(f, "a monomial passes the exponent bound {max}")
            }
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
