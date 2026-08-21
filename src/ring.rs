//! The polynomial ring and every construction that goes through it.

pub(crate) mod field;

use std::fmt;
use std::sync::Arc;

use crate::ideal::Ideal;
use crate::poly::{Monomial, Polynomial, Term};
use field::{Felt, is_prime};

/// The largest number of variables a ring accepts.
///
/// The certificate contract caps `nvars` at 256, so a ring past that bound
/// could never certify. `PolynomialRing::prime_field` shares the cap.
const MAX_VARIABLES: usize = 256;

/// The largest modulus a ring accepts.
///
/// The certificate contract caps `p` at 2^31 - 1, so a ring past that bound
/// could never certify. `PolynomialRing::prime_field` shares the cap. Below
/// this bound the product of two representatives fits a `u64`, so field
/// multiplication needs no wider integer.
const MAX_MODULUS: u64 = (1 << 31) - 1;

#[derive(Debug, PartialEq, Eq)]
struct RingInner {
    modulus: u64,
    variables: Vec<String>,
}

/// A polynomial ring over a prime field, under the grevlex order.
///
/// The ring owns the modulus, the variable names, and the variable order.
/// Every polynomial is built through the ring, so every polynomial stores
/// its ring, holds reduced coefficients, and holds one exponent per
/// variable. Cloning a ring shares one allocation.
///
/// Two rings are equal when they name the same prime and the same
/// variables in the same order. Two separate constructions of one ring are
/// therefore equal, and polynomials of either belong to both.
///
/// The monomial order is grevlex over the variable order given at
/// construction, with the first variable the largest. Certificates name it
/// `grevlex-v1`. There is no other order in this release.
///
/// ```
/// use sylvester::PolynomialRing;
///
/// let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
/// assert_eq!(ring.modulus(), 32003);
/// assert_eq!(ring.nvars(), 3);
/// # Ok::<(), sylvester::RingError>(())
/// ```
#[derive(Clone, Debug)]
pub struct PolynomialRing {
    inner: Arc<RingInner>,
}

impl PartialEq for PolynomialRing {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner) || self.inner == other.inner
    }
}

impl Eq for PolynomialRing {}

impl PolynomialRing {
    /// Build the ring `F_p[variables]`.
    ///
    /// `modulus` must be a prime of at most 2^31 - 1, and the ring accepts
    /// at most 256 variables. Both bounds match the `sylv-gb-cert-v1`
    /// certificate contract, so every ring this constructor builds can
    /// certify. Variable names must be distinct. Each name starts with an
    /// ASCII letter or an underscore and holds only ASCII letters, digits,
    /// and underscores. The order of the names is the variable order.
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::prime_field(7, ["x", "y"])?;
    /// assert_eq!(ring.variables(), ["x", "y"]);
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    pub fn prime_field<I, S>(modulus: u64, variables: I) -> Result<Self, RingError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        if modulus > MAX_MODULUS {
            return Err(RingError::ModulusTooLarge { modulus });
        }
        if !is_prime(modulus) {
            return Err(RingError::ModulusNotPrime { modulus });
        }

        let variables: Vec<String> = variables
            .into_iter()
            .map(|name| name.as_ref().to_string())
            .collect();
        if variables.len() > MAX_VARIABLES {
            return Err(RingError::TooManyVariables {
                count: variables.len(),
            });
        }
        for (index, name) in variables.iter().enumerate() {
            if !is_variable_name(name) {
                return Err(RingError::InvalidVariableName { name: name.clone() });
            }
            if variables[..index].contains(name) {
                return Err(RingError::DuplicateVariable { name: name.clone() });
            }
        }

        Ok(PolynomialRing {
            inner: Arc::new(RingInner { modulus, variables }),
        })
    }

    /// The prime the coefficients live in.
    pub fn modulus(&self) -> u64 {
        self.inner.modulus
    }

    /// The variable names, largest variable first.
    pub fn variables(&self) -> &[String] {
        &self.inner.variables
    }

    /// The number of variables.
    pub fn nvars(&self) -> usize {
        self.inner.variables.len()
    }

    /// The position of a variable, or `None` when the ring has no such
    /// variable.
    pub fn variable_index(&self, name: &str) -> Option<usize> {
        self.inner.variables.iter().position(|held| held == name)
    }

    /// The zero polynomial of this ring.
    pub fn zero(&self) -> Polynomial {
        Polynomial::from_sorted_terms(self.clone(), Vec::new())
    }

    /// The constant polynomial 1 of this ring.
    pub fn one(&self) -> Polynomial {
        Polynomial::from_sorted_terms(
            self.clone(),
            vec![Term {
                coeff: Felt::one(),
                mono: Monomial::one(self.nvars()),
            }],
        )
    }

    /// Build a polynomial from coefficient and exponent pairs.
    ///
    /// Every exponent vector must hold one entry per variable. Repeated
    /// monomials add up, and a term that reduces to zero drops out.
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::prime_field(7, ["x", "y"])?;
    /// let f = ring.polynomial([(1, [2, 0]), (-1, [0, 0])])?;
    /// assert_eq!(f.to_string(), "x^2 + 6");
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    pub fn polynomial<I, E>(&self, terms: I) -> Result<Polynomial, RingError>
    where
        I: IntoIterator<Item = (i64, E)>,
        E: AsRef<[u16]>,
    {
        let modulus = self.modulus();
        let mut collected = Vec::new();
        for (coeff, exps) in terms {
            let exps = exps.as_ref();
            if exps.len() != self.nvars() {
                return Err(RingError::ExponentCount {
                    found: exps.len(),
                    expected: self.nvars(),
                });
            }
            collected.push(Term {
                coeff: Felt::new(coeff, modulus),
                mono: Monomial::from_exps(exps.iter().copied().collect()),
            });
        }
        Ok(Polynomial::from_terms(self.clone(), collected))
    }

    /// Read a polynomial from text.
    ///
    /// The grammar is
    ///
    /// ```text
    /// polynomial := [ sign ] term { sign term }
    /// sign       := "+" | "-"
    /// term       := integer [ "*" powers ] | powers
    /// powers     := power { "*" power }
    /// power      := variable [ "^" integer ]
    /// ```
    ///
    /// A variable is a name of this ring. Integers are decimal and have no
    /// sign of their own; the sign in front of a term carries it.
    /// Whitespace between tokens is ignored. Coefficients are reduced
    /// modulo the ring's prime, so `-1` and `p - 1` name one polynomial.
    /// [`Polynomial`]'s `Display` writes this syntax back.
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
    /// let f = ring.parse_polynomial("x^2*y - 3*z + 1")?;
    /// assert_eq!(f, ring.polynomial([(1, [2, 1, 0]), (-3, [0, 0, 1]), (1, [0, 0, 0])])?);
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    pub fn parse_polynomial(&self, text: &str) -> Result<Polynomial, RingError> {
        Parser::new(self, text)
            .polynomial()
            .map_err(RingError::Parse)
    }

    /// Build the ideal these polynomials generate.
    ///
    /// Every polynomial must belong to this ring.
    pub fn ideal<I>(&self, generators: I) -> Result<Ideal, RingError>
    where
        I: IntoIterator<Item = Polynomial>,
    {
        let generators: Vec<Polynomial> = generators.into_iter().collect();
        for (index, generator) in generators.iter().enumerate() {
            if generator.ring() != self {
                return Err(RingError::ForeignPolynomial { index });
            }
        }
        Ok(Ideal::new(self.clone(), generators))
    }
}

/// Report whether the parser can read `name` back as a variable.
///
/// The syntax is ASCII, so a name outside ASCII has no text form.
fn is_variable_name(name: &str) -> bool {
    let mut chars = name.chars();
    match chars.next() {
        Some(first) if first.is_ascii_alphabetic() || first == '_' => {}
        _ => return false,
    }
    chars.all(|c| c.is_ascii_alphanumeric() || c == '_')
}

/// Why a construction through the ring stops.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RingError {
    /// The modulus is not prime, so the coefficients have no field.
    ModulusNotPrime {
        /// The modulus the caller gave.
        modulus: u64,
    },
    /// The modulus is above 2^31 - 1, the largest prime the certificate
    /// contract accepts.
    ModulusTooLarge {
        /// The modulus the caller gave.
        modulus: u64,
    },
    /// The ring holds more variables than the crate supports.
    TooManyVariables {
        /// The number of names the caller gave.
        count: usize,
    },
    /// A variable name is empty or holds a character the syntax reserves.
    InvalidVariableName {
        /// The name the caller gave.
        name: String,
    },
    /// Two variables share a name.
    DuplicateVariable {
        /// The repeated name.
        name: String,
    },
    /// An exponent vector does not hold one entry per variable.
    ExponentCount {
        /// The width the caller gave.
        found: usize,
        /// The number of variables in the ring.
        expected: usize,
    },
    /// A polynomial of another ring reached a ring operation.
    ForeignPolynomial {
        /// The position of the polynomial in the argument.
        index: usize,
    },
    /// The text does not parse.
    Parse(ParseError),
}

impl fmt::Display for RingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RingError::ModulusNotPrime { modulus } => {
                write!(f, "the modulus {modulus} is not prime")
            }
            RingError::ModulusTooLarge { modulus } => {
                write!(f, "the modulus {modulus} is above 2^31 - 1")
            }
            RingError::TooManyVariables { count } => write!(
                f,
                "the ring names {count} variables, the largest supported count is {MAX_VARIABLES}"
            ),
            RingError::InvalidVariableName { name } => {
                write!(f, "\"{name}\" is not a variable name")
            }
            RingError::DuplicateVariable { name } => {
                write!(f, "the ring names the variable \"{name}\" twice")
            }
            RingError::ExponentCount { found, expected } => write!(
                f,
                "the exponent vector holds {found} entries, the ring has {expected}"
            ),
            RingError::ForeignPolynomial { index } => {
                write!(f, "polynomial {index} belongs to another ring")
            }
            RingError::Parse(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for RingError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            RingError::Parse(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ParseError> for RingError {
    fn from(error: ParseError) -> Self {
        RingError::Parse(error)
    }
}

/// Why [`PolynomialRing::parse_polynomial`] stops.
///
/// Every position is a byte offset into the text the parser read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ParseError {
    /// The text names something the ring has no variable for.
    UnknownVariable {
        /// The name the text holds.
        name: String,
        /// Where the name starts.
        position: usize,
    },
    /// The syntax has no place for this character here.
    UnexpectedCharacter {
        /// The character the text holds.
        character: char,
        /// Where the character is.
        position: usize,
    },
    /// The text stops in the middle of a term.
    UnexpectedEnd {
        /// The end of the text.
        position: usize,
    },
    /// An exponent is above 65535.
    ExponentTooLarge {
        /// Where the exponent starts.
        position: usize,
    },
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ParseError::UnknownVariable { name, position } => {
                write!(f, "the ring has no variable \"{name}\" at byte {position}")
            }
            ParseError::UnexpectedCharacter {
                character,
                position,
            } => write!(f, "unexpected character '{character}' at byte {position}"),
            ParseError::UnexpectedEnd { position } => {
                write!(f, "the text ends inside a term at byte {position}")
            }
            ParseError::ExponentTooLarge { position } => {
                write!(f, "the exponent at byte {position} is above 65535")
            }
        }
    }
}

impl std::error::Error for ParseError {}

struct Parser<'a> {
    ring: &'a PolynomialRing,
    text: &'a [u8],
    at: usize,
}

impl<'a> Parser<'a> {
    fn new(ring: &'a PolynomialRing, text: &'a str) -> Self {
        Parser {
            ring,
            text: text.as_bytes(),
            at: 0,
        }
    }

    fn polynomial(mut self) -> Result<Polynomial, ParseError> {
        let mut terms = Vec::new();
        let mut negative = self.sign()?;
        loop {
            terms.push(self.term(negative)?);
            self.skip_spaces();
            match self.peek() {
                None => break,
                Some(b'+') => {
                    self.at += 1;
                    negative = false;
                }
                Some(b'-') => {
                    self.at += 1;
                    negative = true;
                }
                Some(byte) => return Err(self.unexpected(byte)),
            }
        }
        Ok(Polynomial::from_terms(self.ring.clone(), terms))
    }

    fn sign(&mut self) -> Result<bool, ParseError> {
        self.skip_spaces();
        match self.peek() {
            Some(b'+') => {
                self.at += 1;
                Ok(false)
            }
            Some(b'-') => {
                self.at += 1;
                Ok(true)
            }
            Some(_) => Ok(false),
            None => Err(ParseError::UnexpectedEnd { position: self.at }),
        }
    }

    fn term(&mut self, negative: bool) -> Result<Term, ParseError> {
        let modulus = self.ring.modulus();
        let mut exps = vec![0u16; self.ring.nvars()];
        self.skip_spaces();

        let mut coeff = 1u64 % modulus;
        let mut expect_power = true;
        if matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            coeff = self.coefficient(modulus);
            self.skip_spaces();
            match self.peek() {
                Some(b'*') => self.at += 1,
                _ => expect_power = false,
            }
        }

        if expect_power {
            loop {
                self.power(&mut exps)?;
                self.skip_spaces();
                if self.peek() == Some(b'*') {
                    self.at += 1;
                    continue;
                }
                break;
            }
        }

        let coeff = if negative && coeff != 0 {
            modulus - coeff
        } else {
            coeff
        };
        Ok(Term {
            coeff: Felt::from_residue(coeff),
            mono: Monomial::from_exps(exps.into_iter().collect()),
        })
    }

    /// Read a decimal integer and reduce it as it goes, so text of any
    /// length parses without overflow.
    fn coefficient(&mut self, modulus: u64) -> u64 {
        let mut value = 0u128;
        while let Some(byte) = self.peek() {
            if !byte.is_ascii_digit() {
                break;
            }
            value = (value * 10 + (byte - b'0') as u128) % modulus as u128;
            self.at += 1;
        }
        value as u64
    }

    fn power(&mut self, exps: &mut [u16]) -> Result<(), ParseError> {
        self.skip_spaces();
        let start = self.at;
        let index = match self.peek() {
            Some(byte) if byte.is_ascii_alphabetic() || byte == b'_' => {
                while matches!(self.peek(), Some(byte) if byte.is_ascii_alphanumeric() || byte == b'_')
                {
                    self.at += 1;
                }
                let name = String::from_utf8_lossy(&self.text[start..self.at]).into_owned();
                match self.ring.variable_index(&name) {
                    Some(index) => index,
                    None => {
                        return Err(ParseError::UnknownVariable {
                            name,
                            position: start,
                        });
                    }
                }
            }
            Some(byte) => return Err(self.unexpected(byte)),
            None => return Err(ParseError::UnexpectedEnd { position: self.at }),
        };

        self.skip_spaces();
        let mut exponent = 1u64;
        if self.peek() == Some(b'^') {
            self.at += 1;
            self.skip_spaces();
            let digits = self.at;
            if !matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
                return match self.peek() {
                    Some(byte) => Err(self.unexpected(byte)),
                    None => Err(ParseError::UnexpectedEnd { position: self.at }),
                };
            }
            exponent = 0;
            while let Some(byte) = self.peek() {
                if !byte.is_ascii_digit() {
                    break;
                }
                exponent = exponent
                    .saturating_mul(10)
                    .saturating_add((byte - b'0') as u64);
                self.at += 1;
            }
            if exponent > u16::MAX as u64 {
                return Err(ParseError::ExponentTooLarge { position: digits });
            }
        }

        let slot = &mut exps[index];
        let sum = *slot as u64 + exponent;
        if sum > u16::MAX as u64 {
            return Err(ParseError::ExponentTooLarge { position: start });
        }
        *slot = sum as u16;
        Ok(())
    }

    fn peek(&self) -> Option<u8> {
        self.text.get(self.at).copied()
    }

    fn skip_spaces(&mut self) {
        while matches!(self.peek(), Some(byte) if byte.is_ascii_whitespace()) {
            self.at += 1;
        }
    }

    fn unexpected(&self, byte: u8) -> ParseError {
        let rest = String::from_utf8_lossy(&self.text[self.at..]).into_owned();
        ParseError::UnexpectedCharacter {
            character: rest.chars().next().unwrap_or(byte as char),
            position: self.at,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y", "z"]).expect("7 is prime")
    }

    #[test]
    fn a_non_prime_modulus_has_no_ring() {
        assert_eq!(
            PolynomialRing::prime_field(4, ["x"]),
            Err(RingError::ModulusNotPrime { modulus: 4 })
        );
        assert_eq!(
            PolynomialRing::prime_field(1, ["x"]),
            Err(RingError::ModulusNotPrime { modulus: 1 })
        );
    }

    #[test]
    fn a_modulus_above_the_bound_has_no_ring() {
        let modulus = (1u64 << 31) - 1 + 2;
        assert_eq!(
            PolynomialRing::prime_field(modulus, ["x"]),
            Err(RingError::ModulusTooLarge { modulus })
        );
    }

    #[test]
    fn a_variable_count_above_the_bound_has_no_ring() {
        let names: Vec<String> = (0..257).map(|i| format!("v{i}")).collect();
        assert_eq!(
            PolynomialRing::prime_field(7, names),
            Err(RingError::TooManyVariables { count: 257 })
        );
    }

    #[test]
    fn variable_names_must_be_distinct_and_well_formed() {
        assert_eq!(
            PolynomialRing::prime_field(7, ["x", "x"]),
            Err(RingError::DuplicateVariable {
                name: "x".to_string()
            })
        );
        assert_eq!(
            PolynomialRing::prime_field(7, ["x y"]),
            Err(RingError::InvalidVariableName {
                name: "x y".to_string()
            })
        );
        assert_eq!(
            PolynomialRing::prime_field(7, [""]),
            Err(RingError::InvalidVariableName {
                name: String::new()
            })
        );
    }

    #[test]
    fn a_ring_with_no_variables_holds_the_field() {
        let ring = PolynomialRing::prime_field(7, Vec::<String>::new()).expect("7 is prime");
        assert_eq!(ring.nvars(), 0);
        assert_eq!(ring.one().to_string(), "1");
    }

    #[test]
    fn rings_with_the_same_shape_are_the_same_ring() {
        assert_eq!(ring(), ring());
        assert_ne!(
            ring(),
            PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
        );
        assert_ne!(
            ring(),
            PolynomialRing::prime_field(11, ["x", "y", "z"]).expect("11 is prime")
        );
    }

    #[test]
    fn an_exponent_vector_must_match_the_variable_count() {
        assert_eq!(
            ring().polynomial([(1, [1, 0])]),
            Err(RingError::ExponentCount {
                found: 2,
                expected: 3
            })
        );
    }

    #[test]
    fn an_ideal_takes_only_polynomials_of_its_ring() {
        let other = PolynomialRing::prime_field(11, ["x", "y", "z"]).expect("11 is prime");
        let foreign = other.one();
        assert_eq!(
            ring().ideal([ring().one(), foreign]),
            Err(RingError::ForeignPolynomial { index: 1 })
        );
    }

    #[test]
    fn parsing_reads_the_documented_grammar() {
        let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"]).expect("32003 is prime");
        let parsed = ring
            .parse_polynomial("x^2*y - 3*z + 1")
            .expect("the text parses");
        let built = ring
            .polynomial([(1, [2, 1, 0]), (-3, [0, 0, 1]), (1, [0, 0, 0])])
            .expect("the terms fit the ring");
        assert_eq!(parsed, built);
    }

    #[test]
    fn parsing_ignores_whitespace_and_reads_a_leading_sign() {
        let ring = ring();
        assert_eq!(
            ring.parse_polynomial("  -  2 * x ^ 2 * y  ")
                .expect("parses"),
            ring.polynomial([(-2, [2, 1, 0])]).expect("fits")
        );
    }

    #[test]
    fn parsing_adds_repeated_variables_and_monomials() {
        let ring = ring();
        assert_eq!(
            ring.parse_polynomial("x*x + x^2").expect("parses"),
            ring.polynomial([(2, [2, 0, 0])]).expect("fits")
        );
    }

    #[test]
    fn a_term_that_reduces_to_zero_drops_out() {
        let ring = ring();
        assert!(ring.parse_polynomial("7*x").expect("parses").is_zero());
        assert!(ring.parse_polynomial("0").expect("parses").is_zero());
        assert!(ring.parse_polynomial("x - x").expect("parses").is_zero());
    }

    #[test]
    fn parsing_reports_where_it_stops() {
        let ring = ring();
        assert_eq!(
            ring.parse_polynomial("x + w"),
            Err(RingError::Parse(ParseError::UnknownVariable {
                name: "w".to_string(),
                position: 4
            }))
        );
        assert_eq!(
            ring.parse_polynomial("x $ y"),
            Err(RingError::Parse(ParseError::UnexpectedCharacter {
                character: '$',
                position: 2
            }))
        );
        assert_eq!(
            ring.parse_polynomial("x +"),
            Err(RingError::Parse(ParseError::UnexpectedEnd { position: 3 }))
        );
        assert_eq!(
            ring.parse_polynomial(""),
            Err(RingError::Parse(ParseError::UnexpectedEnd { position: 0 }))
        );
        assert_eq!(
            ring.parse_polynomial("x^70000"),
            Err(RingError::Parse(ParseError::ExponentTooLarge {
                position: 2
            }))
        );
    }

    #[test]
    fn a_long_coefficient_parses_without_overflow() {
        let ring = ring();
        let digits = "9".repeat(64);
        let expected = digits
            .bytes()
            .fold(0u64, |acc, byte| (acc * 10 + (byte - b'0') as u64) % 7);
        assert_eq!(
            ring.parse_polynomial(&digits).expect("parses"),
            ring.polynomial([(expected as i64, [0, 0, 0])])
                .expect("fits")
        );
    }
}
