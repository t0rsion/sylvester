//! The polynomial ring and every construction that goes through it.

pub(crate) mod field;
mod parser;
pub(crate) mod rational;

use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use num_bigint::BigInt;
use num_integer::Integer;
use num_rational::BigRational;
use num_traits::{One, Signed, Zero};

use crate::ideal::Ideal;
use crate::poly::{Monomial, Polynomial, Term};
use field::is_prime;
pub use field::{Felt, PrimeOps};
pub use parser::ExpressionError;
pub use rational::{Established, ModularLift, RationalMeta, RationalOps};

pub(crate) mod sealed {
    /// The bound that keeps [`super::Domain`] and [`super::DomainOps`]
    /// closed to the domains this crate defines.
    pub trait Sealed {}
}

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

/// The coefficient domain of a ring.
///
/// The trait is sealed to the two domains of this crate, [`PrimeField`]
/// and [`Rationals`]. A ring, a polynomial, an ideal, and a basis all
/// carry the domain as a type parameter, so an operation one domain does
/// not offer is a compile error rather than a runtime one.
pub trait Domain: sealed::Sealed + Clone + fmt::Debug + Eq + Hash + Send + Sync + 'static {
    /// One coefficient of the domain.
    type Coeff: Clone + fmt::Debug + Eq + Hash + Send + Sync;
    /// The arithmetic the shared polynomial code calls.
    type Ops: DomainOps<Coeff = Self::Coeff> + Clone + fmt::Debug + Eq + Hash + Send + Sync;
    /// What a basis of this domain records about its own origin.
    ///
    /// It is `()` over [`PrimeField`], which records nothing, and
    /// [`RationalMeta`] over [`Rationals`].
    type BasisMeta: Clone + fmt::Debug + Eq + Send + Sync;
}

/// The domain of the prime field `F_p`.
///
/// It is the default type parameter of [`PolynomialRing`],
/// [`Polynomial`], [`Ideal`], and [`GroebnerBasis`], so a prime-field
/// signature names no domain.
///
/// [`GroebnerBasis`]: crate::GroebnerBasis
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct PrimeField;

impl sealed::Sealed for PrimeField {}

/// The domain of the rational numbers `Q`.
///
/// A coefficient is exact: no modulus reduces it and no rounding touches
/// it. There is no certified path over `Q`, so
/// [`Ideal::groebner_basis_certified`] is a method of `Ideal<PrimeField>`
/// alone.
///
/// [`Ideal::groebner_basis_certified`]: crate::Ideal::groebner_basis_certified
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Rationals;

impl sealed::Sealed for Rationals {}

/// The coefficient arithmetic of one domain.
///
/// The ring holds one value of this type and hands it to every polynomial
/// operation. No stored coefficient is zero, so [`DomainOps::add`],
/// [`DomainOps::sub`], and [`DomainOps::convert`] report a zero result as
/// `None` instead of returning it.
///
/// The trait is sealed. [`PrimeOps`] and [`RationalOps`] are its two
/// implementations.
pub trait DomainOps: sealed::Sealed {
    /// One coefficient of the domain.
    type Coeff;

    /// The coefficient 1.
    fn one(&self) -> Self::Coeff;
    /// Report whether the coefficient is 1.
    fn is_one(&self, a: &Self::Coeff) -> bool;
    /// Add, or report a sum of zero.
    fn add(&self, a: &Self::Coeff, b: &Self::Coeff) -> Option<Self::Coeff>;
    /// Subtract, or report a difference of zero.
    fn sub(&self, a: &Self::Coeff, b: &Self::Coeff) -> Option<Self::Coeff>;
    /// Negate.
    fn neg(&self, a: &Self::Coeff) -> Self::Coeff;
    /// Multiply. Both factors are nonzero, so the product is nonzero.
    fn mul(&self, a: &Self::Coeff, b: &Self::Coeff) -> Self::Coeff;
    /// Invert. The caller passes a nonzero coefficient.
    fn inv(&self, a: &Self::Coeff) -> Self::Coeff;
    /// Read a caller's value into the domain, or report a value of zero.
    fn convert(&self, value: &Coefficient) -> Result<Option<Self::Coeff>, RingError>;
    /// Write the coefficient in the text syntax the parser reads.
    fn write(&self, a: &Self::Coeff, f: &mut fmt::Formatter<'_>) -> fmt::Result;
    /// Report whether `Display` writes the coefficient behind a minus
    /// sign.
    fn is_negative(&self, a: &Self::Coeff) -> bool;
    /// The heap bytes one coefficient holds.
    fn heap_bytes(&self, a: &Self::Coeff) -> usize;
}

/// A coefficient a caller gives to [`PolynomialRing::polynomial`].
///
/// One type serves every domain, and the ring converts it through
/// [`DomainOps::convert`]. Every ordinary integer type converts into it,
/// so `ring.polynomial([(1, [2, 0])])` needs no annotation.
///
/// [`Coefficient::Fraction`] is the raw pair, unreduced: `From` cannot
/// fail, so a zero denominator is [`RingError::ZeroDenominator`] at
/// conversion and never a panic.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Coefficient {
    /// An integer that fits an `i64`.
    Small(i64),
    /// An integer of any width.
    Integer(BigInt),
    /// A quotient of two integers, in any form.
    Fraction {
        /// The numerator.
        numerator: BigInt,
        /// The denominator. Zero is an error at conversion.
        denominator: BigInt,
    },
    /// A rational number in lowest terms.
    Rational(BigRational),
}

impl Coefficient {
    /// The value as a fraction in lowest terms with a positive
    /// denominator.
    ///
    /// A zero denominator is [`RingError::ZeroDenominator`]. Every domain
    /// reduces the fraction before it reads it, so `3/3` is the value 1 in
    /// every domain that holds 1.
    pub(crate) fn normalized(&self) -> Result<(BigInt, BigInt), RingError> {
        let (numerator, denominator) = match self {
            Coefficient::Small(value) => (BigInt::from(*value), BigInt::one()),
            Coefficient::Integer(value) => (value.clone(), BigInt::one()),
            Coefficient::Fraction {
                numerator,
                denominator,
            } => (numerator.clone(), denominator.clone()),
            Coefficient::Rational(value) => (value.numer().clone(), value.denom().clone()),
        };
        if denominator.is_zero() {
            return Err(RingError::ZeroDenominator);
        }
        let divisor = numerator.gcd(&denominator);
        let mut numerator = numerator / &divisor;
        let mut denominator = denominator / &divisor;
        if denominator.is_negative() {
            numerator = -numerator;
            denominator = -denominator;
        }
        Ok((numerator, denominator))
    }
}

macro_rules! coefficient_from_small {
    ($($type:ty),*) => {
        $(
            impl From<$type> for Coefficient {
                fn from(value: $type) -> Self {
                    Coefficient::Small(value as i64)
                }
            }
        )*
    };
}

coefficient_from_small!(i8, i16, i32, i64, u8, u16, u32);

macro_rules! coefficient_from_wide {
    ($($type:ty),*) => {
        $(
            impl From<$type> for Coefficient {
                fn from(value: $type) -> Self {
                    match i64::try_from(value) {
                        Ok(small) => Coefficient::Small(small),
                        Err(_) => Coefficient::Integer(BigInt::from(value)),
                    }
                }
            }
        )*
    };
}

coefficient_from_wide!(i128, isize, u64, u128, usize);

impl From<BigInt> for Coefficient {
    fn from(value: BigInt) -> Self {
        Coefficient::Integer(value)
    }
}

impl From<BigRational> for Coefficient {
    fn from(value: BigRational) -> Self {
        Coefficient::Rational(value)
    }
}

#[derive(Debug, PartialEq, Eq)]
struct RingInner<D: Domain> {
    ops: D::Ops,
    variables: Vec<String>,
}

/// A polynomial ring over the domain `D`, under the grevlex order.
///
/// The ring owns the domain arithmetic, the variable names, and the
/// variable order. Every polynomial is built through the ring, so every
/// polynomial stores its ring, holds reduced coefficients, and holds one
/// exponent per variable. Cloning a ring shares one allocation.
///
/// The domain parameter defaults to [`PrimeField`], so
/// `PolynomialRing` alone names a ring over `F_p`.
///
/// Two rings are equal when they hold the same domain value (the prime)
/// and the same variables in the same order. Two separate constructions of
/// one ring are therefore equal, hash equally, and polynomials of either
/// belong to both.
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
pub struct PolynomialRing<D: Domain = PrimeField> {
    inner: Arc<RingInner<D>>,
}

impl<D: Domain> PartialEq for PolynomialRing<D> {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.inner, &other.inner) || self.inner == other.inner
    }
}

impl<D: Domain> Eq for PolynomialRing<D> {}

impl<D: Domain> Hash for PolynomialRing<D> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.inner.ops.hash(state);
        self.inner.variables.hash(state);
    }
}

impl PolynomialRing<PrimeField> {
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

        Ok(PolynomialRing {
            inner: Arc::new(RingInner {
                ops: PrimeOps::new(modulus),
                variables: variable_names(variables)?,
            }),
        })
    }

    /// The prime the coefficients live in.
    pub fn modulus(&self) -> u64 {
        self.inner.ops.modulus()
    }
}

impl PolynomialRing<Rationals> {
    /// Build the ring `Q[variables]`.
    ///
    /// The ring accepts at most 256 variables. Variable names must be
    /// distinct. Each name starts with an ASCII letter or an underscore and
    /// holds only ASCII letters, digits, and underscores. The order of the
    /// names is the variable order. The name rules and the variable cap are
    /// the ones [`PolynomialRing::prime_field`] applies.
    ///
    /// A rational ring cannot certify: [`Ideal::groebner_basis_certified`]
    /// is a method of `Ideal<PrimeField>` alone.
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::rationals(["x", "y"])?;
    /// let f = ring.parse_polynomial("1/2*x - y")?;
    /// assert_eq!(f.to_string(), "1/2*x - y");
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    ///
    /// [`Ideal::groebner_basis_certified`]: crate::Ideal::groebner_basis_certified
    pub fn rationals<I, S>(variables: I) -> Result<Self, RingError>
    where
        I: IntoIterator<Item = S>,
        S: AsRef<str>,
    {
        Ok(PolynomialRing {
            inner: Arc::new(RingInner {
                ops: RationalOps::new(),
                variables: variable_names(variables)?,
            }),
        })
    }
}

impl<D: Domain> PolynomialRing<D> {
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

    /// Return the variables as degree-one polynomials in ring order.
    ///
    /// The first item is the polynomial for `variables()[0]`. An empty ring
    /// has no generators. Each call builds a fresh vector of polynomials.
    pub fn generators(&self) -> Vec<Polynomial<D>> {
        (0..self.nvars())
            .map(|index| {
                self.generator(index)
                    .expect("a generated index is inside the ring")
            })
            .collect()
    }

    /// Return the degree-one polynomial for one variable.
    ///
    /// `index` follows [`PolynomialRing::variables`]. An index past the ring
    /// width returns `None`.
    pub fn generator(&self, index: usize) -> Option<Polynomial<D>> {
        if index >= self.nvars() {
            return None;
        }
        let mut exps = vec![0u16; self.nvars()];
        exps[index] = 1;
        Some(Polynomial::from_sorted_terms(
            self.clone(),
            vec![Term {
                coeff: self.ops().one(),
                mono: Monomial::from_exps(exps.into_iter().collect()),
            }],
        ))
    }

    /// The domain arithmetic every polynomial operation of this ring
    /// takes.
    pub(crate) fn ops(&self) -> &D::Ops {
        &self.inner.ops
    }

    /// The zero polynomial of this ring.
    pub fn zero(&self) -> Polynomial<D> {
        Polynomial::from_sorted_terms(self.clone(), Vec::new())
    }

    /// The constant polynomial 1 of this ring.
    pub fn one(&self) -> Polynomial<D> {
        Polynomial::from_sorted_terms(
            self.clone(),
            vec![Term {
                coeff: self.ops().one(),
                mono: Monomial::one(self.nvars()),
            }],
        )
    }

    /// Read a polynomial from text.
    ///
    /// The grammar is
    ///
    /// ```text
    /// polynomial  := [ sign ] term { sign term }
    /// sign        := "+" | "-"
    /// term        := coefficient [ "*" powers ] | powers
    /// coefficient := integer [ "/" integer ]
    /// powers      := power { "*" power }
    /// power       := variable [ "^" integer ]
    /// ```
    ///
    /// A variable is a name of this ring. Integers are decimal and have no
    /// sign of their own; the sign in front of a term carries it. A
    /// coefficient may name a fraction, and the `/` takes no whitespace
    /// around it. Whitespace between the other tokens is ignored.
    /// [`Polynomial`]'s `Display` writes this syntax back.
    ///
    /// The domain reads the coefficient. Over `F_p` a coefficient is
    /// reduced modulo the prime, so `-1` and `p - 1` name one polynomial,
    /// and a denominator the prime divides is
    /// [`RingError::CoefficientNotInvertible`]. Over `Q` the value is kept
    /// exactly.
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::prime_field(32003, ["x", "y", "z"])?;
    /// let f = ring.parse_polynomial("x^2*y - 3*z + 1")?;
    /// assert_eq!(f, ring.polynomial([(1, [2, 1, 0]), (-3, [0, 0, 1]), (1, [0, 0, 0])])?);
    ///
    /// let ring = PolynomialRing::rationals(["x", "y"])?;
    /// let g = ring.parse_polynomial("2/3*x^2 - y")?;
    /// assert_eq!(g.to_string(), "2/3*x^2 - y");
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    pub fn parse_polynomial(&self, text: &str) -> Result<Polynomial<D>, RingError> {
        Parser::new(self, text).polynomial()
    }

    /// Read an expression under one time and memory budget.
    ///
    /// The expression grammar accepts `+`, `-`, `*`, parentheses, unary
    /// signs, fractions of decimal integers, and nonnegative powers with
    /// `^` or `**`. Division applies only inside a numeric fraction.
    ///
    /// The budget covers the whole parse, including every intermediate
    /// polynomial. An exhausted budget returns [`ExpressionError`].
    pub fn parse_polynomial_with_budget(
        &self,
        text: &str,
        budget: crate::Budget,
    ) -> Result<Polynomial<D>, ExpressionError> {
        let limits = crate::compute::ComputeLimits::of_budget(&budget);
        self.parse_polynomial_with_limits(text, &limits)
    }

    /// Read an expression under limits already shared by a caller.
    pub(crate) fn parse_polynomial_with_limits(
        &self,
        text: &str,
        limits: &crate::compute::ComputeLimits,
    ) -> Result<Polynomial<D>, ExpressionError> {
        parser::parse_with_limits(self, text, limits)
    }

    /// Build a polynomial from coefficient and exponent pairs.
    ///
    /// Every exponent vector must hold one entry per variable. Repeated
    /// monomials add up, and a term that reduces to zero drops out. Every
    /// ordinary integer type is a coefficient, through
    /// [`Coefficient`].
    ///
    /// ```
    /// use sylvester::PolynomialRing;
    ///
    /// let ring = PolynomialRing::prime_field(7, ["x", "y"])?;
    /// let f = ring.polynomial([(1, [2, 0]), (-1, [0, 0])])?;
    /// assert_eq!(f.to_string(), "x^2 + 6");
    /// # Ok::<(), sylvester::RingError>(())
    /// ```
    pub fn polynomial<I, C, E>(&self, terms: I) -> Result<Polynomial<D>, RingError>
    where
        I: IntoIterator<Item = (C, E)>,
        C: Into<Coefficient>,
        E: AsRef<[u16]>,
    {
        let ops = self.ops();
        let mut collected = Vec::new();
        for (coeff, exps) in terms {
            let exps = exps.as_ref();
            if exps.len() != self.nvars() {
                return Err(RingError::ExponentCount {
                    found: exps.len(),
                    expected: self.nvars(),
                });
            }
            let Some(coeff) = ops.convert(&coeff.into())? else {
                continue;
            };
            collected.push(Term {
                coeff,
                mono: Monomial::from_exps(exps.iter().copied().collect()),
            });
        }
        Ok(Polynomial::from_terms(self.clone(), collected))
    }

    /// Build the ideal these polynomials generate.
    ///
    /// Every polynomial must belong to this ring.
    pub fn ideal<I>(&self, generators: I) -> Result<Ideal<D>, RingError>
    where
        I: IntoIterator<Item = Polynomial<D>>,
    {
        let generators: Vec<Polynomial<D>> = generators.into_iter().collect();
        for (index, generator) in generators.iter().enumerate() {
            if generator.ring() != self {
                return Err(RingError::ForeignPolynomial { index });
            }
        }
        Ok(Ideal::new(self.clone(), generators))
    }
}

/// Check the variable names of a ring under construction.
///
/// The rules are the same in every domain: at most 256 names, each one
/// readable by the parser, and no name twice.
fn variable_names<I, S>(variables: I) -> Result<Vec<String>, RingError>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
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
    Ok(variables)
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
    /// A coefficient names a fraction whose denominator is zero.
    ZeroDenominator,
    /// A coefficient names a fraction the domain cannot read, because the
    /// denominator has no inverse there.
    ///
    /// Over `F_p` the prime divides the denominator. The fraction is
    /// reduced first, so `3/3` over `F_3` is the value 1 and not this
    /// error.
    CoefficientNotInvertible {
        /// The denominator, reduced against the numerator.
        denominator: BigInt,
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
            RingError::ZeroDenominator => f.write_str("a coefficient has a denominator of zero"),
            RingError::CoefficientNotInvertible { denominator } => write!(
                f,
                "the denominator {denominator} has no inverse in the coefficient domain"
            ),
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
    /// A coefficient names a fraction whose denominator is zero.
    ZeroDenominator {
        /// Where the denominator starts.
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
            ParseError::ZeroDenominator { position } => {
                write!(f, "the denominator at byte {position} is zero")
            }
        }
    }
}

impl std::error::Error for ParseError {}

/// The reader of the text syntax.
///
/// One reader serves both domains. It collects each coefficient as a
/// [`Coefficient`] and the ring's domain converts it, so a syntax error is
/// the same everywhere and only the conversion differs.
struct Parser<'a, D: Domain> {
    ring: &'a PolynomialRing<D>,
    text: &'a [u8],
    at: usize,
}

impl<'a, D: Domain> Parser<'a, D> {
    fn new(ring: &'a PolynomialRing<D>, text: &'a str) -> Self {
        Parser {
            ring,
            text: text.as_bytes(),
            at: 0,
        }
    }

    fn polynomial(mut self) -> Result<Polynomial<D>, RingError> {
        let mut terms = Vec::new();
        let mut negative = self.sign()?;
        loop {
            // A term whose coefficient reads as zero drops out here,
            // because a polynomial holds no zero coefficient.
            if let Some(term) = self.term(negative)? {
                terms.push(term);
            }
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
                Some(byte) => return Err(self.unexpected(byte).into()),
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

    /// Read one term, or report that its coefficient is zero.
    fn term(&mut self, negative: bool) -> Result<Option<Term<D>>, RingError> {
        let mut exps = vec![0u16; self.ring.nvars()];
        self.skip_spaces();

        let mut value = Coefficient::Small(1);
        let mut expect_power = true;
        if matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            value = self.coefficient()?;
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

        if negative {
            value = negated(value);
        }
        let Some(coeff) = self.ring.ops().convert(&value)? else {
            return Ok(None);
        };
        Ok(Some(Term {
            coeff,
            mono: Monomial::from_exps(exps.into_iter().collect()),
        }))
    }

    /// Read an integer, or a fraction of two integers.
    ///
    /// The `/` takes no whitespace around it, so `1 / 2` is one term of
    /// value 1 followed by a character the syntax has no place for.
    fn coefficient(&mut self) -> Result<Coefficient, ParseError> {
        let numerator = self.integer();
        if self.peek() != Some(b'/') {
            return Ok(integer_coefficient(numerator));
        }
        self.at += 1;
        let position = self.at;
        if !matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            return match self.peek() {
                Some(byte) => Err(self.unexpected(byte)),
                None => Err(ParseError::UnexpectedEnd { position: self.at }),
            };
        }
        let denominator = self.integer();
        if denominator.is_zero() {
            return Err(ParseError::ZeroDenominator { position });
        }
        Ok(Coefficient::Fraction {
            numerator,
            denominator,
        })
    }

    /// Read a run of decimal digits. The caller has seen the first one.
    fn integer(&mut self) -> BigInt {
        let start = self.at;
        while matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            self.at += 1;
        }
        BigInt::parse_bytes(&self.text[start..self.at], 10)
            .expect("a run of decimal digits is an integer")
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

/// The integer as the narrowest coefficient that holds it.
///
/// [`Coefficient::Small`] is the common case and the one every domain
/// reads without a big-integer step.
fn integer_coefficient(value: BigInt) -> Coefficient {
    match i64::try_from(&value) {
        Ok(small) => Coefficient::Small(small),
        Err(_) => Coefficient::Integer(value),
    }
}

/// The coefficient with its sign flipped.
///
/// The text syntax has no sign inside a coefficient, so the sign in front
/// of the term lands in the numerator here.
fn negated(value: Coefficient) -> Coefficient {
    match value {
        Coefficient::Small(small) => match small.checked_neg() {
            Some(negated) => Coefficient::Small(negated),
            None => Coefficient::Integer(-BigInt::from(small)),
        },
        Coefficient::Integer(value) => Coefficient::Integer(-value),
        Coefficient::Fraction {
            numerator,
            denominator,
        } => Coefficient::Fraction {
            numerator: -numerator,
            denominator,
        },
        Coefficient::Rational(value) => Coefficient::Rational(-value),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y", "z"]).expect("7 is prime")
    }

    /// The unreduced fraction `numerator / denominator`.
    fn fraction(numerator: i64, denominator: i64) -> Coefficient {
        Coefficient::Fraction {
            numerator: BigInt::from(numerator),
            denominator: BigInt::from(denominator),
        }
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
    fn a_fraction_is_reduced_before_the_field_reads_it() {
        let ring = PolynomialRing::prime_field(3, ["x"]).expect("3 is prime");
        // 3/3 is the value 1, not a denominator the field cannot invert.
        let one = ring
            .polynomial([(fraction(3, 3), [0])])
            .expect("3/3 is the value 1");
        assert_eq!(one, ring.one());
        // 3/6 is 1/2, and 2 inverts modulo 3.
        assert_eq!(
            ring.polynomial([(fraction(3, 6), [0])])
                .expect("3/6 is 1/2"),
            ring.polynomial([(2, [0])]).expect("fits")
        );
        // 0/3 is zero, so the term drops out.
        assert!(
            ring.polynomial([(fraction(0, 3), [1])])
                .expect("0/3 is zero")
                .is_zero()
        );
    }

    #[test]
    fn a_denominator_the_field_cannot_invert_is_an_error() {
        let ring = PolynomialRing::prime_field(3, ["x"]).expect("3 is prime");
        assert_eq!(
            ring.polynomial([(fraction(1, 3), [0])]),
            Err(RingError::CoefficientNotInvertible {
                denominator: BigInt::from(3)
            })
        );
        assert_eq!(
            ring.polynomial([(fraction(1, 0), [0])]),
            Err(RingError::ZeroDenominator)
        );
    }

    #[test]
    fn every_integer_width_is_a_coefficient() {
        let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
        let expected = ring.polynomial([(1, [0])]).expect("fits");
        assert_eq!(ring.polynomial([(8i8, [0])]).expect("fits"), expected);
        assert_eq!(ring.polynomial([(8u32, [0])]).expect("fits"), expected);
        assert_eq!(
            ring.polynomial([(u64::MAX, [0])]).expect("fits"),
            ring.polynomial([(u64::MAX % 7, [0])]).expect("fits")
        );
        assert_eq!(
            ring.polynomial([(BigInt::from(15), [0])]).expect("fits"),
            expected
        );
        assert_eq!(
            ring.polynomial([(BigRational::new(BigInt::from(2), BigInt::from(4)), [0])])
                .expect("1/2 inverts modulo 7"),
            ring.polynomial([(4, [0])]).expect("fits")
        );
    }

    #[test]
    fn a_normalized_fraction_carries_the_sign_in_its_numerator() {
        assert_eq!(
            fraction(1, -2)
                .normalized()
                .expect("the denominator is not zero"),
            (BigInt::from(-1), BigInt::from(2))
        );
        assert_eq!(
            fraction(-4, -2)
                .normalized()
                .expect("the denominator is not zero"),
            (BigInt::from(2), BigInt::from(1))
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
    fn parsing_reads_a_fraction_over_the_rationals() {
        let ring = PolynomialRing::rationals(["x", "y"]).expect("the names are variable names");
        let parsed = ring
            .parse_polynomial("1/2*x^2 - 3*y + 2/4")
            .expect("the text parses");
        let built = ring
            .polynomial([
                (fraction(1, 2), [2, 0]),
                (Coefficient::Small(-3), [0, 1]),
                (fraction(1, 2), [0, 0]),
            ])
            .expect("the terms fit the ring");
        assert_eq!(parsed, built);
        assert_eq!(parsed.to_string(), "1/2*x^2 - 3*y + 1/2");
    }

    #[test]
    fn a_rational_coefficient_holds_every_width() {
        let ring = PolynomialRing::rationals(["x"]).expect("the name is a variable name");
        let text = "123456789012345678901234567891/7*x";
        let parsed = ring.parse_polynomial(text).expect("the text parses");
        assert_eq!(parsed.to_string(), text);
    }

    #[test]
    fn a_zero_denominator_in_text_stops_the_parser() {
        let ring = PolynomialRing::rationals(["x"]).expect("the name is a variable name");
        assert_eq!(
            ring.parse_polynomial("1/0*x"),
            Err(RingError::Parse(ParseError::ZeroDenominator {
                position: 2
            }))
        );
        assert_eq!(
            ring.parse_polynomial("1/*x"),
            Err(RingError::Parse(ParseError::UnexpectedCharacter {
                character: '*',
                position: 2
            }))
        );
    }

    #[test]
    fn the_fraction_bar_takes_no_whitespace() {
        let ring = PolynomialRing::rationals(["x"]).expect("the name is a variable name");
        assert_eq!(
            ring.parse_polynomial("1 / 2"),
            Err(RingError::Parse(ParseError::UnexpectedCharacter {
                character: '/',
                position: 2
            }))
        );
    }

    #[test]
    fn a_denominator_the_field_cannot_invert_stops_a_parse() {
        let ring = PolynomialRing::prime_field(3, ["x"]).expect("3 is prime");
        assert_eq!(
            ring.parse_polynomial("1/3*x"),
            Err(RingError::CoefficientNotInvertible {
                denominator: BigInt::from(3)
            })
        );
        // The fraction is reduced first, so 3/3 is the value 1.
        assert_eq!(
            ring.parse_polynomial("3/3*x").expect("3/3 is the value 1"),
            ring.parse_polynomial("x").expect("the text parses")
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
