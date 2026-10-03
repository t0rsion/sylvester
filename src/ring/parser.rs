//! Parenthesized polynomial expressions with one call-wide budget.

use std::any::TypeId;
use std::fmt;
use std::mem::size_of;

use num_bigint::BigInt;
use num_traits::Zero;

use crate::arithmetic::ArithmeticError;
use crate::compute::{ComputeError, ComputeLimits, RunError};
use crate::poly::{Monomial, Polynomial, Term, heap_exps_bytes};

use super::{
    Coefficient, Domain, DomainOps, ParseError, PolynomialRing, Rationals, RingError,
    integer_coefficient,
};

/// The largest parenthesis depth an expression parser enters.
const MAX_NESTING: usize = 256;

/// Why [`PolynomialRing::parse_polynomial_with_budget`] stops.
///
/// Syntax errors retain the positions used by [`ParseError`]. Arithmetic
/// errors report a coefficient conversion failure or an exponent limit.
/// Timeout and memory errors apply to the whole expression.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ExpressionError {
    /// The expression has invalid syntax.
    Parse(ParseError),
    /// A coefficient could not be read by the ring's domain.
    Ring(RingError),
    /// The expression passed its deadline.
    Timeout,
    /// The expression passed its memory limit.
    MemoryLimitExceeded,
    /// A monomial exponent passed the `u16` width.
    ExponentLimit {
        /// The largest exponent one monomial stores.
        limit: u32,
    },
    /// The expression entered more parenthesis levels than the parser holds.
    NestingLimit {
        /// Where the rejected opening parenthesis starts.
        position: usize,
        /// The largest parenthesis depth the parser accepts.
        limit: usize,
    },
    /// A compute error outside the expression limits reached the parser.
    Compute(ComputeError),
    /// An arithmetic error outside the expression limits reached the parser.
    Arithmetic(ArithmeticError),
}

impl fmt::Display for ExpressionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ExpressionError::Parse(error) => error.fmt(f),
            ExpressionError::Ring(error) => error.fmt(f),
            ExpressionError::Timeout => f.write_str("the expression passed its deadline"),
            ExpressionError::MemoryLimitExceeded => {
                f.write_str("the expression passed its memory limit")
            }
            ExpressionError::ExponentLimit { limit } => {
                write!(f, "an expression exponent is above {limit}")
            }
            ExpressionError::NestingLimit { position, limit } => write!(
                f,
                "expression nesting at byte {position} is above the limit {limit}"
            ),
            ExpressionError::Compute(error) => error.fmt(f),
            ExpressionError::Arithmetic(error) => error.fmt(f),
        }
    }
}

impl std::error::Error for ExpressionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ExpressionError::Parse(error) => Some(error),
            ExpressionError::Ring(error) => Some(error),
            ExpressionError::Timeout
            | ExpressionError::MemoryLimitExceeded
            | ExpressionError::ExponentLimit { .. }
            | ExpressionError::NestingLimit { .. }
            | ExpressionError::Compute(_)
            | ExpressionError::Arithmetic(_) => None,
        }
    }
}

impl From<ParseError> for ExpressionError {
    fn from(error: ParseError) -> Self {
        ExpressionError::Parse(error)
    }
}

impl From<RingError> for ExpressionError {
    fn from(error: RingError) -> Self {
        match error {
            RingError::Parse(error) => ExpressionError::Parse(error),
            other => ExpressionError::Ring(other),
        }
    }
}

impl From<RunError> for ExpressionError {
    fn from(error: RunError) -> Self {
        match error.reported() {
            ComputeError::Timeout => ExpressionError::Timeout,
            ComputeError::MemoryLimitExceeded => ExpressionError::MemoryLimitExceeded,
            ComputeError::DegreeLimit { limit } | ComputeError::ExponentLimit { limit } => {
                ExpressionError::ExponentLimit { limit }
            }
            other => ExpressionError::Compute(other),
        }
    }
}

impl From<ArithmeticError> for ExpressionError {
    fn from(error: ArithmeticError) -> Self {
        match error {
            ArithmeticError::CoefficientConversion(error) => error.into(),
            ArithmeticError::ExponentLimit { limit } => ExpressionError::ExponentLimit { limit },
            ArithmeticError::Timeout => ExpressionError::Timeout,
            ArithmeticError::MemoryLimitExceeded => ExpressionError::MemoryLimitExceeded,
            other => ExpressionError::Arithmetic(other),
        }
    }
}

/// Parse one expression under caller-owned absolute limits.
pub(super) fn parse_with_limits<D: Domain>(
    ring: &PolynomialRing<D>,
    text: &str,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, ExpressionError> {
    let mut parser = ExpressionParser::new(ring, text, limits)?;
    let result = parser.sum()?;
    parser.skip_spaces()?;
    if let Some(byte) = parser.peek() {
        return Err(parser.unexpected(byte).into());
    }
    parser.budget.check()?;
    Ok(result)
}

/// A recursive-descent parser. Parenthesis depth is bounded before descent.
struct ExpressionParser<'a, D: Domain> {
    ring: &'a PolynomialRing<D>,
    text: &'a str,
    at: usize,
    depth: usize,
    budget: ExpressionBudget,
}

impl<'a, D: Domain> ExpressionParser<'a, D> {
    fn new(
        ring: &'a PolynomialRing<D>,
        text: &'a str,
        limits: &ComputeLimits,
    ) -> Result<Self, ExpressionError> {
        let state = ExpressionParser {
            ring,
            text,
            at: 0,
            depth: 0,
            budget: ExpressionBudget::new(limits.clone()),
        };
        state.budget.check()?;
        Ok(state)
    }

    fn sum(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        let mut result = self.product()?;
        loop {
            self.skip_spaces()?;
            let subtract = match self.peek() {
                Some(b'+') => false,
                Some(b'-') => true,
                _ => break,
            };
            self.at += 1;
            let right = self.product()?;
            result = self.combine_sum(&result, &right, subtract)?;
        }
        Ok(result)
    }

    fn product(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        let mut result = self.unary()?;
        loop {
            self.skip_spaces()?;
            if self.peek() != Some(b'*') {
                break;
            }
            self.at += 1;
            let right = self.unary()?;
            result = self.multiply(&result, &right)?;
        }
        Ok(result)
    }

    fn unary(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        self.skip_spaces()?;
        let mut negative = false;
        loop {
            match self.peek() {
                Some(b'+') => {
                    self.at += 1;
                    self.budget.tick()?;
                }
                Some(b'-') => {
                    self.at += 1;
                    negative = !negative;
                    self.budget.tick()?;
                }
                _ => break,
            }
            self.skip_spaces()?;
        }
        let value = self.power()?;
        if negative {
            self.negate(&value)
        } else {
            Ok(value)
        }
    }

    fn power(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        let value = self.atom()?;
        self.skip_spaces()?;
        let operator = match self.peek() {
            Some(b'^') => {
                self.at += 1;
                true
            }
            Some(b'*') if self.text.as_bytes().get(self.at.saturating_add(1)) == Some(&b'*') => {
                self.at += 2;
                true
            }
            _ => false,
        };
        if !operator {
            return Ok(value);
        }
        self.skip_spaces()?;
        let exponent = self.exponent()?;
        if exponent > u64::from(u16::MAX) && value.degree().is_some_and(|degree| degree > 0) {
            return Err(ExpressionError::ExponentLimit {
                limit: u32::from(u16::MAX),
            });
        }
        self.power_value(value, exponent)
    }

    fn atom(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        self.skip_spaces()?;
        let Some(byte) = self.peek() else {
            return Err(ParseError::UnexpectedEnd { position: self.at }.into());
        };
        self.atom_byte(byte)
    }

    fn atom_byte(&mut self, byte: u8) -> Result<Polynomial<D>, ExpressionError> {
        if byte == b'(' {
            return self.parenthesized();
        }
        if byte.is_ascii_digit() {
            let value = self.number()?;
            return self.constant(value);
        }
        if is_identifier_start(byte) {
            return self.variable();
        }
        Err(self.unexpected(byte).into())
    }

    fn parenthesized(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        let position = self.at;
        if self.depth >= MAX_NESTING {
            return Err(ExpressionError::NestingLimit {
                position,
                limit: MAX_NESTING,
            });
        }
        self.at += 1;
        self.depth += 1;
        let result = self.sum()?;
        self.skip_spaces()?;
        if self.peek() != Some(b')') {
            return match self.peek() {
                Some(byte) => Err(self.unexpected(byte).into()),
                None => Err(ParseError::UnexpectedEnd { position: self.at }.into()),
            };
        }
        self.at += 1;
        self.depth -= 1;
        Ok(result)
    }

    fn number(&mut self) -> Result<Coefficient, ExpressionError> {
        let numerator = self.integer(0)?;
        self.skip_spaces()?;
        let Some(b'/') = self.peek() else {
            return Ok(integer_coefficient(numerator));
        };
        self.at += 1;
        self.skip_spaces()?;
        let denominator_position = self.at;
        if !matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            return match self.peek() {
                Some(byte) => Err(self.unexpected(byte).into()),
                None => Err(ParseError::UnexpectedEnd { position: self.at }.into()),
            };
        }
        let denominator = self.integer(bigint_bytes(&numerator))?;
        if denominator.is_zero() {
            return Err(ParseError::ZeroDenominator {
                position: denominator_position,
            }
            .into());
        }
        Ok(Coefficient::Fraction {
            numerator,
            denominator,
        })
    }

    fn integer(&mut self, held: usize) -> Result<BigInt, ExpressionError> {
        let start = self.at;
        while matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            self.at += 1;
            self.budget.tick()?;
            if self.budget.limits.memory.is_some() {
                self.budget.integer(self.at.saturating_sub(start), held)?;
            }
        }
        Ok(
            BigInt::parse_bytes(&self.text.as_bytes()[start..self.at], 10)
                .expect("a run of decimal digits is an integer"),
        )
    }

    fn exponent(&mut self) -> Result<u64, ExpressionError> {
        if !matches!(self.peek(), Some(byte) if byte.is_ascii_digit()) {
            return match self.peek() {
                Some(byte) => Err(self.unexpected(byte).into()),
                None => Err(ParseError::UnexpectedEnd { position: self.at }.into()),
            };
        }
        let mut exponent = 0u64;
        while let Some(byte) = self.peek() {
            if !byte.is_ascii_digit() {
                break;
            }
            exponent = match exponent
                .checked_mul(10)
                .and_then(|value| value.checked_add(u64::from(byte - b'0')))
            {
                Some(value) => value,
                None => {
                    return Err(ExpressionError::ExponentLimit {
                        limit: u32::from(u16::MAX),
                    });
                }
            };
            self.at += 1;
            self.budget.tick()?;
        }
        Ok(exponent)
    }

    fn variable(&mut self) -> Result<Polynomial<D>, ExpressionError> {
        let start = self.at;
        while matches!(self.peek(), Some(byte) if byte.is_ascii_alphanumeric() || byte == b'_') {
            self.at += 1;
            self.budget.tick()?;
        }
        let name = &self.text[start..self.at];
        let Some(index) = self.ring.variable_index(name) else {
            self.budget.ensure_bytes(string_bytes(name.len()))?;
            return Err(ParseError::UnknownVariable {
                name: name.to_string(),
                position: start,
            }
            .into());
        };
        self.budget
            .ensure_bytes(single_term_bound::<D>(self.ring.nvars()))?;
        let mut exps = vec![0u16; self.ring.nvars()];
        exps[index] = 1;
        let polynomial = Polynomial::from_sorted_terms(
            self.ring.clone(),
            vec![Term {
                coeff: self.ring.ops().one(),
                mono: Monomial::from_exps(exps.into_iter().collect()),
            }],
        );
        self.budget.result(&polynomial)?;
        Ok(polynomial)
    }

    fn constant(&mut self, value: Coefficient) -> Result<Polynomial<D>, ExpressionError> {
        let term_bound = single_term_bound::<D>(self.ring.nvars());
        self.budget
            .ensure_bytes(coefficient_bound(&value).saturating_add(term_bound))?;
        let Some(coeff) = self.ring.ops().convert(&value)? else {
            return Ok(self.ring.zero());
        };
        let polynomial = Polynomial::from_sorted_terms(
            self.ring.clone(),
            vec![Term {
                coeff,
                mono: Monomial::one(self.ring.nvars()),
            }],
        );
        self.budget.result(&polynomial)?;
        Ok(polynomial)
    }

    fn combine_sum(
        &mut self,
        left: &Polynomial<D>,
        right: &Polynomial<D>,
        subtract: bool,
    ) -> Result<Polynomial<D>, ExpressionError> {
        let result = if subtract {
            left.try_sub_with_limits(right, &self.budget.limits)?
        } else {
            left.try_add_with_limits(right, &self.budget.limits)?
        };
        Ok(result)
    }

    fn negate(&mut self, value: &Polynomial<D>) -> Result<Polynomial<D>, ExpressionError> {
        Ok(value.try_neg_with_limits(&self.budget.limits)?)
    }

    fn multiply(
        &mut self,
        left: &Polynomial<D>,
        right: &Polynomial<D>,
    ) -> Result<Polynomial<D>, ExpressionError> {
        Ok(left.try_mul_with_limits(right, &self.budget.limits)?)
    }

    fn power_value(
        &mut self,
        value: Polynomial<D>,
        mut exponent: u64,
    ) -> Result<Polynomial<D>, ExpressionError> {
        let one_bound = single_term_bound::<D>(self.ring.nvars());
        self.budget
            .ensure_bytes(value.retained_bytes().saturating_add(one_bound))?;
        let mut result = self.ring.one();
        let mut power = value;
        while exponent != 0 {
            self.budget.tick()?;
            if exponent & 1 == 1 {
                result = self.multiply(&result, &power)?;
            }
            exponent >>= 1;
            if exponent != 0 {
                power = self.multiply(&power, &power)?;
            }
        }
        self.budget.result(&result)?;
        Ok(result)
    }

    fn skip_spaces(&mut self) -> Result<(), ExpressionError> {
        while matches!(self.peek(), Some(byte) if byte.is_ascii_whitespace()) {
            self.at += 1;
            self.budget.tick()?;
        }
        Ok(())
    }

    fn peek(&self) -> Option<u8> {
        self.text.as_bytes().get(self.at).copied()
    }

    fn unexpected(&self, byte: u8) -> ParseError {
        let character = self
            .text
            .get(self.at..)
            .and_then(|rest| rest.chars().next())
            .unwrap_or(byte as char);
        ParseError::UnexpectedCharacter {
            character,
            position: self.at,
        }
    }
}

fn is_identifier_start(byte: u8) -> bool {
    byte.is_ascii_alphabetic() || byte == b'_'
}

/// A call-wide estimate of parser allocations and deadline checks.
struct ExpressionBudget {
    limits: ComputeLimits,
    steps: usize,
}

impl ExpressionBudget {
    fn new(limits: ComputeLimits) -> Self {
        ExpressionBudget { limits, steps: 0 }
    }

    fn check(&self) -> Result<(), ExpressionError> {
        self.limits.stop().map_or(Ok(()), |error| Err(error.into()))
    }

    fn tick(&mut self) -> Result<(), ExpressionError> {
        let index = self.steps;
        self.steps = self.steps.saturating_add(1);
        self.limits
            .stop_every(index)
            .map_or(Ok(()), |error| Err(error.into()))
    }

    fn result<D: Domain>(&mut self, value: &Polynomial<D>) -> Result<(), ExpressionError> {
        self.check()?;
        self.ensure_bytes(value.retained_bytes())
    }

    fn integer(&self, digits: usize, held: usize) -> Result<(), ExpressionError> {
        self.ensure_bytes(held.saturating_add(integer_bound(digits)))
    }

    fn ensure_bytes(&self, bytes: usize) -> Result<(), ExpressionError> {
        self.check()?;
        if self.limits.memory.is_some_and(|limit| bytes > limit) {
            return Err(ExpressionError::MemoryLimitExceeded);
        }
        Ok(())
    }
}

fn bigint_bytes(value: &BigInt) -> usize {
    let bits = usize::try_from(value.bits()).unwrap_or(usize::MAX);
    bits.div_ceil(u64::BITS as usize)
        .saturating_mul(size_of::<u64>())
        .saturating_add(size_of::<BigInt>())
}

fn integer_bound(digits: usize) -> usize {
    digits
        .saturating_mul(4)
        .div_ceil(u64::BITS as usize)
        .saturating_mul(size_of::<u64>())
        .saturating_add(size_of::<BigInt>())
}

fn coefficient_bound(value: &Coefficient) -> usize {
    match value {
        Coefficient::Small(_) => size_of::<Coefficient>(),
        Coefficient::Integer(integer) => size_of::<Coefficient>()
            .saturating_add(bigint_bytes(integer))
            .saturating_mul(4),
        Coefficient::Fraction {
            numerator,
            denominator,
        } => size_of::<Coefficient>()
            .saturating_add(bigint_bytes(numerator))
            .saturating_add(bigint_bytes(denominator))
            .saturating_mul(4),
        Coefficient::Rational(value) => size_of::<Coefficient>()
            .saturating_add(bigint_bytes(value.numer()))
            .saturating_add(bigint_bytes(value.denom()))
            .saturating_mul(4),
    }
}

fn string_bytes(length: usize) -> usize {
    size_of::<String>().saturating_add(length)
}

fn single_term_bound<D: Domain>(nvars: usize) -> usize {
    size_of::<Term<D>>()
        .saturating_add(heap_exps_bytes(nvars))
        .saturating_add(nvars.saturating_mul(size_of::<u16>()))
        .saturating_add(size_of::<Vec<u16>>())
        .saturating_add(domain_coefficient_bound::<D>())
}

fn domain_coefficient_bound<D: Domain>() -> usize {
    if TypeId::of::<D>() == TypeId::of::<Rationals>() {
        2 * size_of::<u64>()
    } else {
        0
    }
}
