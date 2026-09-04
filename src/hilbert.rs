//! The Hilbert series of the quotient by a leading monomial ideal.
//!
//! The series reads the leading monomials of a Gröbner basis and nothing
//! else, so one function serves both domains and no coefficient arithmetic
//! runs here. [`GroebnerBasis::hilbert_series`] is where a caller reaches
//! it.
//!
//! [`GroebnerBasis::hilbert_series`]: crate::GroebnerBasis::hilbert_series

use std::fmt;
use std::mem::size_of;

use num_bigint::BigInt;
use num_traits::{One, Signed, Zero};
use rustc_hash::FxHashMap;

use crate::compute::{ComputeError, ComputeLimits, RunError};
use crate::poly::{Monomial, Polynomial, heap_exps_bytes};
use crate::ring::Domain;
use crate::ring::rational::limb_bytes;

/// The Hilbert series of `k[x] / L` as a rational function.
///
/// `L` is the ideal the leading monomials of a Gröbner basis `G` generate.
/// The value is `N(t) / (1 - t)^n`, with `N` a polynomial over `Z` and `n`
/// the number of variables. The coefficient of `t^d` in the power series
/// is the dimension over `k` of the degree `d` piece of `k[x] / L`, which
/// [`HilbertSeries::coefficient`] returns.
///
/// Three statements relate the value to the ideal `I` the basis
/// generates.
///
/// - If `I` is homogeneous, which [`GroebnerBasis::is_homogeneous`]
///   decides, the series of `k[x] / L` is the Hilbert series of
///   `k[x] / I`.
/// - For any `I`, grevlex is a graded order, so the coefficient of `t^d`
///   is the first difference `H(d) - H(d - 1)` of the affine Hilbert
///   function `H(d) = dim_k R_{<=d} / (I ∩ R_{<=d})`, with `H(-1) = 0`.
///   The affine Hilbert function itself is the running sum of the
///   coefficients. The coefficients are not that function: for the zero
///   ideal in one variable they are 1 in every degree, while
///   `H(d) = d + 1`.
/// - The Krull dimension of `k[x] / I` is the dimension
///   [`HilbertSeries::dimension`] reads off the series. Passing to an
///   initial ideal keeps the dimension, so the value is the dimension of
///   the ideal the caller asked about.
///
/// `Display` writes the numerator in parentheses over the denominator,
/// for example `(1 - 2*t^2 + t^3)/(1 - t)^2`.
///
/// [`GroebnerBasis::is_homogeneous`]: crate::GroebnerBasis::is_homogeneous
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct HilbertSeries {
    numerator: Vec<BigInt>,
    denominator_power: usize,
}

impl HilbertSeries {
    /// The numerator, dense and lowest degree first.
    ///
    /// The last coefficient is nonzero, and the zero polynomial is the
    /// empty slice. That case is the unit ideal, whose quotient ring is
    /// zero.
    pub fn numerator(&self) -> &[BigInt] {
        &self.numerator
    }

    /// The power of `1 - t` in the denominator, which is the number of
    /// variables.
    ///
    /// It is the power before any cancellation. The numerator can hold
    /// factors of `1 - t`, and [`HilbertSeries::dimension`] is what counts
    /// them.
    pub fn denominator_power(&self) -> usize {
        self.denominator_power
    }

    /// The Krull dimension of the quotient ring, or `None` for the zero
    /// ring.
    ///
    /// The value is `n - k`, with `k` the number of times `1 - t` divides
    /// the numerator exactly. A zero numerator is the unit ideal: the
    /// quotient ring is zero and its dimension is -1 by the usual
    /// convention, which `None` reports because the return type counts up
    /// from 0.
    pub fn dimension(&self) -> Option<usize> {
        self.cancelled().map(|(_, cancelled)| {
            debug_assert!(
                cancelled <= self.denominator_power,
                "the numerator of a monomial ideal holds at most n factors of 1 - t"
            );
            self.denominator_power.saturating_sub(cancelled)
        })
    }

    /// The multiplicity of the quotient ring, or `None` for the zero ring.
    ///
    /// It is the value at `t = 1` of the numerator after every factor of
    /// `1 - t` is cancelled. For a quotient of dimension 0 it is the
    /// dimension of the whole quotient over `k`, which is the number of
    /// solutions with multiplicity.
    pub fn multiplicity(&self) -> Option<BigInt> {
        self.cancelled()
            .map(|(numerator, _)| numerator.into_iter().sum())
    }

    /// The dimension over `k` of the piece of degree `degree`.
    ///
    /// It is the coefficient of `t^degree` in the power series, which is
    /// the number of monomials of that degree no leading monomial
    /// divides.
    pub fn coefficient(&self, degree: u32) -> BigInt {
        let degree = degree as usize;
        if self.denominator_power == 0 {
            return self
                .numerator
                .get(degree)
                .cloned()
                .unwrap_or_else(BigInt::zero);
        }
        let rank = self.denominator_power - 1;
        self.numerator
            .iter()
            .take(degree + 1)
            .enumerate()
            .map(|(index, coeff)| coeff * binomial(degree - index + rank, rank))
            .sum()
    }

    /// The numerator with every factor of `1 - t` cancelled, next to the
    /// number cancelled, or `None` for the zero numerator.
    ///
    /// The zero numerator is tested first. Cancelling `1 - t` out of it
    /// would not stop.
    fn cancelled(&self) -> Option<(Vec<BigInt>, usize)> {
        if self.numerator.is_empty() {
            return None;
        }
        let mut numerator = self.numerator.clone();
        let mut cancelled = 0;
        while let Some(quotient) = divide_by_one_minus_t(&numerator) {
            numerator = quotient;
            cancelled += 1;
        }
        Some((numerator, cancelled))
    }
}

impl fmt::Display for HilbertSeries {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("(")?;
        if self.numerator.is_empty() {
            f.write_str("0")?;
        } else {
            write_numerator(f, &self.numerator)?;
        }
        write!(f, ")/(1 - t)^{}", self.denominator_power)
    }
}

fn write_numerator(f: &mut fmt::Formatter<'_>, numerator: &[BigInt]) -> fmt::Result {
    let mut first = true;
    for (degree, coefficient) in numerator.iter().enumerate() {
        if coefficient.is_zero() {
            continue;
        }
        write_sign(f, coefficient.is_negative(), first)?;
        write_magnitude(f, coefficient.abs(), degree)?;
        first = false;
    }
    Ok(())
}

fn write_sign(f: &mut fmt::Formatter<'_>, negative: bool, first: bool) -> fmt::Result {
    match (first, negative) {
        (true, true) => f.write_str("-"),
        (true, false) => Ok(()),
        (false, true) => f.write_str(" - "),
        (false, false) => f.write_str(" + "),
    }
}

fn write_magnitude(f: &mut fmt::Formatter<'_>, magnitude: BigInt, degree: usize) -> fmt::Result {
    let coefficient_written = !magnitude.is_one() || degree == 0;
    if coefficient_written {
        write!(f, "{magnitude}")?;
    }
    if degree == 0 {
        return Ok(());
    }
    if coefficient_written {
        f.write_str("*")?;
    }
    if degree == 1 {
        f.write_str("t")
    } else {
        write!(f, "t^{degree}")
    }
}

/// Why a Hilbert series is not computed.
///
/// Both variants report an exhausted budget. The recursion of the
/// numerator branches, and this release offers no bound on its work, so
/// the partiality is typed like every other in the crate. Neither variant
/// says anything about the ideal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HilbertError {
    /// The deadline passed.
    Timeout,
    /// The live data of the recursion passed the memory limit.
    MemoryLimitExceeded,
}

impl fmt::Display for HilbertError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            HilbertError::Timeout => f.write_str("the Hilbert series passed its deadline"),
            HilbertError::MemoryLimitExceeded => {
                f.write_str("the Hilbert series passed its memory limit")
            }
        }
    }
}

impl std::error::Error for HilbertError {}

impl From<RunError> for HilbertError {
    fn from(stop: RunError) -> Self {
        match stop.reported() {
            ComputeError::MemoryLimitExceeded => HilbertError::MemoryLimitExceeded,
            _ => HilbertError::Timeout,
        }
    }
}

/// The Hilbert series of the quotient by the leading monomial ideal of
/// `basis`.
pub(crate) fn hilbert_series<D: Domain>(
    basis: &[Polynomial<D>],
    nvars: usize,
    limits: &ComputeLimits,
) -> Result<HilbertSeries, HilbertError> {
    let leading = basis.iter().filter_map(|f| f.lm()).cloned().collect();
    let generators = minimal(leading);
    let mut recursion = Recursion {
        limits,
        nvars,
        memo: FxHashMap::default(),
        memo_bytes: 0,
        live_bytes: 0,
        frames: 0,
    };
    let numerator = recursion.numerator(&generators)?;
    Ok(HilbertSeries {
        numerator,
        denominator_power: nvars,
    })
}

/// Report whether every polynomial holds terms of one degree.
///
/// The zero polynomial is homogeneous, so an empty list is too.
pub(crate) fn all_homogeneous<D: Domain>(polynomials: &[Polynomial<D>]) -> bool {
    polynomials.iter().all(|f| {
        let mut degrees = f.terms.iter().map(|term| term.mono.deg);
        let first = degrees.next();
        degrees.all(|degree| Some(degree) == first)
    })
}

/// The recursion of section 6.2 of `docs/rational-design.md`, with its memo
/// table and its meters.
struct Recursion<'a> {
    limits: &'a ComputeLimits,
    nvars: usize,
    memo: FxHashMap<Vec<Monomial>, Vec<BigInt>>,
    memo_bytes: usize,
    live_bytes: usize,
    /// The frames the worklist holds.
    frames: usize,
}

/// One node of the worklist.
struct Frame {
    generators: Vec<Monomial>,
    state: State,
    /// The bytes the two child sets hold, released when the node ends.
    held: usize,
    /// The bytes the numerator of the left child holds.
    left_bytes: usize,
}

/// How far a node has got.
enum State {
    /// Nothing is done yet.
    Fresh,
    /// The left child runs. `colon` is the right child, still to run.
    Left { colon: Vec<Monomial> },
    /// The right child runs. `left` is what the left child returned.
    Right { left: Vec<BigInt> },
}

impl Frame {
    fn fresh(generators: Vec<Monomial>) -> Self {
        Frame {
            generators,
            state: State::Fresh,
            held: 0,
            left_bytes: 0,
        }
    }
}

/// What one step of the worklist does next.
enum Step {
    /// Start a child node on `generators`.
    Down(Vec<Monomial>),
    /// End this node with `numerator`.
    Up(Vec<BigInt>),
}

impl Recursion<'_> {
    /// The numerator of the monomial ideal `generators` generates.
    ///
    /// `generators` must be minimal and sorted ascending, which is what
    /// makes it a memo key: one ideal has one key.
    ///
    /// The worklist is a heap allocated stack, not the calling stack. The
    /// colon branch reaches one node per unit of the degree sum, so the
    /// basis `[x^65535*y, z]` reaches 65,535 of them, and the calling
    /// stack does not hold that many frames.
    fn numerator(&mut self, generators: &[Monomial]) -> Result<Vec<BigInt>, HilbertError> {
        let mut stack = vec![Frame::fresh(generators.to_vec())];
        // What the node that ended last returned.
        let mut value: Option<Vec<BigInt>> = None;
        while !stack.is_empty() {
            self.frames = stack.len();
            self.charge(0)?;
            let frame = stack.last_mut().expect("the loop tested the stack");
            let step = self.advance(frame, &mut value)?;
            match step {
                Step::Down(generators) => stack.push(Frame::fresh(generators)),
                Step::Up(numerator) => {
                    stack.pop();
                    value = Some(numerator);
                }
            }
        }
        Ok(value.expect("the root node ended"))
    }

    fn advance(
        &mut self,
        frame: &mut Frame,
        value: &mut Option<Vec<BigInt>>,
    ) -> Result<Step, HilbertError> {
        match std::mem::replace(&mut frame.state, State::Fresh) {
            State::Fresh => self.expand(frame),
            State::Left { colon } => {
                let left = value.take().expect("the left child ended");
                frame.left_bytes = series_bytes(&left);
                self.live_bytes = self.live_bytes.saturating_add(frame.left_bytes);
                frame.state = State::Right { left };
                Ok(Step::Down(colon))
            }
            State::Right { left } => {
                let right = value.take().expect("the right child ended");
                self.charge(dense_bytes(left.len().max(right.len() + 1)))?;
                self.live_bytes = self
                    .live_bytes
                    .saturating_sub(frame.held.saturating_add(frame.left_bytes));
                let sum = add_shifted(left, &right);
                self.remember(&frame.generators, &sum);
                Ok(Step::Up(sum))
            }
        }
    }

    fn expand(&mut self, frame: &mut Frame) -> Result<Step, HilbertError> {
        self.charge(set_bytes(&frame.generators, self.nvars))?;
        if let Some(known) = self.memo.get(&frame.generators) {
            self.charge(series_bytes(known))?;
            return Ok(Step::Up(known.clone()));
        }
        if let Some(base) = self.base(&frame.generators)? {
            self.remember(&frame.generators, &base);
            return Ok(Step::Up(base));
        }
        let pivot = pivot(&frame.generators, self.nvars);
        let plus = minimal(with_variable(&frame.generators, pivot, self.nvars));
        let colon = minimal(colon_by_variable(&frame.generators, pivot));
        debug_assert!(degree_sum(&plus) < degree_sum(&frame.generators));
        debug_assert!(degree_sum(&colon) < degree_sum(&frame.generators));
        frame.held = set_bytes(&plus, self.nvars).saturating_add(set_bytes(&colon, self.nvars));
        self.live_bytes = self.live_bytes.saturating_add(frame.held);
        frame.state = State::Left { colon };
        Ok(Step::Down(plus))
    }

    /// Hold the numerator of `generators` for a node that reaches the same
    /// set again.
    fn remember(&mut self, generators: &[Monomial], value: &[BigInt]) {
        self.memo_bytes = self
            .memo_bytes
            .saturating_add(set_bytes(generators, self.nvars))
            .saturating_add(series_bytes(value));
        self.memo.insert(generators.to_vec(), value.to_vec());
    }

    /// Report why the worklist stops at this step, or `None` to carry on.
    ///
    /// The count covers the sets and the numerators of the live nodes, the
    /// memo table with its keys, the frames of the worklist, and `extra`,
    /// which is what the caller is about to allocate. The numerators are
    /// the large values, so a count of the sets alone would not bound what
    /// the worklist holds.
    fn charge(&self, extra: usize) -> Result<(), HilbertError> {
        if let Some(stop) = self.limits.stop() {
            return Err(HilbertError::from(stop));
        }
        if let Some(limit) = self.limits.memory {
            let bytes = self
                .live_bytes
                .saturating_add(self.memo_bytes)
                .saturating_add(self.frames.saturating_mul(size_of::<Frame>()))
                .saturating_add(extra);
            if bytes > limit {
                return Err(HilbertError::MemoryLimitExceeded);
            }
        }
        Ok(())
    }

    /// The numerator from a base case, or `None` when the node splits.
    ///
    /// Every vector here is charged before it is allocated. The degree of
    /// one generator sets the length, and it reaches 65,535.
    fn base(&self, generators: &[Monomial]) -> Result<Option<Vec<BigInt>>, HilbertError> {
        if generators.iter().any(|m| m.deg == 0) {
            return Ok(Some(Vec::new()));
        }
        match generators {
            [] => return Ok(Some(vec![BigInt::one()])),
            [only] => return self.one_minus_power(only.deg, 0).map(Some),
            _ => {}
        }
        let Some(degrees) = pure_power_degrees(generators) else {
            return Ok(None);
        };
        let mut product = vec![BigInt::one()];
        for degree in degrees {
            let held = series_bytes(&product);
            let factor = self.one_minus_power(degree, held)?;
            product = self.multiply(
                &product,
                &factor,
                held.saturating_add(series_bytes(&factor)),
            )?;
        }
        Ok(Some(product))
    }

    /// `1 - t^degree` as a dense coefficient vector.
    ///
    /// `held` is the bytes the caller holds that the meters do not count
    /// yet.
    fn one_minus_power(&self, degree: u32, held: usize) -> Result<Vec<BigInt>, HilbertError> {
        let length = degree as usize + 1;
        self.charge(held.saturating_add(dense_bytes(length)))?;
        let mut coefficients = vec![BigInt::zero(); length];
        coefficients[0] = BigInt::one();
        coefficients[degree as usize] -= 1;
        Ok(trim(coefficients))
    }

    /// The product of two dense coefficient vectors.
    ///
    /// `held` is the bytes the caller holds that the meters do not count
    /// yet. The loop reads the deadline once every `TICK` products.
    fn multiply(
        &self,
        left: &[BigInt],
        right: &[BigInt],
        held: usize,
    ) -> Result<Vec<BigInt>, HilbertError> {
        if left.is_empty() || right.is_empty() {
            return Ok(Vec::new());
        }
        let length = left.len() + right.len() - 1;
        self.charge(held.saturating_add(dense_bytes(length)))?;
        let mut product = vec![BigInt::zero(); length];
        let mut index = 0usize;
        for (offset, a) in left.iter().enumerate() {
            for (shift, b) in right.iter().enumerate() {
                if let Some(stop) = self.limits.stop_every(index) {
                    return Err(HilbertError::from(stop));
                }
                index += 1;
                product[offset + shift] += a * b;
            }
        }
        Ok(trim(product))
    }
}

/// The minimal generators of the ideal `generators` generates, sorted
/// ascending under grevlex.
fn minimal(mut generators: Vec<Monomial>) -> Vec<Monomial> {
    generators.sort();
    generators.dedup();
    let mut minimal = Vec::with_capacity(generators.len());
    for (index, candidate) in generators.iter().enumerate() {
        let divided = generators
            .iter()
            .enumerate()
            .any(|(other, m)| other != index && m.divides(candidate));
        if !divided {
            minimal.push(candidate.clone());
        }
    }
    minimal
}

/// The variable to split on.
///
/// It comes from a generator that is not a pure power, and among those
/// variables it is one that appears in the most generators, the smallest
/// index winning a tie. The restriction to a generator that is not a pure
/// power is what makes the recursion terminate: a pivot without it can
/// return the parent problem, as `x` does for `(x, y*z)`. The tie break
/// changes neither the value nor the termination, only which subproblems
/// the recursion walks.
fn pivot(generators: &[Monomial], nvars: usize) -> usize {
    let mut counts = vec![0usize; nvars];
    let mut candidate = vec![false; nvars];
    for m in generators {
        let mixed = pure_power(m).is_none();
        for (index, &exp) in m.exps.iter().enumerate() {
            if exp == 0 {
                continue;
            }
            counts[index] += 1;
            candidate[index] |= mixed;
        }
    }
    let mut best: Option<usize> = None;
    for index in 0..nvars {
        if !candidate[index] {
            continue;
        }
        match best {
            Some(current) if counts[current] >= counts[index] => {}
            _ => best = Some(index),
        }
    }
    best.expect("a generator that is not a pure power holds a variable")
}

/// The variable of a pure power, or `None` for every other monomial.
fn pure_power(m: &Monomial) -> Option<usize> {
    let mut variable = None;
    for (index, &exp) in m.exps.iter().enumerate() {
        if exp == 0 {
            continue;
        }
        if variable.is_some() {
            return None;
        }
        variable = Some(index);
    }
    variable
}

/// The degrees of the generators when every one is a pure power of its own
/// variable, and `None` otherwise.
fn pure_power_degrees(generators: &[Monomial]) -> Option<Vec<u32>> {
    let mut variables = Vec::with_capacity(generators.len());
    let mut degrees = Vec::with_capacity(generators.len());
    for m in generators {
        let variable = pure_power(m)?;
        if variables.contains(&variable) {
            return None;
        }
        variables.push(variable);
        degrees.push(m.deg);
    }
    Some(degrees)
}

/// The generators of `L + (x_j)`.
fn with_variable(generators: &[Monomial], pivot: usize, nvars: usize) -> Vec<Monomial> {
    let mut result: Vec<Monomial> = generators
        .iter()
        .filter(|m| m.exps[pivot] == 0)
        .cloned()
        .collect();
    let mut variable = Monomial::one(nvars);
    variable.exps[pivot] = 1;
    variable.deg = 1;
    result.push(variable);
    result
}

/// The generators of `L : x_j`.
fn colon_by_variable(generators: &[Monomial], pivot: usize) -> Vec<Monomial> {
    generators
        .iter()
        .map(|m| {
            let mut quotient = m.clone();
            if quotient.exps[pivot] > 0 {
                quotient.exps[pivot] -= 1;
                quotient.deg -= 1;
            }
            quotient
        })
        .collect()
}

/// The termination measure: the sum of the degrees of the minimal
/// generators.
fn degree_sum(generators: &[Monomial]) -> u64 {
    generators.iter().map(|m| u64::from(m.deg)).sum()
}

/// `left + t * right`.
fn add_shifted(mut left: Vec<BigInt>, right: &[BigInt]) -> Vec<BigInt> {
    if left.len() < right.len() + 1 {
        left.resize(right.len() + 1, BigInt::zero());
    }
    for (index, coeff) in right.iter().enumerate() {
        left[index + 1] += coeff;
    }
    trim(left)
}

/// The quotient of `numerator` by `1 - t`, or `None` when the division
/// leaves a remainder.
///
/// The quotient is the running sum of the coefficients, and the division
/// is exact exactly when the last running sum, which is the value at
/// `t = 1`, is zero.
fn divide_by_one_minus_t(numerator: &[BigInt]) -> Option<Vec<BigInt>> {
    let mut quotient = Vec::with_capacity(numerator.len());
    let mut running = BigInt::zero();
    for coeff in numerator {
        running += coeff;
        quotient.push(running.clone());
    }
    match quotient.pop() {
        Some(last) if last.is_zero() => Some(trim(quotient)),
        _ => None,
    }
}

/// Drop the trailing zero coefficients, so one polynomial has one vector.
fn trim(mut coefficients: Vec<BigInt>) -> Vec<BigInt> {
    while coefficients.last().is_some_and(BigInt::is_zero) {
        coefficients.pop();
    }
    coefficients
}

/// The binomial coefficient `a` over `b`.
fn binomial(a: usize, b: usize) -> BigInt {
    if b > a {
        return BigInt::zero();
    }
    let b = b.min(a - b);
    let mut value = BigInt::one();
    for step in 1..=b {
        // C(a - b + step, step) is an integer at every step, so the
        // division leaves no remainder.
        value = value * BigInt::from(a - b + step) / BigInt::from(step);
    }
    value
}

/// The bytes a set of minimal generators holds.
fn set_bytes(generators: &[Monomial], nvars: usize) -> usize {
    generators
        .len()
        .saturating_mul(size_of::<Monomial>().saturating_add(heap_exps_bytes(nvars)))
}

/// The bytes a dense coefficient vector of `length` coefficients holds.
///
/// It is the projection [`Recursion::charge`] takes before the vector
/// exists. [`series_bytes`] reads the built vector and adds the limbs of
/// every coefficient that spilled to the heap.
fn dense_bytes(length: usize) -> usize {
    length.saturating_mul(size_of::<BigInt>())
}

/// The bytes one numerator holds.
///
/// Each coefficient counts by the estimate of `DomainOps::heap_bytes`:
/// the used bits rounded up to whole limbs.
fn series_bytes(coefficients: &[BigInt]) -> usize {
    coefficients.iter().fold(
        coefficients.len().saturating_mul(size_of::<BigInt>()),
        |bytes, coeff| bytes.saturating_add(limb_bytes(coeff)),
    )
}
