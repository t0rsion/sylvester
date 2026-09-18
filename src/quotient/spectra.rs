//! Characteristic and minimal polynomials of finite quotient operators.

use std::fmt;
use std::mem::size_of;

use crate::compute::{ComputeError, TICK};
use crate::ring::{Domain, DomainOps, PolynomialRing, PrimeField};

use super::linear::{MultiplicationMatrix, QuotientBudget};

type CoefficientVector<D> = Vec<Option<<D as Domain>::Coeff>>;
type ChargedVector<D> = (CoefficientVector<D>, usize);
type OptionalChargedVector<D> = Option<ChargedVector<D>>;
type ChargedDense<D> = (Vec<CoefficientVector<D>>, usize);

/// A polynomial in one indeterminate over a ring's coefficient field.
///
/// Coefficients are stored from the constant term upward. `None` represents
/// zero. Trailing zero coefficients are removed, and a zero polynomial keeps
/// one `None` entry so its coefficient slice is never empty.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnivariatePolynomial<D: Domain = PrimeField> {
    ring: PolynomialRing<D>,
    coefficients: Vec<Option<D::Coeff>>,
}

impl<D: Domain> UnivariatePolynomial<D> {
    /// Build a polynomial from low-first coefficients.
    pub(crate) fn from_coefficients(
        ring: PolynomialRing<D>,
        mut coefficients: Vec<Option<D::Coeff>>,
    ) -> Self {
        while coefficients.len() > 1 && coefficients.last().is_some_and(Option::is_none) {
            coefficients.pop();
        }
        if coefficients.is_empty() {
            coefficients.push(None);
        }
        UnivariatePolynomial { ring, coefficients }
    }

    /// The coefficient ring used by the polynomial.
    pub fn ring(&self) -> &PolynomialRing<D> {
        &self.ring
    }

    /// The coefficients from the constant term upward.
    pub fn coefficients(&self) -> &[Option<D::Coeff>] {
        &self.coefficients
    }

    /// Read the coefficient of `t^degree`, or `None` for zero.
    pub fn coefficient(&self, degree: usize) -> Option<&D::Coeff> {
        self.coefficients.get(degree).and_then(Option::as_ref)
    }

    /// The largest nonzero degree, or `None` for zero.
    pub fn degree(&self) -> Option<usize> {
        self.coefficients.iter().rposition(Option::is_some)
    }

    /// Report whether every coefficient is zero.
    pub fn is_zero(&self) -> bool {
        self.degree().is_none()
    }

    /// Take the low-first coefficients.
    pub fn into_coefficients(self) -> Vec<Option<D::Coeff>> {
        self.coefficients
    }
}

impl<D: Domain> fmt::Display for UnivariatePolynomial<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let ops = self.ring.ops();
        let mut first = true;
        for (degree, coefficient) in self.coefficients.iter().enumerate().rev() {
            let Some(coefficient) = coefficient else {
                continue;
            };
            write_univariate_term::<D>(f, ops, degree, coefficient, &mut first)?;
        }
        if first {
            f.write_str("0")?;
        }
        Ok(())
    }
}

fn write_univariate_term<D: Domain>(
    f: &mut fmt::Formatter<'_>,
    ops: &D::Ops,
    degree: usize,
    coefficient: &D::Coeff,
    first: &mut bool,
) -> fmt::Result {
    let negative = ops.is_negative(coefficient);
    write_univariate_sign(f, *first, negative)?;
    let magnitude = if negative {
        ops.neg(coefficient)
    } else {
        coefficient.clone()
    };
    let coefficient_written = !ops.is_one(&magnitude) || degree == 0;
    if coefficient_written {
        ops.write(&magnitude, f)?;
    }
    write_univariate_power(f, degree, coefficient_written)?;
    *first = false;
    Ok(())
}

fn write_univariate_sign(f: &mut fmt::Formatter<'_>, first: bool, negative: bool) -> fmt::Result {
    match (first, negative) {
        (true, true) => f.write_str("-"),
        (true, false) => Ok(()),
        (false, true) => f.write_str(" - "),
        (false, false) => f.write_str(" + "),
    }
}

fn write_univariate_power(
    f: &mut fmt::Formatter<'_>,
    degree: usize,
    coefficient_written: bool,
) -> fmt::Result {
    if degree == 0 {
        return Ok(());
    }
    if coefficient_written {
        f.write_str("*")?;
    }
    match degree {
        1 => f.write_str("t"),
        _ => write!(f, "t^{degree}"),
    }
}

/// Compute `det(t I - A)` by the division-free Berkowitz recurrence.
///
/// Berkowitz uses additions and multiplications only. It therefore remains
/// valid when the characteristic divides a Newton identity index.
pub(crate) fn characteristic_polynomial<D: Domain>(
    matrix: &MultiplicationMatrix<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<UnivariatePolynomial<D>, ComputeError> {
    let ring = matrix.ring().clone();
    let (coefficients, _bytes) = berkowitz(matrix, &ring, budget)?;
    Ok(UnivariatePolynomial::from_coefficients(ring, coefficients))
}

/// Compute the minimal polynomial of multiplication by the residue class.
///
/// The first dependence among the vectors `1, a, a^2, ...` is found by rank
/// checks. The leading coefficient of its relation is then normalized to 1.
pub(crate) fn minimal_polynomial<D: Domain>(
    matrix: &MultiplicationMatrix<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<UnivariatePolynomial<D>, ComputeError> {
    let dimension = matrix.dimension();
    let ring = matrix.ring().clone();
    if dimension == 0 {
        let (one, _one_bytes) = charged_one(&ring, budget)?;
        return Ok(UnivariatePolynomial::from_coefficients(
            ring,
            vec![Some(one)],
        ));
    }

    let outer_bytes = (dimension + 1).saturating_mul(size_of::<ChargedVector<D>>());
    budget.reserve(outer_bytes)?;
    let mut columns: Vec<ChargedVector<D>> = Vec::with_capacity(dimension + 1);
    let (mut current, mut current_bytes) = zero_vector(dimension, &ring, budget)?;
    let (one, one_bytes) = charged_one(&ring, budget)?;
    current[0] = Some(one);
    current_bytes = current_bytes.saturating_add(one_bytes);

    if let Some((relation, _relation_bytes)) = first_krylov_relation(
        matrix,
        dimension,
        &ring,
        &mut columns,
        current,
        current_bytes,
        budget,
    )? {
        let result = UnivariatePolynomial::from_coefficients(ring, relation);
        release_columns::<D>(&mut columns, budget);
        drop(columns);
        budget.release(outer_bytes);
        return Ok(result);
    }

    release_columns::<D>(&mut columns, budget);
    drop(columns);
    budget.release(outer_bytes);
    unreachable!("a square matrix has a dependence among d + 1 Krylov vectors")
}

fn first_krylov_relation<D: Domain>(
    matrix: &MultiplicationMatrix<D>,
    dimension: usize,
    ring: &PolynomialRing<D>,
    columns: &mut Vec<ChargedVector<D>>,
    mut current: CoefficientVector<D>,
    mut current_bytes: usize,
    budget: &mut QuotientBudget<'_>,
) -> Result<OptionalChargedVector<D>, ComputeError> {
    for degree in 0..=dimension {
        budget.check()?;
        columns.push((current, current_bytes));
        if let Some(relation) = first_dependence(columns, degree, dimension, ring, budget)? {
            return Ok(Some(relation));
        }
        if degree == dimension {
            return Ok(None);
        }
        let (next, next_bytes) = matrix_vector(matrix, &columns[degree].0, ring, budget)?;
        current = next;
        current_bytes = next_bytes;
    }
    Ok(None)
}

fn first_dependence<D: Domain>(
    columns: &[ChargedVector<D>],
    degree: usize,
    dimension: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<OptionalChargedVector<D>, ComputeError> {
    let rank = rank_columns(columns, dimension, ring, budget)?;
    if rank < degree + 1 {
        debug_assert_eq!(rank, degree);
        return solve_relation(columns, degree, ring, budget).map(Some);
    }
    Ok(None)
}

/// Compute the characteristic polynomial recursively on leading principal
/// blocks using the Berkowitz block recurrence.
fn berkowitz<M, D>(
    matrix: &M,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    if let Some(result) = berkowitz_base(matrix, ring, budget)? {
        return Ok(result);
    }
    let dimension = matrix.dimension();
    let (minor, minor_bytes) = principal_minor(matrix, ring, budget)?;
    let minor_matrix = DenseMatrix {
        dimension: dimension - 1,
        entries: minor,
    };
    let (q, q_bytes) = berkowitz(&minor_matrix, ring, budget)?;
    let (s, s_bytes) = block_power_sums(matrix, &minor_matrix, ring, budget)?;
    let (diagonal, diagonal_bytes) = charged_clone_option(matrix.entry(0, 0), ring, budget)?;
    let result = assemble_characteristic(&q, &s, diagonal.as_ref(), ring, budget)?;
    drop(diagonal);
    budget.release(diagonal_bytes);
    drop(minor_matrix);
    budget.release(minor_bytes);
    drop(s);
    budget.release(s_bytes);
    drop(q);
    budget.release(q_bytes);
    Ok(result)
}

fn berkowitz_base<M, D>(
    matrix: &M,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<OptionalChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let dimension = matrix.dimension();
    if dimension == 0 {
        let (mut result, mut bytes) = coefficient_vector::<D>(1, budget)?;
        let (one, one_bytes) = charged_one(ring, budget)?;
        result.push(Some(one));
        bytes = bytes.saturating_add(one_bytes);
        return Ok(Some((result, bytes)));
    }
    if dimension == 1 {
        let value = matrix.entry(0, 0);
        let (mut result, mut bytes) = coefficient_vector::<D>(2, budget)?;
        let (negative, negative_bytes) = charged_negate(value, ring, budget)?;
        result.push(negative);
        bytes = bytes.saturating_add(negative_bytes);
        let (one, one_bytes) = charged_one(ring, budget)?;
        result.push(Some(one));
        bytes = bytes.saturating_add(one_bytes);
        return Ok(Some((result, bytes)));
    }
    Ok(None)
}

fn block_power_sums<M, D>(
    matrix: &M,
    minor_matrix: &DenseMatrix<D>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let dimension = matrix.dimension();
    let (r, r_bytes) = border_row(matrix, dimension, ring, budget)?;
    let (vector, vector_bytes) = border_column(matrix, dimension, ring, budget)?;
    let result = power_sums(r, vector, vector_bytes, minor_matrix, ring, budget);
    budget.release(r_bytes);
    result
}

fn border_row<M, D>(
    matrix: &M,
    dimension: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let payload = (1..dimension).fold(0usize, |bytes, column| {
        bytes.saturating_add(coefficient_bytes(ring, matrix.entry(0, column)))
    });
    let (mut row, bytes) = coefficient_vector_with_payload::<D>(dimension - 1, payload, budget)?;
    for column in 1..dimension {
        if column.is_multiple_of(TICK) {
            budget.check()?;
        }
        row.push(matrix.entry(0, column).cloned());
    }
    budget.check()?;
    Ok((row, bytes))
}

fn border_column<M, D>(
    matrix: &M,
    dimension: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let payload = (1..dimension).fold(0usize, |bytes, row| {
        bytes.saturating_add(coefficient_bytes(ring, matrix.entry(row, 0)))
    });
    let (mut column, bytes) = coefficient_vector_with_payload::<D>(dimension - 1, payload, budget)?;
    for row in 1..dimension {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        column.push(matrix.entry(row, 0).cloned());
    }
    budget.check()?;
    Ok((column, bytes))
}

fn power_sums<D: Domain>(
    r: CoefficientVector<D>,
    mut vector: CoefficientVector<D>,
    mut vector_bytes: usize,
    minor_matrix: &DenseMatrix<D>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let length = minor_matrix.dimension;
    let (mut sums, sums_structural) = coefficient_vector::<D>(length, budget)?;
    let mut sums_bytes = sums_structural;
    for index in 0..length {
        let sum = dot(&r, &vector, ring, budget)?;
        sums_bytes = sums_bytes.saturating_add(coefficient_bytes(ring, sum.as_ref()));
        sums.push(sum);
        if index + 1 < length {
            let (next, next_bytes) = matrix_vector(minor_matrix, &vector, ring, budget)?;
            drop(vector);
            budget.release(vector_bytes);
            vector = next;
            vector_bytes = next_bytes;
        }
    }
    budget.check()?;
    drop(vector);
    budget.release(vector_bytes);
    Ok((sums, sums_bytes))
}

fn assemble_characteristic<D: Domain>(
    q: &[Option<D::Coeff>],
    sums: &[Option<D::Coeff>],
    diagonal: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let dimension = q.len();
    let (mut result, result_structural) = coefficient_vector::<D>(dimension + 1, budget)?;
    let mut result_payload = 0usize;
    for index in 0..=dimension {
        if index.is_multiple_of(TICK) {
            budget.check()?;
        }
        let (coefficient, actual) =
            assembled_coefficient(index, dimension, q, sums, diagonal, ring, budget)?;
        result.push(coefficient);
        result_payload = result_payload.saturating_add(actual);
    }
    budget.check()?;
    result.reverse();
    let result_bytes = result_structural.saturating_add(result_payload);
    Ok((result, result_bytes))
}

fn assembled_coefficient<D: Domain>(
    index: usize,
    dimension: usize,
    q: &[Option<D::Coeff>],
    sums: &[Option<D::Coeff>],
    diagonal: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(Option<D::Coeff>, usize), ComputeError> {
    let bound = assemble_coefficient_bound(index, dimension, q, sums, diagonal, ring, budget)?;
    budget.reserve(bound)?;
    let value = if index == 0 {
        Ok(Some(ring.ops().one()))
    } else {
        assemble_coefficient(index, dimension, q, sums, diagonal, ring, budget)
    };
    let value = match value {
        Ok(value) => value,
        Err(error) => {
            budget.release(bound);
            return Err(error);
        }
    };
    settle_value(value, bound, ring, budget)
}

fn assemble_coefficient_bound<D: Domain>(
    index: usize,
    dimension: usize,
    q: &[Option<D::Coeff>],
    sums: &[Option<D::Coeff>],
    diagonal: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<usize, ComputeError> {
    if index == 0 {
        return Ok(size_of::<D::Coeff>());
    }
    let minor_degree = dimension - 1;
    let mut input = coefficient_bytes(ring, diagonal).saturating_add(coefficient_bytes(
        ring,
        q.get(minor_degree - (index - 1)).and_then(Option::as_ref),
    ));
    if index <= minor_degree {
        input = input.saturating_add(coefficient_bytes(
            ring,
            q.get(minor_degree - index).and_then(Option::as_ref),
        ));
    }
    if index >= 2 {
        for left in 0..=index - 2 {
            if left.is_multiple_of(TICK) {
                budget.check()?;
            }
            input = input
                .saturating_add(coefficient_bytes(
                    ring,
                    q.get(minor_degree - left).and_then(Option::as_ref),
                ))
                .saturating_add(coefficient_bytes(
                    ring,
                    sums.get(index - 2 - left).and_then(Option::as_ref),
                ));
        }
    }
    Ok(coefficient_growth_bound(
        input.saturating_add(size_of::<D::Coeff>()),
    ))
}

fn assemble_coefficient<D: Domain>(
    index: usize,
    dimension: usize,
    q: &[Option<D::Coeff>],
    sums: &[Option<D::Coeff>],
    diagonal: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<Option<D::Coeff>, ComputeError> {
    let minor_degree = dimension - 1;
    let mut coefficient = (index <= minor_degree)
        .then(|| q[minor_degree - index].clone())
        .flatten();
    let previous = q
        .get(minor_degree - (index - 1))
        .and_then(Option::as_ref)
        .and_then(|value| diagonal.map(|a| ring.ops().mul(a, value)));
    coefficient = subtract(coefficient.as_ref(), previous.as_ref(), ring);
    if index >= 2 {
        for left in 0..=index - 2 {
            if left.is_multiple_of(TICK) {
                budget.check()?;
            }
            let product = multiply(
                q.get(minor_degree - left).and_then(Option::as_ref),
                sums.get(index - 2 - left).and_then(Option::as_ref),
                ring,
            );
            coefficient = subtract(coefficient.as_ref(), product.as_ref(), ring);
        }
    }
    Ok(coefficient)
}

/// A dense minor used by Berkowitz.
struct DenseMatrix<D: Domain> {
    dimension: usize,
    entries: Vec<CoefficientVector<D>>,
}

trait MatrixView<D: Domain> {
    fn dimension(&self) -> usize;

    fn entry(&self, row: usize, column: usize) -> Option<&D::Coeff>;
}

impl<D: Domain> MatrixView<D> for MultiplicationMatrix<D> {
    fn dimension(&self) -> usize {
        self.dimension()
    }

    fn entry(&self, row: usize, column: usize) -> Option<&D::Coeff> {
        self.entry_value(row, column)
    }
}

impl<D: Domain> MatrixView<D> for DenseMatrix<D> {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn entry(&self, row: usize, column: usize) -> Option<&D::Coeff> {
        self.entries[row][column].as_ref()
    }
}

fn principal_minor<M, D>(
    matrix: &M,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedDense<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let dimension = matrix.dimension() - 1;
    let payload = (1..matrix.dimension()).fold(0usize, |bytes, row| {
        bytes.saturating_add((1..matrix.dimension()).fold(0usize, |row_bytes, column| {
            row_bytes.saturating_add(coefficient_bytes(ring, matrix.entry(row, column)))
        }))
    });
    let structural = dense_structural_bytes::<D>(dimension, dimension);
    let charged = structural.saturating_add(payload);
    budget.reserve(charged)?;
    let mut entries = Vec::with_capacity(dimension);
    for row in 1..matrix.dimension() {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        let mut values = Vec::with_capacity(dimension);
        for column in 1..matrix.dimension() {
            values.push(matrix.entry(row, column).cloned());
        }
        entries.push(values);
    }
    budget.check()?;
    Ok((entries, charged))
}

fn zero_vector<D: Domain>(
    dimension: usize,
    _ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let (mut vector, bytes) = coefficient_vector::<D>(dimension, budget)?;
    vector.resize_with(dimension, || None);
    Ok((vector, bytes))
}

fn coefficient_vector<D: Domain>(
    length: usize,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let bytes = length.saturating_mul(size_of::<Option<D::Coeff>>());
    budget.reserve(bytes)?;
    Ok((Vec::with_capacity(length), bytes))
}

fn coefficient_vector_with_payload<D: Domain>(
    length: usize,
    payload: usize,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let structural = length.saturating_mul(size_of::<Option<D::Coeff>>());
    let bytes = structural.saturating_add(payload);
    budget.reserve(bytes)?;
    Ok((Vec::with_capacity(length), bytes))
}

fn coefficient_bytes<D: Domain>(ring: &PolynomialRing<D>, value: Option<&D::Coeff>) -> usize {
    value.map_or(0, |coefficient| ring.ops().heap_bytes(coefficient))
}

fn dense_structural_bytes<D: Domain>(rows: usize, columns: usize) -> usize {
    rows.saturating_mul(size_of::<Vec<Option<D::Coeff>>>())
        .saturating_add(
            rows.saturating_mul(columns)
                .saturating_mul(size_of::<Option<D::Coeff>>()),
        )
}

fn coefficient_payload<D: Domain>(ring: &PolynomialRing<D>, values: &[Option<D::Coeff>]) -> usize {
    values.iter().fold(0usize, |bytes, value| {
        bytes.saturating_add(coefficient_bytes(ring, value.as_ref()))
    })
}

fn coefficient_growth_bound(input_bytes: usize) -> usize {
    // A rational update can retain both operands, an unreduced result, and its reduction.
    input_bytes.saturating_mul(8)
}

fn settle_value<D: Domain>(
    value: Option<D::Coeff>,
    bound: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(Option<D::Coeff>, usize), ComputeError> {
    let actual = coefficient_bytes(ring, value.as_ref());
    if actual > bound {
        if let Err(error) = budget.reserve(actual - bound) {
            budget.release(bound);
            return Err(error);
        }
    } else {
        budget.release(bound - actual);
    }
    Ok((value, actual))
}

fn operation_bound<D: Domain>(
    ring: &PolynomialRing<D>,
    left: Option<&D::Coeff>,
    right: Option<&D::Coeff>,
) -> usize {
    let input = coefficient_bytes(ring, left).saturating_add(coefficient_bytes(ring, right));
    coefficient_growth_bound(input)
}

fn three_value_bound<D: Domain>(
    ring: &PolynomialRing<D>,
    first: Option<&D::Coeff>,
    second: Option<&D::Coeff>,
    third: Option<&D::Coeff>,
) -> usize {
    let input = coefficient_bytes(ring, first)
        .saturating_add(coefficient_bytes(ring, second))
        .saturating_add(coefficient_bytes(ring, third));
    coefficient_growth_bound(input)
}

fn replace_value<D: Domain>(
    slot: &mut Option<D::Coeff>,
    value: Option<D::Coeff>,
    bound: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    let old = coefficient_bytes(ring, slot.as_ref());
    let actual = coefficient_bytes(ring, value.as_ref());
    let available = bound.saturating_add(old);
    if actual > available {
        if let Err(error) = budget.reserve(actual - available) {
            budget.release(bound);
            return Err(error);
        }
    } else {
        budget.release(available - actual);
    }
    *slot = value;
    Ok(())
}

fn charged_clone<D: Domain>(
    value: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(D::Coeff, usize), ComputeError> {
    let bound = coefficient_growth_bound(coefficient_bytes(ring, Some(value)));
    budget.reserve(bound)?;
    let clone = value.clone();
    let (Some(clone), bytes) = settle_value(Some(clone), bound, ring, budget)? else {
        unreachable!("a coefficient clone is nonzero")
    };
    Ok((clone, bytes))
}

fn charged_clone_option<D: Domain>(
    value: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(Option<D::Coeff>, usize), ComputeError> {
    let bound = coefficient_growth_bound(coefficient_bytes(ring, value));
    budget.reserve(bound)?;
    settle_value(value.cloned(), bound, ring, budget)
}

fn charged_one<D: Domain>(
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(D::Coeff, usize), ComputeError> {
    let bound = size_of::<D::Coeff>();
    budget.reserve(bound)?;
    let one = ring.ops().one();
    let (Some(one), bytes) = settle_value(Some(one), bound, ring, budget)? else {
        unreachable!("one is nonzero")
    };
    Ok((one, bytes))
}

fn charged_negate<D: Domain>(
    value: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(Option<D::Coeff>, usize), ComputeError> {
    let bound = coefficient_growth_bound(coefficient_bytes(ring, value));
    budget.reserve(bound)?;
    settle_value(negate(value, ring), bound, ring, budget)
}

fn charged_inverse<D: Domain>(
    value: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(D::Coeff, usize), ComputeError> {
    let bound = coefficient_growth_bound(coefficient_bytes(ring, Some(value)));
    budget.reserve(bound)?;
    let inverse = ring.ops().inv(value);
    let (Some(inverse), bytes) = settle_value(Some(inverse), bound, ring, budget)? else {
        unreachable!("the inverse of a nonzero coefficient is nonzero")
    };
    Ok((inverse, bytes))
}

fn charged_product<D: Domain>(
    left: Option<&D::Coeff>,
    right: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(Option<D::Coeff>, usize), ComputeError> {
    let bound = operation_bound(ring, left, right);
    budget.reserve(bound)?;
    let product = multiply(left, right, ring);
    settle_value(product, bound, ring, budget)
}

fn matrix_vector<M, D>(
    matrix: &M,
    vector: &[Option<D::Coeff>],
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let dimension = matrix.dimension();
    let (mut result, bytes) = coefficient_vector::<D>(dimension, budget)?;
    for row in 0..dimension {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        result.push(dot_row(matrix, row, vector, ring, budget)?);
    }
    budget.check()?;
    let bytes = bytes.saturating_add(coefficient_payload(ring, &result));
    Ok((result, bytes))
}

fn dot_row<M, D>(
    matrix: &M,
    row: usize,
    vector: &[Option<D::Coeff>],
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<Option<D::Coeff>, ComputeError>
where
    M: MatrixView<D>,
    D: Domain,
{
    let values = (0..matrix.dimension()).map(|column| matrix.entry(row, column));
    dot_iter(values, vector, ring, budget)
}

fn dot<D: Domain>(
    row: &[Option<D::Coeff>],
    vector: &[Option<D::Coeff>],
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<Option<D::Coeff>, ComputeError> {
    dot_iter(row.iter().map(Option::as_ref), vector, ring, budget)
}

fn dot_iter<'a, D, I>(
    row: I,
    vector: &'a [Option<D::Coeff>],
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<Option<D::Coeff>, ComputeError>
where
    D: Domain,
    I: Iterator<Item = Option<&'a D::Coeff>> + Clone,
{
    let bound = coefficient_growth_bound(dot_input_bytes(row.clone(), vector, ring));
    budget.reserve(bound)?;
    let mut sum = None;
    for (index, (left, right)) in row.zip(vector).enumerate() {
        if index.is_multiple_of(TICK)
            && let Err(error) = budget.check()
        {
            budget.release(bound);
            return Err(error);
        }
        let product = multiply(left, right.as_ref(), ring);
        sum = add(sum.as_ref(), product.as_ref(), ring);
    }
    let (sum, _bytes) = settle_value(sum, bound, ring, budget)?;
    Ok(sum)
}

fn dot_input_bytes<'a, D, I>(
    row: I,
    vector: &'a [Option<D::Coeff>],
    ring: &PolynomialRing<D>,
) -> usize
where
    D: Domain,
    I: Iterator<Item = Option<&'a D::Coeff>>,
{
    row.zip(vector).fold(0usize, |bytes, (left, right)| {
        bytes
            .saturating_add(coefficient_bytes(ring, left))
            .saturating_add(coefficient_bytes(ring, right.as_ref()))
    })
}

fn rank_columns<D: Domain>(
    columns: &[ChargedVector<D>],
    rows: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<usize, ComputeError> {
    let count = columns.len();
    let (mut work, bytes) = rank_work(columns, rows, ring, budget)?;
    let rank = eliminate_rank(&mut work, rows, count, ring, budget)?;
    drop(work);
    budget.release(bytes);
    Ok(rank)
}

fn rank_work<D: Domain>(
    columns: &[ChargedVector<D>],
    rows: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let count = columns.len();
    let structural = rows
        .saturating_mul(count)
        .saturating_mul(size_of::<Option<D::Coeff>>());
    let payload = columns.iter().take(count).fold(0usize, |bytes, column| {
        bytes.saturating_add(coefficient_payload(ring, &column.0[..rows]))
    });
    let charged = structural.saturating_add(payload);
    budget.reserve(charged)?;
    let mut work = vec![None; rows.saturating_mul(count)];
    for row in 0..rows {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        for column in 0..count {
            work[row * count + column] = columns[column].0[row].clone();
        }
    }
    budget.check()?;
    Ok((work, charged))
}

fn eliminate_rank<D: Domain>(
    work: &mut [Option<D::Coeff>],
    rows: usize,
    count: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<usize, ComputeError> {
    let mut rank = 0;
    let shape = RankShape { rows, count };
    for column in 0..count {
        if column.is_multiple_of(TICK) {
            budget.check()?;
        }
        if let Some(pivot) = (rank..rows).find(|&row| work[row * count + column].is_some()) {
            eliminate_rank_pivot(work, shape, rank, column, pivot, ring, budget)?;
            rank += 1;
        }
    }
    budget.check()?;
    Ok(rank)
}

fn eliminate_rank_pivot<D: Domain>(
    work: &mut [Option<D::Coeff>],
    shape: RankShape,
    rank: usize,
    column: usize,
    pivot: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    if pivot != rank {
        swap_rank_rows::<D>(work, rank, pivot, column, shape.count, budget)?;
    }
    let pivot_value = work[rank * shape.count + column]
        .as_ref()
        .expect("the pivot row has a pivot");
    let (inverse, inverse_bytes) = charged_inverse(pivot_value, ring, budget)?;
    let result = eliminate_rank_below(work, shape, rank, column, &inverse, ring, budget);
    budget.release(inverse_bytes);
    result
}

#[derive(Clone, Copy)]
struct RankShape {
    rows: usize,
    count: usize,
}

struct RankRow {
    row: usize,
    pivot_row: usize,
    column: usize,
}

fn swap_rank_rows<D: Domain>(
    work: &mut [Option<D::Coeff>],
    first: usize,
    second: usize,
    start: usize,
    width: usize,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for right in start..width {
        if right.is_multiple_of(TICK) {
            budget.check()?;
        }
        work.swap(first * width + right, second * width + right);
    }
    budget.check()
}

fn eliminate_rank_below<D: Domain>(
    work: &mut [Option<D::Coeff>],
    shape: RankShape,
    rank: usize,
    column: usize,
    inverse: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for row in rank + 1..shape.rows {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        let Some(value) = work[row * shape.count + column].as_ref() else {
            continue;
        };
        let (Some(factor), factor_bytes) =
            charged_product(Some(value), Some(inverse), ring, budget)?
        else {
            unreachable!("a product of nonzero pivot values is nonzero")
        };
        let layout = RankRow {
            row,
            pivot_row: rank,
            column,
        };
        let result = eliminate_rank_row(work, shape, layout, &factor, ring, budget);
        budget.release(factor_bytes);
        result?;
    }
    budget.check()
}

fn eliminate_rank_row<D: Domain>(
    work: &mut [Option<D::Coeff>],
    shape: RankShape,
    layout: RankRow,
    factor: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for right in layout.column..shape.count {
        if right.is_multiple_of(TICK) {
            budget.check()?;
        }
        let target = layout.row * shape.count + right;
        let pivot = layout.pivot_row * shape.count + right;
        let bound = three_value_bound(
            ring,
            Some(factor),
            work[pivot].as_ref(),
            work[target].as_ref(),
        );
        budget.reserve(bound)?;
        let product = multiply(Some(factor), work[pivot].as_ref(), ring);
        let value = subtract(work[target].as_ref(), product.as_ref(), ring);
        drop(product);
        replace_value(&mut work[target], value, bound, ring, budget)?;
    }
    budget.check()
}

fn solve_relation<D: Domain>(
    columns: &[ChargedVector<D>],
    previous: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let rows = columns[0].0.len();
    let width = previous + 1;
    let (mut work, work_bytes) = relation_work(columns, previous, rows, width, ring, budget)?;
    let pivot_bytes = previous.saturating_mul(size_of::<usize>());
    budget.reserve(pivot_bytes)?;
    let pivots = reduce_relation(&mut work, rows, width, previous, ring, budget)?;
    let relation = extract_relation(&work, &pivots, previous, width, ring, budget)?;
    drop(work);
    drop(pivots);
    budget.release(work_bytes.saturating_add(pivot_bytes));
    Ok(relation)
}

fn relation_work<D: Domain>(
    columns: &[ChargedVector<D>],
    previous: usize,
    rows: usize,
    width: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let structural = rows
        .saturating_mul(width)
        .saturating_mul(size_of::<Option<D::Coeff>>());
    let payload = columns
        .iter()
        .take(previous + 1)
        .fold(0usize, |bytes, column| {
            bytes.saturating_add(coefficient_payload(ring, &column.0[..rows]))
        });
    let charged = structural.saturating_add(payload);
    budget.reserve(charged)?;
    let mut work = vec![None; rows.saturating_mul(width)];
    for row in 0..rows {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        for column in 0..previous {
            work[row * width + column] = columns[column].0[row].clone();
        }
        work[row * width + previous] = columns[previous].0[row]
            .as_ref()
            .map(|value| ring.ops().neg(value));
    }
    budget.check()?;
    Ok((work, charged))
}

fn reduce_relation<D: Domain>(
    work: &mut [Option<D::Coeff>],
    rows: usize,
    width: usize,
    previous: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<Vec<usize>, ComputeError> {
    let mut pivots = Vec::with_capacity(previous);
    for column in 0..previous {
        budget.check()?;
        let Some(pivot) = (pivots.len()..rows).find(|&row| work[row * width + column].is_some())
        else {
            unreachable!("the preceding Krylov columns are independent")
        };
        let row = pivots.len();
        normalize_relation_pivot(work, row, pivot, column, width, ring, budget)?;
        eliminate_relation_rows(work, row, column, rows, width, ring, budget)?;
        pivots.push(row);
    }
    budget.check()?;
    Ok(pivots)
}

fn normalize_relation_pivot<D: Domain>(
    work: &mut [Option<D::Coeff>],
    row: usize,
    pivot: usize,
    column: usize,
    width: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    if pivot != row {
        swap_relation_rows::<D>(work, row, pivot, column, width, budget)?;
    }
    let pivot_value = work[row * width + column]
        .as_ref()
        .expect("the pivot row has a pivot");
    let (inverse, inverse_bytes) = charged_inverse(pivot_value, ring, budget)?;
    let result = scale_relation_row(work, row, column, width, &inverse, ring, budget);
    budget.release(inverse_bytes);
    result
}

fn swap_relation_rows<D: Domain>(
    work: &mut [Option<D::Coeff>],
    first: usize,
    second: usize,
    start: usize,
    width: usize,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for right in start..width {
        if right.is_multiple_of(TICK) {
            budget.check()?;
        }
        work.swap(first * width + right, second * width + right);
    }
    budget.check()
}

fn scale_relation_row<D: Domain>(
    work: &mut [Option<D::Coeff>],
    row: usize,
    column: usize,
    width: usize,
    inverse: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for right in column..width {
        if right.is_multiple_of(TICK) {
            budget.check()?;
        }
        let target = row * width + right;
        let bound = operation_bound(ring, work[target].as_ref(), Some(inverse));
        budget.reserve(bound)?;
        let value = scale(work[target].as_ref(), inverse, ring);
        replace_value(&mut work[target], value, bound, ring, budget)?;
    }
    budget.check()
}

fn eliminate_relation_rows<D: Domain>(
    work: &mut [Option<D::Coeff>],
    pivot_row: usize,
    column: usize,
    rows: usize,
    width: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for row in 0..rows {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        if row == pivot_row {
            continue;
        }
        let Some(factor_value) = work[row * width + column].as_ref() else {
            continue;
        };
        let (factor, factor_bytes) = charged_clone(factor_value, ring, budget)?;
        let layout = RelationRow {
            row,
            pivot_row,
            column,
            width,
        };
        let result = eliminate_relation_row(work, layout, &factor, ring, budget);
        budget.release(factor_bytes);
        result?;
    }
    budget.check()
}

struct RelationRow {
    row: usize,
    pivot_row: usize,
    column: usize,
    width: usize,
}

fn eliminate_relation_row<D: Domain>(
    work: &mut [Option<D::Coeff>],
    layout: RelationRow,
    factor: &D::Coeff,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<(), ComputeError> {
    for right in layout.column..layout.width {
        if right.is_multiple_of(TICK) {
            budget.check()?;
        }
        let target = layout.row * layout.width + right;
        let pivot = layout.pivot_row * layout.width + right;
        let bound = three_value_bound(
            ring,
            Some(factor),
            work[pivot].as_ref(),
            work[target].as_ref(),
        );
        budget.reserve(bound)?;
        let product = multiply(Some(factor), work[pivot].as_ref(), ring);
        let value = subtract(work[target].as_ref(), product.as_ref(), ring);
        drop(product);
        replace_value(&mut work[target], value, bound, ring, budget)?;
    }
    budget.check()
}

fn extract_relation<D: Domain>(
    work: &[Option<D::Coeff>],
    pivots: &[usize],
    previous: usize,
    width: usize,
    ring: &PolynomialRing<D>,
    budget: &mut QuotientBudget<'_>,
) -> Result<ChargedVector<D>, ComputeError> {
    let payload = pivots.iter().fold(0usize, |bytes, &row| {
        bytes.saturating_add(coefficient_bytes(
            ring,
            work[row * width + previous].as_ref(),
        ))
    });
    let (mut relation, mut relation_bytes) =
        coefficient_vector_with_payload::<D>(width, payload, budget)?;
    for &row in pivots {
        if row.is_multiple_of(TICK) {
            budget.check()?;
        }
        relation.push(work[row * width + previous].clone());
    }
    budget.check()?;
    let (one, one_bytes) = charged_one(ring, budget)?;
    relation.push(Some(one));
    relation_bytes = relation_bytes.saturating_add(one_bytes);
    Ok((relation, relation_bytes))
}

fn release_columns<D: Domain>(
    columns: &mut Vec<ChargedVector<D>>,
    budget: &mut QuotientBudget<'_>,
) {
    for (_, bytes) in columns.drain(..) {
        budget.release(bytes);
    }
}

fn negate<D: Domain>(value: Option<&D::Coeff>, ring: &PolynomialRing<D>) -> Option<D::Coeff> {
    value.map(|value| ring.ops().neg(value))
}

fn add<D: Domain>(
    left: Option<&D::Coeff>,
    right: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
) -> Option<D::Coeff> {
    match (left, right) {
        (None, None) => None,
        (Some(value), None) | (None, Some(value)) => Some(value.clone()),
        (Some(left), Some(right)) => ring.ops().add(left, right),
    }
}

fn subtract<D: Domain>(
    left: Option<&D::Coeff>,
    right: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
) -> Option<D::Coeff> {
    match (left, right) {
        (None, None) => None,
        (Some(value), None) => Some(value.clone()),
        (None, Some(value)) => Some(ring.ops().neg(value)),
        (Some(left), Some(right)) => ring.ops().sub(left, right),
    }
}

fn multiply<D: Domain>(
    left: Option<&D::Coeff>,
    right: Option<&D::Coeff>,
    ring: &PolynomialRing<D>,
) -> Option<D::Coeff> {
    match (left, right) {
        (Some(left), Some(right)) => Some(ring.ops().mul(left, right)),
        _ => None,
    }
}

fn scale<D: Domain>(
    value: Option<&D::Coeff>,
    scalar: &D::Coeff,
    ring: &PolynomialRing<D>,
) -> Option<D::Coeff> {
    value.map(|value| ring.ops().mul(value, scalar))
}
