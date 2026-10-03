//! Linear algebra over the coefficient field of a finite quotient.

use std::fmt;
use std::mem::size_of;

use crate::compute::{ComputeError, ComputeLimits, TICK};
use crate::ring::{Domain, DomainOps, PolynomialRing, PrimeField};

/// The resource meter shared by one finite quotient operation.
///
/// The caller charges the quotient basis before it constructs this meter.
/// Temporary vectors and matrices charge themselves through [`reserve`].
/// The count estimates live allocations and is not a process memory cap.
pub(crate) struct QuotientBudget<'a> {
    limits: &'a ComputeLimits,
    held: usize,
}

impl<'a> QuotientBudget<'a> {
    /// Start a meter with bytes the caller already holds.
    pub(crate) fn new(limits: &'a ComputeLimits, held: usize) -> Result<Self, ComputeError> {
        let meter = QuotientBudget { limits, held };
        meter.check()?;
        Ok(meter)
    }

    /// Charge an allocation before it is made.
    pub(crate) fn reserve(&mut self, bytes: usize) -> Result<(), ComputeError> {
        self.check_stop()?;
        if let Some(limit) = self.limits.memory
            && self.held.saturating_add(bytes) > limit
        {
            return Err(ComputeError::MemoryLimitExceeded);
        }
        self.held = self.held.saturating_add(bytes);
        Ok(())
    }

    /// Release bytes after an allocation is dropped.
    pub(crate) fn release(&mut self, bytes: usize) {
        self.held = self.held.saturating_sub(bytes);
    }

    /// Check the deadline and the current memory estimate.
    pub(crate) fn check(&self) -> Result<(), ComputeError> {
        self.check_stop()?;
        if let Some(limit) = self.limits.memory
            && self.held > limit
        {
            return Err(ComputeError::MemoryLimitExceeded);
        }
        Ok(())
    }

    /// Check the deadline, including cancellation used by internal callers.
    fn check_stop(&self) -> Result<(), ComputeError> {
        match self.limits.stop() {
            Some(stop) => Err(stop.reported()),
            None => Ok(()),
        }
    }
}

/// A matrix for multiplication by one residue class.
///
/// The matrix is square. Its columns are the coordinates of multiplication
/// by the corresponding basis element, and its entries use row-major order.
/// `None` represents zero. The ring supplies the coefficient arithmetic and
/// the display format.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MultiplicationMatrix<D: Domain = PrimeField> {
    ring: PolynomialRing<D>,
    dimension: usize,
    entries: Vec<Option<D::Coeff>>,
}

impl<D: Domain> MultiplicationMatrix<D> {
    /// Build a matrix from columns of coordinates.
    ///
    /// The first coordinate is conventionally the residue class of `1`.
    /// The constructor stores the columns in row-major order.
    pub(crate) fn from_columns(
        ring: PolynomialRing<D>,
        columns: Vec<Vec<Option<D::Coeff>>>,
        budget: &mut QuotientBudget<'_>,
    ) -> Result<Self, ComputeError> {
        let dimension = columns.len();
        assert!(
            columns.iter().all(|column| column.len() == dimension),
            "matrix columns must have one entry per row"
        );
        let slots = dimension.saturating_mul(dimension);
        let structural = slots.saturating_mul(size_of::<Option<D::Coeff>>());
        let payload = columns.iter().fold(0usize, |bytes, column| {
            bytes.saturating_add(coefficient_payload(&ring, column))
        });
        budget.reserve(structural.saturating_add(payload))?;
        let mut entries = Vec::with_capacity(slots);
        for row in 0..dimension {
            if row.is_multiple_of(TICK) {
                budget.check()?;
            }
            for column in &columns {
                entries.push(column[row].clone());
            }
        }
        budget.check()?;
        Ok(MultiplicationMatrix {
            ring,
            dimension,
            entries,
        })
    }

    /// The coefficient ring used by the matrix.
    pub fn ring(&self) -> &PolynomialRing<D> {
        &self.ring
    }

    /// The number of rows and columns.
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// The entries in row-major order.
    pub fn entries(&self) -> &[Option<D::Coeff>] {
        &self.entries
    }

    /// Read one entry, or `None` for an out-of-range coordinate or zero.
    pub fn entry(&self, row: usize, column: usize) -> Option<&D::Coeff> {
        (row < self.dimension && column < self.dimension)
            .then(|| self.entries[row * self.dimension + column].as_ref())
            .flatten()
    }

    /// Take the row-major entries.
    pub fn into_entries(self) -> Vec<Option<D::Coeff>> {
        self.entries
    }

    /// Estimate the bytes of the entries allocation and its coefficients.
    pub(crate) fn heap_bytes(&self) -> usize {
        entries_bytes(&self.ring, &self.entries)
    }

    /// Return one entry by row and column for internal arithmetic.
    pub(crate) fn entry_value(&self, row: usize, column: usize) -> Option<&D::Coeff> {
        self.entries[row * self.dimension + column].as_ref()
    }
}

fn coefficient_payload<D: Domain>(ring: &PolynomialRing<D>, values: &[Option<D::Coeff>]) -> usize {
    values.iter().fold(0usize, |bytes, value| {
        bytes.saturating_add(
            value
                .as_ref()
                .map_or(0, |coefficient| ring.ops().heap_bytes(coefficient)),
        )
    })
}

/// Estimate the bytes held by a row-major coefficient vector.
pub(crate) fn entries_bytes<D: Domain>(
    ring: &PolynomialRing<D>,
    entries: &[Option<D::Coeff>],
) -> usize {
    let values = entries.iter().fold(0usize, |bytes, value| {
        bytes.saturating_add(
            value
                .as_ref()
                .map_or(0, |coefficient| ring.ops().heap_bytes(coefficient)),
        )
    });
    entries
        .len()
        .saturating_mul(size_of::<Option<D::Coeff>>())
        .saturating_add(values)
}

impl<D: Domain> fmt::Display for MultiplicationMatrix<D> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        for row in 0..self.dimension {
            if row > 0 {
                f.write_str("; ")?;
            }
            self.write_row(row, f)?;
        }
        f.write_str("]")
    }
}

impl<D: Domain> MultiplicationMatrix<D> {
    fn write_row(&self, row: usize, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        for column in 0..self.dimension {
            if column > 0 {
                f.write_str(", ")?;
            }
            self.write_entry(row, column, f)?;
        }
        f.write_str("]")
    }

    fn write_entry(&self, row: usize, column: usize, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.entry_value(row, column) {
            Some(coefficient) => self.ring.ops().write(coefficient, f),
            None => f.write_str("0"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::compute::{Budget, ComputeLimits};
    use crate::ring::{Domain, DomainOps, Rationals};

    #[test]
    fn rational_matrix_payload_is_reserved_before_cloning() {
        let ring = PolynomialRing::<Rationals>::rationals(["x"]).expect("the rational ring builds");
        let source = ring
            .parse_polynomial("123456789012345678901234567890*x")
            .expect("the coefficient parses");
        let coefficient = source
            .leading_term()
            .expect("the source has a leading term")
            .0
            .clone();
        let payload = ring.ops().heap_bytes(&coefficient);
        assert!(payload > 0);
        let structural = size_of::<Option<<Rationals as Domain>::Coeff>>();
        let limits =
            ComputeLimits::of_budget(&Budget::new().memory_limit(structural + payload - 1));
        let mut budget = QuotientBudget::new(&limits, 0).expect("the empty meter fits");
        let error = MultiplicationMatrix::<Rationals>::from_columns(
            ring,
            vec![vec![Some(coefficient)]],
            &mut budget,
        )
        .expect_err("the payload must be charged before the clone");
        assert_eq!(error, ComputeError::MemoryLimitExceeded);
    }
}
