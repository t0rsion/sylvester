//! Module signatures and the labeled polynomials that carry them.

use std::cmp::Ordering;

use crate::poly::Monomial;

/// The module monomial that labels one polynomial.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Signature {
    pub(crate) index: usize,
    pub(crate) term: Monomial,
}

// POT (Position Over Term).
impl Ord for Signature {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.index.cmp(&other.index) {
            Ordering::Equal => self.term.cmp(&other.term),
            ord => ord,
        }
    }
}

impl PartialOrd for Signature {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Signature {
    /// Report whether `self * (num / den)` is below `target`.
    ///
    /// This is the regular reduction test of F5: a reducer may cancel a
    /// leading term only while its own signature stays below the
    /// signature of the polynomial it reduces. The order is POT, so the
    /// index decides first and the term multiple is never built.
    pub(crate) fn shifted_is_below(&self, num: &Monomial, den: &Monomial, target: &Self) -> bool {
        match self.index.cmp(&target.index) {
            Ordering::Equal => self.term.cmp_shifted(num, den, &target.term) == Ordering::Less,
            ord => ord == Ordering::Less,
        }
    }
}

/// One polynomial with its signature and its place in the basis.
#[derive(Clone, Debug)]
pub(crate) struct LabeledPoly {
    pub(crate) sig: Signature,
    pub(crate) poly: crate::poly::Polynomial,
    pub(crate) index: usize,
}

#[cfg(test)]
mod tests {
    use super::Signature;
    use crate::poly::Monomial;
    use smallvec::smallvec;

    #[test]
    fn signature_ordering_is_pot() {
        let a = Signature {
            index: 0,
            term: Monomial::from_exps(smallvec![1u16, 0]),
        };
        let b = Signature {
            index: 1,
            term: Monomial::from_exps(smallvec![0u16, 0]),
        };
        assert!(a < b);

        let c = Signature {
            index: 0,
            term: Monomial::from_exps(smallvec![2u16, 0]),
        };
        assert!(a < c);
    }
}
