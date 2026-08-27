//! Module signatures and the labeled polynomials that carry them.

use std::cmp::Ordering;

use crate::poly::Monomial;

/// The module monomial that labels one polynomial, ordered position over term.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Signature {
    pub(crate) index: usize,
    pub(crate) term: Monomial,
}

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

/// Report whether a syzygy rule already covers `sig`.
pub(super) fn is_syzygy(sig: &Signature, rules: &[Signature]) -> bool {
    rules
        .iter()
        .any(|rule| rule.index == sig.index && rule.term.divides(&sig.term))
}

/// Record `sig` as a syzygy, keeping the rule set minimal under division.
pub(super) fn add_syzygy_rule(rules: &mut Vec<Signature>, sig: Signature) {
    if is_syzygy(&sig, rules) {
        return;
    }
    rules.retain(|rule| !(rule.index == sig.index && sig.term.divides(&rule.term)));
    rules.push(sig);
}

/// Report whether some basis element `g` makes `p` sig-redundant.
///
/// True when `sig(g) | sig(p)` and `lm(g) | lm(p)`, with the two quotients
/// independent (Arri-Perry; Eder-Faugère survey). Call this only on a
/// regular normal form: the completeness of the drop rests on `p` having
/// finished `f5_reduce`. With a = sig(p)/sig(g) and b = lm(p)/lm(g), a > b
/// is impossible, since b*g would still be a legal regular top-reducer of
/// lm(p), so a <= b. If a < b, subtracting a coefficient-scaled a*g leaves
/// the same lead at a strictly smaller signature. If a = b, the witness
/// multiple already carries this signature and lead, the singular case.
/// Either way the drop loses neither a new lead nor signature coverage.
///
/// Termination: per signature index, the accepted (signature term, lead)
/// pairs form a Dickson-bad sequence in N^(2n), so no earlier pair divides
/// a later one componentwise and the sequence is finite. Finitely many
/// indices keep the basis and the pair queues finite. Accepting the
/// equal-quotient (singular) case would instead create equal-lead,
/// equal-signature duplicates that multiply without bound. Every untouched
/// reducer multiple m*g falls in that singular class.
pub(super) fn is_sig_redundant(p: &LabeledPoly, basis: &[LabeledPoly]) -> bool {
    let Some(lm_p) = p.poly.lm() else {
        return false;
    };
    basis.iter().any(|g| {
        g.sig.index == p.sig.index
            && g.sig.term.divides(&p.sig.term)
            && g.poly.lm().is_some_and(|lm_g| lm_g.divides(lm_p))
    })
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
