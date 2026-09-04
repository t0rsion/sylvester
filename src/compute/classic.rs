//! Classic F5: one critical pair at a time, in signature order.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;

use super::interreduce::{reduced_groebner_basis_checked, reduced_groebner_basis_tracked};
use super::signature::{LabeledPoly, Signature};
use super::{ComputeError, ComputeLimits, RunError};
use crate::cert::origin::{self, Origin};
use crate::poly::{ExponentOverflow, Monomial, Polynomial, Term};
use crate::ring::field::Felt;
use crate::ring::{PolynomialRing, PrimeOps};

/// One term multiple `coeff * mono * g_index` of a held polynomial.
struct Multiple<'a> {
    index: usize,
    coeff: Felt,
    mono: &'a Monomial,
}

/// Whether a run carries origin cofactors.
///
/// The engine body is written once and instantiated twice. [`Untracked`]
/// has a zero-sized origin and empty operations, so an uncertified run
/// carries no per-step cost. [`Tracked`] keeps one cofactor vector per
/// value the engine holds.
trait Track {
    /// The cofactors of one polynomial.
    type Origin: Clone;
    /// The input the origin identity is checked against.
    type Witness;

    fn witness(generators: &[Polynomial]) -> Self::Witness;
    fn unit(ring: &PolynomialRing, index: usize, count: usize) -> Self::Origin;
    /// The bytes the held cofactors occupy.
    fn bytes(origins: &[Self::Origin]) -> usize;
    fn spair(
        origins: &[Self::Origin],
        left: Multiple,
        right: Multiple,
        ops: &PrimeOps,
    ) -> Result<Self::Origin, ExponentOverflow>;
    fn sub_scaled(
        target: &mut Self::Origin,
        sources: &[Self::Origin],
        source: usize,
        coeff: Felt,
        mono: &Monomial,
        ops: &PrimeOps,
    ) -> Result<(), ExponentOverflow>;
    fn monic(target: &mut Self::Origin, poly: &Polynomial, ops: &PrimeOps);
    fn holds(
        origin: &Self::Origin,
        witness: &Self::Witness,
        poly: &Polynomial,
        ops: &PrimeOps,
    ) -> bool;
}

struct Untracked;

impl Track for Untracked {
    type Origin = ();
    type Witness = ();

    fn witness(_generators: &[Polynomial]) {}
    fn unit(_ring: &PolynomialRing, _index: usize, _count: usize) {}
    fn bytes(_origins: &[()]) -> usize {
        0
    }
    fn spair(
        _origins: &[()],
        _left: Multiple,
        _right: Multiple,
        _ops: &PrimeOps,
    ) -> Result<(), ExponentOverflow> {
        Ok(())
    }
    fn sub_scaled(
        _target: &mut (),
        _sources: &[()],
        _source: usize,
        _coeff: Felt,
        _mono: &Monomial,
        _ops: &PrimeOps,
    ) -> Result<(), ExponentOverflow> {
        Ok(())
    }
    fn monic(_target: &mut (), _poly: &Polynomial, _ops: &PrimeOps) {}
    fn holds(_origin: &(), _witness: &(), _poly: &Polynomial, _ops: &PrimeOps) -> bool {
        true
    }
}

struct Tracked;

impl Track for Tracked {
    type Origin = Origin;
    type Witness = Vec<Polynomial>;

    fn witness(generators: &[Polynomial]) -> Vec<Polynomial> {
        generators.to_vec()
    }

    fn unit(ring: &PolynomialRing, index: usize, count: usize) -> Origin {
        origin::unit(ring, index, count)
    }

    fn bytes(origins: &[Origin]) -> usize {
        let mut bytes = std::mem::size_of_val(origins);
        for entry in origins {
            bytes = bytes.saturating_add(std::mem::size_of_val(entry.as_slice()));
            for cofactor in entry {
                bytes = bytes.saturating_add(cofactor.heap_bytes());
            }
        }
        bytes
    }

    fn spair(
        origins: &[Origin],
        left: Multiple,
        right: Multiple,
        ops: &PrimeOps,
    ) -> Result<Origin, ExponentOverflow> {
        origin::combine(
            &origins[left.index],
            left.coeff,
            left.mono,
            &origins[right.index],
            right.coeff,
            right.mono,
            ops,
        )
    }

    fn sub_scaled(
        target: &mut Origin,
        sources: &[Origin],
        source: usize,
        coeff: Felt,
        mono: &Monomial,
        ops: &PrimeOps,
    ) -> Result<(), ExponentOverflow> {
        origin::sub_scaled(target, &sources[source], coeff, mono, ops)
    }

    fn monic(target: &mut Origin, poly: &Polynomial, ops: &PrimeOps) {
        if let Some(lc) = poly.lc() {
            origin::scale(target, lc.inv(ops.modulus()), ops);
        }
    }

    fn holds(
        origin: &Origin,
        witness: &Vec<Polynomial>,
        poly: &Polynomial,
        ops: &PrimeOps,
    ) -> bool {
        origin::holds(origin, witness, poly, ops)
    }
}

fn is_syzygy(sig: &Signature, rules: &[Signature]) -> bool {
    rules
        .iter()
        .any(|rule| rule.index == sig.index && rule.term.divides(&sig.term))
}

fn add_syzygy_rule(rules: &mut Vec<Signature>, sig: Signature) {
    if is_syzygy(&sig, rules) {
        return;
    }

    rules.retain(|rule| !(rule.index == sig.index && sig.term.divides(&rule.term)));
    rules.push(sig);
}

fn f5_reduce<T: Track>(
    mut p: LabeledPoly,
    mut cofactors: T::Origin,
    basis: &[LabeledPoly],
    basis_origins: &[T::Origin],
    ops: &PrimeOps,
) -> Result<(LabeledPoly, T::Origin), ExponentOverflow> {
    // The remainder grows by the largest term left, so it is built
    // descending and turned around once.
    let mut remainder: Vec<Term> = Vec::new();
    let leads: Vec<(u32, u64)> = basis
        .iter()
        .map(|g| {
            g.poly
                .lm()
                .map_or((u32::MAX, 0), |lm| (lm.deg, lm.var_mask()))
        })
        .collect();
    let reducers = Reducers {
        basis,
        origins: basis_origins,
        leads: &leads,
        ops,
    };

    while let Some(lt_p) = p.poly.lt().cloned() {
        if !reduce_lead::<T>(&mut p, &mut cofactors, &lt_p, &reducers)? {
            let term = p.poly.pop_lt().expect("polynomial should not be empty");
            remainder.push(term);
        }
    }

    remainder.reverse();
    p.poly = Polynomial::from_sorted_terms(p.poly.ring().clone(), remainder);
    Ok((p, cofactors))
}

/// True when some basis element `g` has `sig(g) | sig(p)` and
/// `lm(g) | lm(p)` (with independent quotients). Such a `p` is
/// sig-redundant (Arri-Perry; Eder-Faugère survey). The caller must only
/// call this on a regular normal form: completeness of the drop rests on
/// `p` having finished `f5_reduce`. Then, with a = sig(p)/sig(g) and b =
/// lm(p)/lm(g), a > b is impossible (b*g would still be a legal regular
/// top-reducer of lm(p)), so a <= b. If a < b, subtracting a
/// coefficient-scaled a*g leaves the same lead at strictly smaller
/// signature; if a = b, the witness multiple already has the same
/// signature and lead, which is the singular-criterion case. Thus the drop
/// loses neither a new lead nor signature coverage. For termination: per
/// signature index, the chronological sequence of accepted (signature
/// term, lead) pairs contains no earlier pair that componentwise divides
/// a later one, so it is a Dickson-bad sequence in N^(2n) and must be
/// finite. With finitely many signature indices, the basis stays finite,
/// and so does the pair queue. Accepting the equal-quotient (singular)
/// case would create equal-lead, equal-signature duplicates that multiply
/// without bound.
fn is_sig_redundant(p: &LabeledPoly, basis: &[LabeledPoly]) -> bool {
    let Some(lm_p) = p.poly.lm() else {
        return false;
    };
    basis.iter().any(|g| {
        g.sig.index == p.sig.index
            && g.sig.term.divides(&p.sig.term)
            && g.poly.lm().is_some_and(|lm_g| lm_g.divides(lm_p))
    })
}

struct Reducers<'a, T: Track> {
    basis: &'a [LabeledPoly],
    origins: &'a [T::Origin],
    leads: &'a [(u32, u64)],
    ops: &'a PrimeOps,
}

fn reduce_lead<T: Track>(
    target: &mut LabeledPoly,
    cofactors: &mut T::Origin,
    lead: &Term,
    reducers: &Reducers<'_, T>,
) -> Result<bool, ExponentOverflow> {
    for (index, reducer) in reducers.basis.iter().enumerate() {
        let Some(mono) = reducer_multiple(reducer, reducers.leads[index], lead, &target.sig) else {
            continue;
        };
        let reducer_lead = reducer.poly.lt().expect("a selected reducer is nonzero");
        let scale = lead.coeff.div(reducer_lead.coeff, reducers.ops.modulus());
        target.poly = target
            .poly
            .sub_scaled_checked(&reducer.poly, &scale, &mono, reducers.ops)?;
        T::sub_scaled(
            cofactors,
            reducers.origins,
            index,
            scale,
            &mono,
            reducers.ops,
        )?;
        return Ok(true);
    }
    Ok(false)
}

fn reducer_multiple(
    reducer: &LabeledPoly,
    lead: (u32, u64),
    target: &Term,
    signature: &Signature,
) -> Option<Monomial> {
    if lead.0 > target.mono.deg || lead.1 & !target.mono.var_mask() != 0 {
        return None;
    }
    let reducer_lead = reducer.poly.lt()?;
    if !reducer_lead.mono.divides(&target.mono) {
        return None;
    }
    if !reducer
        .sig
        .shifted_is_below(&target.mono, &reducer_lead.mono, signature)
    {
        return None;
    }
    target.mono.quotient(&reducer_lead.mono)
}

fn spolynomial<T: Track>(
    basis: &[LabeledPoly],
    origins: &[T::Origin],
    i: usize,
    j: usize,
    ops: &PrimeOps,
) -> Result<(Polynomial, T::Origin), ExponentOverflow> {
    let (f, g) = (&basis[i], &basis[j]);
    // f and g are non-zero when called from the F5 pipeline.
    let lt_f = f.poly.lt().expect("non-zero polynomial must have lt");
    // f and g are non-zero when called from the F5 pipeline.
    let lt_g = g.poly.lt().expect("non-zero polynomial must have lt");
    let lcm = lt_f.mono.lcm(&lt_g.mono);
    // lcm is a multiple of lm(f).
    let m_f = lcm.quotient(&lt_f.mono).expect("lm(f) divides lcm");
    // lcm is a multiple of lm(g).
    let m_g = lcm.quotient(&lt_g.mono).expect("lm(g) divides lcm");

    let f_scaled = f.poly.scale_monomial_checked(&lt_g.coeff, &m_f, ops)?;
    let g_scaled = g.poly.scale_monomial_checked(&lt_f.coeff, &m_g, ops)?;
    let cofactors = T::spair(
        origins,
        Multiple {
            index: i,
            coeff: lt_g.coeff,
            mono: &m_f,
        },
        Multiple {
            index: j,
            coeff: lt_f.coeff,
            mono: &m_g,
        },
        ops,
    )?;
    Ok((f_scaled.sub(&g_scaled, ops), cofactors))
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct CriticalPair {
    sig: Signature,
    i: usize,
    j: usize,
}

impl Ord for CriticalPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.sig
            .cmp(&other.sig)
            .then_with(|| self.i.cmp(&other.i))
            .then_with(|| self.j.cmp(&other.j))
    }
}

impl PartialOrd for CriticalPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn check_limits<T: Track>(
    limits: &ComputeLimits,
    basis: &[LabeledPoly],
    origins: &[T::Origin],
    syzygy_rules: &[Signature],
    queue_len: usize,
) -> Result<(), RunError> {
    if let Some(stop) = limits.stop() {
        return Err(stop);
    }

    if let Some(limit) = limits.memory {
        let mut bytes = 0usize;
        bytes = bytes.saturating_add(std::mem::size_of_val(basis));
        bytes = bytes.saturating_add(basis.iter().map(|lp| lp.poly.heap_bytes()).sum::<usize>());
        bytes = bytes.saturating_add(std::mem::size_of_val(syzygy_rules));
        bytes = bytes.saturating_add(queue_len * size_of::<CriticalPair>());
        bytes = bytes.saturating_add(T::bytes(origins));
        if bytes > limit {
            return Err(RunError::Compute(ComputeError::MemoryLimitExceeded));
        }
    }

    Ok(())
}

fn push_pair(
    queue: &mut BinaryHeap<std::cmp::Reverse<CriticalPair>>,
    syzygy_rules: &[Signature],
    basis: &[LabeledPoly],
    i: usize,
    j: usize,
) -> Result<(), RunError> {
    let f = &basis[i];
    let g = &basis[j];
    let Some(lm_f) = f.poly.lm() else {
        return Ok(());
    };
    let Some(lm_g) = g.poly.lm() else {
        return Ok(());
    };
    let lcm = lm_f.lcm(lm_g);
    if lcm.deg > super::DEGREE_LIMIT {
        return Err(RunError::Compute(ComputeError::DegreeLimit {
            limit: super::DEGREE_LIMIT,
        }));
    }

    // lcm is a multiple of lm_f.
    let m_f = lcm.quotient(lm_f).expect("lm_f divides lcm");
    // lcm is a multiple of lm_g.
    let m_g = lcm.quotient(lm_g).expect("lm_g divides lcm");

    let sig_f = Signature {
        index: f.sig.index,
        term: f.sig.term.checked_mul(&m_f)?,
    };
    let sig_g = Signature {
        index: g.sig.index,
        term: g.sig.term.checked_mul(&m_g)?,
    };
    let sig = match sig_f.cmp(&sig_g) {
        Ordering::Greater => sig_f,
        Ordering::Less => sig_g,
        // Non-regular pair: the component signatures agree as module
        // monomials (coefficients are not tracked, so this covers both the
        // singular case, where the module leading terms cancel and the
        // S-polynomial's true signature is strictly smaller, and the
        // super-regular case). Either way the common value is not a trusted
        // signature for the S-polynomial: keeping the pair could mislabel it
        // and, on a zero reduction, record an inflated syzygy signature that
        // later discards necessary pairs. The signature-based Buchberger
        // criterion only requires regular S-pairs, so rejecting is sound.
        Ordering::Equal => return Ok(()),
    };

    if is_syzygy(&sig, syzygy_rules) {
        return Ok(());
    }

    queue.push(std::cmp::Reverse(CriticalPair { sig, i, j }));
    Ok(())
}

struct ClassicState<T: Track> {
    basis: Vec<LabeledPoly>,
    origins: Vec<T::Origin>,
    syzygy_rules: Vec<Signature>,
    queue: BinaryHeap<std::cmp::Reverse<CriticalPair>>,
    witness: T::Witness,
    ops: PrimeOps,
    nvars: usize,
    limits: ComputeLimits,
}

impl<T: Track> ClassicState<T> {
    fn new(ring: &PolynomialRing, generators: &[Polynomial], limits: &ComputeLimits) -> Self {
        ClassicState {
            basis: Vec::new(),
            origins: Vec::new(),
            syzygy_rules: Vec::new(),
            queue: BinaryHeap::new(),
            witness: T::witness(generators),
            ops: *ring.ops(),
            nvars: ring.nvars(),
            limits: limits.clone(),
        }
    }

    fn check_limits(&self) -> Result<(), RunError> {
        check_limits::<T>(
            &self.limits,
            &self.basis,
            &self.origins,
            &self.syzygy_rules,
            self.queue.len(),
        )
    }

    fn seed(
        &mut self,
        ring: &PolynomialRing,
        generator: &Polynomial,
        index: usize,
        count: usize,
    ) -> Result<(), RunError> {
        self.check_limits()?;
        let signature = Signature {
            index,
            term: Monomial::one(self.nvars),
        };
        let mut origin = T::unit(ring, index, count);
        T::monic(&mut origin, generator, &self.ops);
        let labeled = LabeledPoly {
            sig: signature,
            poly: generator.make_monic(&self.ops),
            index: self.basis.len(),
        };
        let (reduced, origin) =
            f5_reduce::<T>(labeled, origin, &self.basis, &self.origins, &self.ops)?;
        if reduced.poly.is_zero() {
            add_syzygy_rule(&mut self.syzygy_rules, reduced.sig);
            return Ok(());
        }
        self.insert(reduced, origin)
    }

    fn insert(&mut self, mut labeled: LabeledPoly, mut origin: T::Origin) -> Result<(), RunError> {
        T::monic(&mut origin, &labeled.poly, &self.ops);
        labeled.poly = labeled.poly.make_monic(&self.ops);
        debug_assert!(
            T::holds(&origin, &self.witness, &labeled.poly, &self.ops),
            "a basis element must equal the combination its cofactors name"
        );
        labeled.index = self.basis.len();
        let new_index = labeled.index;
        self.basis.push(labeled);
        self.origins.push(origin);
        for index in 0..new_index {
            push_pair(
                &mut self.queue,
                &self.syzygy_rules,
                &self.basis,
                index,
                new_index,
            )?;
        }
        Ok(())
    }

    fn drain_pairs(&mut self) -> Result<(), RunError> {
        while let Some(std::cmp::Reverse(pair)) = self.queue.pop() {
            self.process_pair(pair)?;
        }
        Ok(())
    }

    fn process_pair(&mut self, pair: CriticalPair) -> Result<(), RunError> {
        self.check_limits()?;
        if is_syzygy(&pair.sig, &self.syzygy_rules) {
            return Ok(());
        }
        let (polynomial, origin) =
            spolynomial::<T>(&self.basis, &self.origins, pair.i, pair.j, &self.ops)?;
        if polynomial.is_zero() {
            add_syzygy_rule(&mut self.syzygy_rules, pair.sig);
            return Ok(());
        }
        let labeled = LabeledPoly {
            sig: pair.sig,
            poly: polynomial,
            index: self.basis.len(),
        };
        let (reduced, origin) =
            f5_reduce::<T>(labeled, origin, &self.basis, &self.origins, &self.ops)?;
        if reduced.poly.is_zero() {
            add_syzygy_rule(&mut self.syzygy_rules, reduced.sig);
            return Ok(());
        }
        if is_sig_redundant(&reduced, &self.basis) {
            return Ok(());
        }
        self.insert(reduced, origin)
    }

    fn finish(self) -> (Vec<Polynomial>, Vec<T::Origin>) {
        (
            self.basis.into_iter().map(|labeled| labeled.poly).collect(),
            self.origins,
        )
    }
}

fn solve_raw_with_limits<T: Track>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<Polynomial>, Vec<T::Origin>), RunError> {
    super::check_input_degrees(generators)?;
    let count = generators.len();
    let mut state = ClassicState::<T>::new(ring, generators, limits);
    for (index, generator) in generators.iter().enumerate() {
        state.seed(ring, generator, index, count)?;
        state.drain_pairs()?;
    }
    state.check_limits()?;
    Ok(state.finish())
}

/// Run the classic backend under the limits of the call.
pub(super) fn solve_checked(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial>, RunError> {
    let (raw, _) = solve_raw_with_limits::<Untracked>(ring, generators, limits)?;
    reduced_groebner_basis_checked(ring, raw, limits)
}

/// Run the classic backend and keep the origin cofactors of the result.
///
/// Returns the reduced basis and one cofactor vector per basis element, so
/// that g_j = sum_i c_ji * f_i over `generators`.
///
/// The tracking runs only here. [`solve_checked`] carries none of it.
pub(super) fn solve_tracked(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    limits: &ComputeLimits,
) -> Result<(Vec<Polynomial>, Vec<Origin>), RunError> {
    let (raw, origins) = solve_raw_with_limits::<Tracked>(ring, generators, limits)?;
    reduced_groebner_basis_tracked(ring, raw, origins, generators, limits)
}

/// Run the classic backend with no budget.
#[cfg(test)]
pub(crate) fn solve(ring: &PolynomialRing, generators: &[Polynomial]) -> Vec<Polynomial> {
    // unlimited runs should only fail on internal invariants.
    solve_checked(ring, generators, &ComputeLimits::default())
        .expect("a run without a budget cannot stop early")
}

#[cfg(test)]
mod tests {
    use super::{Untracked, f5_reduce, solve, solve_checked, spolynomial};
    use crate::compute::signature::{LabeledPoly, Signature};
    use crate::compute::{ComputeError, ComputeLimits, RunError};
    use crate::poly::Monomial;
    use crate::ring::PolynomialRing;

    fn ring() -> PolynomialRing {
        PolynomialRing::prime_field(7, ["x", "y"]).expect("7 is prime")
    }

    #[test]
    fn a_pair_past_the_degree_limit_is_a_typed_error() {
        let ring = ring();
        let f = ring
            .polynomial([(1, [65535, 0]), (1, [0, 65535])])
            .expect("fits");
        let g = ring.polynomial([(1, [0, 65535])]).expect("fits");
        assert_eq!(
            solve_checked(&ring, &[f, g], &ComputeLimits::default()),
            Err(RunError::Compute(ComputeError::DegreeLimit {
                limit: 65535
            }))
        );
    }

    #[test]
    fn spolynomial_cancels_the_lcm_term() {
        let ring = ring();
        let f = LabeledPoly {
            sig: Signature {
                index: 0,
                term: Monomial::one(2),
            },
            poly: ring.polynomial([(2, [2, 0]), (3, [0, 1])]).expect("fits"),
            index: 0,
        };
        let g = LabeledPoly {
            sig: Signature {
                index: 1,
                term: Monomial::one(2),
            },
            poly: ring.polynomial([(5, [1, 1]), (1, [0, 0])]).expect("fits"),
            index: 1,
        };

        let (s, ()) = spolynomial::<Untracked>(&[f, g], &[(), ()], 0, 1, ring.ops())
            .expect("the product fits");
        let lcm = Monomial::from_exps([2u16, 1].into_iter().collect());
        assert!(s.terms.iter().all(|t| t.mono != lcm));
    }

    #[test]
    fn f5_reduction_respects_signature_strictness() {
        let ring = PolynomialRing::prime_field(7, ["x"]).expect("7 is prime");
        let basis = vec![LabeledPoly {
            sig: Signature {
                index: 0,
                term: Monomial::one(1),
            },
            poly: ring.polynomial([(1, [1])]).expect("fits"),
            index: 0,
        }];

        // Same signature: the reduction is forbidden and the polynomial
        // stays as it is.
        let labeled = LabeledPoly {
            sig: Signature {
                index: 0,
                term: Monomial::one(1),
            },
            poly: ring.polynomial([(1, [1])]).expect("fits"),
            index: 0,
        };
        let (out, ()) = f5_reduce::<Untracked>(labeled, (), &basis, &[()], ring.ops())
            .expect("the product fits");
        assert_eq!(out.poly.terms.len(), 1);

        // Larger signature: the reduction runs and the polynomial goes to
        // zero.
        let labeled = LabeledPoly {
            sig: Signature {
                index: 1,
                term: Monomial::one(1),
            },
            poly: ring.polynomial([(1, [1])]).expect("fits"),
            index: 0,
        };
        let (out, ()) = f5_reduce::<Untracked>(labeled, (), &basis, &[()], ring.ops())
            .expect("the product fits");
        assert!(out.poly.is_zero());
    }

    #[test]
    fn solve_handles_a_monomial_ideal() {
        let ring = ring();
        let x2 = ring.polynomial([(1, [2, 0])]).expect("fits");
        let xy = ring.polynomial([(1, [1, 1])]).expect("fits");

        let gb = solve(&ring, &[x2, xy]);
        assert_eq!(gb.len(), 2);
        let leads: Vec<Vec<u16>> = gb
            .iter()
            .map(|f| f.leading_term().expect("non-zero").1.to_vec())
            .collect();
        assert!(leads.contains(&vec![2, 0]));
        assert!(leads.contains(&vec![1, 1]));
    }
}
