//! Classic F5: one critical pair at a time, in signature order.

use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;
use std::time::Instant;

use super::interreduce::{reduced_groebner_basis_checked, reduced_groebner_basis_tracked};
use super::signature::{LabeledPoly, Signature, add_syzygy_rule, is_sig_redundant, is_syzygy};
use super::{ComputeError, REDUCE_DEADLINE_STRIDE, poll_deadline};
use crate::cert::origin::{self, Origin};
use crate::poly::{ExponentOverflow, Monomial, Polynomial, Term, key_divides};
use crate::ring::PolynomialRing;
use crate::ring::field::Felt;

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
    /// The bytes the held cofactors occupy, over a ring of `nvars`
    /// variables.
    fn bytes(origins: &[Self::Origin], nvars: usize) -> usize;
    fn spair(
        origins: &[Self::Origin],
        left: Multiple,
        right: Multiple,
        p: u64,
    ) -> Result<Self::Origin, ExponentOverflow>;
    fn sub_scaled(
        target: &mut Self::Origin,
        sources: &[Self::Origin],
        source: usize,
        coeff: Felt,
        mono: &Monomial,
        p: u64,
    ) -> Result<(), ExponentOverflow>;
    fn monic(target: &mut Self::Origin, poly: &Polynomial, p: u64);
    fn holds(origin: &Self::Origin, witness: &Self::Witness, poly: &Polynomial, p: u64) -> bool;
}

struct Untracked;

impl Track for Untracked {
    type Origin = ();
    type Witness = ();

    fn witness(_generators: &[Polynomial]) {}
    fn unit(_ring: &PolynomialRing, _index: usize, _count: usize) {}
    fn bytes(_origins: &[()], _nvars: usize) -> usize {
        0
    }
    fn spair(
        _origins: &[()],
        _left: Multiple,
        _right: Multiple,
        _p: u64,
    ) -> Result<(), ExponentOverflow> {
        Ok(())
    }
    fn sub_scaled(
        _target: &mut (),
        _sources: &[()],
        _source: usize,
        _coeff: Felt,
        _mono: &Monomial,
        _p: u64,
    ) -> Result<(), ExponentOverflow> {
        Ok(())
    }
    fn monic(_target: &mut (), _poly: &Polynomial, _p: u64) {}
    fn holds(_origin: &(), _witness: &(), _poly: &Polynomial, _p: u64) -> bool {
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

    fn bytes(origins: &[Origin], nvars: usize) -> usize {
        let per_term = super::per_term_bytes(nvars);
        let mut bytes = std::mem::size_of_val(origins);
        for entry in origins {
            bytes = bytes.saturating_add(std::mem::size_of_val(entry.as_slice()));
            for cofactor in entry {
                bytes = bytes.saturating_add(cofactor.terms.len() * per_term);
            }
        }
        bytes
    }

    fn spair(
        origins: &[Origin],
        left: Multiple,
        right: Multiple,
        p: u64,
    ) -> Result<Origin, ExponentOverflow> {
        origin::combine(
            &origins[left.index],
            left.coeff,
            left.mono,
            &origins[right.index],
            right.coeff,
            right.mono,
            p,
        )
    }

    fn sub_scaled(
        target: &mut Origin,
        sources: &[Origin],
        source: usize,
        coeff: Felt,
        mono: &Monomial,
        p: u64,
    ) -> Result<(), ExponentOverflow> {
        origin::sub_scaled(target, &sources[source], coeff, mono, p)
    }

    fn monic(target: &mut Origin, poly: &Polynomial, p: u64) {
        if let Some(lc) = poly.lc() {
            origin::scale(target, lc.inv(p), p);
        }
    }

    fn holds(origin: &Origin, witness: &Vec<Polynomial>, poly: &Polynomial, p: u64) -> bool {
        origin::holds(origin, witness, poly, p)
    }
}

fn f5_reduce<T: Track>(
    mut p: LabeledPoly,
    mut cofactors: T::Origin,
    basis: &[LabeledPoly],
    basis_origins: &[T::Origin],
    modulus: u64,
    deadline: Option<Instant>,
) -> Result<(LabeledPoly, T::Origin), ComputeError> {
    // The remainder grows by the largest term left, so it is built
    // descending and turned around once.
    let mut remainder: Vec<Term> = Vec::new();
    let mut since_check = 0usize;
    let leads: Vec<(u32, u64)> = basis
        .iter()
        .map(|g| {
            g.poly
                .lm()
                .map_or((u32::MAX, u64::MAX), |lm| (lm.deg, lm.divisor_key()))
        })
        .collect();

    while let Some(lt_p) = p.poly.lt().cloned() {
        since_check += p.poly.terms.len();
        if since_check >= REDUCE_DEADLINE_STRIDE {
            since_check = 0;
            poll_deadline(deadline)?;
        }
        let mut reduced = false;
        let lead_key = lt_p.mono.divisor_key();

        for (index, g) in basis.iter().enumerate() {
            let (deg, key) = leads[index];
            if deg > lt_p.mono.deg || !key_divides(key, lead_key) {
                continue;
            }
            let Some(lt_g) = g.poly.lt() else {
                continue;
            };
            if !lt_g.mono.divides(&lt_p.mono) {
                continue;
            }
            if !g.sig.shifted_is_below(&lt_p.mono, &lt_g.mono, &p.sig) {
                continue;
            }
            let m = lt_p
                .mono
                .quotient(&lt_g.mono)
                // divides() implies a quotient exists.
                .expect("divides() implies quotient()");

            let scale = lt_p.coeff.div(lt_g.coeff, modulus);
            p.poly = p.poly.sub_scaled(&g.poly, scale, &m, modulus)?;
            T::sub_scaled(&mut cofactors, basis_origins, index, scale, &m, modulus)?;
            reduced = true;
            break;
        }

        if !reduced {
            let term = p
                .poly
                .pop_lt()
                // lt_p was Some, so the polynomial is non-empty here.
                .expect("polynomial should not be empty");
            remainder.push(term);
        }
    }

    remainder.reverse();
    p.poly = Polynomial::from_sorted_terms(p.poly.ring().clone(), remainder);
    Ok((p, cofactors))
}

fn spolynomial<T: Track>(
    basis: &[LabeledPoly],
    origins: &[T::Origin],
    i: usize,
    j: usize,
    p: u64,
) -> Result<(Polynomial, T::Origin), ComputeError> {
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

    let f_scaled = f.poly.scale_monomial(lt_g.coeff, &m_f, p)?;
    let g_scaled = g.poly.scale_monomial(lt_f.coeff, &m_g, p)?;
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
        p,
    )?;
    Ok((f_scaled.sub(&g_scaled, p), cofactors))
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
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
    basis: &[LabeledPoly],
    origins: &[T::Origin],
    syzygy_rules: &[Signature],
    queue_len: usize,
    nvars: usize,
) -> Result<(), ComputeError> {
    poll_deadline(deadline)?;

    if let Some(limit) = max_memory_bytes {
        let per_term = super::per_term_bytes(nvars);
        let mut bytes = 0usize;
        bytes = bytes.saturating_add(std::mem::size_of_val(basis));
        bytes = bytes.saturating_add(
            basis
                .iter()
                .map(|lp| lp.poly.terms.len() * per_term)
                .sum::<usize>(),
        );
        bytes = bytes.saturating_add(std::mem::size_of_val(syzygy_rules));
        bytes = bytes.saturating_add(queue_len * size_of::<CriticalPair>());
        bytes = bytes.saturating_add(T::bytes(origins, nvars));
        if bytes > limit {
            return Err(ComputeError::MemoryLimitExceeded);
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
) -> Result<(), ComputeError> {
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
        return Err(ComputeError::DegreeLimit {
            limit: super::DEGREE_LIMIT,
        });
    }

    // lcm is a multiple of lm_f.
    let m_f = lcm.quotient(lm_f).expect("lm_f divides lcm");
    // lcm is a multiple of lm_g.
    let m_g = lcm.quotient(lm_g).expect("lm_g divides lcm");

    // A signature's degree is unrelated to the lcm's, so the lcm gate above
    // does not bound this product and it reports the limit itself.
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

fn solve_raw_with_limits<T: Track>(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<(Vec<Polynomial>, Vec<T::Origin>), ComputeError> {
    super::check_input_degrees(generators)?;
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    let count = generators.len();
    let witness = T::witness(generators);

    let mut basis: Vec<LabeledPoly> = vec![];
    let mut origins: Vec<T::Origin> = vec![];
    let mut syzygy_rules: Vec<Signature> = vec![];
    let mut queue: BinaryHeap<std::cmp::Reverse<CriticalPair>> = BinaryHeap::new();

    for (k, generator) in generators.iter().enumerate() {
        check_limits::<T>(
            deadline,
            max_memory_bytes,
            &basis,
            &origins,
            &syzygy_rules,
            queue.len(),
            nvars,
        )?;
        let sig = Signature {
            index: k,
            term: Monomial::one(nvars),
        };
        let mut cofactors = T::unit(ring, k, count);
        T::monic(&mut cofactors, generator, modulus);
        let labeled = LabeledPoly {
            sig,
            poly: generator.make_monic(modulus),
            index: basis.len(),
        };
        let (reduced, cofactors) =
            f5_reduce::<T>(labeled, cofactors, &basis, &origins, modulus, deadline)?;
        if reduced.poly.is_zero() {
            add_syzygy_rule(&mut syzygy_rules, reduced.sig);
            continue;
        }

        let mut reduced = reduced;
        let mut cofactors = cofactors;
        T::monic(&mut cofactors, &reduced.poly, modulus);
        reduced.poly = reduced.poly.make_monic(modulus);
        debug_assert!(
            T::holds(&cofactors, &witness, &reduced.poly, modulus),
            "a basis element must equal the combination its cofactors name"
        );
        reduced.index = basis.len();
        let new_idx = reduced.index;
        basis.push(reduced);
        origins.push(cofactors);

        for i in 0..new_idx {
            push_pair(&mut queue, &syzygy_rules, &basis, i, new_idx)?;
        }

        while let Some(std::cmp::Reverse(pair)) = queue.pop() {
            check_limits::<T>(
                deadline,
                max_memory_bytes,
                &basis,
                &origins,
                &syzygy_rules,
                queue.len(),
                nvars,
            )?;
            if is_syzygy(&pair.sig, &syzygy_rules) {
                continue;
            }

            let (s, cofactors) = spolynomial::<T>(&basis, &origins, pair.i, pair.j, modulus)?;
            if s.is_zero() {
                add_syzygy_rule(&mut syzygy_rules, pair.sig);
                continue;
            }

            let labeled = LabeledPoly {
                sig: pair.sig,
                poly: s,
                index: basis.len(),
            };
            let (reduced, cofactors) =
                f5_reduce::<T>(labeled, cofactors, &basis, &origins, modulus, deadline)?;
            if reduced.poly.is_zero() {
                add_syzygy_rule(&mut syzygy_rules, reduced.sig);
                continue;
            }
            if is_sig_redundant(&reduced, &basis) {
                continue;
            }

            let mut reduced = reduced;
            let mut cofactors = cofactors;
            T::monic(&mut cofactors, &reduced.poly, modulus);
            reduced.poly = reduced.poly.make_monic(modulus);
            debug_assert!(
                T::holds(&cofactors, &witness, &reduced.poly, modulus),
                "a basis element must equal the combination its cofactors name"
            );
            reduced.index = basis.len();
            let new_idx = reduced.index;
            basis.push(reduced);
            origins.push(cofactors);

            for i in 0..new_idx {
                push_pair(&mut queue, &syzygy_rules, &basis, i, new_idx)?;
            }
        }
    }

    // An input that never enters the pair loop still spends time and
    // memory. Without this check a run over a single generator, or over no
    // generator at all, reports success past an exhausted budget.
    check_limits::<T>(
        deadline,
        max_memory_bytes,
        &basis,
        &origins,
        &syzygy_rules,
        queue.len(),
        nvars,
    )?;

    Ok((basis.into_iter().map(|lp| lp.poly).collect(), origins))
}

/// Run the classic backend under a deadline and a memory cap.
pub(super) fn solve_checked(
    ring: &PolynomialRing,
    generators: &[Polynomial],
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<Vec<Polynomial>, ComputeError> {
    let (raw, _) =
        solve_raw_with_limits::<Untracked>(ring, generators, deadline, max_memory_bytes)?;
    reduced_groebner_basis_checked(ring, raw, deadline, max_memory_bytes)
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
    deadline: Option<Instant>,
    max_memory_bytes: Option<usize>,
) -> Result<(Vec<Polynomial>, Vec<Origin>), ComputeError> {
    let (raw, origins) =
        solve_raw_with_limits::<Tracked>(ring, generators, deadline, max_memory_bytes)?;
    reduced_groebner_basis_tracked(ring, raw, origins, generators, deadline, max_memory_bytes)
}

/// Run the classic backend with no budget.
#[cfg(test)]
pub(crate) fn solve(ring: &PolynomialRing, generators: &[Polynomial]) -> Vec<Polynomial> {
    // unlimited runs should only fail on internal invariants.
    solve_checked(ring, generators, None, None).expect("a run without a budget cannot stop early")
}

#[cfg(test)]
mod tests {
    use super::{Untracked, f5_reduce, solve, solve_checked, spolynomial};
    use crate::compute::ComputeError;
    use crate::compute::signature::{LabeledPoly, Signature};
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
            solve_checked(&ring, &[f, g], None, None),
            Err(ComputeError::DegreeLimit { limit: 65535 })
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

        let (s, ()) =
            spolynomial::<Untracked>(&[f, g], &[(), ()], 0, 1, 7).expect("the product fits");
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
        let (out, ()) =
            f5_reduce::<Untracked>(labeled, (), &basis, &[()], 7, None).expect("the product fits");
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
        let (out, ()) =
            f5_reduce::<Untracked>(labeled, (), &basis, &[()], 7, None).expect("the product fits");
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
