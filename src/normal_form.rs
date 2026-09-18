//! Division by a Gröbner basis, and the check that a list is one.
//!
//! One algorithm serves every caller: [`GroebnerBasis::normal_form`],
//! [`GroebnerBasis::contains`], and the checked constructor
//! `GroebnerBasis::from_polynomials`. It is generic over the domain, so
//! `F_p` and `Q` run the same code. The verifiers keep their own division
//! and call nothing here, which `tests/isolation.rs` holds.
//!
//! [`GroebnerBasis::normal_form`]: crate::GroebnerBasis::normal_form
//! [`GroebnerBasis::contains`]: crate::GroebnerBasis::contains

use std::fmt;
use std::mem::size_of;

use crate::compute::{Budget, ComputeError, ComputeLimits, DEGREE_LIMIT, RunError};
use crate::ideal::GroebnerBasis;
use crate::poly::{ExponentOverflow, Polynomial, Term, divide_coefficients, heap_exps_bytes};
use crate::ring::{Domain, DomainOps, PolynomialRing, PrimeField};

/// Why division does not return a result.
///
/// [`NormalFormError::Timeout`] and
/// [`NormalFormError::MemoryLimitExceeded`] report an exhausted budget.
/// [`NormalFormError::ExponentLimit`] reports a monomial the division
/// needs and one exponent cannot hold. None of them says anything about
/// the ideal.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NormalFormError {
    /// The polynomial and the basis belong to different rings.
    RingMismatch,
    /// A monomial the division needs holds an exponent above `limit`.
    ///
    /// Dividing `x^65535*y` by `x^65535 + y^65535` reaches it: the
    /// quotient monomial is `y`, and the tail multiple is `y^65536`.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// The deadline passed.
    Timeout,
    /// The live data of the division passed the memory limit.
    MemoryLimitExceeded,
}

impl fmt::Display for NormalFormError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NormalFormError::RingMismatch => {
                f.write_str("the polynomial and the basis belong to different rings")
            }
            NormalFormError::ExponentLimit { limit } => write!(
                f,
                "the division reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            NormalFormError::Timeout => f.write_str("the division passed its deadline"),
            NormalFormError::MemoryLimitExceeded => {
                f.write_str("the division passed its memory limit")
            }
        }
    }
}

impl std::error::Error for NormalFormError {}

/// The quotients and remainder of division by a Gröbner basis.
///
/// If the basis is `G = [g_0, ..., g_k]`, the result holds one quotient
/// `q_i` for every element of `G` and a remainder `r` such that
/// `f = sum_i q_i * g_i + r`. The identity is against the supplied basis
/// `G`, not an original input list `F`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DivisionResult<D: Domain = PrimeField> {
    /// Quotients in the order of the basis elements.
    pub quotients: Vec<Polynomial<D>>,
    /// The remainder, with no monomial divisible by a basis leading
    /// monomial.
    pub remainder: Polynomial<D>,
}

impl<D: Domain> DivisionResult<D> {
    /// Borrow the quotients in basis order.
    pub fn quotients(&self) -> &[Polynomial<D>] {
        &self.quotients
    }

    /// Borrow the remainder.
    pub fn remainder(&self) -> &Polynomial<D> {
        &self.remainder
    }

    /// Take the quotients and remainder.
    pub fn into_parts(self) -> (Vec<Polynomial<D>>, Polynomial<D>) {
        (self.quotients, self.remainder)
    }
}

impl<D: Domain> GroebnerBasis<D> {
    /// Divide `f` by the basis and return every quotient and the remainder.
    ///
    /// If the basis is `G = [g_0, ..., g_k]`, the result holds one quotient
    /// per basis element and satisfies `f = sum_i q_i * g_i + r`. The
    /// identity is against this basis `G`, not an original input list `F`.
    /// The basis is a Gröbner basis, so `r` is the unique normal form.
    ///
    /// The budget stops division between reduction steps. It also charges
    /// the input, working polynomial, quotients, remainder, and temporary
    /// reduction data before each allocation. A polynomial from another
    /// ring returns [`NormalFormError::RingMismatch`].
    pub fn divide(
        &self,
        f: &Polynomial<D>,
        budget: Budget,
    ) -> Result<DivisionResult<D>, NormalFormError> {
        if f.ring() != self.ring() {
            return Err(NormalFormError::RingMismatch);
        }
        let limits = ComputeLimits::of_budget(&budget);
        divide_with_limits(self, f, &limits).map_err(NormalFormError::from)
    }
}

/// Why a list of polynomials is not accepted as a reduced Gröbner basis.
///
/// `GroebnerBasis::from_polynomials` reports it. The first six variants
/// name a property the list does not have, with the index or the index
/// pair that shows it. The last three report a limit of the check itself
/// and say nothing about the list.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisError {
    /// A polynomial belongs to a different ring than the basis.
    RingMismatch,
    /// The polynomial at `index` is zero.
    ZeroPolynomial {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The polynomial at `index` has a leading coefficient other than 1.
    NotMonic {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The leading monomial at `index` is not smaller than the one before
    /// it.
    ///
    /// A reduced basis runs strictly descending, so two equal leading
    /// monomials are reported here too.
    NotSorted {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// A monomial of the polynomial at `index` is divisible by the
    /// leading monomial of another element.
    NotInterreduced {
        /// The position in the list the caller gave.
        index: usize,
    },
    /// The S-polynomial of the pair does not reduce to zero, so the list
    /// is not a Gröbner basis.
    ///
    /// The pair is the first one in ascending order over the caller's own
    /// indices, so the report is a function of the list and not of the
    /// order the check walks.
    NotGroebner {
        /// The smaller position of the pair.
        left: usize,
        /// The larger position of the pair.
        right: usize,
    },
    /// A monomial the check needs holds an exponent above `limit`.
    ///
    /// Building an S-polynomial multiplies the tail of each side by a
    /// quotient monomial, and that product can pass the width even when
    /// both leading monomials and their least common multiple fit.
    ExponentLimit {
        /// The largest value one exponent holds.
        limit: u32,
    },
    /// The deadline passed.
    Timeout,
    /// The live data of the check passed the memory limit.
    MemoryLimitExceeded,
}

impl fmt::Display for BasisError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            BasisError::RingMismatch => {
                f.write_str("a polynomial belongs to a different ring than the basis")
            }
            BasisError::ZeroPolynomial { index } => {
                write!(f, "the polynomial at index {index} is zero")
            }
            BasisError::NotMonic { index } => write!(
                f,
                "the polynomial at index {index} has a leading coefficient other than 1"
            ),
            BasisError::NotSorted { index } => write!(
                f,
                "the leading monomial at index {index} is not smaller than the one before it"
            ),
            BasisError::NotInterreduced { index } => write!(
                f,
                "a monomial of the polynomial at index {index} is divisible by the leading monomial of another element"
            ),
            BasisError::NotGroebner { left, right } => write!(
                f,
                "the S-polynomial of the elements at indices {left} and {right} does not reduce to zero"
            ),
            BasisError::ExponentLimit { limit } => write!(
                f,
                "the check reached an exponent above {limit}, which is the largest one exponent holds"
            ),
            BasisError::Timeout => f.write_str("the check passed its deadline"),
            BasisError::MemoryLimitExceeded => f.write_str("the check passed its memory limit"),
        }
    }
}

impl std::error::Error for BasisError {}

/// Why a division stops before it has a remainder.
///
/// The public errors name the same cases. The division returns this one
/// and each caller maps it: [`NormalFormError`] for a normal form,
/// [`BasisError`] for a basis check.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum DivisionStop {
    /// The budget ran out, or the cancellation flag was set.
    Run(RunError),
    /// A monomial the division needs holds an exponent past the width.
    ExponentLimit,
}

impl From<ExponentOverflow> for DivisionStop {
    fn from(_: ExponentOverflow) -> Self {
        DivisionStop::ExponentLimit
    }
}

impl From<DivisionStop> for NormalFormError {
    fn from(stop: DivisionStop) -> Self {
        match stop {
            // A cancelled division reports a timeout, as
            // `RunError::reported` does: it stopped before it had a
            // remainder.
            DivisionStop::Run(RunError::Compute(ComputeError::MemoryLimitExceeded)) => {
                NormalFormError::MemoryLimitExceeded
            }
            DivisionStop::Run(_) => NormalFormError::Timeout,
            DivisionStop::ExponentLimit => NormalFormError::ExponentLimit {
                limit: DEGREE_LIMIT,
            },
        }
    }
}

impl From<DivisionStop> for BasisError {
    fn from(stop: DivisionStop) -> Self {
        match NormalFormError::from(stop) {
            NormalFormError::MemoryLimitExceeded => BasisError::MemoryLimitExceeded,
            NormalFormError::ExponentLimit { limit } => BasisError::ExponentLimit { limit },
            _ => BasisError::Timeout,
        }
    }
}

impl From<ExponentOverflow> for BasisError {
    fn from(_: ExponentOverflow) -> Self {
        BasisError::ExponentLimit {
            limit: DEGREE_LIMIT,
        }
    }
}

/// Reduce `f` modulo `basis` and return the remainder.
///
/// The algorithm is full division. While the largest monomial of the
/// working polynomial is divisible by the leading monomial of some
/// element, the matching multiple of that element is subtracted. A
/// monomial no element divides moves to the remainder. `basis` must be a
/// Gröbner basis, and then the remainder is the unique normal form
/// whatever order the divisors are picked in. This takes the first
/// divisor in the order of `basis`, so the steps are deterministic as
/// well.
///
/// The limits stop the division between two reduction steps. One step
/// over `Q` is one exact subtraction, so that step is the granularity of
/// the deadline. The memory limit is read twice per step: once on the
/// working value and the remainder, and once more on the replacement the
/// step is about to allocate next to them.
pub(crate) fn normal_form<D: Domain>(
    basis: &[Polynomial<D>],
    f: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<Polynomial<D>, DivisionStop> {
    normal_form_with_extra(basis, f, limits, 0)
}

struct NormalFormState<'a, D: Domain> {
    basis: &'a [Polynomial<D>],
    ring: &'a PolynomialRing<D>,
    limits: &'a ComputeLimits,
    extra_bytes: usize,
    input_bytes: usize,
    ops: &'a D::Ops,
}

fn normal_form_with_extra<D: Domain>(
    basis: &[Polynomial<D>],
    f: &Polynomial<D>,
    limits: &ComputeLimits,
    extra_bytes: usize,
) -> Result<Polynomial<D>, DivisionStop> {
    let ring = f.ring();
    let state = NormalFormState {
        basis,
        ring,
        limits,
        extra_bytes,
        input_bytes: f.retained_bytes(),
        ops: ring.ops(),
    };
    let mut working = normal_form_working(&state, f)?;
    // The remainder grows by the largest term left, so it is built
    // descending and turned around once.
    let mut remainder: Vec<Term<D>> = Vec::new();
    while working.lt().is_some() {
        normal_form_step(&state, &mut working, &mut remainder)?;
    }
    finish_normal_form(&state, working, remainder)
}

fn normal_form_working<D: Domain>(
    state: &NormalFormState<'_, D>,
    f: &Polynomial<D>,
) -> Result<Polynomial<D>, DivisionStop> {
    check(
        state.limits,
        state.extra_bytes.saturating_add(state.input_bytes),
    )?;
    check(
        state.limits,
        state
            .extra_bytes
            .saturating_add(state.input_bytes.saturating_add(state.input_bytes)),
    )?;
    let working = f.clone();
    check(
        state.limits,
        state
            .extra_bytes
            .saturating_add(state.input_bytes)
            .saturating_add(working.retained_bytes()),
    )?;
    Ok(working)
}

fn normal_form_held<D: Domain>(
    state: &NormalFormState<'_, D>,
    working: &Polynomial<D>,
    remainder: &Vec<Term<D>>,
) -> usize {
    state
        .extra_bytes
        .saturating_add(state.input_bytes)
        .saturating_add(working.retained_bytes())
        .saturating_add(retained_terms_bytes(
            remainder,
            remainder.capacity(),
            working.ring().nvars(),
            state.ops,
        ))
}

fn normal_form_step<D: Domain>(
    state: &NormalFormState<'_, D>,
    working: &mut Polynomial<D>,
    remainder: &mut Vec<Term<D>>,
) -> Result<(), DivisionStop> {
    let held = normal_form_held(state, working, remainder);
    check(state.limits, held)?;
    let lead = working
        .lt()
        .expect("the loop condition found a leading term");
    match find_reducer(state.basis, lead, state.limits)? {
        Some(reducer) => {
            let lead_bytes = term_bytes(lead, state.ops, working.ring().nvars());
            check(state.limits, held.saturating_add(lead_bytes))?;
            let lead = lead.clone();
            reduce_normal_form(
                state,
                working,
                lead,
                reducer,
                held.saturating_add(lead_bytes),
            )
        }
        None => move_to_remainder(working, remainder, held, state.limits),
    }
}

fn reduce_normal_form<D: Domain>(
    state: &NormalFormState<'_, D>,
    working: &mut Polynomial<D>,
    lead: Term<D>,
    reducer: Reducer<'_, D>,
    held: usize,
) -> Result<(), DivisionStop> {
    let ring = working.ring();
    let lead_bytes = term_bytes(&lead, state.ops, ring.nvars());
    let multiple_bytes = monomial_bytes(ring.nvars());
    check(state.limits, held.saturating_add(multiple_bytes))?;
    let multiple = lead
        .mono
        .quotient(&reducer.leading.mono)
        // divides() implies a quotient exists.
        .expect("divides() implies quotient()");
    let scale_bound =
        coefficient_division_bound::<D>(state.ops, &lead.coeff, &reducer.leading.coeff);
    check(
        state.limits,
        held.saturating_add(lead_bytes)
            .saturating_add(multiple_bytes)
            .saturating_add(scale_bound),
    )?;
    let scale = divide_coefficients::<D>(state.ops, &lead.coeff, &reducer.leading.coeff);
    let scale_bytes = state.ops.heap_bytes(&scale);
    check(
        state.limits,
        held.saturating_add(lead_bytes)
            .saturating_add(multiple_bytes)
            .saturating_add(scale_bytes),
    )?;
    // The replacement is built while the working value is still live.
    check(
        state.limits,
        held.saturating_add(lead_bytes)
            .saturating_add(multiple_bytes)
            .saturating_add(scale_bytes)
            .saturating_add(working.sub_scaled_bytes(reducer.polynomial, &scale, state.ops)),
    )?;
    *working = working.sub_scaled_checked(reducer.polynomial, &scale, &multiple, state.ops)?;
    Ok(())
}

fn finish_normal_form<D: Domain>(
    state: &NormalFormState<'_, D>,
    working: Polynomial<D>,
    mut remainder: Vec<Term<D>>,
) -> Result<Polynomial<D>, DivisionStop> {
    check(state.limits, normal_form_held(state, &working, &remainder))?;
    drop(working);
    let held = state
        .extra_bytes
        .saturating_add(state.input_bytes)
        .saturating_add(retained_terms_bytes(
            &remainder,
            remainder.capacity(),
            state.ring.nvars(),
            state.ops,
        ));
    check(state.limits, held)?;
    remainder.reverse();
    check(state.limits, held)?;
    Ok(Polynomial::from_sorted_terms(state.ring.clone(), remainder))
}

/// Divide `f` by `basis` under absolute limits and retain the quotients.
///
/// The helper returns the same identity as [`GroebnerBasis::divide`]. Its
/// absolute limits let one caller cover a sequence of divisions with one
/// deadline and cancellation flag.
pub(crate) fn divide_with_limits<D: Domain>(
    basis: &[Polynomial<D>],
    f: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<DivisionResult<D>, DivisionStop> {
    let ring = f.ring();
    let mut state = DivideState {
        basis,
        ring,
        limits,
        ops: ring.ops(),
        input_bytes: f.retained_bytes(),
        quotient_slots: size_of::<Polynomial<D>>().saturating_mul(basis.len()),
    };
    let mut parts = initialize_division(&mut state, f)?;
    divide_steps(
        &state,
        &mut parts.working,
        &mut parts.quotients,
        &mut parts.remainder,
    )?;
    finish_division(&state, parts.quotients, parts.working, &mut parts.remainder)
}

struct DivideState<'a, D: Domain> {
    basis: &'a [Polynomial<D>],
    ring: &'a PolynomialRing<D>,
    limits: &'a ComputeLimits,
    ops: &'a D::Ops,
    input_bytes: usize,
    quotient_slots: usize,
}

struct DivisionParts<D: Domain> {
    quotients: Vec<Polynomial<D>>,
    working: Polynomial<D>,
    remainder: Vec<Term<D>>,
}

fn initialize_division<D: Domain>(
    state: &mut DivideState<'_, D>,
    f: &Polynomial<D>,
) -> Result<DivisionParts<D>, DivisionStop> {
    check(state.limits, state.input_bytes)?;
    check(
        state.limits,
        state.input_bytes.saturating_add(state.quotient_slots),
    )?;
    let quotients = empty_quotients(state.basis.len(), f, state.limits)?;
    let quotient_slots = size_of::<Polynomial<D>>().saturating_mul(quotients.capacity());
    state.quotient_slots = quotient_slots;
    check(
        state.limits,
        state.input_bytes.saturating_add(quotient_slots),
    )?;
    check(
        state.limits,
        state
            .input_bytes
            .saturating_add(quotient_slots)
            .saturating_add(state.input_bytes),
    )?;
    let working = f.clone();
    check(
        state.limits,
        state
            .input_bytes
            .saturating_add(quotient_slots)
            .saturating_add(working.retained_bytes()),
    )?;
    Ok(DivisionParts {
        quotients,
        working,
        remainder: Vec::new(),
    })
}

fn finish_division<D: Domain>(
    state: &DivideState<'_, D>,
    quotients: Vec<Polynomial<D>>,
    working: Polynomial<D>,
    remainder: &mut Vec<Term<D>>,
) -> Result<DivisionResult<D>, DivisionStop> {
    reverse_remainder(remainder, state.limits)?;
    check(
        state.limits,
        live_bytes(
            state.input_bytes,
            state.quotient_slots,
            &quotients,
            &working,
            remainder,
            state.ops,
        ),
    )?;
    drop(working);
    check(
        state.limits,
        state
            .input_bytes
            .saturating_add(state.quotient_slots)
            .saturating_add(quotients.iter().fold(0usize, |bytes, quotient| {
                bytes.saturating_add(quotient.retained_bytes())
            }))
            .saturating_add(retained_terms_bytes(
                remainder,
                remainder.capacity(),
                state.ring.nvars(),
                state.ops,
            )),
    )?;
    Ok(DivisionResult {
        quotients,
        remainder: Polynomial::from_sorted_terms(state.ring.clone(), std::mem::take(remainder)),
    })
}

fn divide_steps<D: Domain>(
    state: &DivideState<'_, D>,
    working: &mut Polynomial<D>,
    quotients: &mut [Polynomial<D>],
    remainder: &mut Vec<Term<D>>,
) -> Result<(), DivisionStop> {
    while working.lt().is_some() {
        divide_step(state, working, quotients, remainder)?;
    }
    Ok(())
}

fn divide_step<D: Domain>(
    state: &DivideState<'_, D>,
    working: &mut Polynomial<D>,
    quotients: &mut [Polynomial<D>],
    remainder: &mut Vec<Term<D>>,
) -> Result<(), DivisionStop> {
    let held = live_bytes(
        state.input_bytes,
        state.quotient_slots,
        quotients,
        working,
        remainder,
        state.ops,
    );
    check(state.limits, held)?;
    let lead = working
        .lt()
        .expect("the loop condition found a leading term");
    let reducer = find_reducer(state.basis, lead, state.limits)?;
    if let Some(reducer) = reducer {
        let lead_bytes = term_bytes(lead, state.ops, working.ring().nvars());
        check(state.limits, held.saturating_add(lead_bytes))?;
        let lead = lead.clone();
        reduce_with_quotient(
            state,
            working,
            quotients,
            reducer,
            &lead,
            held.saturating_add(lead_bytes),
        )?;
    } else {
        move_to_remainder(working, remainder, held, state.limits)?;
    }
    Ok(())
}

fn empty_quotients<D: Domain>(
    count: usize,
    template: &Polynomial<D>,
    limits: &ComputeLimits,
) -> Result<Vec<Polynomial<D>>, DivisionStop> {
    let mut quotients = Vec::with_capacity(count);
    for index in 0..count {
        if let Some(stop) = limits.stop_every(index) {
            return Err(DivisionStop::Run(stop));
        }
        quotients.push(template.zero_like());
    }
    Ok(quotients)
}

fn reverse_remainder<D: Domain>(
    remainder: &mut [Term<D>],
    limits: &ComputeLimits,
) -> Result<(), DivisionStop> {
    for index in 0..remainder.len() / 2 {
        if let Some(stop) = limits.stop_every(index) {
            return Err(DivisionStop::Run(stop));
        }
        let other = remainder.len() - 1 - index;
        remainder.swap(index, other);
    }
    Ok(())
}

struct Reducer<'a, D: Domain> {
    index: usize,
    polynomial: &'a Polynomial<D>,
    leading: &'a Term<D>,
}

fn find_reducer<'a, D: Domain>(
    basis: &'a [Polynomial<D>],
    lead: &Term<D>,
    limits: &ComputeLimits,
) -> Result<Option<Reducer<'a, D>>, DivisionStop> {
    for (index, g) in basis.iter().enumerate() {
        if let Some(stop) = limits.stop_every(index) {
            return Err(DivisionStop::Run(stop));
        }
        let Some(lead_g) = g.lt() else { continue };
        if lead_g.mono.divides(&lead.mono) {
            return Ok(Some(Reducer {
                index,
                polynomial: g,
                leading: lead_g,
            }));
        }
    }
    Ok(None)
}

fn reduce_with_quotient<D: Domain>(
    state: &DivideState<'_, D>,
    working: &mut Polynomial<D>,
    quotients: &mut [Polynomial<D>],
    reducer: Reducer<'_, D>,
    lead: &Term<D>,
    held: usize,
) -> Result<(), DivisionStop> {
    let ring = working.ring();
    let index = reducer.index;
    let reducer_polynomial = reducer.polynomial;
    let lead_reducer = reducer.leading;
    let multiple_bytes = monomial_bytes(ring.nvars());
    let scale_bound = coefficient_division_bound::<D>(state.ops, &lead.coeff, &lead_reducer.coeff);
    let qterm_bound = size_of::<Term<D>>()
        .saturating_add(heap_exps_bytes(ring.nvars()))
        .saturating_add(scale_bound);
    let replacement_bound = sub_scaled_bound(working, reducer_polynomial, scale_bound, state.ops);
    let quotient_capacity_delta = term_capacity_delta_for_vec(&quotients[index].terms, 1);
    check(
        state.limits,
        held.saturating_add(multiple_bytes)
            .saturating_add(scale_bound)
            .saturating_add(qterm_bound)
            .saturating_add(replacement_bound)
            .saturating_add(quotient_capacity_delta),
    )?;

    let multiple = lead
        .mono
        .quotient(&lead_reducer.mono)
        .expect("divides() implies quotient()");
    let scale = divide_coefficients::<D>(state.ops, &lead.coeff, &lead_reducer.coeff);
    let term = Term {
        coeff: scale,
        mono: multiple,
    };
    let term_bytes = term_bytes(&term, state.ops, ring.nvars());
    let scale_bytes = state.ops.heap_bytes(&term.coeff);
    let replacement_bound = sub_scaled_bound(working, reducer_polynomial, scale_bytes, state.ops);
    let quotient_update_bytes = quotient_update_bound(&quotients[index], &term, state.ops);
    check(
        state.limits,
        held.saturating_add(term_bytes)
            .saturating_add(replacement_bound)
            .saturating_add(quotient_update_bytes),
    )?;
    reserve_term_slot(
        &mut quotients[index].terms,
        held.saturating_add(term_bytes)
            .saturating_add(replacement_bound)
            .saturating_add(quotient_update_bytes),
        state.limits,
    )?;
    *working =
        working.sub_scaled_checked(reducer_polynomial, &term.coeff, &term.mono, state.ops)?;
    quotients[index].push_term(term, state.ops);
    Ok(())
}

fn move_to_remainder<D: Domain>(
    working: &mut Polynomial<D>,
    remainder: &mut Vec<Term<D>>,
    held: usize,
    limits: &ComputeLimits,
) -> Result<(), DivisionStop> {
    reserve_term_slot(remainder, held, limits)?;
    let term = working
        .pop_lt()
        .expect("the leading term is held by the working polynomial");
    remainder.push(term);
    Ok(())
}

fn live_bytes<D: Domain>(
    input_bytes: usize,
    quotient_slots: usize,
    quotients: &[Polynomial<D>],
    working: &Polynomial<D>,
    remainder: &Vec<Term<D>>,
    ops: &D::Ops,
) -> usize {
    let quotient_bytes = quotients.iter().fold(0usize, |bytes, quotient| {
        bytes.saturating_add(quotient.retained_bytes())
    });
    input_bytes
        .saturating_add(quotient_slots)
        .saturating_add(quotient_bytes)
        .saturating_add(working.retained_bytes())
        .saturating_add(retained_terms_bytes(
            remainder,
            remainder.capacity(),
            working.ring().nvars(),
            ops,
        ))
}

fn sub_scaled_bound<D: Domain>(
    working: &Polynomial<D>,
    reducer: &Polynomial<D>,
    scale_bytes: usize,
    ops: &D::Ops,
) -> usize {
    let per_term = size_of::<Term<D>>() + heap_exps_bytes(working.ring().nvars());
    let terms = working
        .terms
        .len()
        .saturating_add(reducer.terms.len())
        .saturating_mul(per_term);
    let scaled = reducer.terms.iter().fold(0usize, |bytes, term| {
        bytes
            .saturating_add(ops.heap_bytes(&term.coeff))
            .saturating_add(scale_bytes)
    });
    terms
        .saturating_add(working.coefficient_bytes())
        .saturating_add(scaled)
        .saturating_add(scaled)
        .saturating_add(heap_exps_bytes(working.ring().nvars()))
}

fn quotient_update_bound<D: Domain>(
    quotient: &Polynomial<D>,
    term: &Term<D>,
    ops: &D::Ops,
) -> usize {
    let Ok(index) = quotient
        .terms
        .binary_search_by(|probe| probe.mono.cmp(&term.mono))
    else {
        return 0;
    };
    let existing = ops.heap_bytes(&quotient.terms[index].coeff);
    let incoming = ops.heap_bytes(&term.coeff);
    existing.saturating_add(incoming).saturating_mul(2)
}

/// Check that `polynomials` is the reduced Gröbner basis of the ideal it
/// generates.
///
/// The checks run in one order, so the report is a function of the list:
/// the ring of every element, then no zero element, then monic, then
/// strictly descending leading monomials, then interreduced, then every
/// S-polynomial reduces to zero. The last check is the Buchberger
/// criterion over the whole list, which is what makes the list a Gröbner
/// basis of its own ideal.
///
/// A list this accepts is a check and not a certificate. It says nothing
/// about the ideal a caller meant, only about the list.
pub(crate) fn check_basis<D: Domain>(
    ring: &PolynomialRing<D>,
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
) -> Result<(), BasisError> {
    let ops = ring.ops();
    stop(limits)?;
    let basis_bytes = basis_input_bytes(polynomials);
    check(limits, basis_bytes).map_err(BasisError::from)?;
    check_elements(ring, polynomials, limits, ops)?;
    check_order(polynomials, limits)?;
    check_interreduced(polynomials, limits)?;
    check_pairs(polynomials, limits, ops, basis_bytes)
}

fn check_elements<D: Domain>(
    ring: &PolynomialRing<D>,
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
    ops: &D::Ops,
) -> Result<(), BasisError> {
    for (index, poly) in polynomials.iter().enumerate() {
        stop(limits)?;
        if poly.ring() != ring {
            return Err(BasisError::RingMismatch);
        }
        if poly.is_zero() {
            return Err(BasisError::ZeroPolynomial { index });
        }
        let lead = poly
            .lt()
            .expect("a nonzero polynomial holds a leading term");
        if !ops.is_one(&lead.coeff) {
            return Err(BasisError::NotMonic { index });
        }
    }
    Ok(())
}

fn check_order<D: Domain>(
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
) -> Result<(), BasisError> {
    for index in 1..polynomials.len() {
        stop_at(limits, index)?;
        let previous = polynomials[index - 1]
            .lm()
            .expect("a nonzero polynomial holds a leading monomial");
        let current = polynomials[index]
            .lm()
            .expect("a nonzero polynomial holds a leading monomial");
        if current >= previous {
            return Err(BasisError::NotSorted { index });
        }
    }
    Ok(())
}

fn check_interreduced<D: Domain>(
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
) -> Result<(), BasisError> {
    for (index, poly) in polynomials.iter().enumerate() {
        stop(limits)?;
        for (other, divisor) in polynomials.iter().enumerate() {
            if other == index {
                continue;
            }
            let lead = divisor
                .lm()
                .expect("a nonzero polynomial holds a leading monomial");
            for (term_index, term) in poly.terms.iter().enumerate() {
                stop_at(limits, term_index)?;
                if lead.divides(&term.mono) {
                    return Err(BasisError::NotInterreduced { index });
                }
            }
        }
    }
    Ok(())
}

fn check_pairs<D: Domain>(
    polynomials: &[Polynomial<D>],
    limits: &ComputeLimits,
    ops: &D::Ops,
    basis_bytes: usize,
) -> Result<(), BasisError> {
    for left in 0..polynomials.len() {
        for right in (left + 1)..polynomials.len() {
            stop(limits)?;
            check(
                limits,
                basis_bytes.saturating_add(s_polynomial_bound(
                    &polynomials[left],
                    &polynomials[right],
                    ops,
                )),
            )
            .map_err(BasisError::from)?;
            let spoly = polynomials[left].s_polynomial_checked(&polynomials[right], ops)?;
            check(limits, basis_bytes.saturating_add(spoly.retained_bytes()))
                .map_err(BasisError::from)?;
            if !normal_form_with_extra(polynomials, &spoly, limits, basis_bytes)?.is_zero() {
                return Err(BasisError::NotGroebner { left, right });
            }
        }
    }
    Ok(())
}

/// Report why the division stops now, or `None` to carry on.
fn check(limits: &ComputeLimits, bytes: usize) -> Result<(), DivisionStop> {
    if let Some(stop) = limits.stop() {
        return Err(DivisionStop::Run(stop));
    }
    if let Some(limit) = limits.memory
        && bytes > limit
    {
        return Err(DivisionStop::Run(RunError::Compute(
            ComputeError::MemoryLimitExceeded,
        )));
    }
    Ok(())
}

/// Report why the check stops now, or `None` to carry on.
fn stop(limits: &ComputeLimits) -> Result<(), BasisError> {
    stop_at(limits, 0)
}

fn stop_at(limits: &ComputeLimits, index: usize) -> Result<(), BasisError> {
    match limits.stop_every(index) {
        Some(stop) => Err(BasisError::from(DivisionStop::Run(stop))),
        None => Ok(()),
    }
}

/// The bytes one term of the remainder holds.
fn term_bytes<D: Domain>(term: &Term<D>, ops: &D::Ops, nvars: usize) -> usize {
    size_of::<Term<D>>()
        .saturating_add(heap_exps_bytes(nvars))
        .saturating_add(ops.heap_bytes(&term.coeff))
}

fn basis_input_bytes<D: Domain>(polynomials: &[Polynomial<D>]) -> usize {
    let slots = size_of::<Polynomial<D>>().saturating_mul(polynomials.len());
    polynomials.iter().fold(slots, |bytes, polynomial| {
        bytes.saturating_add(polynomial.retained_bytes())
    })
}

fn scaled_polynomial_bound<D: Domain>(
    polynomial: &Polynomial<D>,
    scale_bytes: usize,
    ops: &D::Ops,
) -> (usize, usize) {
    let per_term = size_of::<Term<D>>() + heap_exps_bytes(polynomial.ring().nvars());
    let slots = per_term.saturating_mul(polynomial.terms.len());
    let coefficients = polynomial.terms.iter().fold(0usize, |bytes, term| {
        bytes
            .saturating_add(ops.heap_bytes(&term.coeff))
            .saturating_add(scale_bytes)
    });
    (slots.saturating_add(coefficients), coefficients)
}

fn s_polynomial_bound<D: Domain>(
    left: &Polynomial<D>,
    right: &Polynomial<D>,
    ops: &D::Ops,
) -> usize {
    let Some(left_lead) = left.lt() else { return 0 };
    let Some(right_lead) = right.lt() else {
        return 0;
    };
    let left_scale = ops.heap_bytes(&right_lead.coeff);
    let right_scale = ops.heap_bytes(&left_lead.coeff);
    let (left_scaled, left_coefficients) = scaled_polynomial_bound(left, left_scale, ops);
    let (right_scaled, right_coefficients) = scaled_polynomial_bound(right, right_scale, ops);
    let output_slots = size_of::<Term<D>>()
        .saturating_add(heap_exps_bytes(left.ring().nvars()))
        .saturating_mul(left.terms.len().saturating_add(right.terms.len()));
    let quotient_monomials = heap_exps_bytes(left.ring().nvars()).saturating_mul(3);
    left_scaled
        .saturating_add(right_scaled)
        .saturating_add(output_slots)
        .saturating_add(left_coefficients)
        .saturating_add(right_coefficients)
        .saturating_add(left_coefficients)
        .saturating_add(right_coefficients)
        .saturating_add(quotient_monomials)
}

fn monomial_bytes(nvars: usize) -> usize {
    heap_exps_bytes(nvars)
}

fn coefficient_division_bound<D: Domain>(
    ops: &D::Ops,
    numerator: &D::Coeff,
    denominator: &D::Coeff,
) -> usize {
    ops.heap_bytes(numerator)
        .saturating_add(ops.heap_bytes(denominator))
}

fn retained_terms_bytes<D: Domain>(
    terms: &[Term<D>],
    capacity: usize,
    nvars: usize,
    ops: &D::Ops,
) -> usize {
    let vector = size_of::<Term<D>>().saturating_mul(capacity);
    let exponents = heap_exps_bytes(nvars).saturating_mul(terms.len());
    let coefficients = terms.iter().fold(0usize, |bytes, term| {
        bytes.saturating_add(ops.heap_bytes(&term.coeff))
    });
    vector
        .saturating_add(exponents)
        .saturating_add(coefficients)
}

fn reserve_term_slot<D: Domain>(
    terms: &mut Vec<Term<D>>,
    held: usize,
    limits: &ComputeLimits,
) -> Result<(), DivisionStop> {
    if terms.len() == terms.capacity() {
        let delta = term_capacity_delta_for_vec::<D>(terms, 1);
        check(limits, held.saturating_add(delta))?;
        terms.reserve_exact(1);
        let actual =
            size_of::<Term<D>>().saturating_mul(terms.capacity().saturating_sub(terms.len()));
        check(limits, held.saturating_add(actual))?;
    }
    Ok(())
}

fn term_capacity_delta_for_vec<D: Domain>(terms: &Vec<Term<D>>, additional: usize) -> usize {
    if additional == 0 || terms.len() < terms.capacity() {
        return 0;
    }
    let needed = terms.len().saturating_add(additional);
    let current = terms.capacity();
    let next = needed.max(current.saturating_mul(2));
    size_of::<Term<D>>().saturating_mul(next.saturating_sub(current))
}
