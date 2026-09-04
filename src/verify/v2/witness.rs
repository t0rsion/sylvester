//! The membership section and the pair section.
//!
//! A division trace reduces a polynomial to zero over G. The membership
//! section holds one trace per input polynomial. The pair section holds one
//! witness per unordered pair of G, in lexicographic pair order.
//!
//! The verifier checks the structural rules of every witness first, then
//! checks the dependency graph for cycles, then validates the witnesses in
//! reverse topological order. That order is part of the contract: it stops
//! a cycle of `Chain` witnesses from justifying itself.

use super::Ctx;
use super::arith::{self, Mono, Poly};
use super::decode::{Pool, Reader, bounded, capped};
use super::replay::{lead_multiple, pair_count};
use crate::verify::error::{
    BinaryFault, Cap, DivisionFault, DivisionSite, VerifyError, WitnessFault,
};

const KIND_COPRIME: u64 = 0;
const KIND_CHAIN: u64 = 1;
const KIND_REDUCE: u64 = 2;

/// One step of a division trace: a monomial and a basis index.
type Step = (Mono, usize);

/// A witness of one pair.
enum Witness {
    Coprime,
    Chain { k: usize },
    Reduce { steps: Vec<Step> },
}

/// The running total of division steps in the whole certificate.
pub(super) struct Steps {
    total: usize,
}

impl Steps {
    pub(super) fn new() -> Self {
        Steps { total: 0 }
    }
}

/// Decode one division trace.
fn division_trace(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    basis_len: usize,
    steps: &mut Steps,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<Step>, VerifyError> {
    let count = division_count(reader, steps, ctx)?;
    let mut trace: Vec<Step> = Vec::with_capacity(count);
    for _ in 0..count {
        trace.push(division_step(reader, pool, basis_len, ctx)?);
    }
    ctx.meter.hold(trace.len())?;
    Ok(trace)
}

fn division_count(
    reader: &mut Reader<'_>,
    steps: &mut Steps,
    ctx: &mut Ctx<'_>,
) -> Result<usize, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(count, ctx.limits.max_division_steps, Cap::DivisionSteps)?;
    steps.total += count;
    if steps.total > ctx.limits.max_total_division_steps {
        return Err(VerifyError::CapExceeded {
            cap: Cap::TotalDivisionSteps,
            limit: ctx.limits.max_total_division_steps,
        });
    }
    if count > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.offset(),
        });
    }
    Ok(count)
}

fn division_step(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    basis_len: usize,
    ctx: &mut Ctx<'_>,
) -> Result<Step, VerifyError> {
    ctx.meter.charge(1)?;
    let pool_index = reader.varint(&mut ctx.meter)?;
    let mono = pool.take(pool_index, &mut ctx.meter)?;
    let element = reader.varint(&mut ctx.meter)?;
    let element = bounded(element, basis_len, "basis index")?;
    Ok((mono, element))
}

/// Replay a division trace and require a zero residual at the end.
fn replay_division(
    start: &Poly,
    steps: &[Step],
    basis: &[Poly],
    at: DivisionSite,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    ctx.meter.charge((start.term_count() as u64).max(1))?;
    ctx.meter.check_intermediate(start.term_count())?;
    let mut residual = start.clone();
    for (step, (mono, element)) in steps.iter().enumerate() {
        residual = replay_step(residual, mono, *element, basis, at, step, ctx)?;
    }
    if !residual.is_zero() {
        return Err(VerifyError::Division {
            at,
            step: steps.len(),
            fault: DivisionFault::NotZero,
        });
    }
    Ok(())
}

fn replay_step(
    residual: Poly,
    mono: &Mono,
    element: usize,
    basis: &[Poly],
    at: DivisionSite,
    step: usize,
    ctx: &mut Ctx<'_>,
) -> Result<Poly, VerifyError> {
    if residual.is_zero() {
        return Err(VerifyError::Division {
            at,
            step,
            fault: DivisionFault::ResidualZero,
        });
    }
    let reducer = &basis[element];
    let lead = lead_multiple(mono, reducer, &mut ctx.meter)?;
    let residual_lead = residual.lm().expect("the residual is not zero");
    arith::charge_monos(&mut ctx.meter, &lead, residual_lead)?;
    if lead != *residual_lead {
        return Err(VerifyError::Division {
            at,
            step,
            fault: DivisionFault::LeadMismatch,
        });
    }
    let coeff = residual.lc().expect("the residual is not zero");
    let multiple = arith::mono_mul(reducer, mono, &mut ctx.meter)?;
    let scaled = arith::scale(&multiple, coeff, ctx.modulus, &mut ctx.meter)?;
    arith::sub(&residual, &scaled, ctx.modulus, &mut ctx.meter)
}

/// Check the membership section: one division trace per input polynomial,
/// in input order, each ending at zero.
pub(super) fn membership(
    reader: &mut Reader<'_>,
    input: &[Poly],
    basis: &[Poly],
    pool: &mut Pool,
    steps: &mut Steps,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    for (index, polynomial) in input.iter().enumerate() {
        let trace = division_trace(reader, pool, basis.len(), steps, ctx)?;
        replay_division(
            polynomial,
            &trace,
            basis,
            DivisionSite::Membership { input: index },
            ctx,
        )?;
        // Membership traces are not live after replay. Pair Reduce traces
        // stay charged through validation.
        ctx.meter.release(trace.len());
    }
    reader.finish()
}

/// The pair index of `(i, j)` with `i < j`, in lexicographic pair order.
fn pair_index(i: u64, j: u64, len: u64) -> u64 {
    i * len - i * (i + 1) / 2 + (j - i - 1)
}

/// The pair `(i, j)` that pair index `t` names.
///
/// The block of pairs that start with `i` begins at `pair_index(i, i+1)`,
/// which increases with `i`, so a binary search finds `i`.
fn pair_of(t: u64, len: u64) -> (usize, usize) {
    let (mut low, mut high) = (0u64, len - 1);
    while low < high {
        let mid = low + (high - low).div_ceil(2);
        if pair_index(mid, mid + 1, len) <= t {
            low = mid;
        } else {
            high = mid - 1;
        }
    }
    let i = low;
    let j = t - pair_index(i, i + 1, len) + i + 1;
    (i as usize, j as usize)
}

/// Check the pair section by the four steps of section 6.4.
pub(super) fn pairs(
    reader: &mut Reader<'_>,
    basis: &[Poly],
    pool: &mut Pool,
    steps: &mut Steps,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let count = pair_count(basis.len());
    let len = basis.len() as u64;
    if count > reader.remaining() as u64 {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.offset(),
        });
    }
    let witnesses = decode_witnesses(
        count,
        len,
        &mut WitnessDecoder {
            reader,
            basis,
            pool,
            steps,
            ctx,
        },
    )?;
    reader.finish()?;

    let order = topological_order(&witnesses, len, ctx)?;
    validate_witnesses(&witnesses, &order, basis, len, ctx)
}

struct WitnessDecoder<'a, 'data, 'ctx> {
    reader: &'a mut Reader<'data>,
    basis: &'a [Poly],
    pool: &'a mut Pool,
    steps: &'a mut Steps,
    ctx: &'a mut Ctx<'ctx>,
}

fn decode_witnesses(
    count: u64,
    len: u64,
    decoder: &mut WitnessDecoder<'_, '_, '_>,
) -> Result<Vec<Witness>, VerifyError> {
    let mut witnesses: Vec<Witness> = Vec::with_capacity(count as usize);
    for t in 0..count {
        decoder.ctx.meter.charge(1)?;
        let (i, j) = pair_of(t, len);
        witnesses.push(decode_witness(
            decoder.reader,
            decoder.basis,
            decoder.pool,
            decoder.steps,
            i,
            j,
            decoder.ctx,
        )?);
    }
    Ok(witnesses)
}

fn decode_witness(
    reader: &mut Reader<'_>,
    basis: &[Poly],
    pool: &mut Pool,
    steps: &mut Steps,
    i: usize,
    j: usize,
    ctx: &mut Ctx<'_>,
) -> Result<Witness, VerifyError> {
    match reader.varint(&mut ctx.meter)? {
        KIND_COPRIME => Ok(Witness::Coprime),
        KIND_CHAIN => decode_chain(reader, basis, i, j, ctx),
        KIND_REDUCE => Ok(Witness::Reduce {
            steps: division_trace(reader, pool, basis.len(), steps, ctx)?,
        }),
        other => Err(VerifyError::Witness {
            i,
            j,
            fault: WitnessFault::UnknownKind { kind: other },
        }),
    }
}

fn decode_chain(
    reader: &mut Reader<'_>,
    basis: &[Poly],
    i: usize,
    j: usize,
    ctx: &mut Ctx<'_>,
) -> Result<Witness, VerifyError> {
    let k = reader.varint(&mut ctx.meter)?;
    let k = bounded(k, basis.len(), "chain index")?;
    if k == i || k == j {
        return Err(VerifyError::Witness {
            i,
            j,
            fault: WitnessFault::ChainNamesPair { k },
        });
    }
    let multiple = lcm_of(basis, i, j, ctx)?;
    let lead = basis[k].lm().expect("every element is nonzero");
    arith::charge_monos(&mut ctx.meter, lead, &multiple)?;
    if !lead.divides(&multiple) {
        return Err(VerifyError::Witness {
            i,
            j,
            fault: WitnessFault::ChainNotDividing { k },
        });
    }
    Ok(Witness::Chain { k })
}

fn validate_witnesses(
    witnesses: &[Witness],
    order: &[usize],
    basis: &[Poly],
    len: u64,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let mut validated = vec![false; witnesses.len()];
    for &t in order {
        let (i, j) = pair_of(t as u64, len);
        validate_witness(&witnesses[t], i, j, basis, &validated, len, ctx)?;
        validated[t] = true;
    }
    Ok(())
}

fn validate_witness(
    witness: &Witness,
    i: usize,
    j: usize,
    basis: &[Poly],
    validated: &[bool],
    len: u64,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    match witness {
        Witness::Coprime => validate_coprime(basis, i, j, ctx),
        Witness::Chain { k } => validate_chain(i, j, *k, validated, len),
        Witness::Reduce { steps } => {
            let polynomial = spoly(basis, i, j, ctx)?;
            replay_division(&polynomial, steps, basis, DivisionSite::Pair { i, j }, ctx)
        }
    }
}

fn validate_coprime(
    basis: &[Poly],
    i: usize,
    j: usize,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let left = basis[i].lm().expect("every element is nonzero");
    let right = basis[j].lm().expect("every element is nonzero");
    arith::charge_monos(&mut ctx.meter, left, right)?;
    if left.shares_variable(right) {
        return Err(VerifyError::Witness {
            i,
            j,
            fault: WitnessFault::NotCoprime,
        });
    }
    Ok(())
}

fn validate_chain(
    i: usize,
    j: usize,
    k: usize,
    validated: &[bool],
    len: u64,
) -> Result<(), VerifyError> {
    for dependency in dependencies(i, j, k, len) {
        if !validated[dependency] {
            return Err(VerifyError::Witness {
                i,
                j,
                fault: WitnessFault::Cycle,
            });
        }
    }
    Ok(())
}

/// The least common multiple of the leading monomials of a pair.
fn lcm_of(basis: &[Poly], i: usize, j: usize, ctx: &mut Ctx<'_>) -> Result<Mono, VerifyError> {
    let left = basis[i].lm().expect("every element is nonzero");
    let right = basis[j].lm().expect("every element is nonzero");
    arith::charge_monos(&mut ctx.meter, left, right)?;
    Ok(left.lcm(right))
}

/// The S-polynomial of a pair of monic basis elements.
fn spoly(basis: &[Poly], i: usize, j: usize, ctx: &mut Ctx<'_>) -> Result<Poly, VerifyError> {
    let multiple = lcm_of(basis, i, j, ctx)?;
    let left_lead = basis[i].lm().expect("every element is nonzero");
    let right_lead = basis[j].lm().expect("every element is nonzero");
    arith::charge_monos(&mut ctx.meter, &multiple, left_lead)?;
    let left_cofactor = multiple.divide(left_lead);
    arith::charge_monos(&mut ctx.meter, &multiple, right_lead)?;
    let right_cofactor = multiple.divide(right_lead);
    let left = arith::mono_mul(&basis[i], &left_cofactor, &mut ctx.meter)?;
    let right = arith::mono_mul(&basis[j], &right_cofactor, &mut ctx.meter)?;
    arith::sub(&left, &right, ctx.modulus, &mut ctx.meter)
}

/// The two pair indices a `Chain` witness depends on.
fn dependencies(i: usize, j: usize, k: usize, len: u64) -> [usize; 2] {
    let sorted = |a: usize, b: usize| {
        let (low, high) = if a < b { (a, b) } else { (b, a) };
        pair_index(low as u64, high as u64, len) as usize
    };
    [sorted(i, k), sorted(k, j)]
}

/// Order the pairs so that every dependency comes before its dependent.
///
/// The traversal is Kahn's algorithm on the edges from a dependency to its
/// dependent. It is iterative, because a recursive one can overflow the
/// stack on a long chain. A graph that does not empty holds a cycle.
fn topological_order(
    witnesses: &[Witness],
    len: u64,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<usize>, VerifyError> {
    let count = witnesses.len();
    ctx.meter.charge(count as u64)?;
    let (waiting, dependents) = dependency_counts(witnesses, len, ctx)?;
    let offsets = dependency_offsets(&dependents);
    let edges = dependency_edges(witnesses, len, &offsets);
    let (order, waiting) = kahn_order(waiting, &offsets, &edges, ctx)?;
    check_cycle(order, waiting, count, len)
}

fn dependency_counts(
    witnesses: &[Witness],
    len: u64,
    ctx: &mut Ctx<'_>,
) -> Result<(Vec<u32>, Vec<usize>), VerifyError> {
    let count = witnesses.len();
    let mut waiting = vec![0u32; count];
    let mut dependents = vec![0usize; count];
    for (t, witness) in witnesses.iter().enumerate() {
        if let Witness::Chain { k } = witness {
            let (i, j) = pair_of(t as u64, len);
            waiting[t] = 2;
            for dependency in dependencies(i, j, *k, len) {
                ctx.meter.charge(1)?;
                dependents[dependency] += 1;
            }
        }
    }
    Ok((waiting, dependents))
}

fn dependency_offsets(dependents: &[usize]) -> Vec<usize> {
    let count = dependents.len();
    let mut offsets = vec![0usize; count + 1];
    let mut at = 0usize;
    for (slot, edges) in offsets.iter_mut().zip(dependents) {
        *slot = at;
        at += edges;
    }
    offsets[count] = at;
    offsets
}

fn dependency_edges(witnesses: &[Witness], len: u64, offsets: &[usize]) -> Vec<usize> {
    let count = witnesses.len();
    let mut edges = vec![0usize; offsets[count]];
    let mut filled = vec![0usize; count];
    for (t, witness) in witnesses.iter().enumerate() {
        if let Witness::Chain { k } = witness {
            let (i, j) = pair_of(t as u64, len);
            for dependency in dependencies(i, j, *k, len) {
                edges[offsets[dependency] + filled[dependency]] = t;
                filled[dependency] += 1;
            }
        }
    }
    edges
}

fn kahn_order(
    mut waiting: Vec<u32>,
    offsets: &[usize],
    edges: &[usize],
    ctx: &mut Ctx<'_>,
) -> Result<(Vec<usize>, Vec<u32>), VerifyError> {
    let count = waiting.len();
    let mut order: Vec<usize> = Vec::with_capacity(count);
    let mut ready: Vec<usize> = (0..count).filter(|&t| waiting[t] == 0).collect();
    while let Some(t) = ready.pop() {
        ctx.meter.charge(1)?;
        order.push(t);
        for &dependent in &edges[offsets[t]..offsets[t + 1]] {
            waiting[dependent] -= 1;
            if waiting[dependent] == 0 {
                ready.push(dependent);
            }
        }
    }
    Ok((order, waiting))
}

fn check_cycle(
    order: Vec<usize>,
    waiting: Vec<u32>,
    count: usize,
    len: u64,
) -> Result<Vec<usize>, VerifyError> {
    if order.len() != count {
        let stuck = waiting
            .iter()
            .position(|&left| left > 0)
            .expect("a graph that does not empty holds a waiting node");
        let (i, j) = pair_of(stuck as u64, len);
        return Err(VerifyError::Witness {
            i,
            j,
            fault: WitnessFault::Cycle,
        });
    }
    Ok(order)
}
