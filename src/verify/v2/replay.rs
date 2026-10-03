//! The trace replay and the basis section.
//!
//! The verifier computes every node value in one forward pass. Each
//! reference points to a smaller node id, so the graph is acyclic by
//! construction. The verifier decrements the use count of every source it
//! reads, frees a value whose count reaches zero, and keeps the live byte
//! estimate inside `max_live_bytes`.

use super::Ctx;
use super::arith::{self, Mono, Poly};
use super::decode::{Pool, Reader, bounded, capped};
use crate::verify::error::{BasisFault, BinaryFault, Cap, NodeFault, VerifyError};
use crate::verify::limits::Meter;

const KIND_INPUT: u64 = 0;
const KIND_MUL: u64 = 1;
const KIND_COMB: u64 = 2;
const KIND_SCALE: u64 = 3;
const POLL_STEPS: usize = 1 << 20;

/// The node values and the use counts the trace declares.
pub(super) struct Nodes {
    values: Vec<Option<Poly>>,
    uses: Vec<u64>,
}

#[derive(Default)]
struct TraceState {
    steps_total: usize,
    seen_other_kind: bool,
    last_input: Option<u64>,
}

struct DecodedNode {
    value: Poly,
    sources: Vec<usize>,
    use_count: u64,
}

impl Nodes {
    fn value(&self, src: usize) -> Result<&Poly, VerifyError> {
        self.values[src].as_ref().ok_or(VerifyError::Trace {
            node: src,
            fault: NodeFault::UseCountMismatch,
        })
    }

    /// Count one reference to a node, and free the value when the last
    /// reference is gone.
    fn use_once(&mut self, src: usize, meter: &mut Meter) -> Result<(), VerifyError> {
        meter.charge(1)?;
        if self.uses[src] == 0 {
            return Err(VerifyError::Trace {
                node: src,
                fault: NodeFault::UseCountMismatch,
            });
        }
        self.uses[src] -= 1;
        if self.uses[src] == 0
            && let Some(value) = self.values[src].take()
        {
            meter.release(value.term_count());
        }
        Ok(())
    }
}

/// Decode the trace section and compute every node value.
pub(super) fn trace(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    input: &[Poly],
    ctx: &mut Ctx<'_>,
) -> Result<Nodes, VerifyError> {
    let count = trace_count(reader, ctx)?;
    let mut nodes = Nodes {
        values: Vec::with_capacity(count),
        uses: Vec::with_capacity(count),
    };
    let mut state = TraceState::default();

    for id in 0..count {
        let decoded = decode_node(
            id,
            &mut NodeDecoder {
                reader,
                pool,
                input,
                nodes: &nodes,
                state: &mut state,
                ctx,
            },
        )?;
        for src in decoded.sources {
            nodes.use_once(src, &mut ctx.meter)?;
        }
        ctx.meter.hold(decoded.value.term_count())?;
        nodes.values.push(Some(decoded.value));
        nodes.uses.push(decoded.use_count);
    }
    reader.finish()?;
    Ok(nodes)
}

fn trace_count(reader: &mut Reader<'_>, ctx: &mut Ctx<'_>) -> Result<usize, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(count, ctx.limits.max_nodes, Cap::Nodes)?;
    if count > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.offset(),
        });
    }
    Ok(count)
}

struct NodeDecoder<'a, 'data, 'ctx> {
    reader: &'a mut Reader<'data>,
    pool: &'a mut Pool,
    input: &'a [Poly],
    nodes: &'a Nodes,
    state: &'a mut TraceState,
    ctx: &'a mut Ctx<'ctx>,
}

fn decode_node(
    id: usize,
    decoder: &mut NodeDecoder<'_, '_, '_>,
) -> Result<DecodedNode, VerifyError> {
    decoder.ctx.meter.charge(1)?;
    let kind = decoder.reader.varint(&mut decoder.ctx.meter)?;
    let use_count = decoder.reader.varint(&mut decoder.ctx.meter)?;
    if use_count == 0 {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::UseCountZero,
        });
    }
    let (value, sources) = decode_node_value(id, kind, decoder)?;
    Ok(DecodedNode {
        value,
        sources,
        use_count,
    })
}

fn decode_node_value(
    id: usize,
    kind: u64,
    decoder: &mut NodeDecoder<'_, '_, '_>,
) -> Result<(Poly, Vec<usize>), VerifyError> {
    match kind {
        KIND_INPUT => input_node(
            decoder.reader,
            id,
            decoder.input,
            decoder.state,
            decoder.ctx,
        ),
        KIND_MUL => {
            decoder.state.seen_other_kind = true;
            mul_node(decoder.reader, id, decoder.pool, decoder.nodes, decoder.ctx)
        }
        KIND_SCALE => {
            decoder.state.seen_other_kind = true;
            scale_node(decoder.reader, id, decoder.nodes, decoder.ctx)
        }
        KIND_COMB => {
            decoder.state.seen_other_kind = true;
            comb_node(
                decoder.reader,
                id,
                decoder.nodes,
                decoder.state,
                decoder.ctx,
            )
        }
        other => Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::UnknownKind { kind: other },
        }),
    }
}

fn input_node(
    reader: &mut Reader<'_>,
    id: usize,
    input: &[Poly],
    state: &mut TraceState,
    ctx: &mut Ctx<'_>,
) -> Result<(Poly, Vec<usize>), VerifyError> {
    if state.seen_other_kind {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::InputOutOfOrder,
        });
    }
    let index = reader.varint(&mut ctx.meter)?;
    if state.last_input.is_some_and(|last| index <= last) {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::InputOutOfOrder,
        });
    }
    state.last_input = Some(index);
    let index = bounded(index, input.len(), "input index")?;
    let value = &input[index];
    ctx.meter.charge((value.term_count() as u64).max(1))?;
    ctx.meter.check_intermediate(value.term_count())?;
    Ok((value.clone(), Vec::new()))
}

fn mul_node(
    reader: &mut Reader<'_>,
    id: usize,
    pool: &mut Pool,
    nodes: &Nodes,
    ctx: &mut Ctx<'_>,
) -> Result<(Poly, Vec<usize>), VerifyError> {
    let src = source(reader, id, ctx)?;
    let pool_index = reader.varint(&mut ctx.meter)?;
    let mono_index = pool.reference(pool_index, &mut ctx.meter)?;
    let mono = pool.get(mono_index);
    if mono.is_identity() {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::IdentityMultiplier,
        });
    }
    let value = arith::mono_mul(nodes.value(src)?, mono, &mut ctx.meter)?;
    Ok((value, vec![src]))
}

fn scale_node(
    reader: &mut Reader<'_>,
    id: usize,
    nodes: &Nodes,
    ctx: &mut Ctx<'_>,
) -> Result<(Poly, Vec<usize>), VerifyError> {
    let src = source(reader, id, ctx)?;
    let scalar = reader.varint(&mut ctx.meter)?;
    if scalar < 2 || scalar >= ctx.modulus {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::ScalarOutOfRange { scalar },
        });
    }
    let value = arith::scale(nodes.value(src)?, scalar, ctx.modulus, &mut ctx.meter)?;
    Ok((value, vec![src]))
}

fn comb_node(
    reader: &mut Reader<'_>,
    id: usize,
    nodes: &Nodes,
    state: &mut TraceState,
    ctx: &mut Ctx<'_>,
) -> Result<(Poly, Vec<usize>), VerifyError> {
    let steps = comb_steps(reader, id, &mut state.steps_total, ctx)?;
    let mut sources = Vec::with_capacity(steps.len());
    let mut iter = steps.into_iter();
    let (src, scalar) = iter.next().expect("a combination has at least two steps");
    let mut value = arith::scale(nodes.value(src)?, scalar, ctx.modulus, &mut ctx.meter)?;
    sources.push(src);
    for (index, (src, scalar)) in iter.enumerate() {
        if index > 0 && index.is_multiple_of(POLL_STEPS) {
            ctx.meter.poll()?;
        }
        value = arith::add_scaled_owned(
            value,
            nodes.value(src)?,
            scalar,
            ctx.modulus,
            &mut ctx.meter,
        )?;
        sources.push(src);
    }
    Ok((value, sources))
}

/// Read a source id and check that it is below the node's own id.
fn source(reader: &mut Reader<'_>, id: usize, ctx: &mut Ctx<'_>) -> Result<usize, VerifyError> {
    let src = reader.varint(&mut ctx.meter)?;
    if src >= id as u64 {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::SourceNotEarlier { src },
        });
    }
    Ok(src as usize)
}

/// Decode the steps of one `Comb` node.
///
/// The source ids are delta coded and strictly increasing. The decoder adds
/// with checked arithmetic, so a delta near `u64::MAX` cannot wrap into a
/// small source and let one source appear twice.
fn comb_steps(
    reader: &mut Reader<'_>,
    id: usize,
    steps_total: &mut usize,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<(usize, u64)>, VerifyError> {
    let step_count = comb_step_count(reader, id, steps_total, ctx)?;
    let mut steps: Vec<(usize, u64)> = Vec::with_capacity(step_count);
    let mut previous: Option<u64> = None;
    for _ in 0..step_count {
        ctx.meter.charge(1)?;
        let raw = reader.varint(&mut ctx.meter)?;
        let src = comb_source(raw, previous, id, reader.offset())?;
        previous = Some(src);
        let scalar = reader.varint(&mut ctx.meter)?;
        check_comb_scalar(scalar, id, ctx.modulus)?;
        steps.push((src as usize, scalar));
    }
    Ok(steps)
}

fn comb_step_count(
    reader: &mut Reader<'_>,
    id: usize,
    steps_total: &mut usize,
    ctx: &mut Ctx<'_>,
) -> Result<usize, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    if count < 2 {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::CombTooShort { steps: count },
        });
    }
    let count = capped(count, ctx.limits.max_comb_steps, Cap::CombSteps)?;
    *steps_total += count;
    if *steps_total > ctx.limits.max_trace_steps {
        return Err(VerifyError::CapExceeded {
            cap: Cap::TraceSteps,
            limit: ctx.limits.max_trace_steps,
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

fn comb_source(
    raw: u64,
    previous: Option<u64>,
    id: usize,
    offset: usize,
) -> Result<u64, VerifyError> {
    let src = match previous {
        None => raw,
        Some(last) => last
            .checked_add(raw)
            .and_then(|sum| sum.checked_add(1))
            .ok_or(VerifyError::MalformedBinary {
                reason: BinaryFault::Overflow {
                    what: "a combination source id",
                },
                offset,
            })?,
    };
    if src >= id as u64 {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::SourceNotEarlier { src },
        });
    }
    Ok(src)
}

fn check_comb_scalar(scalar: u64, id: usize, modulus: u64) -> Result<(), VerifyError> {
    if scalar == 0 || scalar >= modulus {
        return Err(VerifyError::Trace {
            node: id,
            fault: NodeFault::ScalarOutOfRange { scalar },
        });
    }
    Ok(())
}

/// Decode the basis section and check the shape of section 5.7.
///
/// The number of pairs is computed and checked against `max_pairs` before
/// anything runs over the pairs of G.
pub(super) fn basis(
    reader: &mut Reader<'_>,
    nodes: &mut Nodes,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<Poly>, VerifyError> {
    let count = basis_count(reader, ctx)?;
    let mut basis: Vec<Poly> = Vec::with_capacity(count);
    for _ in 0..count {
        basis.push(take_basis_value(reader, nodes, ctx)?);
    }
    reader.finish()?;
    check_node_uses(nodes, ctx)?;
    check_shape(&basis, ctx)?;
    Ok(basis)
}

fn basis_count(reader: &mut Reader<'_>, ctx: &mut Ctx<'_>) -> Result<usize, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(count, ctx.limits.max_basis, Cap::BasisElements)?;
    if pair_count(count) > ctx.limits.max_pairs as u64 {
        return Err(VerifyError::CapExceeded {
            cap: Cap::Pairs,
            limit: ctx.limits.max_pairs,
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

fn take_basis_value(
    reader: &mut Reader<'_>,
    nodes: &mut Nodes,
    ctx: &mut Ctx<'_>,
) -> Result<Poly, VerifyError> {
    let id = reader.varint(&mut ctx.meter)?;
    let id = bounded(id, nodes.values.len(), "basis node id")?;
    let terms = nodes.value(id)?.term_count();
    ctx.meter.charge((terms as u64).max(1))?;
    let value = nodes.value(id)?.clone();
    nodes.use_once(id, &mut ctx.meter)?;
    ctx.meter.hold(value.term_count())?;
    Ok(value)
}

fn check_node_uses(nodes: &Nodes, ctx: &mut Ctx<'_>) -> Result<(), VerifyError> {
    for (node, uses) in nodes.uses.iter().enumerate() {
        ctx.meter.charge(1)?;
        if *uses != 0 {
            return Err(VerifyError::Trace {
                node,
                fault: NodeFault::UseCountMismatch,
            });
        }
    }
    Ok(())
}

/// The number of unordered pairs of a basis of `count` elements.
pub(super) fn pair_count(count: usize) -> u64 {
    let count = count as u64;
    if count < 2 {
        return 0;
    }
    count * (count - 1) / 2
}

/// Check that G is monic, sorted, minimal, and interreduced.
fn check_shape(basis: &[Poly], ctx: &mut Ctx<'_>) -> Result<(), VerifyError> {
    check_elements(basis, ctx)?;
    check_reduced(basis, ctx)
}

fn check_elements(basis: &[Poly], ctx: &mut Ctx<'_>) -> Result<(), VerifyError> {
    for (index, element) in basis.iter().enumerate() {
        check_element(index, element, basis, ctx)?;
    }
    Ok(())
}

fn check_element(
    index: usize,
    element: &Poly,
    basis: &[Poly],
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let Some(lead) = element.lm() else {
        return Err(VerifyError::Basis {
            index,
            fault: BasisFault::Zero,
        });
    };
    let lc = element.lc().unwrap_or_default();
    if lc != 1 {
        return Err(VerifyError::Basis {
            index,
            fault: BasisFault::NotMonic { lc },
        });
    }
    if index == 0 {
        return Ok(());
    }
    let previous = basis[index - 1]
        .lm()
        .expect("an earlier element is nonzero");
    arith::charge_monos(&mut ctx.meter, previous, lead)?;
    if previous.cmp_grevlex(lead) != std::cmp::Ordering::Greater {
        return Err(VerifyError::Basis {
            index,
            fault: BasisFault::NotDescending,
        });
    }
    Ok(())
}

fn check_reduced(basis: &[Poly], ctx: &mut Ctx<'_>) -> Result<(), VerifyError> {
    for (index, element) in basis.iter().enumerate() {
        for (term, entry) in element.terms().iter().enumerate() {
            check_term_reduced(index, term, &entry.mono, basis, ctx)?;
        }
    }
    Ok(())
}

fn check_term_reduced(
    index: usize,
    term: usize,
    mono: &Mono,
    basis: &[Poly],
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    for (by, other) in basis.iter().enumerate() {
        if by == index {
            continue;
        }
        let divisor = other.lm().expect("every element is nonzero");
        arith::charge_monos(&mut ctx.meter, divisor, mono)?;
        if divisor.divides(mono) {
            let fault = match term {
                0 => BasisFault::LeadDivisible { by },
                _ => BasisFault::TailReducible { term, by },
            };
            return Err(VerifyError::Basis { index, fault });
        }
    }
    Ok(())
}

/// Compare `mono * lm(g)` with a target, or report exponent overflow.
pub(super) fn lead_matches(
    mono: &Mono,
    element: &Poly,
    target: &Mono,
    meter: &mut Meter,
) -> Result<bool, VerifyError> {
    let lead = element.lm().expect("every basis element is nonzero");
    arith::charge_monos(meter, mono, lead)?;
    let (matches, support) = mono.mul_matches(lead, target)?;
    meter.charge((support + target.support()).max(1) as u64)?;
    Ok(matches)
}
