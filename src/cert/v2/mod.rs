//! Certificate writing for the `sylv-gb-cert-v2` contract.
//!
//! This module is untrusted. It reads the trace the F4 engine recorded and
//! writes the canonical bytes of `docs/certificate-v2.md`. The verifier in
//! [`crate::verify`] is the trust boundary: it reads those bytes back and
//! accepts or rejects them. A wrong node here becomes a rejection there,
//! never a silent repair. This module never imports [`crate::verify`]:
//! [`crate::compute::groebner_basis_certified`] runs the writer here, then
//! hands the bytes to the verifier itself.
//!
//! [`Recorder`] records one run. [`assemble`] then writes one certificate
//! from the recording, the input, and the basis the run returned. It holds
//! to the deadline and the memory limit of the run, and it charges the
//! certificate bytes to that budget before it returns them.
//!
//! The writer generates the membership traces and the pair witnesses
//! itself, by division over the final basis.

mod divide;
mod encode;
mod record;
mod witness;

use crate::certificate::{CertificateCap, CertifyError, EmitterFault, Place, TraceFault};
use crate::poly::{Monomial, Polynomial};
use crate::ring::PolynomialRing;

use super::{WriterBudget, check_width};
use divide::{Step, divide};
use encode::{Parts, Pool};
use record::{Node, Recording};
use witness::{Witness, pair_count, witnesses};

pub(crate) use record::Recorder;

struct PrunedTrace {
    nodes: Vec<Node>,
    roots: Vec<u32>,
    uses: Vec<u32>,
}

/// The largest number of monomials in the pool.
const MAX_POOL_MONOMIALS: usize = 4_000_000;
/// The largest number of variable and exponent pairs in the pool.
const MAX_POOL_ENTRIES: usize = 16_000_000;
/// The largest number of polynomials in the input section.
const MAX_INPUT_POLYS: usize = 1_000_000;
/// The largest number of terms in one polynomial.
const MAX_TERMS_PER_POLY: usize = 1 << 20;
/// The largest number of input terms in the certificate.
const MAX_TOTAL_TERMS: usize = 8_000_000;
/// The largest number of nodes in the trace.
const MAX_NODES: usize = 8_000_000;
/// The largest number of steps in one `Comb` node.
const MAX_COMB_STEPS: usize = 1 << 20;
/// The largest number of `Comb` steps in the whole trace.
const MAX_TRACE_STEPS: usize = 64_000_000;
/// The largest number of basis elements.
const MAX_BASIS: usize = 1_000_000;
/// The largest number of pairs of the basis.
const MAX_PAIRS: usize = 4_000_000;
/// The largest number of steps in one division trace.
const MAX_DIVISION_STEPS: usize = 1 << 20;
/// The largest number of division steps in the whole certificate.
const MAX_TOTAL_DIVISION_STEPS: usize = 64_000_000;

/// Return an error when `value` exceeds `limit`.
///
/// The verifier applies the same cap before it allocates, so a certificate
/// above a cap is bytes nobody can check. The writer stops instead.
fn cap(value: usize, limit: usize, cap: CertificateCap) -> Result<(), CertifyError> {
    if value > limit {
        return Err(CertifyError::CapExceeded { cap, limit });
    }
    Ok(())
}

/// Write the certificate of one recorded run.
///
/// `input` is the input F in the caller's order, `basis` the basis the run
/// returned, and `recorder` the trace of that run. The bytes are a
/// function of those three alone. The budget stops the work early; it
/// never changes the bytes.
///
/// The writer reports a defect only where it cannot write a certificate at
/// all: a division that leaves a remainder, a basis element that is zero
/// or not monic, an exponent vector of the wrong width, or a trace that
/// does not describe the basis. It checks no identity of its own: the
/// verifier does that.
///
/// The second value is what the budget holds, which the caller uses to
/// size the verifier's caps.
pub(crate) fn assemble(
    ring: &PolynomialRing,
    input: &[Polynomial],
    basis: &[Polynomial],
    recorder: Recorder,
) -> Result<(Vec<u8>, usize), CertifyError> {
    let modulus = ring.modulus();
    let nvars = ring.nvars();
    let (recording, mut budget) = recorder.finish()?;
    validate_parts(input, basis, &recording, nvars)?;
    hold_parts(input, basis, &mut budget)?;
    let trace = prepare_trace(recording, &mut budget)?;
    let (membership, divisions) = membership_traces(input, basis, ring.ops(), &mut budget)?;
    let pairs = pair_witnesses(basis, ring.ops(), divisions, &mut budget)?;
    write_certificate(
        modulus,
        nvars,
        input,
        &trace,
        &membership,
        &pairs,
        &mut budget,
    )
}

fn validate_parts(
    input: &[Polynomial],
    basis: &[Polynomial],
    recording: &Recording,
    nvars: usize,
) -> Result<(), CertifyError> {
    if recording.basis.len() != basis.len() {
        return Err(TraceFault::BasisCount {
            found: recording.basis.len(),
            expected: basis.len(),
        }
        .into());
    }
    cap(input.len(), MAX_INPUT_POLYS, CertificateCap::InputPolys)?;
    cap(basis.len(), MAX_BASIS, CertificateCap::Basis)?;
    validate_input(input, nvars)?;
    validate_basis(basis, nvars)
}

fn validate_input(input: &[Polynomial], nvars: usize) -> Result<(), CertifyError> {
    let mut terms = 0usize;
    for (index, poly) in input.iter().enumerate() {
        check_width(poly, Place::Input(index), nvars)?;
        cap(
            poly.terms.len(),
            MAX_TERMS_PER_POLY,
            CertificateCap::TermsPerPoly,
        )?;
        terms += poly.terms.len();
    }
    cap(terms, MAX_TOTAL_TERMS, CertificateCap::TotalTerms)
}

fn validate_basis(basis: &[Polynomial], nvars: usize) -> Result<(), CertifyError> {
    for (index, poly) in basis.iter().enumerate() {
        if poly.is_zero() {
            return Err(EmitterFault::BasisElementZero { index }.into());
        }
        if poly.lc() != Some(&crate::ring::field::Felt::one()) {
            return Err(EmitterFault::BasisElementNotMonic { index }.into());
        }
        check_width(poly, Place::Basis(index), nvars)?;
    }
    Ok(())
}

fn hold_parts(
    input: &[Polynomial],
    basis: &[Polynomial],
    budget: &mut WriterBudget,
) -> Result<(), CertifyError> {
    budget.check_stop()?;
    budget.hold_polys(input)?;
    budget.hold_polys(basis)?;
    Ok(())
}

fn prepare_trace(
    recording: Recording,
    budget: &mut WriterBudget,
) -> Result<PrunedTrace, CertifyError> {
    let keep = reachable(&recording);
    let (kept, steps) = trace_size(&recording.nodes, &keep)?;
    cap(kept, MAX_NODES, CertificateCap::Nodes)?;
    cap(steps, MAX_TRACE_STEPS, CertificateCap::TraceSteps)?;
    budget.hold_bytes(pruned_bytes(kept, steps, recording.basis.len()))?;
    let (nodes, roots, uses) = prune(recording, &keep, kept);
    Ok(PrunedTrace { nodes, roots, uses })
}

fn trace_size(nodes: &[Node], keep: &[bool]) -> Result<(usize, usize), CertifyError> {
    let mut kept = 0usize;
    let mut steps = 0usize;
    for (node, &keep) in nodes.iter().zip(keep) {
        if !keep {
            continue;
        }
        kept += 1;
        if let Node::Comb { steps: comb } = node {
            cap(comb.len(), MAX_COMB_STEPS, CertificateCap::CombSteps)?;
            steps += comb.len();
        }
    }
    Ok((kept, steps))
}

fn membership_traces(
    input: &[Polynomial],
    basis: &[Polynomial],
    ops: &crate::ring::PrimeOps,
    budget: &mut WriterBudget,
) -> Result<(Vec<Vec<Step>>, usize), CertifyError> {
    let mut divisions = 0usize;
    let mut membership = Vec::with_capacity(input.len());
    for (index, poly) in input.iter().enumerate() {
        budget.check_stop()?;
        let Some(trace) = divide(poly, basis, ops, budget)? else {
            return Err(EmitterFault::InputHasRemainder { input: index }.into());
        };
        divisions = charge_division(&trace, divisions)?;
        budget.hold(1, trace.len())?;
        membership.push(trace);
    }
    Ok((membership, divisions))
}

fn pair_witnesses(
    basis: &[Polynomial],
    ops: &crate::ring::PrimeOps,
    mut divisions: usize,
    budget: &mut WriterBudget,
) -> Result<Vec<Witness>, CertifyError> {
    let pairs = witnesses(basis, ops, budget)?;
    debug_assert_eq!(
        Some(pairs.len()),
        pair_count(basis.len()),
        "one witness per pair"
    );
    for witness in &pairs {
        if let Witness::Reduce(trace) = witness {
            divisions = charge_division(trace, divisions)?;
        }
    }
    Ok(pairs)
}

fn write_certificate(
    modulus: u64,
    nvars: usize,
    input: &[Polynomial],
    trace: &PrunedTrace,
    membership: &[Vec<Step>],
    pairs: &[Witness],
    budget: &mut WriterBudget,
) -> Result<(Vec<u8>, usize), CertifyError> {
    budget.check_stop()?;
    let pool = Pool::build(monomials(input, &trace.nodes, membership, pairs), budget)?;

    budget.check_stop()?;
    let certificate = encode::write(
        &Parts {
            modulus,
            nvars,
            pool: &pool,
            input,
            nodes: &trace.nodes,
            uses: &trace.uses,
            basis: &trace.roots,
            membership,
            pairs,
        },
        budget,
    )?;
    budget.hold_bytes(certificate.len())?;
    Ok((certificate, budget.held()))
}

/// Add the steps of one division trace to the running total.
fn charge_division(trace: &[Step], total: usize) -> Result<usize, CertifyError> {
    cap(
        trace.len(),
        MAX_DIVISION_STEPS,
        CertificateCap::DivisionSteps,
    )?;
    let total = total + trace.len();
    cap(
        total,
        MAX_TOTAL_DIVISION_STEPS,
        CertificateCap::TotalDivisionSteps,
    )?;
    Ok(total)
}

/// Mark the nodes the basis reaches.
///
/// Every reference points backward, so one backward pass marks every
/// reachable node.
fn reachable(recording: &Recording) -> Vec<bool> {
    let mut keep = vec![false; recording.nodes.len()];
    for &root in &recording.basis {
        keep[root as usize] = true;
    }
    for index in (0..recording.nodes.len()).rev() {
        if !keep[index] {
            continue;
        }
        recording.nodes[index].sources(|src| keep[src as usize] = true);
    }
    keep
}

/// The bytes a pruned trace of `kept` nodes, `steps` combination steps,
/// and `roots` basis roots holds.
fn pruned_bytes(kept: usize, steps: usize, roots: usize) -> usize {
    kept.saturating_mul(size_of::<Node>() + size_of::<u32>())
        .saturating_add(steps.saturating_mul(size_of::<(u32, u64)>()))
        .saturating_add(roots.saturating_mul(size_of::<u32>()))
}

/// Keep the marked nodes, renumbered from 0.
///
/// `keep` comes from [`reachable`] and `kept` is the number of marks in
/// it. The result holds the nodes in recording order, the node of each
/// basis element, and the use count of each node. The use count is the
/// number of references from later nodes and from the basis section, so
/// every kept node has a count of at least 1.
fn prune(recording: Recording, keep: &[bool], kept: usize) -> (Vec<Node>, Vec<u32>, Vec<u32>) {
    let Recording { nodes, basis } = recording;
    let mut renamed = vec![0u32; nodes.len()];
    let mut next = 0u32;
    for (index, &mark) in keep.iter().enumerate() {
        if mark {
            renamed[index] = next;
            next += 1;
        }
    }
    debug_assert_eq!(next as usize, kept, "the marking pass counted the marks");

    let mut uses = vec![0u32; kept];
    let mut out = Vec::with_capacity(kept);
    for (index, node) in nodes.into_iter().enumerate() {
        if !keep[index] {
            continue;
        }
        let node = match node {
            Node::Input { index } => Node::Input { index },
            Node::Mul { src, mono } => Node::Mul {
                src: renamed[src as usize],
                mono,
            },
            Node::Scale { src, scalar } => Node::Scale {
                src: renamed[src as usize],
                scalar,
            },
            Node::Comb { steps } => Node::Comb {
                steps: steps
                    .into_iter()
                    .map(|(src, scalar)| (renamed[src as usize], scalar))
                    .collect(),
            },
        };
        node.sources(|src| uses[src as usize] += 1);
        out.push(node);
    }
    let roots: Vec<u32> = basis.iter().map(|&root| renamed[root as usize]).collect();
    for &root in &roots {
        uses[root as usize] += 1;
    }
    (out, roots, uses)
}

/// Every monomial the certificate references.
///
/// The pool holds these and nothing else: the verifier rejects an entry no
/// section names.
fn monomials<'a>(
    input: &'a [Polynomial],
    nodes: &'a [Node],
    membership: &'a [Vec<Step>],
    pairs: &'a [Witness],
) -> impl Iterator<Item = Monomial> + 'a {
    let terms = input
        .iter()
        .flat_map(|poly| poly.terms.iter().map(|term| term.mono.clone()));
    let multipliers = nodes.iter().filter_map(|node| match node {
        Node::Mul { mono, .. } => Some(mono.clone()),
        _ => None,
    });
    let divisions = membership
        .iter()
        .map(Vec::as_slice)
        .chain(pairs.iter().filter_map(|witness| match witness {
            Witness::Reduce(steps) => Some(steps.as_slice()),
            _ => None,
        }))
        .flat_map(|steps| steps.iter().map(|step| step.mono.clone()));
    terms.chain(multipliers).chain(divisions)
}

#[cfg(test)]
mod tests;
