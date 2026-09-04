//! One witness per pair of the basis (contract sections 6 and 9.5).

use crate::certificate::{CertificateCap, CertifyError, EmitterFault};
use crate::poly::Polynomial;
use crate::ring::PrimeOps;

use super::WriterBudget;
use super::divide::{Step, divide};
use super::{MAX_PAIRS, cap};

/// Why one pair of the basis has an lcm representation.
pub(super) enum Witness {
    /// The two leading monomials share no variable.
    Coprime,
    /// The pairs on `{i, k}` and on `{k, j}` cover this pair.
    Chain(u32),
    /// A division trace of the S-polynomial, ending at zero.
    Reduce(Vec<Step>),
}

struct Search {
    deps: Vec<Option<[u32; 2]>>,
    seen: Vec<u32>,
    stack: Vec<u32>,
    generation: u32,
}

/// The index of pair `(i, j)` with `i < j` (contract section 5.9).
pub(super) fn pair_index(i: usize, j: usize, count: usize) -> usize {
    i * count - i * (i + 1) / 2 + (j - i - 1)
}

/// The number of unordered pairs of a basis of `count` elements.
///
/// `None` when the count does not fit a `usize`, which no basis this
/// writer accepts reaches: the pairs cap stops it first.
pub(super) fn pair_count(count: usize) -> Option<usize> {
    if count < 2 {
        return Some(0);
    }
    count.checked_mul(count - 1).map(|product| product / 2)
}

/// Choose one witness per pair, in lexicographic pair order.
///
/// The choice is the one contract section 9.5 fixes: `Coprime` first, then
/// the first eligible `k` whose two dependency edges keep the graph
/// acyclic, then `Reduce`. The verifier validates the witnesses in reverse
/// topological order, so the writer must leave the graph acyclic.
///
/// Every element of `basis` is monic and nonzero, which the caller checks.
pub(super) fn witnesses(
    basis: &[Polynomial],
    ops: &PrimeOps,
    budget: &mut WriterBudget,
) -> Result<Vec<Witness>, CertifyError> {
    let count = basis.len();
    let Some(pairs) = pair_count(count) else {
        return Err(CertifyError::CapExceeded {
            cap: CertificateCap::Pairs,
            limit: MAX_PAIRS,
        });
    };
    cap(pairs, MAX_PAIRS, CertificateCap::Pairs)?;
    let per_pair = size_of::<Witness>() + size_of::<Option<[u32; 2]>>() + size_of::<u32>();
    budget.hold_bytes(pairs.saturating_mul(per_pair))?;

    let mut out: Vec<Witness> = Vec::with_capacity(pairs);
    let mut search = Search {
        deps: vec![None; pairs],
        seen: vec![0u32; pairs],
        stack: Vec::new(),
        generation: 0,
    };

    for i in 0..count {
        budget.check_stop()?;
        for j in i + 1..count {
            out.push(choose_witness(basis, i, j, ops, budget, &mut search)?);
        }
    }
    Ok(out)
}

fn choose_witness(
    basis: &[Polynomial],
    i: usize,
    j: usize,
    ops: &PrimeOps,
    budget: &WriterBudget,
    search: &mut Search,
) -> Result<Witness, CertifyError> {
    let left = lead(basis, i);
    let right = lead(basis, j);
    if left.is_coprime(right) {
        return Ok(Witness::Coprime);
    }
    let target = pair_index(i, j, basis.len()) as u32;
    if let Some((k, edges)) = find_chain(basis, i, j, target, &left.lcm(right), search) {
        search.deps[target as usize] = Some(edges);
        return Ok(Witness::Chain(k));
    }
    reduction_witness(basis, i, j, ops, budget)
}

fn find_chain(
    basis: &[Polynomial],
    i: usize,
    j: usize,
    target: u32,
    lcm: &crate::poly::Monomial,
    search: &mut Search,
) -> Option<(u32, [u32; 2])> {
    for k in 0..basis.len() {
        if !eligible_chain(basis, i, j, k, lcm) {
            continue;
        }
        let edges = chain_edges(i, j, k, basis.len());
        search.generation += 1;
        if closes_cycle(edges, target, search) {
            continue;
        }
        return Some((k as u32, edges));
    }
    None
}

fn eligible_chain(
    basis: &[Polynomial],
    i: usize,
    j: usize,
    k: usize,
    lcm: &crate::poly::Monomial,
) -> bool {
    k != i && k != j && lead(basis, k).divides(lcm)
}

fn chain_edges(i: usize, j: usize, k: usize, count: usize) -> [u32; 2] {
    [
        pair_index(i.min(k), i.max(k), count) as u32,
        pair_index(k.min(j), k.max(j), count) as u32,
    ]
}

fn closes_cycle(edges: [u32; 2], target: u32, search: &mut Search) -> bool {
    reaches(
        &search.deps,
        edges[0],
        target,
        &mut search.seen,
        &mut search.stack,
        search.generation,
    ) || reaches(
        &search.deps,
        edges[1],
        target,
        &mut search.seen,
        &mut search.stack,
        search.generation,
    )
}

fn reduction_witness(
    basis: &[Polynomial],
    i: usize,
    j: usize,
    ops: &PrimeOps,
    budget: &WriterBudget,
) -> Result<Witness, CertifyError> {
    let spoly = basis[i].s_polynomial(&basis[j], ops);
    let Some(steps) = divide(&spoly, basis, ops, budget)? else {
        return Err(EmitterFault::NotAGroebnerBasis { i, j }.into());
    };
    Ok(Witness::Reduce(steps))
}

/// The leading monomial of a basis element.
fn lead(basis: &[Polynomial], index: usize) -> &crate::poly::Monomial {
    basis[index]
        .lm()
        .expect("the writer checked that no basis element is zero")
}

/// Report whether `target` is reachable from `from` in the graph.
///
/// The traversal is iterative, because a recursive one can overflow the
/// stack on a long chain. `seen` holds the generation each node was
/// visited in, so one array serves every search.
fn reaches(
    deps: &[Option<[u32; 2]>],
    from: u32,
    target: u32,
    seen: &mut [u32],
    stack: &mut Vec<u32>,
    generation: u32,
) -> bool {
    if from == target {
        return true;
    }
    stack.clear();
    stack.push(from);
    seen[from as usize] = generation;
    while let Some(node) = stack.pop() {
        let Some(edges) = deps[node as usize] else {
            continue;
        };
        for edge in edges {
            if edge == target {
                return true;
            }
            if seen[edge as usize] != generation {
                seen[edge as usize] = generation;
                stack.push(edge);
            }
        }
    }
    false
}
