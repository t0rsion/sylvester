//! The recorder that turns the engine's reports into certificate nodes.
//!
//! The engine reports one row at a time through
//! [`crate::compute::f4::trace::Trace`]. This module lowers a row into the
//! nodes of `docs/certificate-v2.md` section 9.3 and keeps the two maps
//! that section names: batch row to node, and global basis index to node.
//!
//! The recorder never stops the engine. It latches the first fault it
//! meets and records nothing after it, and [`Recorder::finish`] reports
//! the fault to the writer, which then writes no bytes.

use std::time::Instant;

use crate::certificate::{CertifyError, TraceFault};
use crate::compute::ComputeError;
use crate::compute::f4::trace::{BatchRows, Insertion, Returned, RowId, Trace};
use crate::poly::{Exps, Monomial, Polynomial};
use crate::ring::PolynomialRing;
use crate::ring::field::Felt;

use super::Budget;

/// One node of the operation DAG (contract section 4.6).
///
/// A source is the recording index of an earlier node, so every reference
/// points backward.
pub(super) enum Node {
    /// The input polynomial at this index.
    Input {
        /// The index in the caller's list.
        index: u32,
    },
    /// The monomial multiple of one node.
    Mul {
        src: u32,
        /// The multiplier, which is not the identity monomial.
        mono: Monomial,
    },
    /// The scalar multiple of one node.
    Scale {
        src: u32,
        /// The scalar, which is at least 2.
        scalar: u64,
    },
    /// The sum of at least two scaled nodes, sorted by node.
    Comb {
        /// The steps, by strictly increasing node.
        steps: Vec<(u32, u64)>,
    },
}

impl Node {
    /// Call `f` on every source of this node.
    pub(super) fn sources(&self, mut f: impl FnMut(u32)) {
        match self {
            Node::Input { .. } => {}
            Node::Mul { src, .. } | Node::Scale { src, .. } => f(*src),
            Node::Comb { steps } => {
                for &(src, _) in steps {
                    f(src);
                }
            }
        }
    }
}

/// What one recorded run gives the writer.
pub(super) struct Recording {
    /// The nodes, in recording order.
    pub(super) nodes: Vec<Node>,
    /// The node of each basis element the run returned, in return order.
    pub(super) basis: Vec<u32>,
}

/// The normalization of the input list (contract section 9.1).
///
/// The engine drops the zero generators, makes each survivor monic, drops
/// a survivor equal to an earlier one, and sorts the rest by leading
/// monomial. The writer repeats that work, because the trace must show the
/// basis over the caller's list and not over the normalized one.
struct Normalization {
    /// The scalar that makes each generator monic, or `None` when the
    /// generator is zero or already monic.
    scales: Vec<Option<u64>>,
    /// The input index of each survivor, in insertion order, so survivor
    /// `j` holds global basis index `j`.
    survivors: Vec<u32>,
}

impl Normalization {
    fn of(input: &[Polynomial], modulus: u64) -> Self {
        let mut scales = vec![None; input.len()];
        let mut kept: Vec<(u32, Polynomial)> = Vec::new();
        for (index, poly) in input.iter().enumerate() {
            let Some(lc) = poly.lc() else { continue };
            let monic = poly.make_monic(modulus);
            if kept.iter().any(|(_, seen)| *seen == monic) {
                continue;
            }
            if lc != Felt::one() {
                scales[index] = Some(lc.inv(modulus).value());
            }
            kept.push((index as u32, monic));
        }
        kept.sort_by(|a, b| a.1.lm().cmp(&b.1.lm()).then(a.0.cmp(&b.0)));
        Normalization {
            scales,
            survivors: kept.into_iter().map(|(index, _)| index).collect(),
        }
    }
}

/// The recorder of one F4 run.
///
/// The caller builds it before the run, hands it to the engine, and calls
/// [`Recorder::finish`] after the run.
pub(crate) struct Recorder {
    modulus: u64,
    nvars: usize,
    normalization: Normalization,
    inputs: usize,
    nodes: Vec<Node>,
    /// The node of the monic form of each generator.
    monic: Vec<u32>,
    /// The node of each global basis index.
    basis: Vec<Option<u32>>,
    /// The node of each pivot slot of the current batch.
    pivots: Vec<Option<u32>>,
    /// The node of each lower slot of the current batch.
    lower: Vec<Option<u32>>,
    /// The summands of the row the kernel is reducing.
    summands: Vec<(u32, u64)>,
    row: Option<RowId>,
    normalizer: u64,
    returned: Vec<u32>,
    fault: Option<CertifyError>,
    budget: Budget,
    /// The bytes the nodes of the current attempt hold.
    charged: usize,
}

/// The bytes one recorded node costs the budget.
///
/// A `Comb` step costs its own slot on top.
const NODE_BYTES: usize = size_of::<Node>();

/// The bytes one `Comb` step costs the budget.
const STEP_BYTES: usize = size_of::<(u32, u64)>();

impl Recorder {
    /// A recorder for one run over `input`.
    ///
    /// The deadline and the memory limit are the ones the run holds to.
    /// They cover the nodes the recorder keeps, as they cover the engine.
    pub(crate) fn new(
        ring: &PolynomialRing,
        input: &[Polynomial],
        deadline: Option<Instant>,
        memory_limit: Option<usize>,
    ) -> Self {
        let modulus = ring.modulus();
        let mut recorder = Recorder {
            modulus,
            nvars: ring.nvars(),
            normalization: Normalization::of(input, modulus),
            inputs: input.len(),
            nodes: Vec::new(),
            monic: Vec::new(),
            basis: Vec::new(),
            pivots: Vec::new(),
            lower: Vec::new(),
            summands: Vec::new(),
            row: None,
            normalizer: 1,
            returned: Vec::new(),
            fault: None,
            budget: Budget::new(deadline, memory_limit, ring.nvars()),
            charged: 0,
        };
        recorder.seed();
        recorder
    }

    /// The recording, or the first fault the recorder met.
    ///
    /// The budget comes back with it, because the writer charges the rest
    /// of its work to the same budget.
    pub(super) fn finish(self) -> Result<(Recording, Budget), CertifyError> {
        if let Some(fault) = self.fault {
            return Err(fault);
        }
        let recording = Recording {
            nodes: self.nodes,
            basis: self.returned,
        };
        Ok((recording, self.budget))
    }

    /// Write the input nodes and the normalization nodes (contract section
    /// 9.2, steps 1 and 2).
    fn seed(&mut self) {
        for index in 0..self.inputs {
            self.push(Node::Input {
                index: index as u32,
            });
        }
        self.monic = (0..self.inputs as u32).collect();
        for index in 0..self.inputs {
            let Some(scalar) = self.normalization.scales[index] else {
                continue;
            };
            let node = self.push(Node::Scale {
                src: index as u32,
                scalar,
            });
            self.monic[index] = node;
        }
        self.basis = self
            .normalization
            .survivors
            .iter()
            .map(|&input| Some(self.monic[input as usize]))
            .collect();
    }

    /// Drop every node and every map, and record from the start.
    ///
    /// The attempt's nodes go with it, so the budget takes their bytes
    /// back. Without that the next attempt starts against a budget the
    /// dropped attempt still holds.
    fn reset(&mut self) {
        self.budget.release_bytes(self.charged);
        self.charged = 0;
        self.nodes.clear();
        self.monic.clear();
        self.basis.clear();
        self.pivots.clear();
        self.lower.clear();
        self.summands.clear();
        self.row = None;
        self.normalizer = 1;
        self.returned.clear();
        self.seed();
    }

    /// Append one node and return its recording index.
    fn push(&mut self, node: Node) -> u32 {
        let steps = match &node {
            Node::Comb { steps } => steps.len(),
            _ => 0,
        };
        let bytes = NODE_BYTES + steps * STEP_BYTES;
        if let Err(error) = self.budget.hold_bytes(bytes) {
            self.stop(error.into());
        }
        self.charged += bytes;
        let index = self.nodes.len() as u32;
        self.nodes.push(node);
        index
    }

    /// Latch the first fault and record nothing after it.
    fn stop(&mut self, fault: CertifyError) {
        if self.fault.is_none() {
            self.fault = Some(fault);
        }
    }

    fn stopped(&self) -> bool {
        self.fault.is_some()
    }

    fn bound(map: &[Option<u32>], slot: u32) -> Option<u32> {
        map.get(slot as usize).copied().flatten()
    }

    /// Bind `slot` of `map` to `node`, growing the map.
    fn bind(map: &mut Vec<Option<u32>>, slot: u32, node: u32) {
        let slot = slot as usize;
        if map.len() <= slot {
            map.resize(slot + 1, None);
        }
        map[slot] = Some(node);
    }

    /// The monomial of one reported multiplier, or `None` for the
    /// identity.
    fn multiplier(&mut self, exps: &[u32]) -> Option<Monomial> {
        if exps.iter().all(|&exp| exp == 0) {
            return None;
        }
        let mut out: Exps = Exps::with_capacity(exps.len());
        for &exp in exps {
            match u16::try_from(exp) {
                Ok(exp) => out.push(exp),
                Err(_) => {
                    self.stop(CertifyError::Engine(ComputeError::ExponentLimit {
                        limit: crate::compute::DEGREE_LIMIT,
                    }));
                    return None;
                }
            }
        }
        Some(Monomial::from_exps(out))
    }

    /// Merge, scale, and sort the summands of the row (contract section
    /// 9.3, steps 3 to 5).
    fn lower_row(&mut self) -> Vec<(u32, u64)> {
        self.summands.sort_by_key(|&(node, _)| node);
        let mut merged: Vec<(u32, u64)> = Vec::with_capacity(self.summands.len());
        for &(node, scalar) in &self.summands {
            match merged.last_mut() {
                Some(last) if last.0 == node => {
                    last.1 = (last.1 + scalar) % self.modulus;
                    if last.1 == 0 {
                        merged.pop();
                    }
                }
                _ => merged.push((node, scalar)),
            }
        }
        if self.normalizer != 1 {
            for step in &mut merged {
                step.1 = Felt::from_residue(step.1)
                    .mul(Felt::from_residue(self.normalizer), self.modulus)
                    .value();
            }
        }
        merged
    }
}

impl Trace for Recorder {
    const RECORDS: bool = true;

    fn held_bytes(&self) -> usize {
        self.budget.held()
    }

    fn restart(&mut self) {
        if self.stopped() {
            return;
        }
        self.reset();
    }

    fn rows(&mut self, rows: &BatchRows<'_>) {
        if self.stopped() {
            return;
        }
        if let Err(error) = self.budget.check_deadline() {
            self.stop(error.into());
            return;
        }
        debug_assert_eq!(
            rows.nvars, self.nvars,
            "the report holds one exponent per variable"
        );
        self.pivots.clear();
        self.lower.clear();
        for (index, row) in rows.rows.iter().enumerate() {
            let Some(src) = Recorder::bound(&self.basis, row.source) else {
                self.stop(TraceFault::UnboundBasis.into());
                return;
            };
            let node = match self.multiplier(rows.mult(index)) {
                Some(mono) => self.push(Node::Mul { src, mono }),
                None if self.stopped() => return,
                None => src,
            };
            match row.place {
                RowId::Pivot(slot) => Recorder::bind(&mut self.pivots, slot, node),
                RowId::Lower(slot) => Recorder::bind(&mut self.lower, slot, node),
            }
        }
    }

    fn pivots_sorted(&mut self, npiv: u32, order: &[u32]) {
        if self.stopped() {
            return;
        }
        let mut sorted: Vec<Option<u32>> = Vec::with_capacity(order.len());
        for &from in order {
            sorted.push(self.pivots.get(from as usize).copied().flatten());
        }
        self.pivots.truncate(npiv as usize);
        self.pivots.append(&mut sorted);
    }

    fn inserted(&mut self, entries: &[Insertion]) {
        if self.stopped() {
            return;
        }
        for entry in entries {
            let Some(index) = entry.basis else { continue };
            let Some(node) = Recorder::bound(&self.pivots, entry.pivot) else {
                self.stop(TraceFault::UnboundRow.into());
                return;
            };
            Recorder::bind(&mut self.basis, index, node);
        }
    }

    fn returned(&mut self, entries: &[Returned]) {
        if self.stopped() {
            return;
        }
        self.returned.clear();
        for entry in entries {
            let node = match *entry {
                Returned::Pivot(slot) => Recorder::bound(&self.pivots, slot),
                Returned::Input(index) => self.monic.get(index as usize).copied(),
            };
            let Some(node) = node else {
                self.stop(TraceFault::UnboundBasis.into());
                return;
            };
            self.returned.push(node);
        }
    }

    fn start(&mut self, row: RowId) {
        if self.stopped() {
            return;
        }
        self.summands.clear();
        self.normalizer = 1;
        self.row = Some(row);
        let node = match row {
            RowId::Pivot(slot) => Recorder::bound(&self.pivots, slot),
            RowId::Lower(slot) => Recorder::bound(&self.lower, slot),
        };
        let Some(node) = node else {
            self.stop(TraceFault::UnboundRow.into());
            return;
        };
        self.summands.push((node, 1));
    }

    fn step(&mut self, pivot: u32, factor: u32) {
        if self.stopped() {
            return;
        }
        let Some(node) = Recorder::bound(&self.pivots, pivot) else {
            self.stop(TraceFault::UnboundRow.into());
            return;
        };
        self.summands.push((node, u64::from(factor)));
    }

    fn scale(&mut self, scalar: u32) {
        if self.stopped() {
            return;
        }
        self.normalizer = u64::from(scalar);
    }

    fn end(&mut self, installed: Option<u32>) {
        if self.stopped() {
            return;
        }
        let Some(row) = self.row.take() else { return };
        let Some(slot) = installed else {
            // A row can cancel to zero while its summands do not cancel
            // structurally. The kernel's report decides, so the node the
            // lowering would give is never written.
            return;
        };
        let merged = self.lower_row();
        let node = match merged.len() {
            0 => {
                self.stop(TraceFault::EmptyRow.into());
                return;
            }
            1 if merged[0].1 == 1 => merged[0].0,
            1 => self.push(Node::Scale {
                src: merged[0].0,
                scalar: merged[0].1,
            }),
            _ => self.push(Node::Comb { steps: merged }),
        };
        Recorder::bind(&mut self.pivots, slot, node);
        if let RowId::Lower(row) = row {
            Recorder::bind(&mut self.lower, row, node);
        }
    }
}
