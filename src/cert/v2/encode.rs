//! The byte format of `sylv-gb-cert-v2` (contract section 5).
//!
//! Every integer is a canonical varint. The sections are written into
//! their own buffers, and the header then declares the byte length of
//! each of them.

use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

use crate::certificate::{CertificateCap, CertifyError};
use crate::poly::{Monomial, Polynomial};

use super::divide::Step;
use super::record::Node;
use super::witness::Witness;
use super::{MAX_POOL_ENTRIES, MAX_POOL_MONOMIALS, WriterBudget, cap};

/// The eight magic bytes: `SYLVGB`, format 2, revision 0.
const MAGIC: [u8; 8] = [0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00];

/// The schema string.
const SCHEMA: &[u8] = b"sylv-gb-cert-v2";

/// The order string.
const ORDER: &[u8] = b"grevlex-v1";

/// The pool of monomials the certificate references.
///
/// The map is ordered by the crate's grevlex comparison, so iteration
/// gives the strictly ascending pool the contract asks for.
pub(super) struct Pool {
    index: BTreeMap<Monomial, u32>,
}

impl Pool {
    /// Number the monomials of `monomials`, ascending.
    ///
    /// The two pool caps and the budget run on each new monomial before it
    /// enters the map, so a pool above a cap costs the writer one monomial
    /// and never the whole pool.
    pub(super) fn build(
        monomials: impl IntoIterator<Item = Monomial>,
        budget: &mut WriterBudget,
    ) -> Result<Self, CertifyError> {
        let mut index: BTreeMap<Monomial, u32> = BTreeMap::new();
        let mut entries = 0usize;
        for mono in monomials {
            let support = mono.exps.iter().filter(|&&exp| exp != 0).count();
            let Entry::Vacant(slot) = index.entry(mono) else {
                continue;
            };
            cap(
                entries.saturating_add(support),
                MAX_POOL_ENTRIES,
                CertificateCap::PoolEntries,
            )?;
            // One pool monomial costs at most what one term costs: the
            // header and one exponent per variable.
            budget.hold(0, 1)?;
            entries += support;
            slot.insert(0);
            cap(
                index.len(),
                MAX_POOL_MONOMIALS,
                CertificateCap::PoolMonomials,
            )?;
        }
        for (position, slot) in index.values_mut().enumerate() {
            *slot = position as u32;
        }
        Ok(Pool { index })
    }

    /// The number of monomials.
    pub(super) fn len(&self) -> usize {
        self.index.len()
    }

    /// The pool index of `mono`.
    fn at(&self, mono: &Monomial) -> u32 {
        *self
            .index
            .get(mono)
            .expect("the writer put every monomial it references into the pool")
    }
}

/// What one certificate holds.
pub(super) struct Parts<'a> {
    /// The prime.
    pub(super) modulus: u64,
    /// The number of variables.
    pub(super) nvars: usize,
    /// The monomial pool.
    pub(super) pool: &'a Pool,
    /// The input list, in the caller's order.
    pub(super) input: &'a [Polynomial],
    /// The nodes, in emission order.
    pub(super) nodes: &'a [Node],
    /// The use count of each node.
    pub(super) uses: &'a [u32],
    /// The node of each basis element.
    pub(super) basis: &'a [u32],
    /// One division trace per input polynomial.
    pub(super) membership: &'a [Vec<Step>],
    /// One witness per pair, in lexicographic pair order.
    pub(super) pairs: &'a [Witness],
}

/// One section writer, in the order the header declares the sections.
type Section = fn(&mut Vec<u8>, &Parts<'_>);

const SECTIONS: [Section; 6] = [
    pool_section,
    input_section,
    trace_section,
    basis_section,
    membership_section,
    pairs_section,
];

/// Write one certificate.
///
/// The header declares the byte length of every section, so the sections
/// go into one buffer and the header goes in front of it. The writer holds
/// that buffer and the certificate at once and nothing else, so its peak
/// is the certificate twice over. The budget takes each section as it is
/// written and the copy before it is made. A section is charged after it
/// is written, so the writer passes the limit by at most one section.
pub(super) fn write(parts: &Parts<'_>, budget: &mut WriterBudget) -> Result<Vec<u8>, CertifyError> {
    let mut body: Vec<u8> = Vec::new();
    let mut lengths = [0usize; SECTIONS.len()];
    for (section, length) in SECTIONS.iter().zip(lengths.iter_mut()) {
        let start = body.len();
        section(&mut body, parts);
        *length = body.len() - start;
        budget.hold_bytes(*length)?;
    }

    let mut out = Vec::new();
    out.extend_from_slice(&MAGIC);
    put_string(&mut out, SCHEMA);
    put_string(&mut out, ORDER);
    put_varint(&mut out, parts.modulus);
    put_varint(&mut out, parts.nvars as u64);
    for length in lengths {
        put_varint(&mut out, length as u64);
    }
    budget.check_bytes(body.len())?;
    out.extend_from_slice(&body);
    budget.release_bytes(body.len());
    Ok(out)
}

/// Append `value` in base-128 little-endian form, canonically.
fn put_varint(out: &mut Vec<u8>, value: u64) {
    let mut rest = value;
    while rest >= 0x80 {
        out.push((rest as u8) | 0x80);
        rest >>= 7;
    }
    out.push(rest as u8);
}

/// Append a length and that many bytes.
fn put_string(out: &mut Vec<u8>, text: &[u8]) {
    put_varint(out, text.len() as u64);
    out.extend_from_slice(text);
}

/// The monomial pool, ascending, each as a sparse support.
fn pool_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    put_varint(out, parts.pool.len() as u64);
    for mono in parts.pool.index.keys() {
        let support = mono.exps.iter().filter(|&&exp| exp != 0).count();
        put_varint(out, support as u64);
        for (variable, &exp) in mono.exps.iter().enumerate() {
            if exp == 0 {
                continue;
            }
            put_varint(out, variable as u64);
            put_varint(out, u64::from(exp));
        }
    }
}

/// The input list, terms strictly descending.
fn input_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    put_varint(out, parts.input.len() as u64);
    for poly in parts.input {
        put_varint(out, poly.terms.len() as u64);
        for term in poly.terms.iter().rev() {
            put_varint(out, term.coeff.value());
            put_varint(out, u64::from(parts.pool.at(&term.mono)));
        }
    }
}

/// The nodes, each with its kind, its use count, and its payload.
fn trace_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    put_varint(out, parts.nodes.len() as u64);
    for (node, uses) in parts.nodes.iter().zip(parts.uses) {
        match node {
            Node::Input { index } => {
                put_varint(out, 0);
                put_varint(out, u64::from(*uses));
                put_varint(out, u64::from(*index));
            }
            Node::Mul { src, mono } => {
                put_varint(out, 1);
                put_varint(out, u64::from(*uses));
                put_varint(out, u64::from(*src));
                put_varint(out, u64::from(parts.pool.at(mono)));
            }
            Node::Comb { steps } => {
                put_varint(out, 2);
                put_varint(out, u64::from(*uses));
                put_varint(out, steps.len() as u64);
                let mut previous: Option<u32> = None;
                for &(src, scalar) in steps {
                    match previous {
                        None => put_varint(out, u64::from(src)),
                        Some(last) => put_varint(out, u64::from(src - last - 1)),
                    }
                    put_varint(out, scalar);
                    previous = Some(src);
                }
            }
            Node::Scale { src, scalar } => {
                put_varint(out, 3);
                put_varint(out, u64::from(*uses));
                put_varint(out, u64::from(*src));
                put_varint(out, *scalar);
            }
        }
    }
}

/// The node of each basis element.
fn basis_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    put_varint(out, parts.basis.len() as u64);
    for &node in parts.basis {
        put_varint(out, u64::from(node));
    }
}

/// One division trace per input polynomial, in input order.
fn membership_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    for steps in parts.membership {
        put_steps(out, parts, steps);
    }
}

/// One witness per pair, in lexicographic pair order.
fn pairs_section(out: &mut Vec<u8>, parts: &Parts<'_>) {
    for witness in parts.pairs {
        match witness {
            Witness::Coprime => put_varint(out, 0),
            Witness::Chain(k) => {
                put_varint(out, 1);
                put_varint(out, u64::from(*k));
            }
            Witness::Reduce(steps) => {
                put_varint(out, 2);
                put_steps(out, parts, steps);
            }
        }
    }
}

/// One division trace: a step count and the steps.
fn put_steps(out: &mut Vec<u8>, parts: &Parts<'_>, steps: &[Step]) {
    put_varint(out, steps.len() as u64);
    for step in steps {
        put_varint(out, u64::from(parts.pool.at(&step.mono)));
        put_varint(out, u64::from(step.element));
    }
}
