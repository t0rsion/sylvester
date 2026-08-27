//! Writer for the canonical certificate encoding.
//!
//! The writer emits the `sylv-gb-cert-v1` object: the nine keys in the
//! contract order, no insignificant whitespace, and decimal integers.
//! Engine terms run ascending under the order and certificate terms run
//! descending, so the writer reverses every term list. The same values
//! always give the same bytes.

use crate::compute::ComputeError;
use crate::poly::Polynomial;

use super::Budget;
use super::represent::SpairEntry;

const SCHEMA: &str = "sylv-gb-cert-v1";
const ORDER: &str = "grevlex-v1";

/// The pieces of one certificate, each already in certificate order.
pub(crate) struct Parts<'a> {
    pub modulus: u64,
    pub nvars: usize,
    pub input: &'a [Polynomial],
    pub basis: &'a [Polynomial],
    pub origin: &'a [Vec<Polynomial>],
    pub membership: &'a [Vec<Polynomial>],
    pub spairs: &'a [SpairEntry],
}

/// Write the canonical bytes, or stop when the buffer would pass the limit.
///
/// The writer charges the growing buffer against `budget` once per
/// polynomial. A certificate too large for the memory limit stops within
/// one polynomial of the limit instead of after the whole buffer is built.
/// The budget stops a write; it never changes the bytes.
pub(crate) fn write(parts: &Parts<'_>, budget: &Budget) -> Result<Vec<u8>, ComputeError> {
    let mut out = Vec::new();
    out.extend_from_slice(br#"{"schema":""#);
    out.extend_from_slice(SCHEMA.as_bytes());
    out.extend_from_slice(br#"","order":""#);
    out.extend_from_slice(ORDER.as_bytes());
    out.extend_from_slice(br#"","modulus":"#);
    uint(&mut out, parts.modulus);
    out.extend_from_slice(br#","nvars":"#);
    uint(&mut out, parts.nvars as u64);
    out.extend_from_slice(br#","input":"#);
    poly_array(&mut out, parts.input, parts.modulus, budget)?;
    out.extend_from_slice(br#","basis":"#);
    poly_array(&mut out, parts.basis, parts.modulus, budget)?;
    out.extend_from_slice(br#","origin":"#);
    entry_array(&mut out, parts.origin, parts.modulus, budget)?;
    out.extend_from_slice(br#","membership":"#);
    entry_array(&mut out, parts.membership, parts.modulus, budget)?;
    out.extend_from_slice(br#","spairs":["#);
    for (index, entry) in parts.spairs.iter().enumerate() {
        if index > 0 {
            out.push(b',');
        }
        budget.check_deadline()?;
        out.push(b'[');
        uint(&mut out, entry.i as u64);
        out.push(b',');
        uint(&mut out, entry.j as u64);
        out.push(b',');
        poly_array(&mut out, &entry.cofactors, parts.modulus, budget)?;
        out.push(b']');
    }
    out.extend_from_slice(b"]}");
    Ok(out)
}

fn uint(out: &mut Vec<u8>, value: u64) {
    let mut digits = [0u8; 20];
    let mut index = digits.len();
    let mut left = value;
    loop {
        index -= 1;
        digits[index] = b'0' + (left % 10) as u8;
        left /= 10;
        if left == 0 {
            break;
        }
    }
    out.extend_from_slice(&digits[index..]);
}

fn entry_array(
    out: &mut Vec<u8>,
    entries: &[Vec<Polynomial>],
    modulus: u64,
    budget: &Budget,
) -> Result<(), ComputeError> {
    out.push(b'[');
    for (index, entry) in entries.iter().enumerate() {
        if index > 0 {
            out.push(b',');
        }
        budget.check_deadline()?;
        poly_array(out, entry, modulus, budget)?;
    }
    out.push(b']');
    Ok(())
}

fn poly_array(
    out: &mut Vec<u8>,
    polys: &[Polynomial],
    modulus: u64,
    budget: &Budget,
) -> Result<(), ComputeError> {
    out.push(b'[');
    for (index, poly) in polys.iter().enumerate() {
        if index > 0 {
            out.push(b',');
        }
        polynomial(out, poly, modulus);
        budget.check_bytes(out.len())?;
    }
    out.push(b']');
    Ok(())
}

/// Write one polynomial, terms strictly descending.
///
/// A coefficient outside [1, p-1] has no canonical form, so the writer
/// reduces it and drops the term it reduces to zero. The polynomial it
/// writes holds the value it was given.
fn polynomial(out: &mut Vec<u8>, poly: &Polynomial, modulus: u64) {
    out.push(b'[');
    let mut written = 0usize;
    for term in poly.terms.iter().rev() {
        let coeff = term.coeff.value() % modulus;
        if coeff == 0 {
            continue;
        }
        if written > 0 {
            out.push(b',');
        }
        written += 1;
        out.push(b'[');
        uint(out, coeff);
        out.extend_from_slice(b",[");
        for (index, exp) in term.mono.exps.iter().enumerate() {
            if index > 0 {
                out.push(b',');
            }
            uint(out, *exp as u64);
        }
        out.extend_from_slice(b"]]");
    }
    out.push(b']');
}
