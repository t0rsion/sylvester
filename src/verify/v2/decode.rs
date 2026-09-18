//! Strict decoder for the canonical v2 binary encoding.
//!
//! The decoder reads the header, the monomial pool, and the input section.
//! It rejects every non-canonical form: a wrong magic, a wrong version, an
//! overlong varint, a section length that does not match the byte count,
//! and an unsorted pool.
//!
//! Every count is compared against its cap and against the bytes left in
//! its section before the vector it sizes is reserved.

use std::cmp::Ordering;

use super::Ctx;
use super::arith::{Exp, MAX_EXP, Mono, Poly, Term, charge_monos, is_prime};
use crate::verify::error::{BinaryFault, Cap, Location, PolyFault, PoolFault, VerifyError};
use crate::verify::limits::Meter;

/// The eight magic bytes: `SYLVGB`, format 2, revision 0.
const MAGIC: [u8; 8] = [0x53, 0x59, 0x4c, 0x56, 0x47, 0x42, 0x02, 0x00];

/// The first byte of the magic. `src/verify/mod.rs` dispatches on it.
pub(crate) const MAGIC_FIRST: u8 = MAGIC[0];

const SCHEMA: &[u8] = b"sylv-gb-cert-v2";
const ORDER: &[u8] = b"grevlex-v1";

/// The largest variable count the contract allows.
const MAX_NVARS: u64 = 256;

/// The largest modulus the contract allows.
const MAX_MODULUS: u64 = (1 << 31) - 1;

/// The number of sections after the header.
pub(super) const SECTIONS: usize = 6;

/// The byte range of one section.
pub(super) type Range = (usize, usize);

/// A cursor over one byte range.
///
/// The cursor never reads past the end of its section, so a varint that
/// would cross a section boundary is a rejection.
pub(super) struct Reader<'a> {
    bytes: &'a [u8],
    pos: usize,
    end: usize,
}

impl<'a> Reader<'a> {
    pub(super) fn new(bytes: &'a [u8], range: Range) -> Self {
        Reader {
            bytes,
            pos: range.0,
            end: range.1,
        }
    }

    fn fail<T>(&self, reason: BinaryFault) -> Result<T, VerifyError> {
        Err(VerifyError::MalformedBinary {
            reason,
            offset: self.pos,
        })
    }

    /// The offset of the cursor in the byte string.
    pub(super) fn offset(&self) -> usize {
        self.pos
    }

    /// The number of bytes left in the section.
    pub(super) fn remaining(&self) -> usize {
        self.end - self.pos
    }

    fn byte(&mut self, meter: &mut Meter) -> Result<u8, VerifyError> {
        if self.pos >= self.end {
            return self.fail(BinaryFault::Truncated);
        }
        meter.charge(1)?;
        let byte = self.bytes[self.pos];
        self.pos += 1;
        Ok(byte)
    }

    /// Read one canonical varint.
    pub(super) fn varint(&mut self, meter: &mut Meter) -> Result<u64, VerifyError> {
        let start = self.pos;
        let mut value: u64 = 0;
        let mut shift = 0u32;
        let mut length = 0u32;
        loop {
            let byte = self.byte(meter)?;
            length += 1;
            if length == 10 {
                // The tenth byte carries the top bit of a u64. A
                // continuation bit on it makes the encoding longer than ten
                // bytes, and any other value does not fit u64.
                if byte & 0x80 != 0 {
                    return Err(VerifyError::MalformedBinary {
                        reason: BinaryFault::VarintTooLong,
                        offset: start,
                    });
                }
                if byte != 0x01 {
                    return Err(VerifyError::MalformedBinary {
                        reason: BinaryFault::VarintOverflow,
                        offset: start,
                    });
                }
            }
            value |= u64::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                if length > 1 && byte == 0 {
                    return Err(VerifyError::MalformedBinary {
                        reason: BinaryFault::NonMinimalVarint,
                        offset: start,
                    });
                }
                return Ok(value);
            }
            shift += 7;
        }
    }

    /// Read a length-prefixed byte string.
    fn string(&mut self, meter: &mut Meter) -> Result<Vec<u8>, VerifyError> {
        let length = self.varint(meter)?;
        if length > self.remaining() as u64 {
            return self.fail(BinaryFault::Truncated);
        }
        let length = length as usize;
        meter.charge(length as u64)?;
        let text = self.bytes[self.pos..self.pos + length].to_vec();
        self.pos += length;
        Ok(text)
    }

    /// Reject a section that leaves bytes unread inside its range.
    pub(super) fn finish(&self) -> Result<(), VerifyError> {
        if self.pos != self.end {
            return self.fail(BinaryFault::SectionNotConsumed);
        }
        Ok(())
    }
}

/// Return an error if a decoded count is above its cap.
pub(super) fn capped(value: u64, limit: usize, cap: Cap) -> Result<usize, VerifyError> {
    if value > limit as u64 {
        return Err(VerifyError::CapExceeded { cap, limit });
    }
    Ok(value as usize)
}

/// Return an error if an index does not address an element.
pub(super) fn bounded(index: u64, bound: usize, what: &'static str) -> Result<usize, VerifyError> {
    if index >= bound as u64 {
        return Err(VerifyError::IndexOutOfRange { what, index, bound });
    }
    Ok(index as usize)
}

/// The header of a v2 certificate.
pub(super) struct Header {
    pub(super) modulus: u64,
    pub(super) nvars: usize,
    pub(super) sections: [Range; SECTIONS],
}

/// Decode the header and split the byte string into its six sections.
pub(super) fn header(bytes: &[u8], meter: &mut Meter) -> Result<Header, VerifyError> {
    check_magic(bytes, meter)?;
    let mut reader = Reader::new(bytes, (MAGIC.len(), bytes.len()));
    check_contract(&mut reader, meter)?;
    let (modulus, nvars) = decode_ring(&mut reader, meter)?;
    let lengths = section_lengths(&mut reader, meter)?;
    let sections = section_ranges(bytes.len(), &reader, lengths)?;
    Ok(Header {
        modulus,
        nvars,
        sections,
    })
}

fn check_magic(bytes: &[u8], meter: &mut Meter) -> Result<(), VerifyError> {
    if bytes.len() < MAGIC.len() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: bytes.len(),
        });
    }
    meter.charge(MAGIC.len() as u64)?;
    if bytes[..MAGIC.len()] != MAGIC {
        if bytes[..6] == MAGIC[..6] {
            return Err(VerifyError::MalformedBinary {
                reason: BinaryFault::Version {
                    format: bytes[6],
                    revision: bytes[7],
                },
                offset: 6,
            });
        }
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Magic,
            offset: 0,
        });
    }
    Ok(())
}

fn check_contract(reader: &mut Reader<'_>, meter: &mut Meter) -> Result<(), VerifyError> {
    let schema = reader.string(meter)?;
    if schema != SCHEMA {
        return Err(VerifyError::Schema {
            found: String::from_utf8_lossy(&schema).into_owned(),
        });
    }
    let order = reader.string(meter)?;
    if order != ORDER {
        return Err(VerifyError::Order {
            found: String::from_utf8_lossy(&order).into_owned(),
        });
    }
    Ok(())
}

fn decode_ring(reader: &mut Reader<'_>, meter: &mut Meter) -> Result<(u64, usize), VerifyError> {
    let modulus = reader.varint(meter)?;
    if !(2..=MAX_MODULUS).contains(&modulus) {
        return Err(VerifyError::Modulus {
            found: modulus,
            composite: false,
        });
    }
    if !is_prime(modulus) {
        return Err(VerifyError::Modulus {
            found: modulus,
            composite: true,
        });
    }

    let nvars = reader.varint(meter)?;
    if nvars > MAX_NVARS {
        return Err(VerifyError::Nvars {
            found: nvars,
            max: MAX_NVARS,
        });
    }
    Ok((modulus, nvars as usize))
}

fn section_lengths(
    reader: &mut Reader<'_>,
    meter: &mut Meter,
) -> Result<[u64; SECTIONS], VerifyError> {
    let mut lengths = [0u64; SECTIONS];
    let mut total: u64 = 0;
    for slot in &mut lengths {
        *slot = reader.varint(meter)?;
        total = match total.checked_add(*slot) {
            Some(sum) => sum,
            None => {
                return Err(VerifyError::MalformedBinary {
                    reason: BinaryFault::Overflow {
                        what: "the sum of the section lengths",
                    },
                    offset: reader.pos,
                });
            }
        };
    }
    Ok(lengths)
}

fn section_ranges(
    bytes_len: usize,
    reader: &Reader<'_>,
    lengths: [u64; SECTIONS],
) -> Result<[Range; SECTIONS], VerifyError> {
    let total = lengths.iter().sum::<u64>();
    let start = reader.pos as u64;
    if start.checked_add(total) != Some(bytes_len as u64) {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::SectionLengths,
            offset: reader.pos,
        });
    }

    let mut sections = [(0usize, 0usize); SECTIONS];
    let mut at = reader.pos;
    for (range, length) in sections.iter_mut().zip(lengths) {
        *range = (at, at + length as usize);
        at = range.1;
    }
    Ok(sections)
}

/// The monomial pool, with one mark per entry.
pub(super) struct Pool {
    monos: Vec<Mono>,
    used: Vec<bool>,
}

impl Pool {
    /// Mark a pool entry as used and return its index.
    ///
    /// The charge covers the sparse support that a caller reads. The index
    /// stays cheap to retain when the caller can revisit the pool later.
    pub(super) fn reference(
        &mut self,
        index: u64,
        meter: &mut Meter,
    ) -> Result<usize, VerifyError> {
        let index = bounded(index, self.monos.len(), "pool index")?;
        let mono = &self.monos[index];
        meter.charge(1 + mono.support() as u64)?;
        self.used[index] = true;
        Ok(index)
    }

    /// Return a pool entry whose index was checked during decoding.
    pub(super) fn get(&self, index: usize) -> &Mono {
        &self.monos[index]
    }

    /// Return the monomial an index names, and mark the entry as used.
    ///
    /// The charge covers the copy: one unit per support entry, plus one.
    pub(super) fn take(&mut self, index: u64, meter: &mut Meter) -> Result<Mono, VerifyError> {
        let index = bounded(index, self.monos.len(), "pool index")?;
        let mono = &self.monos[index];
        meter.charge(1 + mono.support() as u64)?;
        self.used[index] = true;
        Ok(mono.clone())
    }

    /// Return an error if an entry has no reference.
    pub(super) fn check_referenced(&self, meter: &mut Meter) -> Result<(), VerifyError> {
        for (index, used) in self.used.iter().enumerate() {
            meter.charge(1)?;
            if !used {
                return Err(VerifyError::Pool {
                    index,
                    fault: PoolFault::Unreferenced,
                });
            }
        }
        Ok(())
    }
}

/// Decode the monomial pool.
pub(super) fn pool(reader: &mut Reader<'_>, ctx: &mut Ctx<'_>) -> Result<Pool, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(count, ctx.limits.max_pool_monomials, Cap::PoolMonomials)?;
    if count > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.pos,
        });
    }
    let mut monos: Vec<Mono> = Vec::with_capacity(count);
    let mut entries_total: usize = 0;
    for index in 0..count {
        monos.push(pool_monomial(
            reader,
            index,
            &mut entries_total,
            &monos,
            ctx,
        )?);
    }
    reader.finish()?;
    let used = vec![false; monos.len()];
    Ok(Pool { monos, used })
}

fn pool_monomial(
    reader: &mut Reader<'_>,
    index: usize,
    entries_total: &mut usize,
    previous: &[Mono],
    ctx: &mut Ctx<'_>,
) -> Result<Mono, VerifyError> {
    let support = pool_support(reader, index, entries_total, ctx)?;
    let entries = pool_entries(reader, index, support, ctx)?;
    let mono = Mono::from_support(entries);
    check_pool_order(index, previous.last(), &mono, ctx)?;
    Ok(mono)
}

fn pool_support(
    reader: &mut Reader<'_>,
    index: usize,
    entries_total: &mut usize,
    ctx: &mut Ctx<'_>,
) -> Result<usize, VerifyError> {
    let support = reader.varint(&mut ctx.meter)?;
    if support > ctx.nvars as u64 {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::SupportTooLarge {
                found: support,
                max: ctx.nvars,
            },
        });
    }
    let support = support as usize;
    ctx.meter.charge(1 + support as u64)?;
    *entries_total += support;
    if *entries_total > ctx.limits.max_pool_entries {
        return Err(VerifyError::CapExceeded {
            cap: Cap::PoolEntries,
            limit: ctx.limits.max_pool_entries,
        });
    }
    if support > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.pos,
        });
    }
    Ok(support)
}

fn pool_entries(
    reader: &mut Reader<'_>,
    index: usize,
    support: usize,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<(u32, Exp)>, VerifyError> {
    let mut entries = Vec::with_capacity(support);
    for _ in 0..support {
        entries.push(pool_entry(reader, index, entries.last().copied(), ctx)?);
    }
    Ok(entries)
}

fn pool_entry(
    reader: &mut Reader<'_>,
    index: usize,
    previous: Option<(u32, Exp)>,
    ctx: &mut Ctx<'_>,
) -> Result<(u32, Exp), VerifyError> {
    let variable = reader.varint(&mut ctx.meter)?;
    if variable >= ctx.nvars as u64 {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::VariableOutOfRange { variable },
        });
    }
    if previous.is_some_and(|(last, _)| variable <= last as u64) {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::VariableNotIncreasing,
        });
    }
    let exponent = reader.varint(&mut ctx.meter)?;
    if exponent == 0 {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::ExponentZero,
        });
    }
    if exponent > MAX_EXP {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::ExponentTooLarge { exponent },
        });
    }
    Ok((variable as u32, exponent as Exp))
}

fn check_pool_order(
    index: usize,
    previous: Option<&Mono>,
    mono: &Mono,
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let Some(previous) = previous else {
        return Ok(());
    };
    charge_monos(&mut ctx.meter, previous, mono)?;
    if previous.cmp_grevlex(mono) != Ordering::Less {
        return Err(VerifyError::Pool {
            index,
            fault: PoolFault::NotAscending,
        });
    }
    Ok(())
}

/// Decode the input section.
pub(super) fn input(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    ctx: &mut Ctx<'_>,
) -> Result<Vec<Poly>, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(count, ctx.limits.max_input_polys, Cap::InputPolynomials)?;
    if count > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.pos,
        });
    }
    let mut polys: Vec<Poly> = Vec::with_capacity(count);
    let mut terms_total: usize = 0;
    for index in 0..count {
        polys.push(input_poly(reader, pool, index, &mut terms_total, ctx)?);
    }
    reader.finish()?;
    Ok(polys)
}

fn input_poly(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    index: usize,
    terms_total: &mut usize,
    ctx: &mut Ctx<'_>,
) -> Result<Poly, VerifyError> {
    let term_count = input_term_count(reader, terms_total, ctx)?;
    let mut terms = Vec::with_capacity(term_count);
    let mut previous = None;
    for term in 0..term_count {
        let (entry, pool_index) =
            input_term(reader, pool, Location::Input(index), term, previous, ctx)?;
        terms.push(entry);
        previous = Some(pool_index);
    }
    ctx.meter.hold(term_count)?;
    Ok(Poly::new(terms))
}

fn input_term_count(
    reader: &mut Reader<'_>,
    terms_total: &mut usize,
    ctx: &mut Ctx<'_>,
) -> Result<usize, VerifyError> {
    let count = reader.varint(&mut ctx.meter)?;
    let count = capped(
        count,
        ctx.limits.max_terms_per_poly,
        Cap::TermsPerPolynomial,
    )?;
    *terms_total += count;
    if *terms_total > ctx.limits.max_total_terms {
        return Err(VerifyError::CapExceeded {
            cap: Cap::TotalTerms,
            limit: ctx.limits.max_total_terms,
        });
    }
    if count > reader.remaining() {
        return Err(VerifyError::MalformedBinary {
            reason: BinaryFault::Truncated,
            offset: reader.pos,
        });
    }
    Ok(count)
}

fn input_term(
    reader: &mut Reader<'_>,
    pool: &mut Pool,
    at: Location,
    term: usize,
    previous: Option<u64>,
    ctx: &mut Ctx<'_>,
) -> Result<(Term, u64), VerifyError> {
    let coeff = reader.varint(&mut ctx.meter)?;
    if coeff == 0 {
        return Err(VerifyError::Polynomial {
            at,
            fault: PolyFault::CoefficientZero { term },
        });
    }
    if coeff >= ctx.modulus {
        return Err(VerifyError::Polynomial {
            at,
            fault: PolyFault::CoefficientOutOfRange { term, coeff },
        });
    }
    let pool_index = reader.varint(&mut ctx.meter)?;
    if previous.is_some_and(|last| pool_index >= last) {
        return Err(VerifyError::Polynomial {
            at,
            fault: PolyFault::NotDescending { term },
        });
    }
    let mono = pool.take(pool_index, &mut ctx.meter)?;
    Ok((Term { coeff, mono }, pool_index))
}
