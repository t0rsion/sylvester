//! Independent verifier for the `sylv-gb-cert-v2` certificate contract.
//!
//! The contract is `docs/certificate-v2.md`. v2 belongs to the F4 backend:
//! the certificate carries a trace of the operations the engine performed,
//! and the verifier replays that trace instead of reading explicit cofactor
//! polynomials.
//!
//! The module shares no code with the engines and none with the v1
//! verifier. It carries its own binary decoder, sparse monomials, field
//! arithmetic, and grevlex comparison. It shares the [`VerifiedGb`] type
//! and the [`Limits`] type, and nothing else.
//!
//! Acceptance proves four facts about the modulus, the input list, and the
//! basis the certificate carries: every basis element lies in the ideal of
//! the input, every input polynomial lies in the ideal of the basis, every
//! pair of basis elements has an lcm representation, and the basis is
//! monic, sorted, minimal, and interreduced. The reduced Gröbner basis of
//! an ideal is unique, so the four facts together say the basis is that
//! reduced Gröbner basis. Acceptance says nothing about the producer of
//! the bytes.
//!
//! The verifier returns a value for every byte string. It checks every
//! decoded index against the collection it addresses, and it adds decoded
//! integers with checked arithmetic, so a malformed certificate is a typed
//! error. The tamper corpus in `tests/verifier_v2.rs` mutates every byte of
//! nine certificates and cuts each one at every offset. It reaches no
//! panic.
//!
//! Each cap bounds one structure the verifier allocates, as section 8.4 of
//! the contract lists. No single cap bounds the memory of the whole check,
//! and no cap bounds wall time. Bytes from an untrusted peer also need a
//! deadline.

mod arith;
mod decode;
mod replay;
mod witness;

use crate::verify::VerifiedGb;
use crate::verify::error::{Cap, VerifyError};
use crate::verify::limits::{Limits, Meter};
use arith::Poly;

pub(crate) use decode::MAGIC_FIRST;

/// The state every section check shares.
pub(super) struct Ctx<'a> {
    limits: &'a Limits,
    meter: Meter,
    modulus: u64,
    nvars: usize,
}

/// A polynomial as terms of a coefficient and a dense exponent vector.
pub(crate) type DensePoly = Vec<(u64, Vec<u32>)>;

/// The values a v2 certificate carries, after every check held.
pub(crate) struct Accepted {
    pub(crate) modulus: u64,
    pub(crate) nvars: usize,
    pub(crate) input: Vec<DensePoly>,
    pub(crate) basis: Vec<DensePoly>,
}

/// Verify v2 certificate bytes under the default caps.
pub fn verify(bytes: &[u8]) -> Result<VerifiedGb, VerifyError> {
    verify_with_limits(bytes, &Limits::default())
}

/// Verify v2 certificate bytes under caller-supplied caps.
///
/// The caps stop the work before the verifier allocates past them. An
/// exhausted verifier reports [`VerifyError::CapExceeded`] or
/// [`VerifyError::DeadlineExceeded`]. Neither value says the certificate is
/// invalid.
pub fn verify_with_limits(bytes: &[u8], limits: &Limits) -> Result<VerifiedGb, VerifyError> {
    let outcome = check(bytes, limits);
    // A passed deadline outranks acceptance and invalidity: work past the
    // deadline is work the caller did not request. A cap already reported
    // remains a cap because the verifier stopped at that boundary.
    if !matches!(&outcome, Err(error) if error.is_exhaustion()) {
        limits.check_deadline()?;
    }
    crate::verify::from_v2(outcome?, limits)
}

/// Run the nine obligations of section 7, in order.
pub(crate) fn check(bytes: &[u8], limits: &Limits) -> Result<Accepted, VerifyError> {
    if bytes.len() > limits.max_bytes {
        return Err(VerifyError::CapExceeded {
            cap: Cap::Bytes,
            limit: limits.max_bytes,
        });
    }
    let mut meter = Meter::new(limits);
    let header = decode::header(bytes, &mut meter)?;
    meter.set_nvars(header.nvars);
    let mut ctx = Ctx {
        limits,
        meter,
        modulus: header.modulus,
        nvars: header.nvars,
    };
    ctx.meter.poll()?;

    let (mut pool, input) = decode_input(bytes, &header.sections, &mut ctx)?;
    let basis = replay_basis(bytes, &header.sections, &mut pool, &input, &mut ctx)?;
    check_witnesses(bytes, &header.sections, &mut pool, &input, &basis, &mut ctx)?;

    Ok(Accepted {
        modulus: ctx.modulus,
        nvars: ctx.nvars,
        input: dense(input, ctx.nvars),
        basis: dense(basis, ctx.nvars),
    })
}

fn decode_input(
    bytes: &[u8],
    sections: &[decode::Range; decode::SECTIONS],
    ctx: &mut Ctx<'_>,
) -> Result<(decode::Pool, Vec<Poly>), VerifyError> {
    let mut pool = decode::pool(&mut decode::Reader::new(bytes, sections[0]), ctx)?;
    ctx.meter.poll()?;

    let input = decode::input(&mut decode::Reader::new(bytes, sections[1]), &mut pool, ctx)?;
    ctx.meter.poll()?;
    Ok((pool, input))
}

fn replay_basis(
    bytes: &[u8],
    sections: &[decode::Range; decode::SECTIONS],
    pool: &mut decode::Pool,
    input: &[Poly],
    ctx: &mut Ctx<'_>,
) -> Result<Vec<Poly>, VerifyError> {
    let mut nodes = replay::trace(
        &mut decode::Reader::new(bytes, sections[2]),
        pool,
        input,
        ctx,
    )?;
    ctx.meter.poll()?;

    let basis = replay::basis(
        &mut decode::Reader::new(bytes, sections[3]),
        &mut nodes,
        ctx,
    )?;
    ctx.meter.poll()?;
    Ok(basis)
}

fn check_witnesses(
    bytes: &[u8],
    sections: &[decode::Range; decode::SECTIONS],
    pool: &mut decode::Pool,
    input: &[Poly],
    basis: &[Poly],
    ctx: &mut Ctx<'_>,
) -> Result<(), VerifyError> {
    let mut steps = witness::Steps::new();
    witness::membership(
        &mut decode::Reader::new(bytes, sections[4]),
        input,
        basis,
        pool,
        &mut steps,
        ctx,
    )?;
    ctx.meter.poll()?;

    witness::pairs(
        &mut decode::Reader::new(bytes, sections[5]),
        basis,
        pool,
        &mut steps,
        ctx,
    )?;
    ctx.meter.poll()?;

    pool.check_referenced(&mut ctx.meter)?;
    ctx.meter.poll()?;
    Ok(())
}

fn dense(polys: Vec<Poly>, nvars: usize) -> Vec<DensePoly> {
    polys
        .into_iter()
        .map(|poly| poly.into_dense(nvars))
        .collect()
}
