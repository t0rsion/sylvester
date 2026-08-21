//! Strict decoder for the canonical certificate encoding.
//!
//! The decoder accepts one shape only: the object the contract lists, with
//! its keys in the contract order and no insignificant whitespace. It reads
//! bytes and returns typed values. It checks the resource caps while it
//! reads, so a rejected certificate never allocates past a cap. Every loop
//! that runs over the bytes polls the deadline every `WORK_STRIDE` steps,
//! so no single field can be read past it.

use super::algebra::{Mono, Poly, Term};
use super::error::{Cap, Syntax, VerifyError};
use super::limits::{Limits, Ticker};

/// The largest exponent the contract allows.
const MAX_EXP: u64 = 65535;

/// The largest variable count the contract allows.
///
/// The decoder checks `nvars` against this bound as soon as it reads the
/// field, and it caps every term's exponent array at the same width, so
/// neither an out-of-range `nvars` nor an oversized exponent array reaches
/// an allocation past it.
pub(super) const MAX_NVARS: u64 = 256;

const KEYS: [&str; 9] = [
    "schema",
    "order",
    "modulus",
    "nvars",
    "input",
    "basis",
    "origin",
    "membership",
    "spairs",
];

/// One S-pair entry, with indices still unchecked.
#[derive(Clone, Debug)]
pub(crate) struct RawSpair {
    pub i: u64,
    pub j: u64,
    pub cofactors: Vec<Poly>,
}

/// The certificate as decoded. Values carry no meaning yet.
#[derive(Clone, Debug)]
pub(crate) struct RawCert {
    pub schema: String,
    pub order: String,
    pub modulus: u64,
    pub nvars: u64,
    pub input: Vec<Poly>,
    pub basis: Vec<Poly>,
    pub origin: Vec<Vec<Poly>>,
    pub membership: Vec<Vec<Poly>>,
    pub spairs: Vec<RawSpair>,
}

type Decoded<T> = Result<T, VerifyError>;

/// Decode certificate bytes.
pub(crate) fn decode(bytes: &[u8], limits: &Limits) -> Decoded<RawCert> {
    if bytes.len() > limits.max_bytes {
        return Err(VerifyError::CapExceeded {
            cap: Cap::Bytes,
            limit: limits.max_bytes,
        });
    }
    let mut decoder = Decoder {
        bytes,
        pos: 0,
        limits,
        ticker: limits.ticker(),
        polys: 0,
        terms: 0,
        entries: 0,
    };
    decoder.certificate()
}

struct Decoder<'a> {
    bytes: &'a [u8],
    pos: usize,
    limits: &'a Limits,
    ticker: Ticker,
    polys: usize,
    terms: usize,
    entries: usize,
}

fn is_whitespace(byte: u8) -> bool {
    matches!(byte, b' ' | b'\t' | b'\n' | b'\r')
}

impl Decoder<'_> {
    fn fail<T>(&self, reason: Syntax) -> Decoded<T> {
        Err(VerifyError::Malformed {
            reason,
            offset: self.pos,
        })
    }

    fn fail_at<T>(&self, reason: Syntax, offset: usize) -> Decoded<T> {
        Err(VerifyError::Malformed { reason, offset })
    }

    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    fn expect(&mut self, byte: u8, expected: &'static str) -> Decoded<()> {
        match self.peek() {
            Some(found) if found == byte => {
                self.pos += 1;
                Ok(())
            }
            Some(found) if is_whitespace(found) => self.fail(Syntax::Whitespace),
            Some(found) => self.fail(Syntax::UnexpectedByte { found, expected }),
            None => self.fail(Syntax::Truncated),
        }
    }

    fn string(&mut self) -> Decoded<String> {
        self.expect(b'"', "a string")?;
        let start = self.pos;
        loop {
            self.ticker.step()?;
            match self.peek() {
                Some(b'"') => {
                    let text = String::from_utf8_lossy(&self.bytes[start..self.pos]).into_owned();
                    self.pos += 1;
                    return Ok(text);
                }
                // The contract needs no escape and no byte outside printable
                // ASCII, so the decoder accepts neither.
                Some(byte) if (0x20..0x7f).contains(&byte) && byte != b'\\' => self.pos += 1,
                Some(_) => return self.fail(Syntax::BadString),
                None => return self.fail(Syntax::Truncated),
            }
        }
    }

    fn uint(&mut self, what: &'static str, max: u64) -> Decoded<u64> {
        match self.peek() {
            Some(b'-' | b'+') => return self.fail(Syntax::Sign),
            Some(byte) if byte.is_ascii_digit() => {}
            Some(byte) if is_whitespace(byte) => return self.fail(Syntax::Whitespace),
            Some(found) => {
                return self.fail(Syntax::UnexpectedByte {
                    found,
                    expected: "an integer",
                });
            }
            None => return self.fail(Syntax::Truncated),
        }
        let start = self.pos;
        let mut value: u64 = 0;
        let mut overflow = false;
        while let Some(byte) = self.peek() {
            if !byte.is_ascii_digit() {
                break;
            }
            self.ticker.step()?;
            value = match value
                .checked_mul(10)
                .and_then(|v| v.checked_add((byte - b'0') as u64))
            {
                Some(next) => next,
                None => {
                    overflow = true;
                    0
                }
            };
            self.pos += 1;
        }
        if self.bytes[start] == b'0' && self.pos - start > 1 {
            return self.fail_at(Syntax::LeadingZero, start);
        }
        match self.peek() {
            Some(b'.') => return self.fail(Syntax::Float),
            Some(b'e' | b'E') => return self.fail(Syntax::Exponent),
            _ => {}
        }
        if overflow {
            return self.fail_at(Syntax::IntegerOverflow { what }, start);
        }
        if value > max {
            return self.fail_at(Syntax::IntegerOutOfRange { what, value, max }, start);
        }
        Ok(value)
    }

    fn array<T>(&mut self, mut item: impl FnMut(&mut Self) -> Decoded<T>) -> Decoded<Vec<T>> {
        self.expect(b'[', "an array")?;
        let mut out = Vec::new();
        if self.peek() == Some(b']') {
            self.pos += 1;
            return Ok(out);
        }
        loop {
            self.ticker.step()?;
            out.push(item(self)?);
            match self.peek() {
                Some(b',') => self.pos += 1,
                Some(b']') => {
                    self.pos += 1;
                    return Ok(out);
                }
                Some(found) if is_whitespace(found) => return self.fail(Syntax::Whitespace),
                Some(found) => {
                    return self.fail(Syntax::UnexpectedByte {
                        found,
                        expected: "a comma or the end of an array",
                    });
                }
                None => return self.fail(Syntax::Truncated),
            }
        }
    }

    /// Read one term's exponent array, capped at the contract's variable
    /// width.
    ///
    /// The width cap applies as the array is read, so an array past the width
    /// fails before the decoder allocates the entry past it. The per-variable
    /// check against the certificate's own `nvars` runs later, once the whole
    /// polynomial is decoded.
    fn exps(&mut self) -> Decoded<Vec<u32>> {
        let mut count = 0u64;
        self.array(|d| {
            if count >= MAX_NVARS {
                return d.fail(Syntax::ExponentWidth { max: MAX_NVARS });
            }
            count += 1;
            d.uint("an exponent", MAX_EXP).map(|e| e as u32)
        })
    }

    fn term(&mut self) -> Decoded<Term> {
        self.expect(b'[', "a term")?;
        let coeff = self.uint("a coefficient", u64::MAX)?;
        self.expect(b',', "a comma")?;
        let exps = self.exps()?;
        self.expect(b']', "the end of a term")?;
        Ok(Term::new(coeff, Mono::new(exps)))
    }

    fn poly(&mut self) -> Decoded<Poly> {
        self.polys += 1;
        if self.polys > self.limits.max_polys {
            return Err(VerifyError::CapExceeded {
                cap: Cap::Polynomials,
                limit: self.limits.max_polys,
            });
        }
        self.ticker.step()?;
        let mut count = 0usize;
        let terms = self.array(|d| {
            count += 1;
            if count > d.limits.max_terms_per_poly {
                return Err(VerifyError::CapExceeded {
                    cap: Cap::TermsPerPolynomial,
                    limit: d.limits.max_terms_per_poly,
                });
            }
            d.terms += 1;
            if d.terms > d.limits.max_total_terms {
                return Err(VerifyError::CapExceeded {
                    cap: Cap::TotalTerms,
                    limit: d.limits.max_total_terms,
                });
            }
            d.term()
        })?;
        Ok(Poly::new(terms))
    }

    fn poly_array(&mut self) -> Decoded<Vec<Poly>> {
        self.array(|d| d.poly())
    }

    /// Count one entry of the `origin`, `membership`, or `spairs` array.
    ///
    /// An entry holds a list of cofactors. An entry with no cofactor still
    /// costs memory, so the cap covers it.
    fn entry(&mut self) -> Decoded<()> {
        self.entries += 1;
        if self.entries > self.limits.max_entries {
            return Err(VerifyError::CapExceeded {
                cap: Cap::Entries,
                limit: self.limits.max_entries,
            });
        }
        self.ticker.step()
    }

    fn spair(&mut self) -> Decoded<RawSpair> {
        self.expect(b'[', "an S-pair entry")?;
        let i = self.uint("an S-pair index", u64::MAX)?;
        self.expect(b',', "a comma")?;
        let j = self.uint("an S-pair index", u64::MAX)?;
        self.expect(b',', "a comma")?;
        let cofactors = self.poly_array()?;
        self.expect(b']', "the end of an S-pair entry")?;
        Ok(RawSpair { i, j, cofactors })
    }

    /// Read the key at `slot` and the colon after it.
    ///
    /// A key from a later slot and a key that is absent look the same at
    /// this position. The decoder reports the later key it found.
    fn key(&mut self, slot: usize) -> Decoded<()> {
        if self.peek() == Some(b'}') {
            return self.fail(Syntax::MissingKey(KEYS[slot]));
        }
        let start = self.pos;
        let found = self.string()?;
        match KEYS.iter().position(|key| *key == found) {
            None => return self.fail_at(Syntax::UnknownKey(found), start),
            Some(index) if index < slot => {
                return self.fail_at(Syntax::DuplicateKey(found), start);
            }
            Some(index) if index > slot => {
                return self.fail_at(
                    Syntax::KeyOutOfOrder {
                        found,
                        expected: KEYS[slot],
                    },
                    start,
                );
            }
            Some(_) => {}
        }
        self.expect(b':', "a colon")
    }

    fn separator(&mut self, next_slot: usize) -> Decoded<()> {
        if self.peek() == Some(b'}') {
            return self.fail(Syntax::MissingKey(KEYS[next_slot]));
        }
        self.expect(b',', "a comma")
    }

    fn certificate(&mut self) -> Decoded<RawCert> {
        self.expect(b'{', "an object")?;
        self.key(0)?;
        let schema = self.string()?;
        self.separator(1)?;
        self.key(1)?;
        let order = self.string()?;
        self.separator(2)?;
        self.key(2)?;
        let modulus = self.uint("the modulus", u64::MAX)?;
        self.separator(3)?;
        self.key(3)?;
        let nvars = self.uint("nvars", u64::MAX)?;
        // The width check runs here, before any polynomial is decoded,
        // rather than waiting for `checks::check` to run it once the whole
        // certificate, however large, is already in memory.
        if nvars > MAX_NVARS {
            return Err(VerifyError::Nvars {
                found: nvars,
                max: MAX_NVARS,
            });
        }
        self.separator(4)?;
        self.key(4)?;
        let input = self.poly_array()?;
        self.separator(5)?;
        self.key(5)?;
        let basis = self.poly_array()?;
        self.separator(6)?;
        self.key(6)?;
        let origin = self.array(|d| {
            d.entry()?;
            d.poly_array()
        })?;
        self.separator(7)?;
        self.key(7)?;
        let membership = self.array(|d| {
            d.entry()?;
            d.poly_array()
        })?;
        self.separator(8)?;
        self.key(8)?;
        let spairs = self.array(|d| {
            d.entry()?;
            d.spair()
        })?;
        self.expect(b'}', "the end of the object")?;
        if self.pos != self.bytes.len() {
            return self.fail(Syntax::TrailingBytes);
        }
        Ok(RawCert {
            schema,
            order,
            modulus,
            nvars,
            input,
            basis,
            origin,
            membership,
            spairs,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SMALL: &str = concat!(
        r#"{"schema":"sylv-gb-cert-v1","order":"grevlex-v1","modulus":7,"nvars":2,"#,
        r#""input":[[[1,[1,0]]]],"basis":[[[1,[1,0]]]],"origin":[[[[1,[0,0]]]]],"#,
        r#""membership":[[[[1,[0,0]]]]],"spairs":[]}"#
    );

    fn reason(bytes: &[u8]) -> Syntax {
        match decode(bytes, &Limits::default()) {
            Err(VerifyError::Malformed { reason, .. }) => reason,
            other => panic!("expected a malformed certificate, got {other:?}"),
        }
    }

    #[test]
    fn decodes_a_canonical_object() {
        let cert = decode(SMALL.as_bytes(), &Limits::default()).expect("decode");
        assert_eq!(cert.schema, "sylv-gb-cert-v1");
        assert_eq!(cert.modulus, 7);
        assert_eq!(cert.nvars, 2);
        assert_eq!(cert.input.len(), 1);
        assert_eq!(cert.input[0].terms().len(), 1);
        assert_eq!(cert.input[0].terms()[0].coeff(), 1);
        assert_eq!(cert.input[0].terms()[0].mono().exps(), &[1, 0]);
        assert!(cert.spairs.is_empty());
    }

    #[test]
    fn rejects_whitespace() {
        assert_eq!(reason(b"{ \"schema\":\"x\"}"), Syntax::Whitespace);
    }

    #[test]
    fn rejects_a_leading_zero() {
        let bytes = SMALL.replace("\"modulus\":7", "\"modulus\":07");
        assert_eq!(reason(bytes.as_bytes()), Syntax::LeadingZero);
    }

    #[test]
    fn rejects_a_float() {
        let bytes = SMALL.replace("\"modulus\":7", "\"modulus\":7.5");
        assert_eq!(reason(bytes.as_bytes()), Syntax::Float);
    }

    #[test]
    fn rejects_an_exponent_part() {
        let bytes = SMALL.replace("\"modulus\":7", "\"modulus\":7e2");
        assert_eq!(reason(bytes.as_bytes()), Syntax::Exponent);
    }

    #[test]
    fn rejects_a_sign() {
        let bytes = SMALL.replace("\"modulus\":7", "\"modulus\":-7");
        assert_eq!(reason(bytes.as_bytes()), Syntax::Sign);
    }

    #[test]
    fn rejects_an_exponent_above_the_contract_range() {
        let bytes = SMALL.replace("[1,[1,0]]", "[1,[65536,0]]");
        assert_eq!(
            reason(bytes.as_bytes()),
            Syntax::IntegerOutOfRange {
                what: "an exponent",
                value: 65536,
                max: 65535
            }
        );
    }

    #[test]
    fn rejects_an_integer_above_64_bits() {
        let bytes = SMALL.replace("\"modulus\":7", "\"modulus\":99999999999999999999999");
        assert_eq!(
            reason(bytes.as_bytes()),
            Syntax::IntegerOverflow {
                what: "the modulus"
            }
        );
    }

    #[test]
    fn rejects_empty_bytes() {
        assert_eq!(reason(b""), Syntax::Truncated);
    }

    #[test]
    fn reports_the_offset_of_the_fault() {
        match decode(b"{\"schema\": \"x\"}", &Limits::default()) {
            Err(VerifyError::Malformed { offset, .. }) => assert_eq!(offset, 10),
            other => panic!("expected a malformed certificate, got {other:?}"),
        }
    }
}
