//! The `syl` format, which `benchmarks/gb-comparison/inputs` holds.
//!
//! ```text
//! 4
//! 1,0,0,0,1;1,0,0,1,0;1,0,1,0,0;1,1,0,0,0
//! -1,0,0,0,0;1,1,1,1,1
//! ```
//!
//! The first line is the variable count. Each further line is one
//! polynomial: terms separated by `;`, each term a coefficient and one
//! exponent per variable, separated by `,`. The format carries no variable
//! names, so a command that reads one names the variables `x1 .. xn` unless
//! another source names them.
//!
//! A term whose coefficient holds a `/` is a rational number. The `.sylq`
//! files of the harness are this same grammar under another extension, and
//! they carry no ring either. The generator writes integer coefficients
//! only, so reading is a superset of what it writes.

use std::io::{self, Write};

use super::{Body, Header, Reading, Rendered, Term};

/// The largest variable count the reader accepts.
///
/// It is the count a ring accepts (`RingError::TooManyVariables`). The
/// number comes from the file, so it is checked before the reader
/// reserves anything for it.
const MAX_VARIABLES: usize = 256;

/// Read a `syl` source.
pub fn read(origin: &str, text: &str) -> Result<Reading, String> {
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty());
    let Some((index, count)) = lines.next() else {
        return Err(format!("{origin} is empty"));
    };
    let count = count.trim();
    let nvars = count.parse::<usize>().map_err(|_| {
        format!(
            "{origin} line {}: \"{count}\" is not a variable count",
            index + 1
        )
    })?;
    if nvars == 0 {
        return Err(format!(
            "{origin} line {}: the ring has no variable",
            index + 1
        ));
    }
    if nvars > MAX_VARIABLES {
        return Err(format!(
            "{origin} line {}: the file names {nvars} variables, the largest supported count is {MAX_VARIABLES}",
            index + 1
        ));
    }
    let mut polynomials = Vec::new();
    for (index, line) in lines {
        let at = format!("{origin} line {}", index + 1);
        let mut terms = Vec::new();
        for field in line.trim().split(';') {
            let field = field.trim();
            if field.is_empty() {
                continue;
            }
            terms.push(term(&at, field, nvars)?);
        }
        polynomials.push(terms);
    }
    Ok(Reading {
        names: None,
        domain: None,
        nvars: Some(nvars),
        body: Body::Terms(polynomials),
    })
}

fn term(at: &str, field: &str, nvars: usize) -> Result<Term, String> {
    let mut parts = field.split(',');
    let coefficient = parts.next().unwrap_or("").trim();
    if !is_coefficient(coefficient) {
        return Err(format!("{at}: \"{coefficient}\" is not a coefficient"));
    }
    let mut exponents = Vec::with_capacity(nvars);
    for part in parts {
        let part = part.trim();
        exponents.push(
            part.parse::<u16>()
                .map_err(|_| format!("{at}: \"{part}\" is not an exponent"))?,
        );
    }
    if exponents.len() != nvars {
        return Err(format!(
            "{at}: a term holds {} exponents, the file names {nvars} variables",
            exponents.len()
        ));
    }
    Ok(Term {
        coefficient: coefficient.to_string(),
        exponents,
    })
}

/// Report whether the text is a coefficient: digits, an optional leading
/// `-`, and an optional `/` denominator.
fn is_coefficient(text: &str) -> bool {
    let text = text.strip_prefix('-').unwrap_or(text);
    let (numerator, denominator) = match text.split_once('/') {
        Some((numerator, denominator)) => (numerator, Some(denominator)),
        None => (text, None),
    };
    let digits = |part: &str| !part.is_empty() && part.bytes().all(|b| b.is_ascii_digit());
    digits(numerator) && denominator.is_none_or(digits)
}

/// Write a `syl` source.
pub fn write(header: &Header<'_>, polynomials: &[Rendered], out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "{}", header.names.len())?;
    for polynomial in polynomials {
        let terms: Vec<String> = polynomial
            .terms
            .iter()
            .map(|(coefficient, exponents)| {
                let mut term = coefficient.clone();
                for exponent in exponents {
                    term.push(',');
                    term.push_str(&exponent.to_string());
                }
                term
            })
            .collect();
        writeln!(out, "{}", terms.join(";"))?;
    }
    Ok(())
}
