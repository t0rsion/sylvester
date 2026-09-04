//! The `ms` format, which msolve reads.
//!
//! ```text
//! x1,x2,x3
//! 32003
//! x1+x2+x3,
//! x1*x2+x2*x3+x1*x3,
//! x1*x2*x3-1
//! ```
//!
//! The first line names the variables, the second the characteristic, and
//! `0` there is the rational numbers. Commas separate the polynomials, and
//! one polynomial may span several lines.

use std::io::{self, Write};

use super::{Body, DomainClaim, Header, Reading, Rendered, compact_expression};

/// Read an `ms` source.
pub fn read(origin: &str, text: &str) -> Result<Reading, String> {
    let mut lines = text
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty());
    let Some((_, variables)) = lines.next() else {
        return Err(format!("{origin} is empty"));
    };
    let names: Vec<String> = variables
        .split(',')
        .map(|name| name.trim().to_string())
        .filter(|name| !name.is_empty())
        .collect();
    if names.is_empty() {
        return Err(format!("{origin} line 1: the variable list is empty"));
    }
    let Some((index, characteristic)) = lines.next() else {
        return Err(format!("{origin} names no characteristic"));
    };
    let characteristic = characteristic.trim();
    let domain = match characteristic.parse::<u64>() {
        Ok(0) => DomainClaim::Rationals,
        Ok(modulus) => DomainClaim::Prime(modulus),
        Err(_) => {
            return Err(format!(
                "{origin} line {}: \"{characteristic}\" is not a characteristic",
                index + 1
            ));
        }
    };
    let body: String = lines.map(|(_, line)| line).collect::<Vec<_>>().join("\n");
    let expressions: Vec<String> = body
        .split(',')
        .map(str::trim)
        .filter(|piece| !piece.is_empty())
        .map(str::to_string)
        .collect();
    let nvars = names.len();
    Ok(Reading {
        names: Some(names),
        domain: Some(domain),
        nvars: Some(nvars),
        body: Body::Expressions(expressions),
    })
}

/// Write an `ms` source.
pub fn write(header: &Header<'_>, polynomials: &[Rendered], out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "{}", header.names.join(","))?;
    match header.domain {
        DomainClaim::Prime(modulus) => writeln!(out, "{modulus}")?,
        DomainClaim::Rationals => writeln!(out, "0")?,
    }
    for (index, polynomial) in polynomials.iter().enumerate() {
        let separator = if index + 1 == polynomials.len() {
            ""
        } else {
            ","
        };
        writeln!(
            out,
            "{}{separator}",
            compact_expression(polynomial, header.names)
        )?;
    }
    Ok(())
}
