//! The `text` format: the crate's own syntax with a ring header.
//!
//! ```text
//! # vars: x, y, z
//! # modulus: 32003
//! x^2*y - 3*z + 1
//! x*z - y
//! ```
//!
//! `# coefficients: rationals` takes the place of `# modulus:` over `Q`.
//! Any other `#` line is a comment. One polynomial goes on one line, in the
//! syntax `parse_polynomial` reads and `Display` writes.

use std::io::{self, Write};

use super::{Body, DomainClaim, Header, Reading, Rendered, expression};

/// Read a `text` source.
pub fn read(origin: &str, text: &str) -> Result<Reading, String> {
    let mut names: Option<Vec<String>> = None;
    let mut modulus: Option<u64> = None;
    let mut rationals = false;
    let mut expressions = Vec::new();
    for (index, line) in text.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let Some(rest) = line.strip_prefix('#') else {
            expressions.push(line.to_string());
            continue;
        };
        let Some((key, value)) = rest.split_once(':') else {
            continue;
        };
        let value = value.trim();
        let at = format!("{origin} line {}", index + 1);
        match key.trim() {
            "vars" => {
                if names.is_some() {
                    return Err(format!("{at}: the source names its variables twice"));
                }
                let read: Vec<String> = value
                    .split([',', ' ', '\t'])
                    .filter(|name| !name.is_empty())
                    .map(str::to_string)
                    .collect();
                if read.is_empty() {
                    return Err(format!("{at}: the variable list is empty"));
                }
                names = Some(read);
            }
            "modulus" => {
                if modulus.is_some() {
                    return Err(format!("{at}: the source names its modulus twice"));
                }
                modulus = Some(
                    value
                        .parse::<u64>()
                        .map_err(|_| format!("{at}: \"{value}\" is not a modulus"))?,
                );
            }
            "coefficients" => {
                if value != "rationals" {
                    return Err(format!(
                        "{at}: \"{value}\" is not a coefficient domain, the only one this line names is rationals"
                    ));
                }
                rationals = true;
            }
            _ => {}
        }
    }
    let domain = match (modulus, rationals) {
        (Some(_), true) => {
            return Err(format!(
                "{origin} names both a modulus and rational coefficients"
            ));
        }
        (Some(modulus), false) => Some(DomainClaim::Prime(modulus)),
        (None, true) => Some(DomainClaim::Rationals),
        (None, false) => None,
    };
    let nvars = names.as_ref().map(Vec::len);
    Ok(Reading {
        names,
        domain,
        nvars,
        body: Body::Expressions(expressions),
        record: None,
    })
}

/// Write a `text` source.
pub fn write(header: &Header<'_>, polynomials: &[Rendered], out: &mut dyn Write) -> io::Result<()> {
    writeln!(out, "# vars: {}", header.names.join(", "))?;
    match header.domain {
        DomainClaim::Prime(modulus) => writeln!(out, "# modulus: {modulus}")?,
        DomainClaim::Rationals => writeln!(out, "# coefficients: rationals")?,
    }
    for polynomial in polynomials {
        writeln!(out, "{}", expression(polynomial, header.names))?;
    }
    Ok(())
}
