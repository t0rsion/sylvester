//! The input formats and the writers that produce them.
//!
//! A reader turns bytes into a [`Reading`]: what the source says about the
//! ring, and the polynomials it holds. A polynomial arrives either as an
//! expression in the syntax `sylvester::PolynomialRing::parse_polynomial`
//! reads, or as the term lists the `syl` format writes. Term lists become
//! expressions once the command has resolved its ring, because the names
//! come from that ring and not from the file.
//!
//! A writer takes [`Rendered`] polynomials, which hold coefficient text and
//! exponents. The same rendering covers both coefficient domains and the
//! polynomials a certificate carries.

mod json;
mod ms;
mod syl;
mod text;

use std::borrow::Cow;
use std::fmt;
use std::io::{self, Write};
use std::path::Path;

use clap::ValueEnum;
use sylvester::{Budget, Domain, EnvelopeError, Polynomial, ResultEnvelope, verify};

/// A format `sylv` reads. Writers support formats whose output schema fits.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum Format {
    /// The msolve input format.
    ///
    /// A line of comma-separated variable names, a line with the
    /// characteristic (`0` for `Q`), then comma-separated polynomials.
    Ms,
    /// The benchmark format, under the `syl` and `sylq` extensions.
    ///
    /// A line with the variable count, then one polynomial per line as
    /// `coefficient,exponents` terms separated by `;`. It carries no ring.
    Syl,
    /// The crate's own syntax.
    ///
    /// A `# vars:` line, a `# modulus:` or `# coefficients: rationals`
    /// line, then one polynomial per line.
    Text,
    /// A `sylv-result-v1` computation record.
    Json,
}

/// Why a source could not be decoded.
#[derive(Debug)]
pub enum ReadError {
    /// The source does not follow its selected text format.
    Format(String),
    /// The source is a computation record with an invalid or exhausted read.
    Envelope(EnvelopeError),
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Format(message) => f.write_str(message),
            Self::Envelope(error) => error.fmt(f),
        }
    }
}

impl Format {
    /// The format a file name extension names.
    ///
    /// `ms`, `syl`, `sylq`, `txt`, `text`, and `json` are the extensions. Any other
    /// one is `None`. The caller reports that as a usage error.
    pub fn of_path(path: &Path) -> Option<Format> {
        match path.extension()?.to_str()? {
            "ms" => Some(Format::Ms),
            "syl" | "sylq" => Some(Format::Syl),
            "txt" | "text" => Some(Format::Text),
            "json" => Some(Format::Json),
            _ => None,
        }
    }
}

impl fmt::Display for Format {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Format::Ms => "ms",
            Format::Syl => "syl",
            Format::Text => "text",
            Format::Json => "json",
        })
    }
}

/// The coefficient domain a source names.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DomainClaim {
    /// The prime field of this modulus.
    Prime(u64),
    /// The rational numbers.
    Rationals,
}

impl fmt::Display for DomainClaim {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            DomainClaim::Prime(modulus) => write!(f, "the prime field of {modulus}"),
            DomainClaim::Rationals => f.write_str("the rational numbers"),
        }
    }
}

/// What one source holds.
///
/// `names` and `domain` are the ring the source names. The `ms` and `text`
/// formats name both, or neither. The `syl` format names no ring and
/// reports the variable count in `nvars`.
#[derive(Debug)]
pub struct Reading {
    /// The variable names, largest variable first.
    pub names: Option<Vec<String>>,
    /// The coefficient domain.
    pub domain: Option<DomainClaim>,
    /// The variable count, when the source states one.
    pub nvars: Option<usize>,
    /// The polynomials.
    pub body: Body,
    /// The untrusted computation record, when the source is JSON.
    pub record: Option<ResultEnvelope>,
}

/// The polynomials of a source, in the shape the format writes them.
#[derive(Debug)]
pub enum Body {
    /// Expressions in the syntax `parse_polynomial` reads.
    Expressions(Vec<String>),
    /// Term lists, which need the variable names to become expressions.
    Terms(Vec<Vec<Term>>),
}

/// One term of a `syl` polynomial.
#[derive(Debug)]
pub struct Term {
    /// The coefficient, as the file writes it: digits, an optional leading
    /// `-`, and an optional `/` denominator.
    pub coefficient: String,
    /// One exponent per variable.
    pub exponents: Vec<u16>,
}

impl Body {
    /// Return one polynomial as an expression without cloning expression bodies.
    pub fn expression(
        &self,
        origin: &str,
        index: usize,
        names: &[String],
    ) -> Result<Cow<'_, str>, String> {
        match self {
            Body::Expressions(expressions) => expressions
                .get(index)
                .map(|expression| Cow::Borrowed(expression.as_str()))
                .ok_or_else(|| format!("{origin} has no polynomial {}", index + 1)),
            Body::Terms(polynomials) => polynomials
                .get(index)
                .ok_or_else(|| format!("{origin} has no polynomial {}", index + 1))
                .and_then(|terms| expression_of_terms(origin, index, terms, names))
                .map(Cow::Owned),
        }
    }

    /// The polynomials as expressions over `names`.
    ///
    /// A term list whose width is not the number of names is an error.
    #[cfg(test)]
    pub fn expressions(&self, origin: &str, names: &[String]) -> Result<Vec<String>, String> {
        match self {
            Body::Expressions(expressions) => Ok(expressions.clone()),
            Body::Terms(polynomials) => polynomials
                .iter()
                .enumerate()
                .map(|(index, terms)| expression_of_terms(origin, index, terms, names))
                .collect(),
        }
    }

    /// The number of polynomials.
    pub fn len(&self) -> usize {
        match self {
            Body::Expressions(expressions) => expressions.len(),
            Body::Terms(polynomials) => polynomials.len(),
        }
    }
}

/// Write one `syl` polynomial as an expression.
fn expression_of_terms(
    origin: &str,
    index: usize,
    terms: &[Term],
    names: &[String],
) -> Result<String, String> {
    let mut expression = String::new();
    for term in terms {
        if term.exponents.len() != names.len() {
            return Err(format!(
                "{origin} polynomial {} holds {} exponents, the ring has {} variables",
                index + 1,
                term.exponents.len(),
                names.len()
            ));
        }
        let (sign, magnitude) = match term.coefficient.strip_prefix('-') {
            Some(magnitude) => ("-", magnitude),
            None => ("+", term.coefficient.as_str()),
        };
        if !expression.is_empty() || sign == "-" {
            expression.push_str(sign);
        }
        let constant = term.exponents.iter().all(|&e| e == 0);
        let mut written = false;
        if magnitude != "1" || constant {
            expression.push_str(magnitude);
            written = true;
        }
        for (name, &exponent) in names.iter().zip(&term.exponents) {
            if exponent == 0 {
                continue;
            }
            if written {
                expression.push('*');
            }
            written = true;
            expression.push_str(name);
            if exponent > 1 {
                expression.push_str(&format!("^{exponent}"));
            }
        }
    }
    if expression.is_empty() {
        expression.push('0');
    }
    Ok(expression)
}

/// Read a source.
///
/// `origin` names the source in every message the reader returns.
#[cfg(test)]
pub fn read(format: Format, origin: &str, text: &str) -> Result<Reading, String> {
    read_with_budget(format, origin, text, Budget::new()).map_err(|error| error.to_string())
}

/// Read a source under one parse budget.
pub fn read_with_budget(
    format: Format,
    origin: &str,
    text: &str,
    budget: Budget,
) -> Result<Reading, ReadError> {
    match format {
        Format::Ms => ms::read(origin, text).map_err(ReadError::Format),
        Format::Syl => syl::read(origin, text).map_err(ReadError::Format),
        Format::Text => text::read(origin, text).map_err(ReadError::Format),
        Format::Json => json::read(origin, text, budget).map_err(ReadError::Envelope),
    }
}

/// A polynomial as coefficient text and exponents, largest monomial first.
///
/// The coefficient carries its own sign, as `Display` writes it: `6` over a
/// prime field, `-2/3` over the rationals.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rendered {
    /// The terms.
    pub terms: Vec<(String, Vec<u32>)>,
}

/// Render a polynomial of either domain.
pub fn render<D: Domain>(polynomial: &Polynomial<D>) -> Rendered
where
    D::Coeff: fmt::Display,
{
    Rendered {
        terms: polynomial
            .terms()
            .map(|(coeff, exps)| {
                (
                    coeff.to_string(),
                    exps.iter().map(|&e| u32::from(e)).collect(),
                )
            })
            .collect(),
    }
}

/// Render a polynomial a verifier accepted.
///
/// The coefficients are residues, so no term carries a sign.
pub fn render_verified(polynomial: &verify::Poly) -> Rendered {
    Rendered {
        terms: polynomial
            .terms()
            .iter()
            .map(|term| (term.coeff().to_string(), term.mono().exps().to_vec()))
            .collect(),
    }
}

/// The polynomial as an expression in the syntax `parse_polynomial` reads.
///
/// It is what `Display` writes for the same polynomial.
pub fn expression(polynomial: &Rendered, names: &[String]) -> String {
    write_expression(polynomial, names, true)
}

/// The expression with no space around a sign.
///
/// It is the form `benchmarks/gb-comparison/gen.py` writes, which every
/// tool of the comparison reads.
pub fn compact_expression(polynomial: &Rendered, names: &[String]) -> String {
    write_expression(polynomial, names, false)
}

fn write_expression(polynomial: &Rendered, names: &[String], spaced: bool) -> String {
    if polynomial.terms.is_empty() {
        return "0".to_string();
    }
    let mut text = String::new();
    for (index, (coefficient, exponents)) in polynomial.terms.iter().enumerate() {
        let (negative, magnitude) = match coefficient.strip_prefix('-') {
            Some(magnitude) => (true, magnitude),
            None => (false, coefficient.as_str()),
        };
        if index > 0 {
            text.push_str(match (negative, spaced) {
                (true, true) => " - ",
                (true, false) => "-",
                (false, true) => " + ",
                (false, false) => "+",
            });
        } else if negative {
            text.push('-');
        }
        let constant = exponents.iter().all(|&e| e == 0);
        let mut written = false;
        if magnitude != "1" || constant {
            text.push_str(magnitude);
            written = true;
        }
        for (name, &exponent) in names.iter().zip(exponents) {
            if exponent == 0 {
                continue;
            }
            if written {
                text.push('*');
            }
            written = true;
            text.push_str(name);
            if exponent > 1 {
                text.push_str(&format!("^{exponent}"));
            }
        }
    }
    text
}

/// The ring a writer states.
pub struct Header<'a> {
    /// The variable names, largest variable first.
    pub names: &'a [String],
    /// The coefficient domain.
    pub domain: DomainClaim,
}

/// Write a polynomial system.
pub fn write(
    format: Format,
    header: &Header<'_>,
    polynomials: &[Rendered],
    out: &mut dyn Write,
) -> io::Result<()> {
    match format {
        Format::Ms => ms::write(header, polynomials, out),
        Format::Syl => syl::write(header, polynomials, out),
        Format::Text => text::write(header, polynomials, out),
        Format::Json => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "JSON output requires a computation record",
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sylvester::PolynomialRing;

    #[test]
    fn an_expression_is_what_display_writes() {
        let prime = PolynomialRing::prime_field(32003, ["x", "y"]).expect("the ring");
        let rationals = PolynomialRing::rationals(["x", "y"]).expect("the ring");
        for text in ["x^2*y - 3*y + 1", "x*y", "-x + y^3", "7"] {
            let f = prime.parse_polynomial(text).expect("the polynomial");
            assert_eq!(expression(&render(&f), prime.variables()), f.to_string());
        }
        for text in ["1/2*x^2 - y", "-2/3*x*y + 1", "x - 1"] {
            let f = rationals.parse_polynomial(text).expect("the polynomial");
            assert_eq!(
                expression(&render(&f), rationals.variables()),
                f.to_string()
            );
        }
    }

    #[test]
    fn a_compact_expression_holds_no_space() {
        let ring = PolynomialRing::prime_field(7, ["x", "y"]).expect("the ring");
        let f = ring.parse_polynomial("x^2 + 3*y").expect("the polynomial");
        assert_eq!(compact_expression(&render(&f), ring.variables()), "x^2+3*y");
    }

    #[test]
    fn the_ms_reader_takes_the_ring_from_the_header() {
        let reading = read(Format::Ms, "input", "x,y\n7\nx+y,\nx*y-1\n").expect("the reading");
        assert_eq!(
            reading.names.as_deref(),
            Some(["x".to_string(), "y".to_string()].as_slice())
        );
        assert_eq!(reading.domain, Some(DomainClaim::Prime(7)));
        assert_eq!(reading.body.len(), 2);
    }

    #[test]
    fn an_ms_characteristic_of_zero_is_the_rationals() {
        let reading = read(Format::Ms, "input", "x\n0\nx-1\n").expect("the reading");
        assert_eq!(reading.domain, Some(DomainClaim::Rationals));
    }

    #[test]
    fn the_syl_reader_names_no_ring() {
        let reading = read(Format::Syl, "input", "2\n1,1,0;-1,0,1\n").expect("the reading");
        assert!(reading.names.is_none());
        assert!(reading.domain.is_none());
        assert_eq!(reading.nvars, Some(2));
        let names = ["u".to_string(), "v".to_string()];
        assert_eq!(
            reading
                .body
                .expressions("input", &names)
                .expect("names of the right width"),
            ["u-v"]
        );
    }

    #[test]
    fn a_syl_term_of_the_wrong_width_is_an_error() {
        let reading = read(Format::Syl, "input", "2\n1,1,0,0\n");
        assert!(reading.is_err());
    }

    #[test]
    fn a_syl_coefficient_may_be_a_fraction() {
        let reading = read(Format::Syl, "input", "1\n-1/3,2\n").expect("the reading");
        let names = ["x".to_string()];
        assert_eq!(
            reading
                .body
                .expressions("input", &names)
                .expect("one variable"),
            ["-1/3*x^2"]
        );
    }

    #[test]
    fn a_text_source_may_name_its_variables_alone() {
        let reading = read(Format::Text, "input", "# vars: a, b\na - b\n").expect("the reading");
        assert_eq!(reading.nvars, Some(2));
        assert!(reading.domain.is_none());
    }

    #[test]
    fn a_text_source_names_one_coefficient_domain() {
        let reading = read(
            Format::Text,
            "input",
            "# vars: a\n# modulus: 7\n# coefficients: rationals\na\n",
        );
        assert!(reading.is_err());
    }
}
