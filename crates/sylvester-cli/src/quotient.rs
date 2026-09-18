//! The finite quotient algebra command.

use clap::Args;
use sylvester::{FiniteQuotient, QuotientError};

use super::{
    AnyRing, Budget, Deadline, Domain, Established, Fault, Format, GroebnerBasis, Header,
    Polynomial, PolynomialRing, PrimeField, Rationals, Rendered, Resolved, Source, SystemArgs,
    basis_memory_estimate, budget_after_sources_with_held, build_ring,
    compute_options_after_sources_with_held, format, ideal_memory_estimate, load, of_basis,
    of_compute, of_expression, of_ring, parse, polynomial_memory_estimate,
    polynomials_memory_estimate, qualification, rational_options_after_sources_with_held,
    write_text,
};

/// Arguments for the finite quotient algebra command.
#[derive(Args)]
pub(super) struct QuotientArgs {
    #[command(flatten)]
    pub(super) system: SystemArgs,
    /// Read a checked basis instead of computing one from generators.
    #[arg(long)]
    pub(super) from_basis: bool,
    /// Run the exact rational ideal equality check after computing a basis.
    #[arg(long)]
    pub(super) check_equality: bool,
    #[command(flatten)]
    pub(super) out: super::OutputArgs,
    /// Print the standard monomial staircase.
    #[arg(long)]
    pub(super) standard_monomials: bool,
    /// Print the vector-space dimension.
    #[arg(long)]
    pub(super) dimension: bool,
    /// Print the multiplication matrix of a polynomial.
    #[arg(long, visible_alias = "multiplication-matrix", value_name = "POLY")]
    pub(super) matrix: Option<String>,
    /// Print the characteristic polynomial of multiplication by a polynomial.
    #[arg(long, visible_alias = "characteristic-polynomial", value_name = "POLY")]
    pub(super) characteristic: Option<String>,
    /// Print the minimal polynomial of multiplication by a polynomial.
    #[arg(long, visible_alias = "minimal-polynomial", value_name = "POLY")]
    pub(super) minimal: Option<String>,
}

pub(super) fn command(started: std::time::Instant, args: QuotientArgs) -> Result<(), Fault> {
    let reporter = args.system.progress.reporter();
    reporter.phase("read");
    let deadline = Deadline::new(started, args.system.limits.timeout)?;
    let mut stdin_taken = false;
    let source = load(
        &args.system.input,
        args.system.in_format,
        &mut stdin_taken,
        deadline.budget(&args.system.limits)?,
        args.system.limits.memory,
        &deadline,
        0,
    )?;
    let resolved = super::resolve(&[&source], &args.system.ring)?;
    if args.from_basis {
        super::reject_engine_options(
            &args.system.engine,
            "describes a computation, and --from-basis runs none",
        )?;
    }
    let header = Header {
        names: &resolved.names,
        domain: resolved.domain,
    };
    let command = QuotientCommand {
        source: &source,
        resolved: &resolved,
        args: &args,
        deadline: &deadline,
        header: &header,
        reporter,
    };
    match build_ring(&resolved)? {
        AnyRing::Prime(ring) => quotient_prime(ring, command),
        AnyRing::Rational(ring) => quotient_rational(ring, command),
    }
}

struct QuotientCommand<'a> {
    source: &'a Source,
    resolved: &'a Resolved,
    args: &'a QuotientArgs,
    deadline: &'a Deadline,
    header: &'a Header<'a>,
    reporter: super::progress::Reporter,
}

fn quotient_prime(
    ring: PolynomialRing<PrimeField>,
    command: QuotientCommand<'_>,
) -> Result<(), Fault> {
    super::reject_rational_options(&command.args.system.engine)?;
    if command.args.check_equality {
        return Err(Fault::usage(
            "--check-equality applies only to the rational numbers",
        ));
    }
    command.reporter.phase("parse");
    let polynomials = parse(
        &ring,
        command.source,
        &command.resolved.names,
        !command.args.from_basis,
        command.deadline,
        &command.args.system.limits,
    )?;
    command.reporter.phase("compute");
    let basis = if command.args.from_basis {
        prime_supplied_basis(&ring, polynomials, &command)?
    } else {
        prime_computed_basis(&ring, polynomials, &command)?
    };
    let quotient = basis
        .finite_quotient(budget_after_sources_with_held(
            command.deadline,
            &command.args.system.limits,
            &[command.source],
            basis_memory_estimate(&basis),
        )?)
        .map_err(of_quotient)?;
    let output = QuotientOutput {
        args: command.args,
        deadline: command.deadline,
        header: command.header,
        source: command.source,
        basis_held: basis_memory_estimate(&basis),
    };
    command.reporter.phase("output");
    quotient_output(&quotient, output, None, false)
}

fn prime_supplied_basis(
    ring: &PolynomialRing<PrimeField>,
    polynomials: Vec<Polynomial<PrimeField>>,
    command: &QuotientCommand<'_>,
) -> Result<GroebnerBasis<PrimeField>, Fault> {
    let parsed_held = polynomials_memory_estimate(&polynomials);
    GroebnerBasis::<PrimeField>::from_polynomials(
        ring,
        polynomials,
        budget_after_sources_with_held(
            command.deadline,
            &command.args.system.limits,
            &[command.source],
            parsed_held,
        )?,
    )
    .map_err(of_basis)
}

fn prime_computed_basis(
    ring: &PolynomialRing<PrimeField>,
    polynomials: Vec<Polynomial<PrimeField>>,
    command: &QuotientCommand<'_>,
) -> Result<GroebnerBasis<PrimeField>, Fault> {
    let parsed_held = polynomials_memory_estimate(&polynomials);
    let options = compute_options_after_sources_with_held(
        &command.args.system.engine,
        &command.args.system.limits,
        command.deadline,
        &[command.source],
        parsed_held,
    )?;
    ring.ideal(polynomials)
        .map_err(of_ring)?
        .groebner_basis(options)
        .map_err(of_compute)
}

fn quotient_rational(
    ring: PolynomialRing<Rationals>,
    command: QuotientCommand<'_>,
) -> Result<(), Fault> {
    require_rational_input(&command)?;
    command.reporter.phase("parse");
    let polynomials = parse(
        &ring,
        command.source,
        &command.resolved.names,
        !command.args.from_basis,
        command.deadline,
        &command.args.system.limits,
    )?;
    command.reporter.phase("compute");
    let mut state = rational_quotient_basis(&ring, polynomials, &command)?;
    let basis_held = basis_memory_estimate(&state.basis);
    let equality_checked = rational_quotient_equality(&state, &command)?.is_some();
    state.input = None;
    let quotient = state
        .basis
        .finite_quotient(budget_after_sources_with_held(
            command.deadline,
            &command.args.system.limits,
            &[command.source],
            basis_held,
        )?)
        .map_err(of_quotient)?;
    let output = QuotientOutput {
        args: command.args,
        deadline: command.deadline,
        header: command.header,
        source: command.source,
        basis_held,
    };
    command.reporter.phase("output");
    quotient_output(&quotient, output, state.established, equality_checked)
}

fn require_rational_input(command: &QuotientCommand<'_>) -> Result<(), Fault> {
    if command.args.check_equality
        && command.args.from_basis
        && command.source.reading.record.is_none()
    {
        return Err(Fault::usage(
            "--check-equality with --from-basis requires a JSON record carrying the original input",
        ));
    }
    Ok(())
}

struct RationalQuotientBasis {
    basis: GroebnerBasis<Rationals>,
    established: Option<Established>,
    input: Option<sylvester::Ideal<Rationals>>,
}

fn rational_quotient_basis(
    ring: &PolynomialRing<Rationals>,
    polynomials: Vec<Polynomial<Rationals>>,
    command: &QuotientCommand<'_>,
) -> Result<RationalQuotientBasis, Fault> {
    if command.args.from_basis {
        let parsed_held = polynomials_memory_estimate(&polynomials);
        let basis = GroebnerBasis::<Rationals>::from_polynomials(
            ring,
            polynomials,
            budget_after_sources_with_held(
                command.deadline,
                &command.args.system.limits,
                &[command.source],
                parsed_held,
            )?,
        )
        .map_err(of_basis)?;
        Ok(RationalQuotientBasis {
            basis,
            established: None,
            input: None,
        })
    } else {
        let parsed_held = polynomials_memory_estimate(&polynomials);
        let options = rational_options_after_sources_with_held(
            &command.args.system.engine,
            &command.args.system.limits,
            command.deadline,
            &[command.source],
            parsed_held,
        )?;
        let ideal = ring.ideal(polynomials).map_err(of_ring)?;
        let basis = ideal.groebner_basis(options).map_err(of_compute)?;
        let established = basis.lift().map(|lift| lift.established);
        Ok(RationalQuotientBasis {
            basis,
            established,
            input: Some(ideal),
        })
    }
}

fn rational_quotient_equality(
    state: &RationalQuotientBasis,
    command: &QuotientCommand<'_>,
) -> Result<Option<sylvester::RationalEqualityCheck>, Fault> {
    if !command.args.check_equality {
        return Ok(None);
    }
    command.reporter.phase("equality-check");
    if command.args.from_basis {
        let record = command.source.reading.record.as_ref().ok_or_else(|| {
            Fault::usage("--check-equality with --from-basis requires a JSON record carrying the original input")
        })?;
        return record
            .check_rational(rational_equality_budget(state, command)?)
            .map(Some)
            .map_err(super::of_envelope);
    }
    let input = state.input.as_ref().ok_or_else(|| {
        Fault::usage("--check-equality requires the command's original input generators")
    })?;
    input
        .check_basis_equality(&state.basis, rational_equality_budget(state, command)?)
        .map(Some)
        .map_err(super::of_equality)
}

fn rational_equality_budget(
    state: &RationalQuotientBasis,
    command: &QuotientCommand<'_>,
) -> Result<Budget, Fault> {
    let input_held = state
        .input
        .as_ref()
        .map(ideal_memory_estimate)
        .unwrap_or_default();
    budget_after_sources_with_held(
        command.deadline,
        &command.args.system.limits,
        &[command.source],
        input_held.saturating_add(basis_memory_estimate(&state.basis)),
    )
}

struct QuotientOutput<'a> {
    args: &'a QuotientArgs,
    deadline: &'a Deadline,
    header: &'a Header<'a>,
    source: &'a Source,
    basis_held: usize,
}

fn quotient_output<D: Domain>(
    quotient: &FiniteQuotient<D>,
    output: QuotientOutput<'_>,
    established: Option<Established>,
    equality_checked: bool,
) -> Result<(), Fault>
where
    D::Coeff: std::fmt::Display,
{
    if output.args.out.out_format != Format::Text {
        return Err(Fault::usage(
            "finite quotient results use labeled text output, so --out-format text is required",
        ));
    }
    let default_output = !has_operation(output.args);
    let mut text = qualification(established);
    if equality_checked {
        text.insert_str(0, "# equality_check: passed\n");
    }
    append_basic_values(
        &mut text,
        quotient,
        output.header.names,
        output.args,
        default_output,
    );
    append_matrix(&mut text, quotient, &output)?;
    append_characteristic(&mut text, quotient, &output)?;
    append_minimal(&mut text, quotient, &output)?;
    write_text(output.args.out.output.as_deref(), &text)
}

fn has_operation(args: &QuotientArgs) -> bool {
    args.standard_monomials
        || args.dimension
        || args.matrix.is_some()
        || args.characteristic.is_some()
        || args.minimal.is_some()
}

fn append_basic_values<D: Domain>(
    text: &mut String,
    quotient: &FiniteQuotient<D>,
    names: &[String],
    args: &QuotientArgs,
    default_output: bool,
) {
    if args.dimension || default_output {
        text.push_str(&format!(
            "dimension: {}\n",
            quotient.vector_space_dimension()
        ));
    }
    if args.standard_monomials || default_output {
        let monomials = quotient
            .standard_monomials()
            .iter()
            .map(|exponents| standard_monomial(exponents, names))
            .collect::<Vec<_>>();
        if monomials.is_empty() {
            text.push_str("standard_monomials: []\n");
        } else {
            text.push_str(&format!("standard_monomials: {}\n", monomials.join(", ")));
        }
    }
}

fn append_matrix<D: Domain>(
    text: &mut String,
    quotient: &FiniteQuotient<D>,
    output: &QuotientOutput<'_>,
) -> Result<(), Fault> {
    let Some(element) = &output.args.matrix else {
        return Ok(());
    };
    let polynomial = quotient_element(quotient, element, output, "--matrix")?;
    let matrix = quotient
        .multiplication_matrix(&polynomial, operation_budget(output, element, &polynomial)?)
        .map_err(of_quotient)?;
    text.push_str(&format!("matrix: {matrix}\n"));
    Ok(())
}

fn append_characteristic<D: Domain>(
    text: &mut String,
    quotient: &FiniteQuotient<D>,
    output: &QuotientOutput<'_>,
) -> Result<(), Fault> {
    let Some(element) = &output.args.characteristic else {
        return Ok(());
    };
    let polynomial = quotient_element(quotient, element, output, "--characteristic")?;
    let characteristic = quotient
        .characteristic_polynomial(&polynomial, operation_budget(output, element, &polynomial)?)
        .map_err(of_quotient)?;
    text.push_str(&format!("characteristic_polynomial: {characteristic}\n"));
    Ok(())
}

fn append_minimal<D: Domain>(
    text: &mut String,
    quotient: &FiniteQuotient<D>,
    output: &QuotientOutput<'_>,
) -> Result<(), Fault> {
    let Some(element) = &output.args.minimal else {
        return Ok(());
    };
    let polynomial = quotient_element(quotient, element, output, "--minimal")?;
    let minimal = quotient
        .minimal_polynomial(&polynomial, operation_budget(output, element, &polynomial)?)
        .map_err(of_quotient)?;
    text.push_str(&format!("minimal_polynomial: {minimal}\n"));
    Ok(())
}

fn quotient_element<D: Domain>(
    quotient: &FiniteQuotient<D>,
    text: &str,
    output: &QuotientOutput<'_>,
    option: &str,
) -> Result<Polynomial<D>, Fault> {
    quotient
        .ring()
        .parse_polynomial_with_budget(
            text,
            budget_after_sources_with_held(
                output.deadline,
                &output.args.system.limits,
                &[output.source],
                output.basis_held.saturating_add(text.len()),
            )?,
        )
        .map_err(|error| {
            let fault = of_expression(error);
            Fault {
                message: format!("{option}: {}", fault.message),
                ..fault
            }
        })
}

fn operation_budget<D: Domain>(
    output: &QuotientOutput<'_>,
    element: &str,
    polynomial: &Polynomial<D>,
) -> Result<Budget, Fault> {
    budget_after_sources_with_held(
        output.deadline,
        &output.args.system.limits,
        &[output.source],
        output
            .basis_held
            .saturating_add(element.len())
            .saturating_add(polynomial_memory_estimate(polynomial)),
    )
}

fn standard_monomial(exponents: &[u16], names: &[String]) -> String {
    let rendered = Rendered {
        terms: vec![(
            "1".to_string(),
            exponents
                .iter()
                .map(|&exponent| u32::from(exponent))
                .collect(),
        )],
    };
    format::expression(&rendered, names)
}

fn of_quotient(error: QuotientError) -> Fault {
    match error {
        QuotientError::RingMismatch => Fault::usage(error.to_string()),
        QuotientError::NotFinite
        | QuotientError::ExponentLimit { .. }
        | QuotientError::Timeout
        | QuotientError::MemoryLimitExceeded => Fault::limit(error.to_string()),
    }
}
