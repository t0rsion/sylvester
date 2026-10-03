//! `sylv`, the command line interface of the sylvester crate.
//!
//! The binary reads a polynomial system in one of four formats and
//! resolves one ring for the command. The domain is a type parameter in
//! the library, so every subcommand matches [`AnyRing`] once and branches
//! no further.
//!
//! `docs/rational-design.md` section 8 fixes the subcommands, the formats, the
//! ring resolution rules, and the exit codes.

#[path = "io.rs"]
mod bounded_io;
mod format;
mod parse;
mod progress;
mod quotient;

use std::io::{self, Read, Write};
use std::mem::size_of;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use clap::{Args, Parser, Subcommand, ValueEnum};

use format::{DomainClaim, Format, Header, Reading, Rendered};
use sylvester::{
    Backend, BasisError, Budget, CertifyError, ComputeError, ComputeOptions, ComputeReport, Domain,
    EnvelopeError, EqualityCheckError, Established, ExpressionError, GroebnerBasis, HilbertError,
    HilbertSeries, NormalFormError, Polynomial, PolynomialRing, PrimeField, RationalEqualityCheck,
    RationalOptions, RationalStop, Rationals, ResultEnvelope, RingError, verify,
};

fn main() -> ExitCode {
    let started = Instant::now();
    let cli = Cli::parse();
    match run(started, cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(fault) => {
            if !fault.silent {
                eprintln!("sylv: {}", fault.message);
            }
            ExitCode::from(fault.code)
        }
    }
}

/// Why a command stopped, and the exit code that reports it.
///
/// The codes are the ones `docs/rational-design.md` section 8.4 lists: 1 for
/// certificate bytes the verifier rejected, 2 for usage and input, 3 for a
/// resource or structural limit, 4 for a defect in this crate, and 5 for a
/// file that cannot be read or written.
struct Fault {
    code: u8,
    message: String,
    /// A closed pipe gets no message. Every other fault writes one to
    /// standard error.
    silent: bool,
}

impl Fault {
    /// The verifier rejected untrusted bytes.
    fn rejected(message: impl Into<String>) -> Fault {
        Fault {
            code: 1,
            message: message.into(),
            silent: false,
        }
    }

    /// The command line, the input, or an operation the domain does not
    /// offer.
    fn usage(message: impl Into<String>) -> Fault {
        Fault {
            code: 2,
            message: message.into(),
            silent: false,
        }
    }

    /// A resource or structural limit.
    fn limit(message: impl Into<String>) -> Fault {
        Fault {
            code: 3,
            message: message.into(),
            silent: false,
        }
    }

    /// A defect in this crate.
    fn internal(message: impl Into<String>) -> Fault {
        Fault {
            code: 4,
            message: message.into(),
            silent: false,
        }
    }

    /// A file that cannot be read or written.
    fn io(message: impl Into<String>) -> Fault {
        Fault {
            code: 5,
            message: message.into(),
            silent: false,
        }
    }

    /// A reader of the output stopped. The code reports it and nothing is
    /// printed, so `sylv gb system.ms | head` stays quiet.
    fn broken_pipe() -> Fault {
        Fault {
            code: 5,
            message: "the reader of the output closed the pipe".to_string(),
            silent: true,
        }
    }

    /// Map an I/O error. A broken pipe is silent.
    fn of_io(context: &str, error: &io::Error) -> Fault {
        if error.kind() == io::ErrorKind::BrokenPipe {
            return Fault::broken_pipe();
        }
        Fault::io(format!("{context}: {error}"))
    }
}

#[derive(Parser)]
#[command(
    name = "sylv",
    version,
    about = "Gröbner bases and finite quotient algebras over prime fields and rationals",
    long_about = "Compute Gröbner bases and inspect finite quotient algebras under the grevlex order.\n\
        Certify prime-field bases and verify certificates.\n\n\
        A system is read from a file or from standard input in one of four formats:\n\
        ms (msolve), syl (the benchmark format), text (the crate's own syntax),\n\
        and json (a sylv-result-v1 computation record).\n\
        The ring comes from the input when the format carries one, and from\n\
        --modulus or --rationals when it does not."
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand)]
enum Command {
    /// Compute the reduced Gröbner basis of an ideal.
    #[command(alias = "basis")]
    Gb(GbArgs),
    /// Compute the basis, certify it, and write the certificate.
    Certify(CertifyArgs),
    /// Reduce a polynomial modulo a basis and print the remainder.
    #[command(name = "normal-form", alias = "nf")]
    NormalForm(NormalFormArgs),
    /// Report whether a polynomial belongs to the ideal a basis generates.
    Member(MemberArgs),
    /// Print the Hilbert series of the quotient by the leading monomial
    /// ideal.
    Hilbert(HilbertArgs),
    /// Print the Krull dimension of the quotient by an ideal.
    Dim(HilbertArgs),
    /// Inspect a finite-dimensional quotient algebra.
    #[command(alias = "standard-monomials")]
    Quotient(quotient::QuotientArgs),
    /// Check certificate bytes with the independent verifier.
    Verify(VerifyArgs),
}

/// The input system, the ring flags, and the limits of one computation.
#[derive(Args)]
struct SystemArgs {
    /// The file to read, or `-` for standard input.
    #[arg(value_name = "FILE", default_value = "-")]
    input: String,
    /// The format of FILE. The default is the extension of FILE, or text
    /// on standard input.
    #[arg(long, value_name = "F")]
    in_format: Option<Format>,
    #[command(flatten)]
    ring: RingArgs,
    #[command(flatten)]
    engine: EngineArgs,
    #[command(flatten)]
    limits: LimitArgs,
    #[command(flatten)]
    progress: progress::ProgressArgs,
}

/// Where a polynomial system goes.
#[derive(Args)]
struct OutputArgs {
    /// The format of the output.
    #[arg(long, value_name = "F", default_value_t = Format::Text)]
    out_format: Format,
    /// Write to this file instead of standard output.
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,
}

/// The coefficient domain, when no source names one.
#[derive(Args)]
struct RingArgs {
    /// The prime field modulus.
    #[arg(long, value_name = "P", conflicts_with = "rationals")]
    modulus: Option<u64>,
    /// Compute over the rational numbers.
    #[arg(long)]
    rationals: bool,
}

/// What the engine does.
#[derive(Args)]
struct EngineArgs {
    /// The backend of the computation. The default is f4.
    #[arg(long)]
    backend: Option<BackendArg>,
    /// Compute on this many threads. The default is the rayon global pool.
    #[arg(long, value_name = "N")]
    threads: Option<usize>,
    /// The stopping rule of the rational engine. The default is contains-input.
    #[arg(long)]
    stop: Option<StopArg>,
    /// The number of further primes the stopping rule observes. The
    /// default is 2.
    #[arg(long, value_name = "N")]
    extra_primes: Option<NonZeroUsize>,
}

/// The deadline and the memory limit of the whole command.
#[derive(Args)]
struct LimitArgs {
    /// Stop the command after this many seconds.
    #[arg(long, value_name = "SECS")]
    timeout: Option<f64>,
    /// Stop the command once the live data passes this many bytes.
    #[arg(
        long,
        alias = "memory-limit",
        value_name = "BYTES",
        value_parser = parse::bytes
    )]
    memory: Option<usize>,
}

#[derive(Args)]
struct GbArgs {
    #[command(flatten)]
    system: SystemArgs,
    #[command(flatten)]
    out: OutputArgs,
    /// Run the certified path: an isolated verifier accepts the basis
    /// before it is printed.
    #[arg(long)]
    certified: bool,
    /// Write the accepted certificate to this file. It implies
    /// --certified.
    #[arg(long, value_name = "PATH")]
    certificate: Option<PathBuf>,
    /// Write the counters of the run to standard error.
    #[arg(long)]
    report: bool,
    /// Run the exact rational ideal equality check after computing a basis.
    #[arg(long)]
    check_equality: bool,
}

#[derive(Args)]
struct CertifyArgs {
    #[command(flatten)]
    system: SystemArgs,
    #[command(flatten)]
    out: OutputArgs,
    /// Write the accepted certificate to this file.
    #[arg(long, value_name = "PATH")]
    certificate: Option<PathBuf>,
}

/// A basis and one polynomial, from two sources.
#[derive(Args)]
struct DivideArgs {
    /// The file the basis comes from.
    #[arg(long, value_name = "BASIS_FILE")]
    basis: String,
    /// The format of the basis file. The default is the extension, or text
    /// on standard input.
    #[arg(long, value_name = "F")]
    basis_format: Option<Format>,
    /// The polynomial, in the syntax of the ring.
    #[arg(
        long,
        value_name = "TEXT",
        conflicts_with = "poly_file",
        required_unless_present = "poly_file"
    )]
    poly: Option<String>,
    /// The file the polynomial comes from. It holds exactly one
    /// polynomial.
    #[arg(long, value_name = "FILE")]
    poly_file: Option<String>,
    /// The format of the polynomial file. The default is the extension, or
    /// text on standard input.
    #[arg(long, value_name = "F")]
    in_format: Option<Format>,
    #[command(flatten)]
    ring: RingArgs,
    #[command(flatten)]
    limits: LimitArgs,
    #[command(flatten)]
    progress: progress::ProgressArgs,
}

#[derive(Args)]
struct NormalFormArgs {
    #[command(flatten)]
    divide: DivideArgs,
    #[command(flatten)]
    out: OutputArgs,
    /// Print each division quotient together with the remainder.
    #[arg(long)]
    quotients: bool,
}

#[derive(Args)]
struct MemberArgs {
    #[command(flatten)]
    divide: DivideArgs,
    /// Write to this file instead of standard output.
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,
}

#[derive(Args)]
struct HilbertArgs {
    #[command(flatten)]
    system: SystemArgs,
    /// Read a basis instead of generators. The basis is checked, and no
    /// computation runs, so --backend, --threads, --stop, and
    /// --extra-primes do not apply.
    #[arg(long)]
    from_basis: bool,
    /// Write to this file instead of standard output.
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,
}

#[derive(Args)]
struct VerifyArgs {
    /// The certificate file, or `-` for standard input.
    #[arg(value_name = "CERT")]
    certificate: String,
    /// Reject a certificate longer than this many bytes, before reading
    /// it.
    #[arg(long, value_name = "N", value_parser = parse::bytes)]
    max_bytes: Option<usize>,
    /// Stop the verifier after this many seconds.
    #[arg(long, value_name = "SECS")]
    timeout: Option<f64>,
    /// The format of the printed basis.
    #[arg(long, value_name = "F", default_value_t = Format::Text)]
    out_format: Format,
    /// Write to this file instead of standard output.
    #[arg(short = 'o', long, value_name = "PATH")]
    output: Option<PathBuf>,
    /// Print nothing. The exit code carries the verdict.
    #[arg(long)]
    quiet: bool,
    #[command(flatten)]
    progress: progress::ProgressArgs,
}

#[derive(Clone, Copy, ValueEnum)]
enum BackendArg {
    /// The F4 engine, which writes `sylv-gb-cert-v2`.
    F4,
    /// The classic F5 oracle, which writes `sylv-gb-cert-v1`.
    Classic,
}

impl From<BackendArg> for Backend {
    fn from(backend: BackendArg) -> Backend {
        match backend {
            BackendArg::F4 => Backend::F4,
            BackendArg::Classic => Backend::Classic,
        }
    }
}

#[derive(Clone, Copy, ValueEnum)]
enum StopArg {
    /// Stop when the lift is unchanged over the extra primes.
    Unchanged,
    /// Stop as unchanged does, and then check the basis over `Q`.
    ContainsInput,
}

fn run(started: Instant, cli: Cli) -> Result<(), Fault> {
    match cli.command {
        Command::Gb(args) => {
            let certified = args.certified || args.certificate.is_some();
            if args.report && certified {
                return Err(Fault::usage(
                    "the certified path reports no counters, so --report does not apply with --certified",
                ));
            }
            basis_command(
                started,
                args.system,
                args.out,
                certified,
                args.certificate,
                args.report,
                args.check_equality,
            )
        }
        Command::Certify(args) => basis_command(
            started,
            args.system,
            args.out,
            true,
            args.certificate,
            false,
            false,
        ),
        Command::NormalForm(args) => divide_command(
            started,
            args.divide,
            if args.quotients {
                Divide::Quotients(args.out)
            } else {
                Divide::Remainder(args.out)
            },
        ),
        Command::Member(args) => divide_command(started, args.divide, Divide::Member(args.output)),
        Command::Hilbert(args) => hilbert_command(started, args, Report::Series),
        Command::Dim(args) => hilbert_command(started, args, Report::Dimension),
        Command::Quotient(args) => quotient::command(started, args),
        Command::Verify(args) => verify_command(started, args),
    }
}

/// The deadline of the whole command.
///
/// The binary takes one instant at start and hands the remaining time to
/// every library call.
struct Deadline {
    started: Instant,
    timeout: Option<Duration>,
}

impl Deadline {
    fn new(started: Instant, seconds: Option<f64>) -> Result<Deadline, Fault> {
        let timeout = match seconds {
            None => None,
            // try_from_secs_f64 rejects a negative, an infinite, a NaN,
            // and a value past what a duration holds. from_secs_f64
            // panics on each of them.
            Some(seconds) => match Duration::try_from_secs_f64(seconds) {
                Ok(timeout) => Some(timeout),
                Err(_) => {
                    return Err(Fault::usage(
                        "--timeout takes a count of seconds that is zero or more and fits a duration",
                    ));
                }
            },
        };
        Ok(Deadline { started, timeout })
    }

    /// The time the command has left, or an exhausted deadline.
    fn remaining(&self) -> Result<Option<Duration>, Fault> {
        match self.timeout {
            None => Ok(None),
            Some(timeout) => timeout
                .checked_sub(self.started.elapsed())
                .filter(|left| !left.is_zero())
                .map(Some)
                .ok_or_else(|| Fault::limit("the command passed its deadline")),
        }
    }

    /// The budget of one library call.
    fn budget(&self, limits: &LimitArgs) -> Result<Budget, Fault> {
        let mut budget = Budget::new();
        if let Some(remaining) = self.remaining()? {
            budget = budget.timeout(remaining);
        }
        if let Some(bytes) = limits.memory {
            budget = budget.memory_limit(bytes);
        }
        Ok(budget)
    }

    fn timed_budget(&self) -> Result<Budget, Fault> {
        let mut budget = Budget::new();
        if let Some(remaining) = self.remaining()? {
            budget = budget.timeout(remaining);
        }
        Ok(budget)
    }

    fn absolute(&self) -> Option<Instant> {
        self.timeout
            .and_then(|timeout| self.started.checked_add(timeout))
    }
}

/// One ring over the two domains.
///
/// The binary learns the domain at run time, so every subcommand matches
/// this enum once.
enum AnyRing {
    Prime(PolynomialRing<PrimeField>),
    Rational(PolynomialRing<Rationals>),
}

trait CliDomain: Domain {
    fn checked_basis(
        ring: &PolynomialRing<Self>,
        polynomials: Vec<Polynomial<Self>>,
        budget: Budget,
    ) -> Result<GroebnerBasis<Self>, BasisError>;
}

impl CliDomain for PrimeField {
    fn checked_basis(
        ring: &PolynomialRing<Self>,
        polynomials: Vec<Polynomial<Self>>,
        budget: Budget,
    ) -> Result<GroebnerBasis<Self>, BasisError> {
        GroebnerBasis::<PrimeField>::from_polynomials(ring, polynomials, budget)
    }
}

impl CliDomain for Rationals {
    fn checked_basis(
        ring: &PolynomialRing<Self>,
        polynomials: Vec<Polynomial<Self>>,
        budget: Budget,
    ) -> Result<GroebnerBasis<Self>, BasisError> {
        GroebnerBasis::<Rationals>::from_polynomials(ring, polynomials, budget)
    }
}

/// A source and what it holds.
struct Source {
    origin: String,
    reading: Reading,
}

/// The ring one command computes in.
struct Resolved {
    names: Vec<String>,
    domain: DomainClaim,
}

/// Read one source.
///
/// `stdin_taken` holds whether an earlier source already read standard
/// input. One stream cannot carry two independently parsed inputs.
fn load(
    input: &str,
    declared: Option<Format>,
    stdin_taken: &mut bool,
    budget: Budget,
    memory: Option<usize>,
    deadline: &Deadline,
    held: usize,
) -> Result<Source, Fault> {
    let read_memory = memory.map(|limit| limit.saturating_sub(held));
    let (origin, text, format) = if input == "-" {
        if *stdin_taken {
            return Err(Fault::usage(
                "two sources read standard input, one stream carries one input",
            ));
        }
        *stdin_taken = true;
        let text = read_source_text(
            io::stdin(),
            "standard input",
            read_memory,
            deadline.absolute(),
        )?;
        (
            "standard input".to_string(),
            text,
            declared.unwrap_or(Format::Text),
        )
    } else {
        let path = Path::new(input);
        let format = match declared {
            Some(format) => format,
            None => Format::of_path(path).ok_or_else(|| {
                Fault::usage(format!(
                    "{input} has no format extension, name the format with a format option"
                ))
            })?,
        };
        let text = bounded_io::read_path_text(
            path,
            bounded_io::ReadLimits::new(deadline.absolute(), read_memory, None),
        )
        .map_err(|error| of_source_read(input, error))?;
        (input.to_string(), text, format)
    };
    let reading_budget =
        budget_for_held_source(budget, held.saturating_add(text.capacity()), memory);
    let reading = format::read_with_budget(format, &origin, &text, reading_budget)
        .map_err(|error| of_format_read(&origin, error))?;
    Ok(Source { origin, reading })
}

fn budget_for_held_source(budget: Budget, held: usize, memory: Option<usize>) -> Budget {
    match memory {
        Some(limit) => budget.memory_limit(limit.saturating_sub(held)),
        None => budget,
    }
}

fn read_source_text<R: Read>(
    reader: R,
    origin: &str,
    memory: Option<usize>,
    deadline: Option<Instant>,
) -> Result<String, Fault> {
    bounded_io::read_text(reader, bounded_io::ReadLimits::new(deadline, memory, None))
        .map_err(|error| of_source_read(origin, error))
}

/// Resolve the one ring of a command.
///
/// Three rules, in order. Every source that names a part of the ring must
/// name the same one. If a source names the coefficient domain, `--modulus`
/// and `--rationals` do not apply. If no source names one, exactly one of
/// the two is required.
///
/// A source that names no variables takes the names of a source that does,
/// and `x1 .. xn` when no source names any.
fn resolve(sources: &[&Source], ring: &RingArgs) -> Result<Resolved, Fault> {
    let names = resolve_names(sources)?;
    let claimed = resolve_claimed_domain(sources)?;
    let domain = select_domain(claimed, ring)?;
    let names = match names {
        Some(names) => names,
        None => resolve_counted_names(sources)?,
    };
    validate_variable_counts(sources, names.len())?;
    Ok(Resolved { names, domain })
}

fn resolve_names(sources: &[&Source]) -> Result<Option<Vec<String>>, Fault> {
    let mut named: Option<(&str, &Vec<String>)> = None;
    for source in sources {
        let Some(names) = &source.reading.names else {
            continue;
        };
        match named {
            None => named = Some((&source.origin, names)),
            Some((first, first_names)) if first_names != names => {
                return Err(Fault::usage(format!(
                    "{first} and {} name different variables",
                    source.origin
                )));
            }
            Some(_) => {}
        }
    }
    Ok(named.map(|(_, names)| names.clone()))
}

fn resolve_claimed_domain<'a>(
    sources: &'a [&'a Source],
) -> Result<Option<(&'a str, DomainClaim)>, Fault> {
    let mut claimed: Option<(&str, DomainClaim)> = None;
    for source in sources {
        let Some(domain) = source.reading.domain else {
            continue;
        };
        match claimed {
            None => claimed = Some((&source.origin, domain)),
            Some((first, first_domain)) if first_domain != domain => {
                return Err(Fault::usage(format!(
                    "{first} names {first_domain}, and {} names {domain}",
                    source.origin
                )));
            }
            Some(_) => {}
        }
    }
    Ok(claimed)
}

fn select_domain(
    claimed: Option<(&str, DomainClaim)>,
    ring: &RingArgs,
) -> Result<DomainClaim, Fault> {
    let flag = match (ring.modulus, ring.rationals) {
        (Some(modulus), _) => Some(DomainClaim::Prime(modulus)),
        (None, true) => Some(DomainClaim::Rationals),
        (None, false) => None,
    };
    match (claimed, flag) {
        (Some((origin, _)), Some(_)) => Err(Fault::usage(format!(
            "{origin} names the coefficient domain, so --modulus and --rationals do not apply"
        ))),
        (Some((_, domain)), None) => Ok(domain),
        (None, Some(flag)) => Ok(flag),
        (None, None) => Err(Fault::usage(
            "no source names a coefficient domain, so --modulus or --rationals is required",
        )),
    }
}

fn resolve_counted_names(sources: &[&Source]) -> Result<Vec<String>, Fault> {
    let mut counted: Option<(&str, usize)> = None;
    for source in sources {
        let Some(count) = source.reading.nvars else {
            continue;
        };
        match counted {
            None => counted = Some((&source.origin, count)),
            Some((first, first_count)) if first_count != count => {
                return Err(Fault::usage(format!(
                    "{first} names {first_count} variables, and {} names {count}",
                    source.origin
                )));
            }
            Some(_) => {}
        }
    }
    let (_, count) = counted.ok_or_else(|| {
        Fault::usage(
            "no source names the variables, a text source names them in a \"# vars:\" line",
        )
    })?;
    Ok((1..=count).map(|index| format!("x{index}")).collect())
}

fn validate_variable_counts(sources: &[&Source], ring_count: usize) -> Result<(), Fault> {
    for source in sources {
        if let Some(count) = source.reading.nvars
            && count != ring_count
        {
            return Err(Fault::usage(format!(
                "{} names {count} variables, the ring has {}",
                source.origin, ring_count
            )));
        }
    }
    Ok(())
}

/// Build the ring the command computes in.
fn build_ring(resolved: &Resolved) -> Result<AnyRing, Fault> {
    match resolved.domain {
        DomainClaim::Prime(modulus) => PolynomialRing::prime_field(modulus, &resolved.names)
            .map(AnyRing::Prime)
            .map_err(of_ring),
        DomainClaim::Rationals => PolynomialRing::rationals(&resolved.names)
            .map(AnyRing::Rational)
            .map_err(of_ring),
    }
}

/// Parse the polynomials of one source in the command's ring.
fn parse<D: Domain>(
    ring: &PolynomialRing<D>,
    source: &Source,
    names: &[String],
    input: bool,
    deadline: &Deadline,
    limits: &LimitArgs,
) -> Result<Vec<Polynomial<D>>, Fault> {
    parse_with_held(ring, source, names, input, deadline, limits, 0)
}

fn parse_with_held<D: Domain>(
    ring: &PolynomialRing<D>,
    source: &Source,
    names: &[String],
    input: bool,
    deadline: &Deadline,
    limits: &LimitArgs,
    extra_held: usize,
) -> Result<Vec<Polynomial<D>>, Fault> {
    let count = parse_count(source, input);
    let mut retained = source_memory_estimate(source).saturating_add(extra_held);
    retained = retained.saturating_add(count.saturating_mul(size_of::<Polynomial<D>>()));
    check_parse_memory(limits.memory, retained)?;
    let mut polynomials = Vec::with_capacity(count);
    for index in 0..count {
        let (expression, transient_bytes) = source_expression(source, names, input, index)?;
        let expression_budget = parse_budget(deadline, limits, retained, transient_bytes)?;
        let polynomial = ring
            .parse_polynomial_with_budget(&expression, expression_budget)
            .map_err(|error| {
                let fault = of_expression(error);
                Fault {
                    message: format!(
                        "{} polynomial {}: {}",
                        source.origin,
                        index + 1,
                        fault.message
                    ),
                    ..fault
                }
            })?;
        retained = retained.saturating_add(polynomial_memory_estimate(&polynomial));
        check_parse_memory(limits.memory, retained)?;
        polynomials.push(polynomial);
    }
    Ok(polynomials)
}

fn parse_count(source: &Source, input: bool) -> usize {
    if input && let Some(record) = &source.reading.record {
        return record.input().len();
    }
    source.reading.body.len()
}

fn source_expression<'a>(
    source: &'a Source,
    names: &[String],
    input: bool,
    index: usize,
) -> Result<(std::borrow::Cow<'a, str>, usize), Fault> {
    if input && let Some(record) = &source.reading.record {
        let expression = record.input().get(index).ok_or_else(|| {
            Fault::usage(format!("{} has no polynomial {}", source.origin, index + 1))
        })?;
        return Ok((std::borrow::Cow::Borrowed(expression.as_str()), 0));
    }
    let expression = source
        .reading
        .body
        .expression(&source.origin, index, names)
        .map_err(Fault::usage)?;
    let transient_bytes = match &expression {
        std::borrow::Cow::Borrowed(_) => 0,
        std::borrow::Cow::Owned(value) => value.len(),
    };
    Ok((expression, transient_bytes))
}

fn parse_budget(
    deadline: &Deadline,
    limits: &LimitArgs,
    retained: usize,
    expression_bytes: usize,
) -> Result<Budget, Fault> {
    let mut budget = deadline.budget(limits)?;
    if let Some(limit) = limits.memory {
        let held = retained.saturating_add(expression_bytes);
        let available = limit.checked_sub(held).ok_or_else(parse_memory_fault)?;
        budget = budget.memory_limit(available);
    }
    Ok(budget)
}

fn check_parse_memory(limit: Option<usize>, retained: usize) -> Result<(), Fault> {
    if limit.is_some_and(|limit| retained > limit) {
        Err(parse_memory_fault())
    } else {
        Ok(())
    }
}

fn parse_memory_fault() -> Fault {
    Fault::limit("the parsed input passed its memory limit")
}

fn source_memory_estimate(source: &Source) -> usize {
    let mut bytes = size_of::<Source>().saturating_add(source.origin.capacity());
    if let Some(names) = &source.reading.names {
        bytes = bytes
            .saturating_add(names.capacity().saturating_mul(size_of::<String>()))
            .saturating_add(
                names
                    .iter()
                    .fold(0, |bytes, name| bytes.saturating_add(name.capacity())),
            );
    }
    bytes = match &source.reading.body {
        format::Body::Expressions(expressions) => bytes
            .saturating_add(expressions.capacity().saturating_mul(size_of::<String>()))
            .saturating_add(expressions.iter().fold(0, |bytes, expression| {
                bytes.saturating_add(expression.capacity())
            })),
        format::Body::Terms(polynomials) => bytes
            .saturating_add(
                polynomials
                    .capacity()
                    .saturating_mul(size_of::<Vec<format::Term>>()),
            )
            .saturating_add(polynomials.iter().fold(0, |bytes, terms| {
                bytes
                    .saturating_add(terms.capacity().saturating_mul(size_of::<format::Term>()))
                    .saturating_add(terms.iter().fold(0, |bytes, term| {
                        bytes
                            .saturating_add(term.coefficient.capacity())
                            .saturating_add(
                                term.exponents.capacity().saturating_mul(size_of::<u16>()),
                            )
                    }))
            })),
    };
    if let Some(record) = &source.reading.record {
        bytes = bytes
            .saturating_add(size_of::<ResultEnvelope>())
            .saturating_add(record.variables().len().saturating_mul(size_of::<String>()))
            .saturating_add(record.input().len().saturating_mul(size_of::<String>()))
            .saturating_add(record.basis().len().saturating_mul(size_of::<String>()))
            .saturating_add(
                record
                    .variables()
                    .iter()
                    .chain(record.input())
                    .chain(record.basis())
                    .fold(0, |bytes, value| {
                        bytes.saturating_add(value.capacity().saturating_add(size_of::<String>()))
                    }),
            )
            .saturating_add(record.certificate().map_or(0, <[u8]>::len));
    }
    bytes
}

fn polynomial_memory_estimate<D: Domain>(polynomial: &Polynomial<D>) -> usize {
    polynomial.estimated_heap_bytes()
}

fn polynomials_memory_estimate<D: Domain>(polynomials: &[Polynomial<D>]) -> usize {
    size_of::<Vec<Polynomial<D>>>()
        .saturating_add(polynomials.len().saturating_mul(size_of::<Polynomial<D>>()))
        .saturating_add(polynomials.iter().fold(0, |bytes, polynomial| {
            bytes.saturating_add(polynomial_memory_estimate(polynomial))
        }))
}

fn ideal_memory_estimate<D: Domain>(ideal: &sylvester::Ideal<D>) -> usize {
    size_of::<sylvester::Ideal<D>>().saturating_add(polynomials_memory_estimate(ideal.generators()))
}

fn basis_memory_estimate<D: Domain>(basis: &GroebnerBasis<D>) -> usize {
    size_of::<GroebnerBasis<D>>().saturating_add(polynomials_memory_estimate(basis))
}

fn budget_after_sources_with_held(
    deadline: &Deadline,
    limits: &LimitArgs,
    sources: &[&Source],
    extra_held: usize,
) -> Result<Budget, Fault> {
    let mut budget = deadline.budget(limits)?;
    if let Some(limit) = limits.memory {
        let held = sources
            .iter()
            .fold(0usize, |bytes, source| {
                bytes.saturating_add(source_memory_estimate(source))
            })
            .saturating_add(extra_held);
        budget = budget.memory_limit(limit.checked_sub(held).ok_or_else(parse_memory_fault)?);
    }
    Ok(budget)
}

fn compute_options_after_sources_with_held(
    engine: &EngineArgs,
    limits: &LimitArgs,
    deadline: &Deadline,
    sources: &[&Source],
    extra_held: usize,
) -> Result<ComputeOptions, Fault> {
    compute_options_with_budget(
        engine,
        budget_after_sources_with_held(deadline, limits, sources, extra_held)?,
    )
}

fn compute_options_with_budget(
    engine: &EngineArgs,
    budget: Budget,
) -> Result<ComputeOptions, Fault> {
    let mut options = ComputeOptions::new().budget(budget);
    if let Some(backend) = engine.backend {
        options = options.backend(backend.into());
    }
    if let Some(threads) = engine.threads {
        options = options.threads(threads);
    }
    Ok(options)
}

fn rational_options_after_sources_with_held(
    engine: &EngineArgs,
    limits: &LimitArgs,
    deadline: &Deadline,
    sources: &[&Source],
    extra_held: usize,
) -> Result<RationalOptions, Fault> {
    rational_options_with_compute(
        engine,
        compute_options_after_sources_with_held(engine, limits, deadline, sources, extra_held)?,
    )
}

fn rational_options_with_compute(
    engine: &EngineArgs,
    compute: ComputeOptions,
) -> Result<RationalOptions, Fault> {
    let extra = engine
        .extra_primes
        .unwrap_or(NonZeroUsize::new(2).expect("2 is not zero"));
    let stop = match engine.stop {
        Some(StopArg::Unchanged) => RationalStop::Unchanged { extra },
        None | Some(StopArg::ContainsInput) => RationalStop::ContainsInput { extra },
    };
    Ok(RationalOptions::new().compute(compute).stop(stop))
}

/// Report the options that describe a computation the command does not run.
fn reject_engine_options(engine: &EngineArgs, why: &str) -> Result<(), Fault> {
    let named = [
        ("--backend", engine.backend.is_some()),
        ("--threads", engine.threads.is_some()),
        ("--stop", engine.stop.is_some()),
        ("--extra-primes", engine.extra_primes.is_some()),
    ];
    for (option, given) in named {
        if given {
            return Err(Fault::usage(format!("{option} {why}")));
        }
    }
    Ok(())
}

/// Report the rational options on a prime-field ring.
fn reject_rational_options(engine: &EngineArgs) -> Result<(), Fault> {
    for (option, given) in [
        ("--stop", engine.stop.is_some()),
        ("--extra-primes", engine.extra_primes.is_some()),
    ] {
        if given {
            return Err(Fault::usage(format!(
                "{option} describes the rational engine, and the ring is a prime field"
            )));
        }
    }
    Ok(())
}

fn basis_command(
    started: Instant,
    system: SystemArgs,
    out: OutputArgs,
    certified: bool,
    certificate: Option<PathBuf>,
    report: bool,
    check_equality: bool,
) -> Result<(), Fault> {
    let reporter = system.progress.reporter();
    reject_output_collision(out.output.as_deref(), certificate.as_deref())?;
    reporter.phase("read");
    let deadline = Deadline::new(started, system.limits.timeout)?;
    let mut stdin_taken = false;
    let source = load(
        &system.input,
        system.in_format,
        &mut stdin_taken,
        deadline.budget(&system.limits)?,
        system.limits.memory,
        &deadline,
        0,
    )?;
    let resolved = resolve(&[&source], &system.ring)?;
    let header = Header {
        names: &resolved.names,
        domain: resolved.domain,
    };
    let command = BasisCommand {
        source: &source,
        resolved: &resolved,
        system: &system,
        deadline: &deadline,
        out: &out,
        header: &header,
        certified,
        certificate: certificate.as_deref(),
        report,
        check_equality,
        reporter,
    };
    match build_ring(&resolved)? {
        AnyRing::Prime(ring) => prime_basis_command(ring, command)?,
        AnyRing::Rational(ring) => rational_basis_command(ring, command)?,
    }
    Ok(())
}

struct BasisCommand<'a> {
    source: &'a Source,
    resolved: &'a Resolved,
    system: &'a SystemArgs,
    deadline: &'a Deadline,
    out: &'a OutputArgs,
    header: &'a Header<'a>,
    certified: bool,
    certificate: Option<&'a Path>,
    report: bool,
    check_equality: bool,
    reporter: progress::Reporter,
}

fn prime_basis_command(
    ring: PolynomialRing<PrimeField>,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    command.reporter.phase("parse");
    if command.check_equality {
        return Err(Fault::usage(
            "--check-equality applies only to the rational numbers",
        ));
    }
    reject_rational_options(&command.system.engine)?;
    let generators = parse(
        &ring,
        command.source,
        &command.resolved.names,
        true,
        command.deadline,
        &command.system.limits,
    )?;
    let parsed_held = polynomials_memory_estimate(&generators);
    let ideal = ring.ideal(generators).map_err(of_ring)?;
    let options = compute_options_after_sources_with_held(
        &command.system.engine,
        &command.system.limits,
        command.deadline,
        &[command.source],
        parsed_held,
    )?;
    if command.certified {
        command.reporter.phase("compute");
        return write_certified_basis(ideal, options, command);
    }
    if command.report {
        command.reporter.phase("compute");
        return write_reported_basis(ideal, options, command);
    }
    command.reporter.phase("compute");
    let basis = ideal.groebner_basis(options).map_err(of_compute)?;
    command.reporter.phase("output");
    write_prime_result(&ideal, &basis, &command, None)
}

fn write_certified_basis(
    ideal: sylvester::Ideal<PrimeField>,
    options: ComputeOptions,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    let accepted = ideal
        .groebner_basis_certified(options)
        .map_err(of_certify)?;
    command.reporter.phase("output");
    if let Some(path) = command.certificate {
        write_bytes(path, accepted.certificate())?;
    }
    write_prime_result(&ideal, accepted.basis(), &command, Some(&accepted))
}

fn write_reported_basis(
    ideal: sylvester::Ideal<PrimeField>,
    options: ComputeOptions,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    let (basis, counters) = ideal
        .groebner_basis_with_report(options)
        .map_err(of_compute)?;
    command.reporter.phase("output");
    write_prime_result(&ideal, &basis, &command, None)?;
    write_report(&counters, false)
}

fn rational_basis_command(
    ring: PolynomialRing<Rationals>,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    command.reporter.phase("parse");
    if command.certified {
        return Err(Fault::usage(
            "the rational numbers have no certified path, so --certified and --certificate do not apply",
        ));
    }
    let generators = parse(
        &ring,
        command.source,
        &command.resolved.names,
        true,
        command.deadline,
        &command.system.limits,
    )?;
    let parsed_held = polynomials_memory_estimate(&generators);
    let ideal = ring.ideal(generators).map_err(of_ring)?;
    let options = rational_options_after_sources_with_held(
        &command.system.engine,
        &command.system.limits,
        command.deadline,
        &[command.source],
        parsed_held,
    )?;
    command.reporter.phase("compute");
    let (basis, report) = rational_computation(&ideal, options, command.report)?;
    let equality = rational_equality(&ideal, &basis, &command)?;
    command.reporter.phase("output");
    write_rational_result(&ideal, &basis, &command, equality.as_ref())?;
    if let Some(report) = report {
        write_report(&report, equality.is_some())
    } else {
        Ok(())
    }
}

fn rational_computation(
    ideal: &sylvester::Ideal<Rationals>,
    options: RationalOptions,
    report: bool,
) -> Result<(GroebnerBasis<Rationals>, Option<ComputeReport>), Fault> {
    if report {
        let (basis, report) = ideal
            .groebner_basis_with_report(options)
            .map_err(of_compute)?;
        Ok((basis, Some(report)))
    } else {
        Ok((ideal.groebner_basis(options).map_err(of_compute)?, None))
    }
}

fn rational_equality(
    ideal: &sylvester::Ideal<Rationals>,
    basis: &GroebnerBasis<Rationals>,
    command: &BasisCommand<'_>,
) -> Result<Option<RationalEqualityCheck>, Fault> {
    if !command.check_equality {
        return Ok(None);
    }
    command.reporter.phase("equality-check");
    ideal
        .check_basis_equality(
            basis,
            budget_after_sources_with_held(
                command.deadline,
                &command.system.limits,
                &[command.source],
                ideal_memory_estimate(ideal).saturating_add(basis_memory_estimate(basis)),
            )?,
        )
        .map(Some)
        .map_err(of_equality)
}

/// What a division command prints.
enum Divide {
    /// The remainder, as a polynomial system.
    Remainder(OutputArgs),
    /// Every quotient and the remainder, as labeled text lines.
    Quotients(OutputArgs),
    /// The membership verdict, as `true` or `false`.
    Member(Option<PathBuf>),
}

fn divide_command(started: Instant, divide: DivideArgs, print: Divide) -> Result<(), Fault> {
    if divide.poly.is_some() && divide.in_format.is_some() {
        return Err(Fault::usage(
            "--in-format applies to --poly-file, not to --poly",
        ));
    }
    let reporter = divide.progress.reporter();
    reporter.phase("read");
    let deadline = Deadline::new(started, divide.limits.timeout)?;
    let (basis_source, poly_source, resolved) = divide_sources(&divide, &deadline)?;
    let header = Header {
        names: &resolved.names,
        domain: resolved.domain,
    };
    let command = DivideCommand {
        basis_source: &basis_source,
        poly_source: poly_source.as_ref(),
        resolved: &resolved,
        deadline: &deadline,
        args: &divide,
        header: &header,
        print,
        reporter,
    };
    match build_ring(&resolved)? {
        AnyRing::Prime(ring) => divide_in_ring(ring, command),
        AnyRing::Rational(ring) => divide_in_ring(ring, command),
    }
}

fn divide_sources(
    divide: &DivideArgs,
    deadline: &Deadline,
) -> Result<(Source, Option<Source>, Resolved), Fault> {
    let mut stdin_taken = false;
    let basis_source = load(
        &divide.basis,
        divide.basis_format,
        &mut stdin_taken,
        deadline.budget(&divide.limits)?,
        divide.limits.memory,
        deadline,
        0,
    )?;
    let poly_source = load_polynomial_source(
        divide,
        deadline,
        &mut stdin_taken,
        source_memory_estimate(&basis_source),
    )?;
    let mut sources = vec![&basis_source];
    if let Some(source) = &poly_source {
        sources.push(source);
    }
    let resolved = resolve(&sources, &divide.ring)?;
    Ok((basis_source, poly_source, resolved))
}

fn load_polynomial_source(
    divide: &DivideArgs,
    deadline: &Deadline,
    stdin_taken: &mut bool,
    held: usize,
) -> Result<Option<Source>, Fault> {
    let Some(path) = &divide.poly_file else {
        return Ok(None);
    };
    let source = load(
        path,
        divide.in_format,
        stdin_taken,
        deadline.budget(&divide.limits)?,
        divide.limits.memory,
        deadline,
        held,
    )?;
    if source.reading.body.len() != 1 {
        return Err(Fault::usage(format!(
            "{} holds {} polynomials, --poly-file holds exactly one",
            source.origin,
            source.reading.body.len()
        )));
    }
    Ok(Some(source))
}

struct DivideCommand<'a> {
    basis_source: &'a Source,
    poly_source: Option<&'a Source>,
    resolved: &'a Resolved,
    deadline: &'a Deadline,
    args: &'a DivideArgs,
    header: &'a Header<'a>,
    print: Divide,
    reporter: progress::Reporter,
}

fn divide_budget_with_held(
    command: &DivideCommand<'_>,
    extra_held: usize,
) -> Result<Budget, Fault> {
    let mut sources = vec![command.basis_source];
    if let Some(source) = command.poly_source {
        sources.push(source);
    }
    budget_after_sources_with_held(command.deadline, &command.args.limits, &sources, extra_held)
}

fn divide_in_ring<D: CliDomain>(
    ring: PolynomialRing<D>,
    command: DivideCommand<'_>,
) -> Result<(), Fault>
where
    D::Coeff: std::fmt::Display,
{
    command.reporter.phase("parse");
    let basis_extra = command
        .poly_source
        .map(source_memory_estimate)
        .unwrap_or_default();
    let polynomials = parse_with_held(
        &ring,
        command.basis_source,
        &command.resolved.names,
        false,
        command.deadline,
        &command.args.limits,
        basis_extra,
    )?;
    let parsed_basis_held = polynomials_memory_estimate(&polynomials);
    let basis = D::checked_basis(
        &ring,
        polynomials,
        divide_budget_with_held(&command, parsed_basis_held)?,
    )
    .map_err(of_basis)?;
    command.reporter.phase("compute");
    let basis_held = basis_memory_estimate(&basis);
    let direct_text_held = command
        .args
        .poly
        .as_ref()
        .map(String::len)
        .unwrap_or_default();
    let polynomial = one_polynomial(
        &ring,
        command.args,
        command.poly_source,
        &command.resolved.names,
        command.deadline,
        basis_held,
        divide_budget_with_held(&command, basis_held.saturating_add(direct_text_held))?,
    )?;
    command.reporter.phase("output");
    let output_budget = divide_budget_with_held(
        &command,
        basis_held
            .saturating_add(polynomial_memory_estimate(&polynomial))
            .saturating_add(direct_text_held),
    )?;
    divide_output(
        &basis,
        &polynomial,
        command.header,
        command.print,
        output_budget,
    )
}

/// The polynomial `normal-form` and `member` divide.
fn one_polynomial<D: Domain>(
    ring: &PolynomialRing<D>,
    divide: &DivideArgs,
    source: Option<&Source>,
    names: &[String],
    deadline: &Deadline,
    extra_held: usize,
    budget: Budget,
) -> Result<Polynomial<D>, Fault> {
    match (source, &divide.poly) {
        (Some(source), _) => {
            let mut polynomials = parse_with_held(
                ring,
                source,
                names,
                false,
                deadline,
                &divide.limits,
                extra_held,
            )?;
            Ok(polynomials.remove(0))
        }
        (None, Some(text)) => ring
            .parse_polynomial_with_budget(text, budget)
            .map_err(|error| {
                let fault = of_expression(error);
                Fault {
                    message: format!("--poly: {}", fault.message),
                    ..fault
                }
            }),
        (None, None) => Err(Fault::usage(
            "name the polynomial with --poly or --poly-file",
        )),
    }
}

/// Print the remainder, or the membership verdict.
fn divide_output<D: Domain>(
    basis: &GroebnerBasis<D>,
    f: &Polynomial<D>,
    header: &Header<'_>,
    print: Divide,
    budget: Budget,
) -> Result<(), Fault>
where
    D::Coeff: std::fmt::Display,
{
    match print {
        Divide::Remainder(out) => {
            let remainder = basis.normal_form(f, budget).map_err(of_normal_form)?;
            write_system(&out, header, &[format::render(&remainder)])
        }
        Divide::Quotients(out) => {
            if out.out_format != Format::Text {
                return Err(Fault::usage(
                    "--quotients uses the labeled text output, so --out-format text is required",
                ));
            }
            let result = basis.divide(f, budget).map_err(of_normal_form)?;
            let mut text = String::new();
            for (index, quotient) in result.quotients().iter().enumerate() {
                text.push_str(&format!(
                    "quotient[{index}]: {}\n",
                    format::expression(&format::render(quotient), header.names)
                ));
            }
            text.push_str(&format!(
                "remainder: {}\n",
                format::expression(&format::render(result.remainder()), header.names)
            ));
            write_text(out.output.as_deref(), &text)
        }
        Divide::Member(path) => {
            let member = basis.contains(f, budget).map_err(of_normal_form)?;
            write_text(path.as_deref(), &format!("{member}\n"))
        }
    }
}

/// What a Hilbert command prints.
enum Report {
    Series,
    Dimension,
}

fn hilbert_command(started: Instant, args: HilbertArgs, report: Report) -> Result<(), Fault> {
    let reporter = args.system.progress.reporter();
    reporter.phase("read");
    let (deadline, source, resolved) = prepare_hilbert(started, &args)?;
    reporter.phase("parse");
    reporter.phase("compute");
    let (series, established) = match build_ring(&resolved)? {
        AnyRing::Prime(ring) => (
            prime_hilbert(ring, &source, &resolved, &args, &deadline)?,
            None,
        ),
        AnyRing::Rational(ring) => rational_hilbert(ring, &source, &resolved, &args, &deadline)?,
    };
    let text = match report {
        Report::Series => hilbert_text(&series),
        Report::Dimension => dimension_line(&series),
    };
    reporter.phase("output");
    write_text(
        args.output.as_deref(),
        &format!("{}{text}", qualification(established)),
    )
}

fn prepare_hilbert(
    started: Instant,
    args: &HilbertArgs,
) -> Result<(Deadline, Source, Resolved), Fault> {
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
    let resolved = resolve(&[&source], &args.system.ring)?;
    if args.from_basis {
        reject_engine_options(
            &args.system.engine,
            "describes a computation, and --from-basis runs none",
        )?;
    }
    Ok((deadline, source, resolved))
}

fn prime_hilbert(
    ring: PolynomialRing<PrimeField>,
    source: &Source,
    resolved: &Resolved,
    args: &HilbertArgs,
    deadline: &Deadline,
) -> Result<HilbertSeries, Fault> {
    reject_rational_options(&args.system.engine)?;
    let polynomials = parse(
        &ring,
        source,
        &resolved.names,
        !args.from_basis,
        deadline,
        &args.system.limits,
    )?;
    let parsed_held = polynomials_memory_estimate(&polynomials);
    let basis = if args.from_basis {
        GroebnerBasis::<PrimeField>::from_polynomials(
            &ring,
            polynomials,
            budget_after_sources_with_held(deadline, &args.system.limits, &[source], parsed_held)?,
        )
        .map_err(of_basis)?
    } else {
        let options = compute_options_after_sources_with_held(
            &args.system.engine,
            &args.system.limits,
            deadline,
            &[source],
            parsed_held,
        )?;
        ring.ideal(polynomials)
            .map_err(of_ring)?
            .groebner_basis(options)
            .map_err(of_compute)?
    };
    let basis_held = basis_memory_estimate(&basis);
    basis
        .hilbert_series(budget_after_sources_with_held(
            deadline,
            &args.system.limits,
            &[source],
            basis_held,
        )?)
        .map_err(of_hilbert)
}

fn rational_hilbert(
    ring: PolynomialRing<Rationals>,
    source: &Source,
    resolved: &Resolved,
    args: &HilbertArgs,
    deadline: &Deadline,
) -> Result<(HilbertSeries, Option<Established>), Fault> {
    let polynomials = parse(
        &ring,
        source,
        &resolved.names,
        !args.from_basis,
        deadline,
        &args.system.limits,
    )?;
    let parsed_held = polynomials_memory_estimate(&polynomials);
    let basis = if args.from_basis {
        GroebnerBasis::<Rationals>::from_polynomials(
            &ring,
            polynomials,
            budget_after_sources_with_held(deadline, &args.system.limits, &[source], parsed_held)?,
        )
        .map_err(of_basis)?
    } else {
        let options = rational_options_after_sources_with_held(
            &args.system.engine,
            &args.system.limits,
            deadline,
            &[source],
            parsed_held,
        )?;
        ring.ideal(polynomials)
            .map_err(of_ring)?
            .groebner_basis(options)
            .map_err(of_compute)?
    };
    let established = basis.lift().map(|lift| lift.established);
    let basis_held = basis_memory_estimate(&basis);
    let series = basis
        .hilbert_series(budget_after_sources_with_held(
            deadline,
            &args.system.limits,
            &[source],
            basis_held,
        )?)
        .map_err(of_hilbert)?;
    Ok((series, established))
}

/// The lines that qualify a value read off a lifted rational basis.
///
/// The multimodular driver is a heuristic. The value describes the ideal
/// the lifted basis generates. `Established` is what relates that ideal
/// to the input. A basis read from a file carries no lift and no
/// qualification.
fn qualification(established: Option<Established>) -> String {
    match established {
        None => String::new(),
        Some(Established::Unchanged) => concat!(
            "# established: unchanged\n",
            "# the value describes the ideal the lifted basis generates, not the input ideal\n",
        )
        .to_string(),
        Some(Established::ContainsInput) => concat!(
            "# established: contains-input\n",
            "# the value describes an ideal that contains the input ideal\n",
        )
        .to_string(),
    }
}

/// The `hilbert` output schema.
fn hilbert_text(series: &HilbertSeries) -> String {
    let numerator: Vec<String> = if series.numerator().is_empty() {
        // The unit ideal has the zero numerator, which holds no
        // coefficient. The line states the zero polynomial instead of
        // nothing.
        vec!["0".to_string()]
    } else {
        series.numerator().iter().map(ToString::to_string).collect()
    };
    let multiplicity = match series.multiplicity() {
        Some(value) => value.to_string(),
        None => "none".to_string(),
    };
    format!(
        "series: {}\ndenominator_power: {}\n{}multiplicity: {multiplicity}\n",
        numerator.join(", "),
        series.denominator_power(),
        dimension_line(series),
    )
}

/// The `dimension:` line, which is `none` for the unit ideal.
fn dimension_line(series: &HilbertSeries) -> String {
    match series.dimension() {
        Some(dimension) => format!("dimension: {dimension}\n"),
        None => "dimension: none\n".to_string(),
    }
}

fn verify_command(started: Instant, args: VerifyArgs) -> Result<(), Fault> {
    let reporter = if args.quiet {
        progress::ProgressArgs::default().reporter()
    } else {
        args.progress.reporter()
    };
    reporter.phase("read");
    let deadline = Deadline::new(started, args.timeout)?;
    let mut limits = verify::Limits::default();
    if let Some(max_bytes) = args.max_bytes {
        limits.max_bytes = max_bytes;
    }
    let bytes = read_capped(&args.certificate, limits.max_bytes, &deadline)?;
    reporter.phase("verify");
    let (verified, names, synthetic_names) =
        verify_record_or_certificate(&bytes, &mut limits, &deadline)?;
    if args.quiet {
        return Ok(());
    }
    let header = Header {
        names: &names,
        domain: DomainClaim::Prime(verified.modulus()),
    };
    let rendered: Vec<Rendered> = verified
        .basis()
        .iter()
        .map(format::render_verified)
        .collect();
    let out = OutputArgs {
        out_format: args.out_format,
        output: args.output,
    };
    reporter.phase("output");
    write_system(&out, &header, &rendered)?;
    if synthetic_names {
        eprintln!(
            "sylv: the certificate carries no variable names, so the basis prints under the synthetic names x1 to x{}",
            verified.nvars()
        );
    }
    Ok(())
}

fn verify_record_or_certificate(
    bytes: &[u8],
    limits: &mut verify::Limits,
    deadline: &Deadline,
) -> Result<(verify::VerifiedGb, Vec<String>, bool), Fault> {
    match ResultEnvelope::from_json(bytes, deadline.timed_budget()?) {
        Ok(record) => {
            let verified = record
                .verify_prime(limits, deadline.timed_budget()?)
                .map_err(of_record_verification)?;
            Ok((verified, record.variables().to_vec(), false))
        }
        Err(EnvelopeError::Format(_)) => {
            if let Some(remaining) = deadline.remaining()? {
                limits.deadline = Instant::now().checked_add(remaining);
            }
            let verified = verify::verify_with_limits(bytes, limits).map_err(|error| {
                if error.is_exhaustion() {
                    Fault::limit(format!("the verifier stopped: {error}"))
                } else {
                    Fault::rejected(format!("the verifier rejected the certificate: {error}"))
                }
            })?;
            let names = (1..=verified.nvars())
                .map(|index| format!("x{index}"))
                .collect();
            Ok((verified, names, true))
        }
        Err(error) => Err(of_record_verification(error)),
    }
}

/// Read a certificate under its byte cap and command deadline.
fn read_capped(input: &str, max_bytes: usize, deadline: &Deadline) -> Result<Vec<u8>, Fault> {
    let limits = bounded_io::ReadLimits::new(deadline.absolute(), None, Some(max_bytes));
    let result = if input == "-" {
        bounded_io::read_bytes(io::stdin(), limits)
    } else {
        bounded_io::read_path_bytes(Path::new(input), limits)
    };
    result.map_err(|error| of_bounded_read(input, error))
}

fn render_all<D: Domain>(polynomials: &[Polynomial<D>]) -> Vec<Rendered>
where
    D::Coeff: std::fmt::Display,
{
    polynomials.iter().map(format::render).collect()
}

fn write_prime_result(
    ideal: &sylvester::Ideal<PrimeField>,
    basis: &GroebnerBasis<PrimeField>,
    command: &BasisCommand<'_>,
    certified: Option<&sylvester::CertifiedGroebnerBasis>,
) -> Result<(), Fault> {
    if command.out.out_format == Format::Json {
        let certificate_bytes = certified.map_or(0, |value| value.certificate().len());
        let construction_budget = result_budget(command, ideal, basis, certificate_bytes)?;
        let record = match certified {
            Some(certified) => {
                ResultEnvelope::from_certified(ideal, certified, construction_budget)
            }
            None => ResultEnvelope::from_prime(ideal, basis, None, construction_budget),
        }
        .map_err(of_envelope)?;
        return write_result_json(
            command.out,
            &record,
            result_budget(command, ideal, basis, certificate_bytes)?,
        );
    }
    write_system(command.out, command.header, &render_all(basis))
}

fn write_rational_result(
    ideal: &sylvester::Ideal<Rationals>,
    basis: &GroebnerBasis<Rationals>,
    command: &BasisCommand<'_>,
    equality: Option<&RationalEqualityCheck>,
) -> Result<(), Fault> {
    if command.out.out_format == Format::Json {
        let equality_bytes = equality.map_or(0, equality_memory_estimate);
        let construction_budget = result_budget(command, ideal, basis, equality_bytes)?;
        let record = match equality {
            Some(equality) => ResultEnvelope::from_checked_rational(equality, construction_budget),
            None => ResultEnvelope::from_rational(ideal, basis, construction_budget),
        }
        .map_err(of_envelope)?;
        return write_result_json(
            command.out,
            &record,
            result_budget(command, ideal, basis, equality_bytes)?,
        );
    }
    let rendered = render_all(basis);
    if equality.is_some() && command.out.out_format == Format::Text {
        let mut bytes = b"# equality_check: passed\n".to_vec();
        format::write(Format::Text, command.header, &rendered, &mut bytes)
            .map_err(|error| Fault::of_io("cannot render output", &error))?;
        return write_output_bytes(command.out.output.as_deref(), &bytes);
    }
    write_system(command.out, command.header, &rendered)
}

fn equality_memory_estimate(check: &RationalEqualityCheck) -> usize {
    ideal_memory_estimate(check.input())
        .saturating_add(basis_memory_estimate(check.basis()))
        .saturating_add(check.origins().iter().fold(0usize, |bytes, row| {
            bytes
                .saturating_add(polynomials_memory_estimate(row))
                .saturating_add(
                    row.capacity()
                        .saturating_sub(row.len())
                        .saturating_mul(size_of::<Polynomial<Rationals>>()),
                )
        }))
}

fn result_budget<D: Domain>(
    command: &BasisCommand<'_>,
    ideal: &sylvester::Ideal<D>,
    basis: &GroebnerBasis<D>,
    extra_held: usize,
) -> Result<Budget, Fault> {
    budget_after_sources_with_held(
        command.deadline,
        &command.system.limits,
        &[command.source],
        ideal_memory_estimate(ideal)
            .saturating_add(basis_memory_estimate(basis))
            .saturating_add(extra_held),
    )
}

fn write_result_json(
    out: &OutputArgs,
    record: &ResultEnvelope,
    budget: Budget,
) -> Result<(), Fault> {
    let mut bytes = record.to_json(budget).map_err(of_envelope)?;
    bytes.push(b'\n');
    write_output_bytes(out.output.as_deref(), &bytes)
}

/// Write a polynomial system to a file or to standard output.
fn write_system(
    out: &OutputArgs,
    header: &Header<'_>,
    polynomials: &[Rendered],
) -> Result<(), Fault> {
    if out.out_format == Format::Json {
        return Err(Fault::usage(
            "JSON output is available for Gröbner basis results only",
        ));
    }
    match &out.output {
        Some(path) => {
            let mut bytes = Vec::new();
            format::write(out.out_format, header, polynomials, &mut bytes)
                .map_err(|error| Fault::of_io("cannot render output", &error))?;
            write_output_bytes(Some(path), &bytes)
        }
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            format::write(out.out_format, header, polynomials, &mut handle)
                .and_then(|()| handle.flush())
                .map_err(|error| Fault::of_io("cannot write standard output", &error))
        }
    }
}

fn write_output_bytes(path: Option<&Path>, bytes: &[u8]) -> Result<(), Fault> {
    match path {
        Some(path) => write_file_atomic(path, bytes)
            .map_err(|error| Fault::io(format!("cannot write {}: {error}", path.display()))),
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            handle
                .write_all(bytes)
                .and_then(|()| handle.flush())
                .map_err(|error| Fault::of_io("cannot write standard output", &error))
        }
    }
}

/// Write plain text to a file or to standard output.
fn write_text(path: Option<&Path>, text: &str) -> Result<(), Fault> {
    match path {
        Some(path) => write_file_atomic(path, text.as_bytes())
            .map_err(|error| Fault::io(format!("cannot write {}: {error}", path.display()))),
        None => {
            let stdout = io::stdout();
            let mut handle = stdout.lock();
            handle
                .write_all(text.as_bytes())
                .and_then(|()| handle.flush())
                .map_err(|error| Fault::of_io("cannot write standard output", &error))
        }
    }
}

fn write_bytes(path: &Path, bytes: &[u8]) -> Result<(), Fault> {
    write_file_atomic(path, bytes)
        .map_err(|error| Fault::io(format!("cannot write {}: {error}", path.display())))
}

fn write_file_atomic(path: &Path, bytes: &[u8]) -> io::Result<()> {
    bounded_io::write_atomic(path, bytes)
}

fn reject_output_collision(output: Option<&Path>, certificate: Option<&Path>) -> Result<(), Fault> {
    bounded_io::reject_same_targets(certificate, output)
        .map_err(|_| Fault::usage("the certificate and result output must use different paths"))
}

/// Write the counters of a run to standard error, one per line.
fn write_report(report: &ComputeReport, equality_checked: bool) -> Result<(), Fault> {
    let mut lines = vec![
        format!(
            "backend: {}",
            match report.backend {
                Backend::F4 => "f4",
                Backend::Classic => "classic",
            }
        ),
        format!("elapsed_seconds: {}", report.elapsed.as_secs_f64()),
        format!("threads_used: {}", report.threads_used),
    ];
    if let Some(counters) = report.counters {
        lines.extend([
            format!("batches: {}", counters.batches),
            format!("matrix_rows: {}", counters.matrix_rows),
            format!("matrix_columns: {}", counters.matrix_columns),
            format!("matrix_nonzeros: {}", counters.matrix_nonzeros),
            format!("zero_rows: {}", counters.zero_rows),
            format!("new_pivots: {}", counters.new_pivots),
            format!("lane_restarts: {}", counters.lane_restarts),
            format!("batch_retries: {}", counters.batch_retries),
            format!("basis_monomials: {}", counters.basis_monomials),
            format!("pairs_generated: {}", counters.pairs_generated),
            format!(
                "pairs_discarded_product: {}",
                counters.pairs_discarded_product
            ),
            format!("pairs_discarded_b: {}", counters.pairs_discarded_b),
            format!("pairs_discarded_m: {}", counters.pairs_discarded_m),
            format!("pairs_discarded_f: {}", counters.pairs_discarded_f),
        ]);
    }
    if let Some(concurrency) = report.modular_concurrency {
        lines.push(format!("modular_concurrency: {concurrency}"));
    }
    if let Some(lift) = report.modular {
        lines.extend([
            format!("primes_consumed: {}", lift.primes_consumed),
            format!("primes_skipped: {}", lift.primes_skipped),
            format!("primes_folded: {}", lift.primes_folded),
            format!("primes_discarded: {}", lift.primes_discarded),
            format!("confirming_primes: {}", lift.confirming_primes),
            format!("modulus_bits: {}", lift.modulus_bits),
            format!(
                "established: {}",
                match lift.established {
                    sylvester::Established::Unchanged => "unchanged",
                    sylvester::Established::ContainsInput => "contains-input",
                }
            ),
        ]);
    }
    if equality_checked {
        lines.push("equality_check: passed".to_string());
    }
    let stderr = io::stderr();
    let mut handle = stderr.lock();
    for line in lines {
        writeln!(handle, "{line}")
            .map_err(|error| Fault::of_io("cannot write standard error", &error))?;
    }
    Ok(())
}

fn of_expression(error: ExpressionError) -> Fault {
    match error {
        ExpressionError::Timeout
        | ExpressionError::MemoryLimitExceeded
        | ExpressionError::ExponentLimit { .. }
        | ExpressionError::Compute(ComputeError::Timeout)
        | ExpressionError::Compute(ComputeError::MemoryLimitExceeded)
        | ExpressionError::Compute(ComputeError::DegreeLimit { .. })
        | ExpressionError::Compute(ComputeError::ExponentLimit { .. })
        | ExpressionError::Compute(ComputeError::TableFull)
        | ExpressionError::Compute(ComputeError::PrimesExhausted) => {
            Fault::limit(error.to_string())
        }
        ExpressionError::Parse(_)
        | ExpressionError::Ring(_)
        | ExpressionError::NestingLimit { .. }
        | ExpressionError::Arithmetic(_) => Fault::usage(error.to_string()),
    }
}

fn of_format_read(origin: &str, error: format::ReadError) -> Fault {
    match error {
        format::ReadError::Format(message) => Fault::usage(message),
        format::ReadError::Envelope(error) => {
            let fault = match error {
                EnvelopeError::Timeout | EnvelopeError::MemoryLimitExceeded => {
                    Fault::limit(error.to_string())
                }
                EnvelopeError::Ring(_) | EnvelopeError::RingMismatch | EnvelopeError::Format(_) => {
                    Fault::usage(error.to_string())
                }
                EnvelopeError::Expression(error) => of_expression(error),
                EnvelopeError::EqualityCheck(error) => of_equality(error),
                EnvelopeError::Verification(error) => {
                    if error.is_exhaustion() {
                        Fault::limit(error.to_string())
                    } else {
                        Fault::rejected(error.to_string())
                    }
                }
                EnvelopeError::CertificateMismatch => Fault::rejected(error.to_string()),
            };
            Fault {
                message: format!("{origin}: {}", fault.message),
                ..fault
            }
        }
    }
}

fn of_source_read(origin: &str, error: bounded_io::ReadError) -> Fault {
    of_bounded_read(origin, error)
}

fn of_bounded_read(origin: &str, error: bounded_io::ReadError) -> Fault {
    match error {
        bounded_io::ReadError::Io(error) => Fault::of_io(&format!("cannot read {origin}"), &error),
        bounded_io::ReadError::Timeout | bounded_io::ReadError::MemoryLimitExceeded => {
            Fault::limit(format!("{origin}: {error}"))
        }
        bounded_io::ReadError::TooLong { limit } => {
            Fault::limit(format!("{origin} is longer than {limit} bytes"))
        }
        bounded_io::ReadError::InvalidUtf8 => Fault::usage(format!("{origin} is not valid UTF-8")),
    }
}

fn of_envelope(error: EnvelopeError) -> Fault {
    match error {
        EnvelopeError::Timeout | EnvelopeError::MemoryLimitExceeded => {
            Fault::limit(error.to_string())
        }
        EnvelopeError::Expression(error) => of_expression(error),
        EnvelopeError::EqualityCheck(error) => of_equality(error),
        EnvelopeError::Verification(error) => {
            if error.is_exhaustion() {
                Fault::limit(error.to_string())
            } else {
                Fault::rejected(error.to_string())
            }
        }
        EnvelopeError::CertificateMismatch => Fault::rejected(error.to_string()),
        EnvelopeError::RingMismatch | EnvelopeError::Format(_) | EnvelopeError::Ring(_) => {
            Fault::internal(error.to_string())
        }
    }
}

fn of_equality(error: EqualityCheckError) -> Fault {
    match error {
        EqualityCheckError::Candidate(error) => {
            let fault = of_basis(error);
            Fault {
                message: format!("the equality check failed: {}", fault.message),
                ..fault
            }
        }
        EqualityCheckError::RingMismatch
        | EqualityCheckError::ReverseMembership { .. }
        | EqualityCheckError::ForwardMembership { .. } => {
            Fault::usage(format!("the equality check failed: {error}"))
        }
        EqualityCheckError::ExponentLimit { .. }
        | EqualityCheckError::Timeout
        | EqualityCheckError::MemoryLimitExceeded => {
            Fault::limit(format!("the equality check stopped: {error}"))
        }
    }
}

fn of_record_verification(error: EnvelopeError) -> Fault {
    match error {
        EnvelopeError::Timeout | EnvelopeError::MemoryLimitExceeded => {
            Fault::limit(error.to_string())
        }
        EnvelopeError::Verification(error) if error.is_exhaustion() => {
            Fault::limit(format!("the verifier stopped: {error}"))
        }
        EnvelopeError::Verification(error) => {
            Fault::rejected(format!("the verifier rejected the certificate: {error}"))
        }
        EnvelopeError::CertificateMismatch => Fault::rejected(error.to_string()),
        EnvelopeError::Expression(_)
        | EnvelopeError::EqualityCheck(_)
        | EnvelopeError::Format(_)
        | EnvelopeError::RingMismatch
        | EnvelopeError::Ring(_) => {
            Fault::rejected(format!("the result record is not verified: {error}"))
        }
    }
}

fn of_ring(error: RingError) -> Fault {
    Fault::usage(error.to_string())
}

fn of_compute(error: ComputeError) -> Fault {
    Fault::limit(error.to_string())
}

fn of_hilbert(error: HilbertError) -> Fault {
    Fault::limit(error.to_string())
}

fn of_normal_form(error: NormalFormError) -> Fault {
    match error {
        NormalFormError::RingMismatch => Fault::usage(error.to_string()),
        _ => Fault::limit(error.to_string()),
    }
}

fn of_basis(error: BasisError) -> Fault {
    match error {
        BasisError::ExponentLimit { .. }
        | BasisError::Timeout
        | BasisError::MemoryLimitExceeded => Fault::limit(error.to_string()),
        _ => Fault::usage(error.to_string()),
    }
}

fn of_certify(error: CertifyError) -> Fault {
    match error {
        CertifyError::Engine(_)
        | CertifyError::WriterExhausted(_)
        | CertifyError::VerifierExhausted(_)
        | CertifyError::CapExceeded { .. } => Fault::limit(error.to_string()),
        CertifyError::Emitter(_) | CertifyError::InputMismatch | CertifyError::Rejected(_) => {
            Fault::internal(error.to_string())
        }
    }
}
