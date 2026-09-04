//! `sylv`, the command line interface of the sylvester crate.
//!
//! The binary reads a polynomial system in one of three formats and
//! resolves one ring for the command. The domain is a type parameter in
//! the library, so every subcommand matches [`AnyRing`] once and branches
//! no further.
//!
//! `docs/rational-design.md` section 8 fixes the subcommands, the formats, the
//! ring resolution rules, and the exit codes.

mod format;

use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::time::{Duration, Instant};

use clap::{Args, Parser, Subcommand, ValueEnum};

use format::{DomainClaim, Format, Header, Reading, Rendered};
use sylvester::{
    Backend, BasisError, Budget, CertifyError, ComputeError, ComputeOptions, ComputeReport, Domain,
    Established, GroebnerBasis, HilbertError, HilbertSeries, NormalFormError, Polynomial,
    PolynomialRing, PrimeField, RationalOptions, RationalStop, Rationals, RingError, verify,
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
    about = "Gröbner bases over prime fields and the rational numbers",
    long_about = "Compute, certify, and check Gröbner bases under the grevlex order.\n\n\
        A system is read from a file or from standard input in one of three formats:\n\
        ms (msolve), syl (the benchmark format), and text (the crate's own syntax).\n\
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
    #[arg(long, alias = "memory-limit", value_name = "BYTES")]
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
}

#[derive(Args)]
struct NormalFormArgs {
    #[command(flatten)]
    divide: DivideArgs,
    #[command(flatten)]
    out: OutputArgs,
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
    #[arg(long, value_name = "N")]
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
            )
        }
        Command::Certify(args) => basis_command(
            started,
            args.system,
            args.out,
            true,
            args.certificate,
            false,
        ),
        Command::NormalForm(args) => {
            divide_command(started, args.divide, Divide::Remainder(args.out))
        }
        Command::Member(args) => divide_command(started, args.divide, Divide::Member(args.output)),
        Command::Hilbert(args) => hilbert_command(started, args, Report::Series),
        Command::Dim(args) => hilbert_command(started, args, Report::Dimension),
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
fn load(input: &str, declared: Option<Format>, stdin_taken: &mut bool) -> Result<Source, Fault> {
    let (origin, text, format) = if input == "-" {
        if *stdin_taken {
            return Err(Fault::usage(
                "two sources read standard input, one stream carries one input",
            ));
        }
        *stdin_taken = true;
        let mut text = String::new();
        io::stdin()
            .read_to_string(&mut text)
            .map_err(|error| Fault::io(format!("cannot read standard input: {error}")))?;
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
        let text = fs::read_to_string(path)
            .map_err(|error| Fault::io(format!("cannot read {input}: {error}")))?;
        (input.to_string(), text, format)
    };
    let reading = format::read(format, &origin, &text).map_err(Fault::usage)?;
    Ok(Source { origin, reading })
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
) -> Result<Vec<Polynomial<D>>, Fault> {
    let expressions = source
        .reading
        .body
        .expressions(&source.origin, names)
        .map_err(Fault::usage)?;
    let mut polynomials = Vec::with_capacity(expressions.len());
    for (index, expression) in expressions.iter().enumerate() {
        polynomials.push(ring.parse_polynomial(expression).map_err(|error| {
            Fault::usage(format!(
                "{} polynomial {}: {error}",
                source.origin,
                index + 1
            ))
        })?);
    }
    Ok(polynomials)
}

/// The options of one prime-field computation.
fn compute_options(
    engine: &EngineArgs,
    limits: &LimitArgs,
    deadline: &Deadline,
) -> Result<ComputeOptions, Fault> {
    let mut options = ComputeOptions::new().budget(deadline.budget(limits)?);
    if let Some(backend) = engine.backend {
        options = options.backend(backend.into());
    }
    if let Some(threads) = engine.threads {
        options = options.threads(threads);
    }
    Ok(options)
}

/// The options of one rational computation.
fn rational_options(
    engine: &EngineArgs,
    limits: &LimitArgs,
    deadline: &Deadline,
) -> Result<RationalOptions, Fault> {
    let extra = engine
        .extra_primes
        .unwrap_or(NonZeroUsize::new(2).expect("2 is not zero"));
    let stop = match engine.stop {
        Some(StopArg::Unchanged) => RationalStop::Unchanged { extra },
        None | Some(StopArg::ContainsInput) => RationalStop::ContainsInput { extra },
    };
    Ok(RationalOptions::new()
        .compute(compute_options(engine, limits, deadline)?)
        .stop(stop))
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
) -> Result<(), Fault> {
    let deadline = Deadline::new(started, system.limits.timeout)?;
    let mut stdin_taken = false;
    let source = load(&system.input, system.in_format, &mut stdin_taken)?;
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
}

fn prime_basis_command(
    ring: PolynomialRing<PrimeField>,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    reject_rational_options(&command.system.engine)?;
    let generators = parse(&ring, command.source, &command.resolved.names)?;
    let ideal = ring.ideal(generators).map_err(of_ring)?;
    let options = compute_options(
        &command.system.engine,
        &command.system.limits,
        command.deadline,
    )?;
    if command.certified {
        return write_certified_basis(
            ideal,
            options,
            command.out,
            command.header,
            command.certificate,
        );
    }
    if command.report {
        return write_reported_basis(ideal, options, command.out, command.header);
    }
    let basis = ideal.groebner_basis(options).map_err(of_compute)?;
    write_system(command.out, command.header, &render_all(&basis))
}

fn write_certified_basis(
    ideal: sylvester::Ideal<PrimeField>,
    options: ComputeOptions,
    out: &OutputArgs,
    header: &Header<'_>,
    certificate: Option<&Path>,
) -> Result<(), Fault> {
    let accepted = ideal
        .groebner_basis_certified(options)
        .map_err(of_certify)?;
    if let Some(path) = certificate {
        write_bytes(path, accepted.certificate())?;
    }
    write_system(out, header, &render_all(accepted.basis()))
}

fn write_reported_basis(
    ideal: sylvester::Ideal<PrimeField>,
    options: ComputeOptions,
    out: &OutputArgs,
    header: &Header<'_>,
) -> Result<(), Fault> {
    let (basis, counters) = ideal
        .groebner_basis_with_report(options)
        .map_err(of_compute)?;
    write_system(out, header, &render_all(&basis))?;
    write_report(&counters)
}

fn rational_basis_command(
    ring: PolynomialRing<Rationals>,
    command: BasisCommand<'_>,
) -> Result<(), Fault> {
    if command.certified {
        return Err(Fault::usage(
            "the rational numbers have no certified path, so --certified and --certificate do not apply",
        ));
    }
    let generators = parse(&ring, command.source, &command.resolved.names)?;
    let ideal = ring.ideal(generators).map_err(of_ring)?;
    let options = rational_options(
        &command.system.engine,
        &command.system.limits,
        command.deadline,
    )?;
    if command.report {
        let (basis, counters) = ideal
            .groebner_basis_with_report(options)
            .map_err(of_compute)?;
        write_system(command.out, command.header, &render_all(&basis))?;
        return write_report(&counters);
    }
    let basis = ideal.groebner_basis(options).map_err(of_compute)?;
    write_system(command.out, command.header, &render_all(&basis))
}

/// What a division command prints.
enum Divide {
    /// The remainder, as a polynomial system.
    Remainder(OutputArgs),
    /// The membership verdict, as `true` or `false`.
    Member(Option<PathBuf>),
}

fn divide_command(started: Instant, divide: DivideArgs, print: Divide) -> Result<(), Fault> {
    let deadline = Deadline::new(started, divide.limits.timeout)?;
    let mut stdin_taken = false;
    let basis_source = load(&divide.basis, divide.basis_format, &mut stdin_taken)?;
    let poly_source = load_polynomial_source(&divide, &mut stdin_taken)?;
    let mut sources = vec![&basis_source];
    if let Some(source) = &poly_source {
        sources.push(source);
    }
    let resolved = resolve(&sources, &divide.ring)?;
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
    };
    match build_ring(&resolved)? {
        AnyRing::Prime(ring) => divide_in_ring(ring, command),
        AnyRing::Rational(ring) => divide_in_ring(ring, command),
    }
}

fn load_polynomial_source(
    divide: &DivideArgs,
    stdin_taken: &mut bool,
) -> Result<Option<Source>, Fault> {
    let Some(path) = &divide.poly_file else {
        return Ok(None);
    };
    let source = load(path, divide.in_format, stdin_taken)?;
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
}

fn divide_in_ring<D: CliDomain>(
    ring: PolynomialRing<D>,
    command: DivideCommand<'_>,
) -> Result<(), Fault>
where
    D::Coeff: std::fmt::Display,
{
    let polynomials = parse(&ring, command.basis_source, &command.resolved.names)?;
    let basis = D::checked_basis(
        &ring,
        polynomials,
        command.deadline.budget(&command.args.limits)?,
    )
    .map_err(of_basis)?;
    let polynomial = one_polynomial(
        &ring,
        command.args,
        command.poly_source,
        &command.resolved.names,
    )?;
    divide_output(
        &basis,
        &polynomial,
        command.deadline,
        command.args,
        command.header,
        command.print,
    )
}

/// The polynomial `normal-form` and `member` divide.
fn one_polynomial<D: Domain>(
    ring: &PolynomialRing<D>,
    divide: &DivideArgs,
    source: Option<&Source>,
    names: &[String],
) -> Result<Polynomial<D>, Fault> {
    match (source, &divide.poly) {
        (Some(source), _) => {
            let mut polynomials = parse(ring, source, names)?;
            Ok(polynomials.remove(0))
        }
        (None, Some(text)) => ring
            .parse_polynomial(text)
            .map_err(|error| Fault::usage(format!("--poly: {error}"))),
        (None, None) => Err(Fault::usage(
            "name the polynomial with --poly or --poly-file",
        )),
    }
}

/// Print the remainder, or the membership verdict.
fn divide_output<D: Domain>(
    basis: &GroebnerBasis<D>,
    f: &Polynomial<D>,
    deadline: &Deadline,
    divide: &DivideArgs,
    header: &Header<'_>,
    print: Divide,
) -> Result<(), Fault>
where
    D::Coeff: std::fmt::Display,
{
    let budget = deadline.budget(&divide.limits)?;
    match print {
        Divide::Remainder(out) => {
            let remainder = basis.normal_form(f, budget).map_err(of_normal_form)?;
            write_system(&out, header, &[format::render(&remainder)])
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
    let (deadline, source, resolved) = prepare_hilbert(started, &args)?;
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
    let source = load(&args.system.input, args.system.in_format, &mut stdin_taken)?;
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
    let polynomials = parse(&ring, source, &resolved.names)?;
    let basis = if args.from_basis {
        GroebnerBasis::<PrimeField>::from_polynomials(
            &ring,
            polynomials,
            deadline.budget(&args.system.limits)?,
        )
        .map_err(of_basis)?
    } else {
        let options = compute_options(&args.system.engine, &args.system.limits, deadline)?;
        ring.ideal(polynomials)
            .map_err(of_ring)?
            .groebner_basis(options)
            .map_err(of_compute)?
    };
    basis
        .hilbert_series(deadline.budget(&args.system.limits)?)
        .map_err(of_hilbert)
}

fn rational_hilbert(
    ring: PolynomialRing<Rationals>,
    source: &Source,
    resolved: &Resolved,
    args: &HilbertArgs,
    deadline: &Deadline,
) -> Result<(HilbertSeries, Option<Established>), Fault> {
    let polynomials = parse(&ring, source, &resolved.names)?;
    let basis = if args.from_basis {
        GroebnerBasis::<Rationals>::from_polynomials(
            &ring,
            polynomials,
            deadline.budget(&args.system.limits)?,
        )
        .map_err(of_basis)?
    } else {
        let options = rational_options(&args.system.engine, &args.system.limits, deadline)?;
        ring.ideal(polynomials)
            .map_err(of_ring)?
            .groebner_basis(options)
            .map_err(of_compute)?
    };
    let established = basis.lift().map(|lift| lift.established);
    let series = basis
        .hilbert_series(deadline.budget(&args.system.limits)?)
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
    let deadline = Deadline::new(started, args.timeout)?;
    let mut limits = verify::Limits::default();
    if let Some(max_bytes) = args.max_bytes {
        limits.max_bytes = max_bytes;
    }
    let bytes = read_capped(&args.certificate, limits.max_bytes)?;
    if let Some(remaining) = deadline.remaining()? {
        limits.deadline = Instant::now().checked_add(remaining);
    }
    let verified = match verify::verify_with_limits(&bytes, &limits) {
        Ok(verified) => verified,
        Err(error) if error.is_exhaustion() => {
            return Err(Fault::limit(format!("the verifier stopped: {error}")));
        }
        Err(error) => {
            return Err(Fault::rejected(format!(
                "the verifier rejected the certificate: {error}"
            )));
        }
    };
    if args.quiet {
        return Ok(());
    }
    let names: Vec<String> = (1..=verified.nvars())
        .map(|index| format!("x{index}"))
        .collect();
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
    write_system(&out, &header, &rendered)?;
    eprintln!(
        "sylv: the certificate carries no variable names, so the basis prints under the synthetic names x1 to x{}",
        verified.nvars()
    );
    Ok(())
}

/// Read at most `max_bytes + 1` bytes, and reject anything longer.
fn read_capped(input: &str, max_bytes: usize) -> Result<Vec<u8>, Fault> {
    let cap = max_bytes.saturating_add(1) as u64;
    let mut bytes = Vec::new();
    if input == "-" {
        io::stdin()
            .take(cap)
            .read_to_end(&mut bytes)
            .map_err(|error| Fault::io(format!("cannot read standard input: {error}")))?;
    } else {
        let file = File::open(input)
            .map_err(|error| Fault::io(format!("cannot read {input}: {error}")))?;
        file.take(cap)
            .read_to_end(&mut bytes)
            .map_err(|error| Fault::io(format!("cannot read {input}: {error}")))?;
    }
    if bytes.len() > max_bytes {
        return Err(Fault::limit(format!(
            "the certificate is longer than {max_bytes} bytes"
        )));
    }
    Ok(bytes)
}

fn render_all<D: Domain>(polynomials: &[Polynomial<D>]) -> Vec<Rendered>
where
    D::Coeff: std::fmt::Display,
{
    polynomials.iter().map(format::render).collect()
}

/// Write a polynomial system to a file or to standard output.
fn write_system(
    out: &OutputArgs,
    header: &Header<'_>,
    polynomials: &[Rendered],
) -> Result<(), Fault> {
    match &out.output {
        Some(path) => {
            let mut file = File::create(path)
                .map_err(|error| Fault::io(format!("cannot write {}: {error}", path.display())))?;
            format::write(out.out_format, header, polynomials, &mut file)
                .map_err(|error| Fault::of_io(&format!("cannot write {}", path.display()), &error))
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

/// Write plain text to a file or to standard output.
fn write_text(path: Option<&Path>, text: &str) -> Result<(), Fault> {
    match path {
        Some(path) => fs::write(path, text)
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
    fs::write(path, bytes)
        .map_err(|error| Fault::io(format!("cannot write {}: {error}", path.display())))
}

/// Write the counters of a run to standard error, one per line.
fn write_report(report: &ComputeReport) -> Result<(), Fault> {
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
    let stderr = io::stderr();
    let mut handle = stderr.lock();
    for line in lines {
        writeln!(handle, "{line}")
            .map_err(|error| Fault::of_io("cannot write standard error", &error))?;
    }
    Ok(())
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
