//! Bounded reads and atomic file writes for the command line.

use std::ffi::OsString;
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::{Component, Path, PathBuf};
use std::process;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const READ_CHUNK: usize = 16 * 1024;
const TEMP_ATTEMPTS: usize = 100;
static NEXT_TEMP: AtomicU64 = AtomicU64::new(0);

/// The limits a CLI read checks between chunks.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct ReadLimits {
    /// The absolute command deadline.
    pub(crate) deadline: Option<Instant>,
    /// The remaining live byte allowance for the returned buffer.
    pub(crate) memory: Option<usize>,
    /// The maximum source size, including one byte used to detect overflow.
    pub(crate) max_bytes: Option<usize>,
}

impl ReadLimits {
    /// Build limits for one bounded read.
    pub(crate) const fn new(
        deadline: Option<Instant>,
        memory: Option<usize>,
        max_bytes: Option<usize>,
    ) -> Self {
        ReadLimits {
            deadline,
            memory,
            max_bytes,
        }
    }
}

/// Why a bounded read stopped.
#[derive(Debug)]
pub(crate) enum ReadError {
    /// The operating system rejected a read.
    Io(io::Error),
    /// The absolute deadline passed.
    Timeout,
    /// The returned buffer would pass its byte allowance.
    MemoryLimitExceeded,
    /// The source contains more than the requested maximum.
    TooLong { limit: usize },
    /// The source is not UTF-8.
    InvalidUtf8,
}

impl std::fmt::Display for ReadError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadError::Io(error) => error.fmt(f),
            ReadError::Timeout => f.write_str("the read passed its deadline"),
            ReadError::MemoryLimitExceeded => f.write_str("the read passed its memory limit"),
            ReadError::TooLong { limit } => write!(f, "the source is longer than {limit} bytes"),
            ReadError::InvalidUtf8 => f.write_str("the source is not valid UTF-8"),
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReadError::Io(error) => Some(error),
            ReadError::Timeout
            | ReadError::MemoryLimitExceeded
            | ReadError::TooLong { .. }
            | ReadError::InvalidUtf8 => None,
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(error: io::Error) -> Self {
        ReadError::Io(error)
    }
}

/// Read bytes while checking the same limits throughout the read.
///
/// A blocking reader can still wait inside one operating system read. The
/// deadline is checked before and after each read without a platform API.
pub(crate) fn read_bytes<R: Read>(mut reader: R, limits: ReadLimits) -> Result<Vec<u8>, ReadError> {
    let mut bytes = Vec::new();
    let mut chunk = [0u8; READ_CHUNK];
    loop {
        check_limits(&bytes, &limits)?;
        let size = read_size(bytes.len(), &limits);
        let count = reader.read(&mut chunk[..size])?;
        check_deadline(&limits)?;
        if count == 0 {
            return Ok(bytes);
        }
        append_chunk(&mut bytes, &chunk[..count], &limits)?;
    }
}

/// Read a UTF-8 source while checking the same limits throughout the read.
pub(crate) fn read_text<R: Read>(reader: R, limits: ReadLimits) -> Result<String, ReadError> {
    let bytes = read_bytes(reader, limits)?;
    check_deadline(&limits)?;
    String::from_utf8(bytes).map_err(|_| ReadError::InvalidUtf8)
}

/// Open and read a bounded byte source from a path.
pub(crate) fn read_path_bytes(path: &Path, limits: ReadLimits) -> Result<Vec<u8>, ReadError> {
    let file = File::open(path)?;
    check_deadline(&limits)?;
    reject_metadata_size(&file, &limits)?;
    read_bytes(file, limits)
}

/// Open and read a bounded UTF-8 source from a path.
pub(crate) fn read_path_text(path: &Path, limits: ReadLimits) -> Result<String, ReadError> {
    let bytes = read_path_bytes(path, limits)?;
    check_deadline(&limits)?;
    String::from_utf8(bytes).map_err(|_| ReadError::InvalidUtf8)
}

/// Reject certificate and result paths that resolve to one target.
pub(crate) fn reject_same_targets(
    certificate: Option<&Path>,
    result: Option<&Path>,
) -> Result<(), TargetConflict> {
    let (Some(certificate), Some(result)) = (certificate, result) else {
        return Ok(());
    };
    if same_target(certificate, result) {
        return Err(TargetConflict {
            certificate: certificate.to_path_buf(),
            result: result.to_path_buf(),
        });
    }
    Ok(())
}

/// Two output paths that name one file.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct TargetConflict {
    /// The certificate path supplied by the caller.
    pub(crate) certificate: PathBuf,
    /// The result path supplied by the caller.
    pub(crate) result: PathBuf,
}

impl std::fmt::Display for TargetConflict {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "certificate path {} and result path {} name the same file",
            self.certificate.display(),
            self.result.display()
        )
    }
}

impl std::error::Error for TargetConflict {}

/// Write a file through a same-directory temporary file.
pub(crate) fn write_atomic(path: &Path, bytes: &[u8]) -> io::Result<()> {
    write_atomic_with(path, |file| file.write_all(bytes))
}

/// Render a file through a same-directory temporary file.
///
/// The final rename follows the platform's replacement semantics. The helper
/// never removes the destination first, so a failed replacement keeps it.
pub(crate) fn write_atomic_with<F>(path: &Path, write: F) -> io::Result<()>
where
    F: FnOnce(&mut File) -> io::Result<()>,
{
    let (temporary_path, mut temporary) = create_temporary(path)?;
    let result = write_temporary(&mut temporary, write);
    drop(temporary);
    if let Err(error) = result {
        remove_temporary(&temporary_path);
        return Err(error);
    }
    if let Err(error) = rename_replacing(&temporary_path, path) {
        remove_temporary(&temporary_path);
        return Err(error);
    }
    Ok(())
}

fn check_limits(bytes: &[u8], limits: &ReadLimits) -> Result<(), ReadError> {
    check_deadline(limits)?;
    check_size(bytes.len(), limits)
}

fn check_deadline(limits: &ReadLimits) -> Result<(), ReadError> {
    if limits
        .deadline
        .is_some_and(|deadline| Instant::now() >= deadline)
    {
        return Err(ReadError::Timeout);
    }
    Ok(())
}

fn check_size(length: usize, limits: &ReadLimits) -> Result<(), ReadError> {
    if let Some(limit) = limits.max_bytes
        && length > limit
    {
        return Err(ReadError::TooLong { limit });
    }
    if limits.memory.is_some_and(|limit| length > limit) {
        return Err(ReadError::MemoryLimitExceeded);
    }
    Ok(())
}

fn read_size(length: usize, limits: &ReadLimits) -> usize {
    let mut size = READ_CHUNK;
    if let Some(limit) = limits.max_bytes {
        size = size.min(limit.saturating_sub(length).saturating_add(1));
    }
    if let Some(limit) = limits.memory {
        size = size.min(limit.saturating_sub(length).saturating_add(1));
    }
    size.max(1)
}

fn append_chunk(bytes: &mut Vec<u8>, chunk: &[u8], limits: &ReadLimits) -> Result<(), ReadError> {
    let length = bytes
        .len()
        .checked_add(chunk.len())
        .ok_or(ReadError::MemoryLimitExceeded)?;
    check_size(length, limits)?;
    bytes
        .try_reserve_exact(chunk.len())
        .map_err(|_| ReadError::MemoryLimitExceeded)?;
    if bytes.capacity() > limits.memory.unwrap_or(usize::MAX) {
        return Err(ReadError::MemoryLimitExceeded);
    }
    bytes.extend_from_slice(chunk);
    Ok(())
}

fn reject_metadata_size(file: &File, limits: &ReadLimits) -> Result<(), ReadError> {
    let length = usize::try_from(file.metadata()?.len()).unwrap_or(usize::MAX);
    check_size(length, limits)
}

fn same_target(first: &Path, second: &Path) -> bool {
    let first_lexical = lexical_path(first);
    let second_lexical = lexical_path(second);
    if first_lexical == second_lexical {
        return true;
    }
    #[cfg(windows)]
    if first_lexical
        .to_string_lossy()
        .eq_ignore_ascii_case(&second_lexical.to_string_lossy())
    {
        return true;
    }
    #[cfg(unix)]
    if same_unix_file(first, second) {
        return true;
    }
    match (canonical_path(first), canonical_path(second)) {
        (Ok(first), Ok(second)) => first == second,
        _ => false,
    }
}

#[cfg(unix)]
fn same_unix_file(first: &Path, second: &Path) -> bool {
    use std::os::unix::fs::MetadataExt;

    match (fs::metadata(first), fs::metadata(second)) {
        (Ok(first), Ok(second)) => {
            first.is_file()
                && second.is_file()
                && first.dev() == second.dev()
                && first.ino() == second.ino()
        }
        _ => false,
    }
}

fn canonical_path(path: &Path) -> io::Result<PathBuf> {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|directory| directory.join(path))
            .unwrap_or_else(|_| path.to_path_buf())
    };
    let mut current = absolute.as_path();
    let mut tail = Vec::new();
    loop {
        match fs::canonicalize(current) {
            Ok(mut base) => {
                for component in tail.iter().rev() {
                    base.push(component);
                }
                return Ok(base);
            }
            Err(error) if !current.as_os_str().is_empty() => {
                let name = current.file_name().ok_or(error)?.to_os_string();
                tail.push(name);
                current = current.parent().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::NotFound, "the path has no existing parent")
                })?;
            }
            Err(error) => return Err(error),
        }
    }
}

fn lexical_path(path: &Path) -> PathBuf {
    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map(|directory| directory.join(path))
            .unwrap_or_else(|_| path.to_path_buf())
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            Component::CurDir => {}
            Component::ParentDir => {
                normalized.pop();
            }
            Component::Prefix(_) | Component::RootDir | Component::Normal(_) => {
                normalized.push(component.as_os_str());
            }
        }
    }
    normalized
}

fn create_temporary(path: &Path) -> io::Result<(PathBuf, File)> {
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let name = path.file_name().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "the output path has no file name",
        )
    })?;
    for _ in 0..TEMP_ATTEMPTS {
        let mut temporary_name = OsString::from(".");
        temporary_name.push(name);
        temporary_name.push(format!(
            ".sylv-{}-{}.tmp",
            process::id(),
            NEXT_TEMP.fetch_add(1, Ordering::Relaxed)
        ));
        let temporary_path = parent.join(temporary_name);
        let result = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary_path);
        match result {
            Ok(file) => match copy_permissions(path, &file) {
                Ok(()) => return Ok((temporary_path, file)),
                Err(error) => {
                    drop(file);
                    remove_temporary(&temporary_path);
                    return Err(error);
                }
            },
            Err(error) if error.kind() == io::ErrorKind::AlreadyExists => continue,
            Err(error) => return Err(error),
        }
    }
    Err(io::Error::new(
        io::ErrorKind::AlreadyExists,
        "could not create a unique temporary output",
    ))
}

fn copy_permissions(path: &Path, temporary: &File) -> io::Result<()> {
    let Ok(metadata) = fs::metadata(path) else {
        return Ok(());
    };
    temporary.set_permissions(metadata.permissions())
}

fn write_temporary<F>(temporary: &mut File, write: F) -> io::Result<()>
where
    F: FnOnce(&mut File) -> io::Result<()>,
{
    write(temporary)?;
    temporary.flush()?;
    temporary.sync_all()
}

fn remove_temporary(path: &Path) {
    let _ = fs::remove_file(path);
}

fn rename_replacing(temporary: &Path, destination: &Path) -> io::Result<()> {
    fs::rename(temporary, destination)
}
