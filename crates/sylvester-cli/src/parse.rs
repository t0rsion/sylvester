//! Parsers for command-line values that need more than `clap`'s defaults.

use std::num::IntErrorKind;

/// Parse a byte count with an optional binary or decimal suffix.
///
/// Bare decimal integers remain valid. Binary suffixes use `KiB`, `MiB`,
/// `GiB`, or `TiB`; decimal suffixes use `KB`, `MB`, `GB`, or `TB`.
pub(crate) fn bytes(text: &str) -> Result<usize, String> {
    let text = text.trim();
    let split = text
        .find(|character: char| !character.is_ascii_digit())
        .unwrap_or(text.len());
    let (digits, suffix) = text.split_at(split);
    if digits.is_empty() {
        return Err(byte_error(text));
    }
    let number = digits.parse::<u128>().map_err(|error| match error.kind() {
        IntErrorKind::PosOverflow => format!("{text} is larger than the supported byte count"),
        _ => byte_error(text),
    })?;
    let multiplier = match suffix.to_ascii_lowercase().as_str() {
        "" | "b" => 1,
        "kib" => 1 << 10,
        "mib" => 1 << 20,
        "gib" => 1 << 30,
        "tib" => 1 << 40,
        "kb" => 1_000,
        "mb" => 1_000_000,
        "gb" => 1_000_000_000,
        "tb" => 1_000_000_000_000,
        _ => return Err(byte_error(text)),
    };
    let value = number
        .checked_mul(multiplier)
        .and_then(|value| usize::try_from(value).ok())
        .ok_or_else(|| format!("{text} is larger than the supported byte count"))?;
    Ok(value)
}

fn byte_error(text: &str) -> String {
    format!("{text:?} is not a byte count; use an integer or a suffix such as 512MiB")
}

#[cfg(test)]
mod tests {
    use super::bytes;

    #[test]
    fn accepts_bare_integer_and_binary_suffixes() {
        assert_eq!(bytes("67108864"), Ok(67_108_864));
        assert_eq!(bytes("512MiB"), Ok(512 * 1024 * 1024));
        assert_eq!(bytes("2gib"), Ok(2 * 1024 * 1024 * 1024));
    }

    #[test]
    fn accepts_decimal_suffixes() {
        assert_eq!(bytes("2MB"), Ok(2_000_000));
        assert_eq!(bytes("4kb"), Ok(4_000));
    }

    #[test]
    fn rejects_unknown_and_overflowing_values() {
        assert!(bytes("4blocks").is_err());
        assert!(bytes("18446744073709551616").is_err());
        assert!(bytes("-1").is_err());
    }
}
