#!/usr/bin/env python3
"""Remap build paths and scan release binaries for home directories."""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from zipfile import BadZipFile, ZipFile

UNIT_SEPARATOR = "\x1f"
UNIX_HOME_NAMES = ("home", "Users", "root")
ASCII_HOME_PATTERNS = tuple(
    re.compile(("/" + name + "/").encode("ascii")) for name in UNIX_HOME_NAMES
)
ASCII_HOME_PATTERNS += (re.compile(rb"[A-Za-z]:[\\/]+Users[\\/]"),)
TEXT_HOME_PATTERNS = tuple(re.compile("/" + name + "/") for name in UNIX_HOME_NAMES)
TEXT_HOME_PATTERNS += (re.compile(r"[A-Za-z]:[\\/]+Users[\\/]"),)


def native_path(raw: str) -> Path | None:
    """Resolve one environment path, including a Git Bash drive path."""
    if not raw:
        return None
    if (
        os.name == "nt"
        and len(raw) >= 3
        and raw[0] == "/"
        and raw[1].isalpha()
        and raw[2] in "/\\"
    ):
        raw = f"{raw[1]}:{raw[2:]}"
    candidate = Path(raw).expanduser()
    try:
        resolved = candidate.resolve()
    except OSError:
        resolved = candidate.absolute()
    if not resolved.is_absolute() or resolved == Path(resolved.anchor):
        return None
    return resolved


def remap_roots() -> list[tuple[Path, str]]:
    """Return the build roots that need stable non-user prefixes."""
    script_root = Path(__file__).resolve().parents[1]
    workspace = native_path(os.environ.get("GITHUB_WORKSPACE", ""))
    if workspace is None:
        workspace = native_path(str(script_root))

    cargo_home = native_path(os.environ.get("CARGO_HOME", ""))
    if cargo_home is None:
        cargo_home = native_path(str(Path.home() / ".cargo"))

    build_home = native_path(str(Path.home()))
    roots: list[tuple[Path, str]] = []
    seen: set[str] = set()
    for root, replacement in (
        (build_home, "/build-home"),
        (cargo_home, "/cargo-home"),
        (workspace, "/workspace"),
    ):
        if root is None:
            continue
        prefix = str(root)
        if not prefix or prefix in seen:
            continue
        if any(char in prefix for char in "\x00\x1f\r\n"):
            raise ValueError("a path prefix contains a control character")
        seen.add(prefix)
        roots.append((root, replacement))
    if not roots:
        raise ValueError("no absolute build path is available for remapping")
    return roots


def existing_rust_flags() -> list[str]:
    """Read the flags Cargo would use before adding release remaps."""
    encoded = os.environ.get("CARGO_ENCODED_RUSTFLAGS")
    if encoded is not None:
        if "\x00" in encoded:
            raise ValueError("CARGO_ENCODED_RUSTFLAGS contains a NUL")
        return encoded.split(UNIT_SEPARATOR) if encoded else []

    raw = os.environ.get("RUSTFLAGS", "")
    if not raw:
        return []
    return raw.split()


def encoded_flags() -> str:
    """Build encoded Rust flags without losing the caller's flags."""
    flags = existing_rust_flags()
    for root, replacement in remap_roots():
        flags.append(f"--remap-path-prefix={root}={replacement}")
    encoded = UNIT_SEPARATOR.join(flags)
    if any(char in encoded for char in "\x00\r\n"):
        raise ValueError("Rust flags contain a control character")
    return encoded


def emit_flags() -> int:
    """Write encoded flags for a caller that exports them."""
    sys.stdout.write(encoded_flags())
    return 0


def emit_github_environment() -> int:
    """Write the encoded flags in GitHub Actions environment-file syntax."""
    sys.stdout.write(f"CARGO_ENCODED_RUSTFLAGS={encoded_flags()}\n")
    return 0


def text_home_match(text: str) -> tuple[int, str] | None:
    """Find a Unix or Windows home path in decoded binary text."""
    for pattern in TEXT_HOME_PATTERNS:
        match = pattern.search(text)
        if match is not None:
            return match.start(), match.group()
    return None


def binary_home_matches(data: bytes) -> list[str]:
    """Find ASCII and UTF-16 Unix or Windows home paths in bytes."""
    matches: list[str] = []
    for pattern in ASCII_HOME_PATTERNS:
        match = pattern.search(data)
        if match is not None:
            value = match.group().decode("ascii", errors="replace")
            matches.append(f"ASCII offset {match.start()}: {value}")

    for encoding in ("utf-16-le", "utf-16-be"):
        for byte_offset in (0, 1):
            text = data[byte_offset:].decode(encoding, errors="ignore")
            match = text_home_match(text)
            if match is not None:
                char_offset, value = match
                offset = byte_offset + (char_offset * 2)
                matches.append(f"{encoding} offset {offset}: {value}")
    return matches


def check_binary(path: Path) -> int:
    """Reject one binary that retains a Unix or Windows home path."""
    data = path.read_bytes()
    matches = binary_home_matches(data)
    if matches:
        for match in matches:
            print(f"{path}: {match}", file=sys.stderr)
        return 1
    print(f"{path}: no home paths")
    return 0


def check_wheel(path: Path) -> int:
    """Reject wheel members that retain a Unix or Windows home path."""
    found = False
    with ZipFile(path) as archive:
        for member in archive.infolist():
            matches = binary_home_matches(archive.read(member))
            for match in matches:
                print(f"{path}!{member.filename}: {match}", file=sys.stderr)
                found = True
    if found:
        return 1
    print(f"{path}: no home paths")
    return 0


def parser() -> argparse.ArgumentParser:
    """Create the command line parser."""
    command_parser = argparse.ArgumentParser(description=__doc__)
    subparsers = command_parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("flags", help="print encoded Rust flags")
    subparsers.add_parser("github-env", help="print a GitHub environment entry")
    binary_parser = subparsers.add_parser("check-binary", help="scan binary files")
    binary_parser.add_argument("paths", nargs="+")
    wheel_parser = subparsers.add_parser("check-wheel", help="scan wheel members")
    wheel_parser.add_argument("paths", nargs="+")
    return command_parser


def main(argv: list[str] | None = None) -> int:
    """Run one path remap or binary scan command."""
    args = parser().parse_args(argv)
    if args.command == "flags":
        return emit_flags()
    if args.command == "github-env":
        return emit_github_environment()
    try:
        if args.command == "check-binary":
            return max(check_binary(Path(path)) for path in args.paths)
        return max(check_wheel(Path(path)) for path in args.paths)
    except (BadZipFile, OSError, ValueError) as error:
        print(f"release-paths: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
