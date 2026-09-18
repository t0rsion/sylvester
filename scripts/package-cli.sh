#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "usage: $0 TAG TARGET OUTPUT_DIR" >&2
    exit 2
fi

tag=$1
target=$2
output_dir=$3

if [[ ! $tag =~ ^v[0-9]+\.[0-9]+\.[0-9]+([.-][0-9A-Za-z.-]+)?$ ]]; then
    echo "invalid release tag: $tag" >&2
    exit 2
fi
if [[ ! $target =~ ^[a-z0-9_]+(-[a-z0-9_]+)+$ ]]; then
    echo "invalid Rust target: $target" >&2
    exit 2
fi

case "$target" in
    *-pc-windows-*)
        binary_name=sylv.exe
        archive_suffix=zip
        ;;
    *)
        binary_name=sylv
        archive_suffix=tar.gz
        ;;
esac

binary_path="target/$target/release/$binary_name"
if [[ ! -f $binary_path ]]; then
    echo "missing release binary: $binary_path" >&2
    exit 1
fi

mkdir -p "$output_dir"
output_dir=$(cd "$output_dir" && pwd)
archive_stem="sylv-$tag-$target"
archive_path="$output_dir/$archive_stem.$archive_suffix"
staging_root=$(mktemp -d)
trap 'rm -rf "$staging_root"' EXIT
staging_dir="$staging_root/$archive_stem"
mkdir -p "$staging_dir"
cp "$binary_path" "$staging_dir/$binary_name"
cp crates/sylvester-cli/README.md LICENSE-APACHE LICENSE-MIT "$staging_dir/"

if [[ $archive_suffix == zip ]]; then
    command -v 7z >/dev/null 2>&1 || {
        echo "7z is required to create Windows archives" >&2
        exit 1
    }
    (cd "$staging_root" && 7z a -tzip "$archive_path" "$archive_stem" >/dev/null)
else
    tar -C "$staging_root" -czf "$archive_path" "$archive_stem"
fi

printf '%s\n' "$archive_path"
