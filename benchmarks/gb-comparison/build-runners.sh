#!/usr/bin/env bash
# Build the sylvester runner and the optional msolve timer.
set -euo pipefail
cd "$(dirname "$0")"

CARGO_TARGET_DIR="runner/target-default" cargo +1.92 build --release \
    --manifest-path runner/Cargo.toml
mkdir -p runner/bin
cp "runner/target-default/release/sylv-runner" "runner/bin/sylv-runner"

if pkg-config --exists msolve; then
    make -C runner/msolve
    mkdir -p runner/bin
    cp runner/msolve/msolve-inproc runner/bin/msolve-inproc
else
    echo "build-runners.sh: msolve not found via pkg-config, skipping msolve-inproc" >&2
fi
