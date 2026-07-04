#!/usr/bin/env bash
set -euo pipefail

echo "=== no-default-features ==="
cargo check --no-default-features

echo "=== swisseph-files only ==="
cargo check --no-default-features --features swisseph-files

echo "=== jpl only ==="
cargo check --no-default-features --features jpl

echo "=== serde ==="
cargo check --features serde

echo "=== all-features ==="
cargo check --all-features

echo "=== lib tests (no-default-features) ==="
cargo test --lib --no-default-features

echo "=== full test suite (default features) ==="
cargo test

echo "=== clippy (all-features) ==="
cargo clippy --all-features -- -D warnings

echo "All feature combinations passed."
