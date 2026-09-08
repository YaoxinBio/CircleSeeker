#!/usr/bin/env bash
# Compatibility entry point: conda build performs package tests by default.
set -euo pipefail
RECIPE_DIR="$(cd "$(dirname "$0")" && pwd)"
exec bash "$RECIPE_DIR/build_and_upload.sh" local
