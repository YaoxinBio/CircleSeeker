#!/usr/bin/env bash
# Build and test locally. Bioconda publication proceeds through a recipe PR.
set -euo pipefail
RECIPE_DIR="$(cd "$(dirname "$0")" && pwd)"
if [[ "${1:-local}" != "local" ]]; then
  echo "This script builds locally only; publish to Bioconda through a recipe PR." >&2
  exit 2
fi
if ! command -v conda-build >/dev/null 2>&1; then
  echo "Install conda-build in a packaging environment before running this script." >&2
  exit 1
fi
conda build "$RECIPE_DIR" -c conda-forge -c bioconda --python=3.12 --no-anaconda-upload
conda build "$RECIPE_DIR" --output -c conda-forge -c bioconda --python=3.12
