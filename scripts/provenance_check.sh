#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
terms=$(grep -v '^#' scripts/provenance_terms.txt | grep -v '^\s*$' | paste -sd'|' -)
if grep -rn -i -E "($terms)" --exclude-dir=.git --exclude-dir=.venv --exclude=provenance_terms.txt . ; then echo "provenance_check: FAIL"; exit 1; fi
echo "provenance_check: OK"
