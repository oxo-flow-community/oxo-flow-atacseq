#!/usr/bin/env bash
# Acceptance test for the oxo-flow-atacseq port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
# The plan goes to stderr (stdout is reserved for machine output), so both
# streams are captured here. --samples is not a dry-run flag; sample
# selection comes from the workflow's [[sample_groups]] declaration.
"$OXO" dry-run main.oxoflow > /tmp/oxo-dryrun-$$.txt 2>&1
grep -Eq "would execute|To execute:" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards}"
# debug prints the expanded plan to stderr too (stdout is machine output)
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1 || true
if grep -q '{sample}' /tmp/oxo-debug-$$.txt || grep -q '{config\.' /tmp/oxo-debug-$$.txt; then
    echo "unexpanded wildcards in debug output"
    exit 1
fi

echo "PASS"
