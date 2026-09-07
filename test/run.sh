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

echo "==> dry-run with replicate samples instantiates the merged-replicate chain"
# sample_groups replacement + merged_samples config (see README
# 'Merged-replicate analysis'); the JSON plan lists instance names.
# Instance naming differs by engine: input_groups engines name instances
# `merge_replicates_S1`, released engines (sample-group fan-out) name them
# `merge_replicates_cohort_S1_REP1` — match both spellings.
"$OXO" dry-run main.oxoflow --json --samples @test/fixtures/samples_replicates.csv --arg merged_samples=S1,S2 > /tmp/oxo-dryrun-reps-$$.json 2>&1
grep -Eq '"name": "merge_replicates(_cohort)?_(S1|S1_REP1)"' /tmp/oxo-dryrun-reps-$$.json
grep -Eq '"name": "merge_replicates(_cohort)?_(S2|S2_REP1)"' /tmp/oxo-dryrun-reps-$$.json
grep -Eq '"name": "macs2_callpeak(_cohort)?_(S1|S1_REP1)"' /tmp/oxo-dryrun-reps-$$.json
grep -Eq '"name": "ucsc_bedgraphtobigwig(_cohort)?_(S1|S1_REP1)"' /tmp/oxo-dryrun-reps-$$.json

echo "==> debug: expanded commands contain no literal {wildcards}"
# debug prints the expanded plan to stderr too (stdout is machine output).
# WARN lines are tracing log lines, not commands — an input_groups no-op
# legitimately quotes its pattern (e.g. merge_replicates on a default
# config) — so strip ANSI and the timestamp prefix before checking.
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1 || true
STRIPPED=/tmp/oxo-debug-stripped-$$.txt
sed 's/\x1b\[[0-9;]*m//g' /tmp/oxo-debug-$$.txt | grep -v '^[0-9]\{4\}-[0-9]\{2\}-[0-9]\{2\}T' > "$STRIPPED"
if grep -q '{sample}' "$STRIPPED" || grep -q '{config\.' "$STRIPPED"; then
    echo "unexpanded wildcards in debug output"
    exit 1
fi

echo "==> dry-run with macs_gsize=\"\" activates the khmer auto-estimate branch"
# Upstream KHMER_UNIQUEKMERS: an empty params.macs_gsize routes MACS2_CALLPEAK's
# --gsize through a khmer estimate at params.read_length. The port gates
# ref::khmer_uniquekmers on config.macs_gsize == '' — override it via a sed'd
# temp copy (the engine cannot pass empty-string --arg overrides; the copy
# must live in this directory because [[include]] paths resolve relative to
# the workflow file).
sed 's/^macs_gsize = "2.7e9"/macs_gsize = ""/' main.oxoflow > .khmer-test-tmp.oxoflow
grep -q '^macs_gsize = ""' .khmer-test-tmp.oxoflow  # the sed must have matched
"$OXO" dry-run .khmer-test-tmp.oxoflow > /tmp/oxo-dryrun-khmer-$$.txt 2>&1
grep -Eq 'khmer_uniquekmers  \[run:' /tmp/oxo-dryrun-khmer-$$.txt
grep -Eq 'macs2_callpeak_S1  \[run:' /tmp/oxo-dryrun-khmer-$$.txt
grep -Eq 'GSIZE=\$\(cat results/genome/kmers\.txt\)' /tmp/oxo-dryrun-khmer-$$.txt
rm -f .khmer-test-tmp.oxoflow

echo "PASS"
