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

# --- Upstream PREPARE_GENOME convenience branches -------------------------------------
# nf-core/atacseq 2.1.2 PREPARE_GENOME auto-decompresses .gz references and
# auto-unpacks .tar.gz index archives before consumers touch them (GUNZIP,
# GFFREAD, UNTAR, GTF2BED, TSS_EXTRACT). The port implements these as producer
# rules whose outputs are the canonical uncompressed paths the consumers read;
# consumers wire them via depends_on (a closed gate counts as satisfied, so
# the default pre-built path is unchanged) and resolve the produced file at
# shell level (khmer-style). Each branch-flip below activates one branch via a
# sed'd temp copy and asserts that (a) the convenience rule schedules and
# (b) the consumers still schedule.

echo "==> dry-run with a .gz GTF source activates the gtf gunzip branch"
# Upstream GUNZIP_GTF (prepare_genome.nf): when params.gtf ends with ".gz",
# decompress before HOMER (and any other GTF consumer) reads it. The gunzip
# output is the archive minus the ".gz" suffix — the canonical {config.gtf}
# path the port's consumers already use. The port keeps the canonical consumer
# key (config.gtf) separate from the compressed source slot (config.gtf_src):
# set gtf_src to your .gz, and ref::gunzip_gtf unzips to config.gtf.
sed 's|^gtf_src = ""|gtf_src = "test/fixtures/genome/genes.gtf.gz"|' main.oxoflow > .refgunzip-test-tmp.oxoflow
grep -q '^gtf_src = "test/fixtures/genome/genes.gtf.gz"' .refgunzip-test-tmp.oxoflow  # the sed must have matched
"$OXO" dry-run .refgunzip-test-tmp.oxoflow > /tmp/oxo-dryrun-refgunzip-$$.txt 2>&1
grep -Eq 'gunzip_gtf  \[run:' /tmp/oxo-dryrun-refgunzip-$$.txt
grep -Eq 'homer_annotatepeaks_S1  \[run:' /tmp/oxo-dryrun-refgunzip-$$.txt
# homer's GTF selection picks the gunzipped canonical path (khmer-style
# fallback: config.gtf prebuilt vs. gffread-derived ref.gtf).
grep -Eq 'GTF="test/fixtures/genome/genes\.gtf"' /tmp/oxo-dryrun-refgunzip-$$.txt
rm -f .refgunzip-test-tmp.oxoflow

echo "==> dry-run with GFF3 annotation and no GTF activates the gffread branch"
# Upstream GFFREAD (prepare_genome.nf): the else-if after `if (params.gtf)` —
# when params.gtf is unset but params.gff is set, convert the GFF3 to GTF
# (--keep-exon-attrs -F -T, gffread 0.12.1). The port gates ref::gffread on
# config.gff != "" && config.gtf == "" — flip config.gtf to "" (the fixture
# always ships a pre-built GTF, which routes the else-if the other way).
sed 's|^gtf = "test/fixtures/genome/genes.gtf"|gtf = ""|; s|^gff = ""|gff = "test/fixtures/genome/genes.gff"|' main.oxoflow > .refgffread-test-tmp.oxoflow
grep -q '^gff = "test/fixtures/genome/genes.gff"' .refgffread-test-tmp.oxoflow  # the sed must have matched
grep -q '^gtf = ""' .refgffread-test-tmp.oxoflow
"$OXO" dry-run .refgffread-test-tmp.oxoflow > /tmp/oxo-dryrun-refgffread-$$.txt 2>&1
grep -Eq 'gffread  \[run:' /tmp/oxo-dryrun-refgffread-$$.txt
grep -Eq 'homer_annotatepeaks_S1  \[run:' /tmp/oxo-dryrun-refgffread-$$.txt
rm -f .refgffread-test-tmp.oxoflow

echo "==> dry-run with a .gz gene/tss BED source activates the bed gunzip branches"
# Upstream GUNZIP_GENE_BED / GUNZIP_TSS_BED (prepare_genome.nf): decompress
# before deepTools reads the regions files. Same source-slot pattern as above.
sed 's|^gene_bed_src = ""|gene_bed_src = "test/fixtures/genome/gene.bed.gz"|; s|^tss_bed_src = ""|tss_bed_src = "test/fixtures/genome/tss.bed.gz"|' main.oxoflow > .refbedgunzip-test-tmp.oxoflow
grep -q '^gene_bed_src = "test/fixtures/genome/gene.bed.gz"' .refbedgunzip-test-tmp.oxoflow  # the sed must have matched
grep -q '^tss_bed_src = "test/fixtures/genome/tss.bed.gz"' .refbedgunzip-test-tmp.oxoflow
"$OXO" dry-run .refbedgunzip-test-tmp.oxoflow > /tmp/oxo-dryrun-refbedgunzip-$$.txt 2>&1
grep -Eq 'gunzip_gene_bed  \[run:' /tmp/oxo-dryrun-refbedgunzip-$$.txt
grep -Eq 'gunzip_tss_bed  \[run:' /tmp/oxo-dryrun-refbedgunzip-$$.txt
grep -Eq 'deeptools_plots_S1  \[run:' /tmp/oxo-dryrun-refbedgunzip-$$.txt
rm -f .refbedgunzip-test-tmp.oxoflow

echo "==> dry-run with a .tar.gz BWA index source activates the untar branch"
# Upstream UNTAR_BWA_INDEX (prepare_genome.nf): when params.bwa_index ends with
# ".tar.gz", unpack it before BWA-MEM reads the index; the archive minus its
# ".tar.gz" suffix is the index directory. The untar output (the index prefix
# written at shell level via the find -name "*.amb" idiom) feeds bwa_mem:
# the assert below pins that bwa_mem still schedules.
sed 's|^bwa_index_src = ""|bwa_index_src = "test/fixtures/genome/genome-bwa-index.tar.gz"|' main.oxoflow > .refuntar-test-tmp.oxoflow
grep -q '^bwa_index_src = "test/fixtures/genome/genome-bwa-index.tar.gz"' .refuntar-test-tmp.oxoflow  # the sed must have matched
"$OXO" dry-run .refuntar-test-tmp.oxoflow > /tmp/oxo-dryrun-refuntar-$$.txt 2>&1
grep -Eq 'untar_bwa_index  \[run:' /tmp/oxo-dryrun-refuntar-$$.txt
grep -Eq 'bwa_mem_cohort_S1  \[run:' /tmp/oxo-dryrun-refuntar-$$.txt
rm -f .refuntar-test-tmp.oxoflow

echo "PASS"
