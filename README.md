# oxo-flow-atacseq — ATAC-seq: peak calling and QC

[![CI](https://github.com/oxo-flow-community/oxo-flow-atacseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-atacseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`nf-core/atacseq`](https://github.com/nf-core/atacseq) @ `2.1.2` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

ATAC-seq (Assay for Transposase-Accessible Chromatin) analysis for a
cohort of samples: raw-read QC (FastQC), adapter trimming (Trim Galore),
BWA-MEM alignment, Picard mark-duplicates, BAMTools filtering, MACS2
broad-peak calling, HOMER peak annotation, FRiP scoring, normalised
bigWig tracks, deepTools QC plots (profile, heatmap and fingerprint) and
a single MultiQC report. Input is a directory of `<sample>.fastq.gz`
reads plus pre-built reference files; output lands in `results/` with
per-step subdirectories and a combined MultiQC HTML report.

The default plan is the upstream single-end `aligner=bwa` main path
(15 rules, 29 rule instances for the fixture cohort). Upstream branches
that are off the main path — paired-end, alternative aligners
(Bowtie2/Chromap/STAR), reference preparation, mitochondrial filtering,
consensus peaks / DESeq2, preseq, Picard metrics, ataqv, IGV, R QC
plots — are ported as **when-gated branches** (see
[Gated branches](#gated-branches)); each activates only when its config
key is set, so the default plan is unchanged.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Recommended — prebuilt release binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternative via conda:

```bash
conda install -c bioconda oxo-flow-cli
```

Note the bioconda package may lag the latest release; other platform
binaries are on the [releases page](https://github.com/Traitome/oxo-flow/releases).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-atacseq.git
cd oxo-flow-atacseq
```

### 3. Requirements

**Reference data** (pre-built files, declared under `[config]` in `main.oxoflow`):

| config key | file |
|---|---|
| `reference` | genome FASTA (with `.fai` beside it) |
| `bwa_index` | BWA index prefix (`<prefix>.amb/.ann/.bwt/.pac/.sa`) |
| `chrom_sizes` | chrom sizes file (e.g. `<genome>.sizes`) |
| `gtf` / `gene_bed` / `tss_bed` | annotation for HOMER / deepTools |
| `blacklist` | optional include-regions BED (upstream ENCODE blacklist + chrM complement) — empty disables `-L` filtering |
| `bowtie2_index` / `chromap_index` / `star_index` | index files for the alternative aligners (only when `aligner` is set to them) |

**Input data**: single-end `raw/<sample>.fastq.gz` reads, named in
`[[sample_groups]]` (see `test/fixtures/` for a tiny working set); the
paired-end branch (`config.paired=true`) reads
`raw/<sample>_1.fastq.gz` + `raw/<sample>_2.fastq.gz`.

**Compute**: up to **12 CPUs / 72 GB per rule** (trimming, alignment and
deepTools rules are the heaviest); a few rules need as little as 1 CPU / 6 GB.

**Tools**: a mixed delivery — 41 of 42 rules run in **pinned Docker images**
(`biocontainers/*` tags, e.g. `biocontainers/macs2:2.2.7.1--py38h4a8c8d9_3`),
executed by oxo-flow via Docker or Singularity; the remaining rule
(`picard_markduplicates`) uses a **pinned conda env** at
`envs/picard-samtools.yaml` (picard 3.0.0, samtools 1.17), which requires
conda/mamba at runtime.

## Usage

```bash
# 1. install oxo-flow (see Installation)
# 2. prepare data: <raw_dir>/<sample>.fastq.gz (single-end), see fixtures/
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 8
# 5. run a subset
oxo-flow run main.oxoflow -t multiqc --samples first:1
```

Sample names are declared in `[[sample_groups]]` (`S1`, `S2` in the
fixture set) — add your own names there or point `raw_dir` at your data.
Pipeline behaviour is tuned through `[config]`: `macs_gsize` (default
`2.7e9`), `narrow_peak`, `broad_cutoff`, `fragment_size`,
`min_trimmed_reads`, `out_dir`, and the `skip_*` toggles (`skip_fastqc`,
`skip_qc`, `skip_trimming`, `skip_plot_profile`, `skip_plot_fingerprint`,
`skip_multiqc`, `skip_peak_annotation`). All can be overridden on the CLI.

## Gated branches

Each non-default branch below is `[[include]]`d from `modules/` and gated
by a `when` condition on a `[config]` key (defaults in parentheses).
Nothing from a gated branch runs unless its key is set, and a
`dry-run`/`run` shows exactly which branch a toggle activates.

| config key(s) | branch (rules) | upstream equivalent |
|---|---|---|
| `paired = true` | paired-end: `pe::fastqc_pe`, `pe::trimgalore_pe`, `pe::bwa_mem_pe`, `pe::bamtools_filter_pe`, `pe::pe_name_sort_remove_orphans`, `pe::bedtools_genomecov_pe`, `pe::plotfingerprint_pe`, `pe::multiqc_pe` (8) | PE input path + BAM_SORT_STATS_ORPHANS (name sort + orphan removal) |
| `aligner = "bowtie2"` / `"chromap"` / `"star"` | `alt::bowtie2_align` / `alt::chromap_align` / `alt::star_align` (3, SE only) | BOWTIE2_ALIGN / CHROMAP_CHROMAP / STAR_ALIGN |
| `prepare_reference = true` | `ref::bwa_index`, `ref::custom_getchromsizes` (2) | BWA_INDEX, CUSTOM_GETCHROMSIZES (prepare_genome) |
| `mito_name = "chrM"` or `raw_blacklist = "<bed>"` | `mito::genome_blacklist_regions` (1) | GENOME_BLACKLIST_REGIONS (mitochondrial filtering) |
| `skip_preseq = false` | `qce::preseq_lcextrap` (1) | PRESEQ_LCEXTRAP |
| `skip_picard_metrics = false` | `qce::picard_collectmultiplemetrics` (1) | PICARD_COLLECTMULTIPLEMETRICS |
| `skip_ataqv = false` | `qce::get_autosomes`, `qce::ataqv`, `qce::mkarv` (3) | GET_AUTOSOMES, ATAQV_ATAQV, ATAQV_MKARV |
| `skip_peak_qc = false` | `qce::plot_macs2_qc`, `qce::plot_homer_annotatepeaks` (2) | PLOT_MACS2_QC, PLOT_HOMER_ANNOTATEPEAKS |
| `skip_igv = false` | `qce::igv` (1) | IGV |
| `multiqc_custom_peaks = true` | `qce::multiqc_custom_peaks` (1) | MULTIQC_CUSTOM_PEAKS |
| `skip_consensus_peaks = false` | `cons::macs2_consensus`, `cons::homer_annotatepeaks_consensus`, `cons::subread_featurecounts`, `cons::deseq2_qc` (4, needs ≥ 2 samples) | MACS2_CONSENSUS_PEAKS, HOMER_ANNOTATEPEAKS, SUBREAD_FEATURECOUNTS, DESEQ2_QC |
| `skip_peak_annotation = true` | turns the default HOMER annotation rules off | `params.skip_peak_annotation` |

Where an upstream branch runs by default (`--save_reference`, ataqv,
Picard metrics, IGV, consensus/DESeq2, R QC plots), this port ships it
**off** so the default plan stays the minimal main path; enabling it is
one config flag. See [Known divergences](#known-divergences).

## Source

Upstream: **[nf-core/atacseq](https://github.com/nf-core/atacseq)** @
`2.1.2` (commit `1a1dbe52ffbd82256c941a032b0e22abbd925b8a`), MIT license.
Created 2026-08-15; this workflow may lag behind upstream releases.
Upstream attribution in [NOTICE.md](NOTICE.md).

## Fidelity

Ported with upstream defaults: `aligner=bwa`, single-end, `narrow_peak=false`
(broad peaks), no control. One row per upstream process; steps not ported are
listed with reasons. `when`-gated rules carry the gate in the Notes column.

| Upstream process | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| FASTQC | `fastqc` / `pe::fastqc_pe` | fastqc 0.11.9 | identical command (`--quiet --threads`); PE variant in the paired branch |
| TRIMGALORE | `trimgalore` / `pe::trimgalore_pe` | trim-galore 0.6.7 | identical command (`--fastqc --cores 8 --gzip`); PE variant uses `--paired` + `_val_1/2` outputs |
| FASTQ_FASTQC_UMITOOLS_TRIMGALORE (UMITOOLS_EXTRACT) | — | — | **not ported** — dead code at 2.1.2: the workflow hardcodes `with_umi=false`, the branch can never fire |
| BWA_MEM | `bwa_mem` / `pe::bwa_mem_pe` | bwa 0.7.17, samtools 1.17 | identical (`-M -R '@RG...'`, secondary-alignment filter `-F 0x0100`), mulled container; PE variant merges read groups |
| BOWTIE2_ALIGN | `alt::bowtie2_align` | bowtie2 2.5.1, samtools 1.17 | identical (end-to-end, `--very-sensitive`, `-k 4`, SAM→BAM filter); when `aligner = "bowtie2"` |
| CHROMAP_CHROMAP | `alt::chromap_align` | chromap 0.2.5, samtools 1.17 | identical (`-l 2000 --Tn5-shift --low-mem`); when `aligner = "chromap"` |
| STAR_ALIGN | `alt::star_align` | star 2.6.1d | identical (EndToEnd, `--alignIntronMax 1`, unsorted BAM); when `aligner = "star"` |
| BAM_SORT_STATS_SAMTOOLS (SAMTOOLS_SORT + SAMTOOLS_INDEX + SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `samtools_sort_stats` / `pe::pe_name_sort_remove_orphans` | samtools 1.17 | identical commands, folded into one rule (same env); PE adds name sort → `bampe_rm_orphan.py --only_fr_pairs` → coordinate sort (upstream BAM_SORT_STATS_ORPHANS) |
| PICARD_MERGESAMFILES_LIBRARY | `picard_mergesamfiles` / `pe::bwa_mem_pe` | picard 3.0.0 | upstream symlink branch for single-library samples replicated; multi-library merge (actual MergeSamFiles) not expressible — no library source in oxo-flow |
| BAM_MARKDUPLICATES_PICARD (PICARD_MARKDUPLICATES + SAMTOOLS_INDEX + SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `picard_markduplicates` | picard 3.0.0, samtools 1.17 | identical commands (`--ASSUME_SORTED --REMOVE_DUPLICATES false`, `XMX` heap sizing); combined conda env |
| BAMTOOLS_FILTER + SAMTOOLS_INDEX + BAM_STATS_SAMTOOLS (SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `bamtools_filter` / `pe::bamtools_filter_pe` | bamtools 2.5.2, samtools 1.17 | identical (`-F 0x004 -F 0x0400 -q 1`, optional `-L blacklist`, `assets/bamtools_filter_{se,pe}.json`); PE adds `-f 0x001 -F 0x0008` |
| MACS2_CALLPEAK | `macs2_callpeak` | macs2 2.2.7.1 | identical (`--keep-dup all --nomodel --broad --broad-cutoff 0.1`, `gsize` from config, `--format BAM`) |
| HOMER_ANNOTATEPEAKS | `homer_annotatepeaks` / `cons::homer_annotatepeaks_consensus` | homer 4.11 | identical (`-gid -gtf`); consensus variant annotates the consensus BED, when `skip_consensus_peaks = false` |
| FRIP_SCORE | `frip_score` | bedtools 2.30.0, samtools 1.17 | identical (intersectBed `-f 0.20`, flagstat `mapped` fraction) |
| BEDTOOLS_GENOMECOV | `bedtools_genomecov` / `pe::bedtools_genomecov_pe` | bedtools 2.30.0 | identical (`-bg -scale 1e6/reads -fs fragment_size`, sort); PE uses `-pc` instead of `-fs` (upstream) |
| UCSC_BEDGRAPHTOBIGWIG | `ucsc_bedgraphtobigwig` | ucsc-bedgraphtobigwig 445 | identical |
| BIGWIG_PLOT_DEEPTOOLS (COMPUTEMATRIX scale-regions + reference-point, PLOTPROFILE, PLOTHEATMAP) | `deeptools_plots` | deeptools 3.5.1 | identical args (regionBodyLength 1000, ±3000, `--missingDataAsZero --skipZeros --smartLabels`) |
| MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT | `plotfingerprint` / `pe::plotfingerprint_pe` | deeptools 3.5.1 | identical (`--extendReads fragment_size`); PE omits `--extendReads` (upstream) |
| MULTIQC | `multiqc` / `pe::multiqc_pe` | multiqc 1.13 | upstream mechanism replicated: `multiqc_config.yml` staged in cwd, `multiqc -f .`; config `path_filters` adapted to this port's `results/` layout; PE adds the PE fastqc/trimgalore patterns |
| BWA_INDEX | `ref::bwa_index` | bwa 0.7.17 | identical (`bwa index -p`); when `prepare_reference = true` (default false — port takes pre-built indexes) |
| CUSTOM_GETCHROMSIZES / CUSTOM_GENOME_FASTA_INDEX (prepare_genome) | `ref::custom_getchromsizes` | samtools 1.16.1 | `samtools faidx` + `cut -f 1,2` into `config.chrom_sizes`; when `prepare_reference = true` |
| GENOME_BLACKLIST_REGIONS (mitochondrial filtering) | `mito::genome_blacklist_regions` | bedtools 2.30.0 | sortBed/complementBed over `config.raw_blacklist`, then optional chrM drop (`config.mito_name`, `keep_mito`); when `mito_name`/`raw_blacklist` set; consumed via `-L` by `bamtools_filter` |
| PRESEQ_LCEXTRAP | `qce::preseq_lcextrap` | preseq 3.1.2 | identical (`-verbose -bam -seed 1`; `-pe` when paired); when `skip_preseq = false` |
| PICARD_COLLECTMULTIPLEMETRICS | `qce::picard_collectmultiplemetrics` | picard 3.0.0 | identical (5 metrics + PDFs into `picard_metrics/{,pdf}/`); when `skip_picard_metrics = false` |
| GET_AUTOSOMES + ATAQV_ATAQV + ATAQV_MKARV | `qce::get_autosomes`, `qce::ataqv`, `qce::mkarv` | python 3.8.3, ataqv 1.3.1 | identical commands (`--ignore-read-groups`, mitochondrial-reference-name when `mito_name`; mkarv HTML index); when `skip_ataqv = false` |
| PLOT_MACS2_QC / PLOT_HOMER_ANNOTATEPEAKS | `qce::plot_macs2_qc` / `qce::plot_homer_annotatepeaks` | mulled R image | scripts copied verbatim from upstream `bin/`; when `skip_peak_qc = false` |
| MULTIQC_CUSTOM_PEAKS | `qce::multiqc_custom_peaks` | multiqc headers | identical count/FRiP TSVs (`assets/multiqc/` headers copied verbatim); when `multiqc_custom_peaks = true` (port-only switch; upstream always emits) |
| IGV | `qce::igv` | python 3.8.3 | `igv_files_to_session.py` copied verbatim, `--path_prefix '../../'`; when `skip_igv = false`; per-library tracks only (no merged-replicate sets, see below) |
| MACS2_CONSENSUS_PEAKS | `cons::macs2_consensus` | mulled (macs2 + bedtools + R) | `sort + mergeBed -c 2,3,4,5,6,7,8,9 -o collapse...` → `macs2_merged_expand.py --min_replicates` → BED/SAF/UpSet plot (bin scripts verbatim); when `skip_consensus_peaks = false`, needs ≥ 2 samples |
| SUBREAD_FEATURECOUNTS | `cons::subread_featurecounts` | subread 2.0.1 | identical (`-F SAF -O --fracOverlap 0.2 -s 0`, `-p` when paired); when `skip_consensus_peaks = false` |
| DESEQ2_QC | `cons::deseq2_qc` | mulled (R + DESeq2) | `deseq2_qc.r` verbatim (`--id_col 1 --count_col 7`, `--vst TRUE` when `deseq2_vst`); when `skip_consensus_peaks = false` and `skip_deseq2_qc = false` |
| SAMTOOLS_MERGE / BAM_MERGED_REPLICATE_PICARD / BAM_FILTER_MERGED_REPLICATE / BAM_MERGE_REPLICATES_AND_PEAKS_BEDTOOLS | — | — | **not ported** — merged-replicate analysis is a structural Nextflow pattern: `groupTuple(by: [0])` folds per-replicate BAMs keyed `_REP\d+` into sets, then the merged set drives replicate-level filtering/peaks/QC. oxo-flow has no replicate-grouping primitive over globbed inputs, so the `_REP` merge (and the replicate tracks in IGV/FRiP) cannot be expressed |
| INPUT_CHECK (samplesheet_check) | — | — | **not ported** — pipeline plumbing; oxo-flow provides native `[[sample_groups]]` declaration + `validate` |
| DUMP_SOFTWARE_VERSIONS | — | — | **not ported** — pipeline plumbing; oxo-flow has native version/audit mechanisms |

### Known divergences

- **Branches that upstream runs by default ship OFF here**: ataqv,
  Picard metrics, IGV, consensus peaks/DESeq2 and the R QC plots all run on
  the upstream default path; this port gates them behind
  `skip_ataqv`/`skip_picard_metrics`/`skip_igv`/`skip_consensus_peaks`/
  `skip_peak_qc` (default `true` = off) so the default plan stays the
  minimal main path. `skip_preseq = true` matches the upstream default
  (preseq is off there too). Enabling a branch is a single config flag.
- **Reference inputs**: upstream builds/derives the reference from
  `--genome` (iGenomes) at runtime; this port consumes pre-built files
  (`reference`, `bwa_index` prefix, `gtf`, `gene_bed`, `tss_bed`,
  `chrom_sizes`, optional `blacklist`). `prepare_reference = true` instead
  generates the BWA index and chrom sizes from the FASTA (the fixture
  already ships them, so the rules report up-to-date there).
- **ataqv and Picard metrics are SE-only**: both rules consume the
  single-end filtered `mLb.clN` BAM; the paired-end variants (over the
  `mLb.flT` orphan-removed BAM) are not ported — with `paired=true` the
  rules are skipped even when their `skip_*` gates are off.
- **Alternative aligners are single-end only** (`when` adds
  `!config.paired`): the paired branch is bwa-only, matching upstream's
  paired alignment options. Misconfigurations (e.g. `paired=true` with
  `aligner="star"`) surface as validation warnings via missing inputs.
- **PE requires `bwa_index` and produces `bwa/library/` outputs**: the
  paired branch's `bwa_mem_pe` merges read groups (`-R '@RG\tID:{sample}\tSM:{sample}'`)
  as upstream does; the default SE path keeps the original `@RG` handling.
- **Broad peaks are hardcoded**: the port's `macs2_callpeak` and all
  consumers use `--broad`; `narrow_peak = true` is honoured by
  `macs2_callpeak` but the downstream rules in this port read
  `*_peaks.broadPeak` only.
- **MultiQC config**: `path_filters`/`module_order` were trimmed to the
  ported steps; preseq/featureCounts/ataqv/DESeq2 sections appear when
  their branches are enabled (the PE MultiQC adds the PE fastqc/trimgalore
  logs). Report comment points at the upstream pipeline.
- **`macs_gsize`**: upstream derives it from the read length keyed genome
  block; the port exposes it as `config.macs_gsize` (default `2.7e9`, the
  upstream GRCh37/38 @ 50 bp value).
- **IGV tracks are per-library (`mLb_*`) only**: upstream additionally
  emits merged-replicate tracks; replicate merging is not ported (above).
- **featureCounts is unstranded**: `-s 0` is passed explicitly (the
  upstream default); a `[SCI-FEATURECOUNTS-STRAND]` preflight hint may
  appear for the gated rule — informational.
- **`size_factors/` stays in the workdir**: upstream DESeq2_QC publishes
  the `size_factors` directory; the port leaves it in the rule workdir
  (not declared as an output) — documented, not lost.
- **Consensus branch needs ≥ 2 samples**: upstream filters the peak
  channel to `size() > 1`; with one sample the consensus rules would
  produce degenerate output.

## Test

```bash
bash test/run.sh
```

Runs `oxo-flow validate` + `lint` + `dry-run` (sample selection comes
from the workflow's `[[sample_groups]]` declaration) against the fixture
data; a debug pass additionally asserts that no literal `{wildcards}` leak
into expanded commands. See `test/run.sh` for details.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream
(nf-core/atacseq) is MIT — see [NOTICE.md](NOTICE.md) and `LICENSE.upstream`.
