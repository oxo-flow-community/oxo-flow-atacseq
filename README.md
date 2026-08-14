# oxo-flow-atacseq — ATAC-seq peak calling and QC

[![CI](https://github.com/oxo-flow-community/oxo-flow-atacseq/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-atacseq/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

ATAC-seq (Assay for Transposase-Accessible Chromatin) analysis pipeline:
raw-read QC (FastQC), adapter trimming (Trim Galore), BWA-MEM alignment,
Picard mark-duplicates, BAMTools filtering, MACS2 broad-peak calling, HOMER
peak annotation, FRiP scoring, normalised bigWig tracks, deepTools QC plots
and a MultiQC report. Ported to oxo-flow from nf-core/atacseq, single-end,
BWA aligner, broad-peak default configuration.

## Source

Ported from **[nf-core/atacseq](https://github.com/nf-core/atacseq)**, version
`2.1.2` (MIT). This port is maintained independently and **may lag the
upstream** — check the version above and the fidelity table below for the
exact ported state.

## Fidelity

Ported with upstream defaults: `aligner=bwa`, single-end, `narrow_peak=false`
(broad peaks), no control. One row per upstream process; steps not ported are
listed with reasons.

| Upstream process | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| FASTQC | `fastqc` | fastqc 0.11.9 | identical command (`--quiet --threads`), single-end |
| TRIMGALORE | `trimgalore` | trim-galore 0.6.7 | identical command (`--fastqc --cores 8 --gzip`), single-end branch |
| BWA_MEM | `bwa_mem` | bwa 0.7.17, samtools 1.17 | identical (`-M -R '@RG...'`, secondary-alignment filter `-F 0x0100`), mulled container |
| BAM_SORT_STATS_SAMTOOLS (SAMTOOLS_SORT + SAMTOOLS_INDEX + SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `samtools_sort_stats` | samtools 1.17 | identical commands, folded into one rule (same env) |
| PICARD_MERGESAMFILES_LIBRARY | `picard_mergesamfiles` | picard 3.0.0 | upstream symlink branch for single-library samples replicated; multi-library merge (actual MergeSamFiles) not expressible — no library source in oxo-flow |
| BAM_MARKDUPLICATES_PICARD (PICARD_MARKDUPLICATES + SAMTOOLS_INDEX + SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `picard_markduplicates` | picard 3.0.0, samtools 1.17 | identical commands (`--ASSUME_SORTED --REMOVE_DUPLICATES false`, `XMX` heap sizing); combined conda env |
| BAMTOOLS_FILTER + SAMTOOLS_INDEX + BAM_STATS_SAMTOOLS (SAMTOOLS_STATS + SAMTOOLS_FLAGSTAT + SAMTOOLS_IDXSTATS) | `bamtools_filter` | bamtools 2.5.2, samtools 1.17 | identical (`-F 0x004 -F 0x0400 -q 1`, optional `-L blacklist`, `assets/bamtools_filter_se.json`); single-end branch — no name-sort/orphan removal |
| MACS2_CALLPEAK | `macs2_callpeak` | macs2 2.2.7.1 | identical (`--keep-dup all --nomodel --broad --broad-cutoff 0.1`, `gsize` from config, `--format BAM`) |
| HOMER_ANNOTATEPEAKS | `homer_annotatepeaks` | homer 4.11 | identical (`-gid -gtf`) |
| FRIP_SCORE | `frip_score` | bedtools 2.30.0, samtools 1.17 | identical (intersectBed `-f 0.20`, flagstat `mapped` fraction) |
| BEDTOOLS_GENOMECOV | `bedtools_genomecov` | bedtools 2.30.0 | identical (`-bg -scale 1e6/reads -fs fragment_size`, sort) |
| UCSC_BEDGRAPHTOBIGWIG | `ucsc_bedgraphtobigwig` | ucsc-bedgraphtobigwig 445 | identical |
| BIGWIG_PLOT_DEEPTOOLS (COMPUTEMATRIX scale-regions + reference-point, PLOTPROFILE, PLOTHEATMAP) | `deeptools_plots` | deeptools 3.5.1 | identical args (regionBodyLength 1000, ±3000, `--missingDataAsZero --skipZeros --smartLabels`) |
| MERGED_LIBRARY_DEEPTOOLS_PLOTFINGERPRINT | `plotfingerprint` | deeptools 3.5.1 | identical (`--extendReads fragment_size`, single-end) |
| MULTIQC | `multiqc` | multiqc 1.13 | upstream mechanism replicated: `multiqc_config.yml` staged in cwd, `multiqc -f .`; config `path_filters` adapted to this port's `results/` layout |
| BWA_INDEX | — | bwa 0.7.17 | **not ported** — reference preparation delegated to the user (pre-built index files; see References) |
| CUSTOM_GETCHROMSIZES / CUSTOM_GENOME_FASTA_INDEX (prepare_genome) | — | — | **not ported** — pre-built `.sizes`/`.fai` files are inputs |
| GENOME_BLACKLIST_REGIONS (mitochondrial filtering) | — | — | **not ported** — blacklist complement step needs prepare_genome; an optional `-L` blacklist file is honored by `bamtools_filter` when `config.blacklist` is set |
| SAMTOOLS_MERGE / BAM_MERGED_REPLICATE_PICARD / BAM_FILTER_MERGED_REPLICATE | — | — | **not ported** — merged-replicate analysis is off the default main path |
| BAM_MERGE_REPLICATES_AND_PEAKS_BEDTOOLS / MACS2_CONSENSUS_PEAKS / DESEQ2_* | — | — | **not ported** — consensus peaks / differential expression (needs replicate/group metadata) |
| ATAQV / PRESEQ_LCEXTRAP / PICARD_METRICS / IGV | — | — | **not ported** — optional QC branches |
| PLOT_MACS2_QC / other R QC plots | — | — | **not ported** — R plotting utilities, off the default main path |
| MULTIQC_CUSTOM_PEAKS | — | — | **not ported** — custom peak-count/FRiP report sections; FRiP scores still produced by `frip_score` |
| INPUT_CHECK (samplesheet_check) / DUMP_SOFTWARE_VERSIONS | — | — | **not ported** — pipeline plumbing (oxo-flow provides native sample/validate handling) |
| UMI handling (umitools branch), STAR/Bowtie2 aligners | — | — | **not ported** — `--umitools`/`aligner` alternatives to the default path |

### Known divergences

- **Reference inputs**: upstream builds/derives the reference from
  `--genome` (iGenomes) at runtime; this port consumes pre-built files
  (`reference`, `bwa_index` prefix, `gtf`, `gene_bed`, `tss_bed`,
  `chrom_sizes`, optional `blacklist`).
- **MultiQC config**: `path_filters`/`module_order` were trimmed to the
  ported steps (preseq, featureCounts, merged-replicate sections removed);
  report comment points at the upstream pipeline.
- **`macs_gsize`**: upstream derives it from the read length keyed genome
  block; the port exposes it as `config.macs_gsize` (default `2.7e9`, the
  upstream GRCh37/38 @ 50 bp value).

## Quickstart

```bash
# 1. install oxo-flow (see Requirements)
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

## Requirements

- **oxo-flow ≥ 0.11.0** — install the prebuilt binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

- Conda users may alternatively `conda install -c bioconda oxo-flow-cli`
  (note: the bioconda package currently lags the release binary at 0.10.2 —
  some 0.11.0 format features may not validate).
- Docker/Singularity/conda at runtime, per the environments declared in
  `main.oxoflow` (`picard_markduplicates` uses `envs/picard-samtools.yaml`;
  all other rules use the upstream container tags).

## References

The port expects pre-built reference files (no `prepare_genome` step):

| config key | file |
|---|---|
| `reference` | genome FASTA (with `.fai` beside it) |
| `bwa_index` | BWA index prefix (`<prefix>.amb/.ann/.bwt/.pac/.sa`) |
| `chrom_sizes` | chrom sizes file (e.g. `<genome>.sizes`) |
| `gtf` / `gene_bed` / `tss_bed` | annotation for HOMER / deepTools |
| `blacklist` | optional include-regions BED (upstream ENCODE blacklist + chrM complement) — empty disables `-L` filtering |

See `test/fixtures/` for a tiny working set.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md).

## Community

https://oxo-flow-community.github.io/
