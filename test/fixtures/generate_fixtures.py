#!/usr/bin/env python3
"""Deterministic fixture generator for oxo-flow-atacseq live tests.

One seed drives every file (failure catalog: "Fixture incoherence" — reads
and annotations must come from one source, or alignment/annotation yields
nothing):

  genome.fa      chr1/chr2 200 kb each + chrM 20 kb, fixed-width lines
  genome.fa.fai  computed from the byte offsets used while writing
  genome.fa.sizes
  genes.gtf      one gene per peak zone (exon = zone)
  gene.bed       peak zones (deepTools computeMatrix regions)
  tss.bed        1 bp TSS entries at zone starts (deepTools + ataqv)
  raw/S*.fastq.gz  single-end reads
  raw/S*_1/_2.fastq.gz  paired-end reads (R2 = reverse complement of R1)

Fixture design (each choice is load-bearing, see the failure catalog):

- ~150 peak zones per autosome, reads concentrated on zones: MACS2 needs
  real peak structure, and the DESeq2 dispersion fit needs HUNDREDS of
  differential features (a handful of peaks makes locfit die with
  "newsplit: out of vertex space").
- Sample-differential zone weights: S1 reads concentrate on even zones,
  S2 on odd zones, so the consensus featureCounts matrix carries real
  per-feature variance for DESeq2.
- ~30 % of reads duplicated at 2x/4x/8x/16x: preseq lc_extrap needs a
  proper duplicate-count curve ("too many defects in the approximation"
  otherwise).
- R2 is the reverse COMPLEMENT of R1: a plain reversal maps to the same
  strand, every pair fails the -f 0x001 proper-pair filter and the PE
  path ends with zero mapped reads.
- chrM contig + chrM-mapped reads: exercises the mito-filter gate
  (mito_name=chrM).

Run:  python3 test/fixtures/generate_fixtures.py [outdir]
(outdir defaults to the fixtures directory; the bwa index beside
genome.fa is built with the SAME mulled bwa+samtools image as bwa_mem —
see the live-test notes in README — not by this script.)
"""

from __future__ import annotations

import gzip
import pathlib
import random
import sys

SEED = 42
READ_LEN = 50
LINE_WIDTH = 60

CONTIGS = {"chr1": 200_000, "chr2": 200_000, "chrM": 20_000}
ZONES_PER_AUTOSOME = 150
ZONE_LEN = 400
ZONE_STRIDE = 1300  # zone start distance
ZONE_OFFSET = 500

# sample -> per-zone read yield: high/low/extra-noise (odd/even zones)
WEIGHTS = {
    "S1": {"high": 100, "low": 30, "noise": 15},
    "S2": {"high": 100, "low": 30, "noise": 15},
}
CHRM_READS = 800
# duplicates: pool 150 unique reads, add copies at these multiplicities
DUP_MULTS = [2] * 60 + [4] * 40 + [8] * 30 + [16] * 20


def zones(contig: str) -> list[tuple[int, int]]:
    return [
        (ZONE_OFFSET + i * ZONE_STRIDE, ZONE_OFFSET + i * ZONE_STRIDE + ZONE_LEN)
        for i in range(ZONES_PER_AUTOSOME)
    ]


def main() -> int:
    outdir = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).parent
    genomedir = outdir / "genome"
    rawdir = outdir / "raw"
    genomedir.mkdir(parents=True, exist_ok=True)
    rawdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    # --- genome + fai (offsets tracked while writing) + sizes ---
    fasta = genomedir / "genome.fa"
    sizes = []
    fai = []
    with open(fasta, "w") as fh:
        for contig, length in CONTIGS.items():
            seq = "".join(rng.choice("ACGT") for _ in range(length))
            offset = fh.tell()
            fh.write(f">{contig}\n")
            for i in range(0, length, LINE_WIDTH):
                fh.write(seq[i : i + LINE_WIDTH] + "\n")
            line_bases = LINE_WIDTH
            n_lines = (length + line_bases - 1) // line_bases
            sizes.append((contig, length))
            fai.append((contig, length, offset, line_bases, line_bases * n_lines + len(f">{contig}\n")))
    (genomedir / "genome.fa.sizes").write_text(
        "".join(f"{c}\t{l}\n" for c, l in sizes)
    )
    (genomedir / "genome.fa.fai").write_text(
        "".join(f"{c}\t{l}\t{o}\t{b}\t{lb}\n" for c, l, o, b, lb in fai)
    )
    print(f"genome: {fasta} ({sum(l for _, l in sizes):,} bp)")

    # --- annotations: one gene per peak zone ---
    gtf_lines = []
    bed_lines = []
    tss_lines = []
    for contig in ("chr1", "chr2"):
        for idx, (start, end) in enumerate(zones(contig)):
            gid = f"g{contig[3:]}_{idx}"
            gtf_lines.append(
                f'{contig}\tfixture\texon\t{start + 1}\t{end}\t.\t+\t.\t'
                f'gene_id "{gid}"; transcript_id "t{gid}";'
            )
            bed_lines.append(f"{contig}\t{start}\t{end}\t{gid}")
            tss_lines.append(f"{contig}\t{start}\t{start + 1}\t{gid}")
    (genomedir / "genes.gtf").write_text("\n".join(gtf_lines) + "\n")
    (genomedir / "gene.bed").write_text("\n".join(bed_lines) + "\n")
    (genomedir / "tss.bed").write_text("\n".join(tss_lines) + "\n")
    print(f"annotations: {len(gtf_lines)} genes (genes.gtf / gene.bed / tss.bed)")

    # --- sequences for read synthesis ---
    seqs = {}
    name = None
    buf = []
    for line in open(fasta):
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name = line[1:].strip()
            buf = []
        else:
            buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")

    def read_from(contig: str, start: int) -> str:
        return seqs[contig][start : start + READ_LEN]

    for sample, weights in WEIGHTS.items():
        rng2 = random.Random(SEED ^ hash(sample))
        reads = []
        for contig in ("chr1", "chr2"):
            zs = zones(contig)
            for idx, (start, end) in enumerate(zs):
                is_high = (idx % 2 == 0) if sample == "S1" else (idx % 2 == 1)
                n = (weights["high"] if is_high else weights["low"]) + rng2.randrange(
                    weights["noise"] + 1
                )
                for k in range(n):
                    pos = rng2.randrange(start, end - READ_LEN)
                    reads.append((f"{sample}_{contig}_{idx}_{k}", read_from(contig, pos)))
        # chrM-mapped reads (mito gate)
        for k in range(CHRM_READS):
            pos = rng2.randrange(0, CONTIGS["chrM"] - READ_LEN)
            reads.append((f"{sample}_chrM_{k}", read_from("chrM", pos)))
        rng2.shuffle(reads)

        # ~30 % duplicates at 2x/4x/8x/16x (preseq duplicate curve)
        pool = random.Random(SEED).sample(reads, min(150, len(reads)))
        dups = []
        for j, mult in enumerate(DUP_MULTS):
            dups.extend((f"dup{j}_{k}", pool[j % len(pool)][1]) for k in range(mult))
        all_reads = reads + dups
        rng2.shuffle(all_reads)

        with gzip.open(rawdir / f"{sample}.fastq.gz", "wt") as fh:
            for h, r in all_reads:
                fh.write(f"@{h}\n{r}\n+\n{'I' * READ_LEN}\n")
        with gzip.open(rawdir / f"{sample}_1.fastq.gz", "wt") as fh1, gzip.open(
            rawdir / f"{sample}_2.fastq.gz", "wt"
        ) as fh2:
            for h, r in all_reads:
                rc = r.translate(comp)[::-1]
                fh1.write(f"@{h}/1\n{r}\n+\n{'I' * READ_LEN}\n")
                fh2.write(f"@{h}/2\n{rc}\n+\n{'I' * READ_LEN}\n")
        print(f"reads: {sample} — {len(all_reads):,} total (~30 % duplicated)")

    print("done. Next: build the bwa index with the bwa_mem mulled image:")
    print("  docker run --rm -v $PWD/test/fixtures:/f -w /f \\")
    print("    quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:219b6c272b25e7e642ae3ff0bf0c5c81a5135ab4-0 \\")
    print("    bwa index -p /f/genome/genome.fa /f/genome/genome.fa")
    return 0


if __name__ == "__main__":
    sys.exit(main())
