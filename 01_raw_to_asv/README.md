# 01 — Raw reads to ASVs (QIIME2)

Followed the pipeline in
[vanearranz/Metabarcoding_CO1_kelpholdfast](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast),
using **QIIME2 2022.11** and the latest release of the
[MARES](https://github.com/wpearman1996/MARES_database_pipeline) reference
database at the time (`MARES_NOBAR_16Jan2023`).

## Input
- Raw paired-end FASTQ files, two lanes per
  sample, concatenated before import (see `01_import_and_concatenate.sh`).

## Steps
1. Concatenate lanes per sample, import into QIIME2.
2. Primer removal with `cutadapt` (COI primers).
3. Denoising, chimera removal, and ASV inference with `DADA2`.
4. Export the ASV table and representative sequences.

All commands: `01_import_and_concatenate.sh`, `02_cutadapt_trim.sh`,
`03_dada2_denoise.sh`.

## Key parameters
| Step | Parameter | Value |
|---|---|---|
| cutadapt | forward primer | `GGWACWRGWTGRACWNTNTAYCCYCC` |
| cutadapt | reverse primer | `TANACYTCNGGRTGNCCRAARAAYCA` |
| cutadapt | error rate | 0.1 |
| DADA2 | trim-left (F/R) | 0 / 0 |
| DADA2 | trunc-len (F/R) | 210 / 180 |
| DADA2 | chimera method | consensus |

Truncation lengths were chosen from the per-base quality plots (R1 stayed
high quality throughout; R2 quality dropped after ~180 bp).

## Output
- `ASVtable/ASV-frequency-table.tsv`, `ASVseqs/rep-seq-ASV.fasta`
  → inputs for `02_taxonomic_assignment/`

