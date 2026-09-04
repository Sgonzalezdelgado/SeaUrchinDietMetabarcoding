#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 01.3 — Denoise, remove chimeras, infer ASVs (DADA2), and export
# ══════════════════════════════════════════════════════════════════════════
conda activate qiime2-2022.11
mkdir -p 03.dada2
cd 02.qiime_ready

# Note: --p-trim-left set to 0/0 because primers were already removed by
# cutadapt in the previous step. Truncation lengths chosen from the R1/R2
# quality plots (R2 quality drops after ~180 bp).
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed-seqs.qza \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --p-trunc-len-f 210 \
  --p-trunc-len-r 180 \
  --p-n-threads 0 \
  --p-chimera-method consensus \
  --o-table ../03.dada2/table-dada2.qza \
  --o-representative-sequences ../03.dada2/rep-seqs-dada2.qza \
  --o-denoising-stats ../03.dada2/denoising-stats-dada2.qza

cd ../03.dada2

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table-dada2.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs-dada2.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats-dada2.qza \
  --o-visualization denoising-stats-dada2.qzv

# ---- Export ASV table and representative sequences ----
mkdir -p ../04.MARES/ASVtable
mkdir -p ../04.MARES/ASVseqs

qiime tools export \
  --input-path table-dada2.qza \
  --output-path ../04.MARES/ASVtable

qiime tools export \
  --input-path rep-seqs-dada2.qza \
  --output-path ../04.MARES/ASVseqs

mv ../04.MARES/ASVseqs/dna-sequences.fasta \
   ../04.MARES/ASVseqs/rep-seq-ASV.fasta

cd ..
biom convert \
  -i 04.MARES/ASVtable/feature-table.biom \
  -o 04.MARES/ASVtable/ASV-frequency-table.tsv \
  --to-tsv
