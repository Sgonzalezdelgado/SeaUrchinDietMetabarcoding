#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 01.2 — Primer removal (cutadapt, via QIIME2)
# ══════════════════════════════════════════════════════════════════════════
conda activate qiime2-2022.11
cd 02.qiime_ready

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences demux-paired-end.qza \
  --p-front-f GGWACWRGWTGRACWNTNTAYCCYCC \
  --p-front-r TANACYTCNGGRTGNCCRAARAAYCA \
  --p-error-rate 0.1 \
  --o-trimmed-sequences trimmed-seqs.qza \
  --verbose

qiime demux summarize \
  --i-data trimmed-seqs.qza \
  --o-visualization trimmed-seqs.qzv
