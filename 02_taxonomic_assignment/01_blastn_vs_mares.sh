#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 02.1 — BLAST ASVs against the MARES reference database
# ══════════════════════════════════════════════════════════════════════════
makeblastdb \
  -in MARES_NOBAR_16012023.fasta \
  -dbtype nucl \
  -parse_seqids

blastn \
  -db MARES_NOBAR_16Jan2023/MARES_NOBAR_16012023.fasta \
  -query 04.MARES/ASVseqs/rep-seq-ASV.fasta \
  -evalue 1e-20 \
  -outfmt 5 \
  -max_target_seqs 500 \
  -num_threads 35 \
  -out 04.MARES/MARES_MEGAN_e20.xml

# Output feeds into MEGAN6 (manual step — see README.md)
