#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 03.1 — OTU clustering at 97% similarity
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  ASV_repseqs_tsc_e20.fasta  (from step 02)
# OUTPUT: OTU97_centroids_e20.fasta, OTU97_mapping_e20.uc

vsearch \
  --cluster_fast ASV_repseqs_tsc_e20.fasta \
  --id 0.97 \
  --centroids OTU97_centroids_e20.fasta \
  --uc OTU97_mapping_e20.uc \
  --threads 35

# Build a simple ASV → OTU97 map from the .uc file
awk 'BEGIN{FS="\t"; OFS="\t"}
$1=="H" {print $9,$10}
$1=="S" {print $9,$9}
' OTU97_mapping_e20.uc > ASV_to_OTU97_e20.tsv
