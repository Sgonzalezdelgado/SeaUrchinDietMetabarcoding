#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 03.3 — Self-match at 84% identity, for LULU's match list
# ══════════════════════════════════════════════════════════════════════════
# INPUT:  OTU97_centroids_e20.fasta
# OUTPUT: match_list_OTU97_3col_e20.txt  (query, target, identity 0-1)

vsearch \
  --usearch_global OTU97_centroids_e20.fasta \
  --db OTU97_centroids_e20.fasta \
  --id 0.84 \
  --self \
  --iddef 1 \
  --maxaccepts 0 \
  --maxrejects 0 \
  --uc match_list_OTU97_e20.uc \
  --threads 35

awk 'BEGIN{FS="\t"; OFS="\t"}
$1=="H" {
  pid=$4
  if (pid > 1) pid = pid/100
  print $9,$10,pid
}' match_list_OTU97_e20.uc > match_list_OTU97_3col_e20.txt
