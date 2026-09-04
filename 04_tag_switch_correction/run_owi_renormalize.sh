#!/bin/bash
# ══════════════════════════════════════════════════════════════════════════
# 02 — Tag-switching correction
# Corrects tag-jump/tag-switch artefacts using owi_renormalize
# (Wangensteen, Project Metabarpark)
# ══════════════════════════════════════════════════════════════════════════
#
# INPUT:  exported OTU/ASV table with sample columns in a known range
# OUTPUT: renormalized table, feeding into 03_otu_clustering_curation
#
# -c : cumulative frequency cutoff (the only scientifically meaningful
#      parameter here — report this one in Methods)
# -s / -e : start/end column of sample data in the input table (indexing
#      detail, not a method parameter). Confirmed here as 2-171, matching
#      the 170 samples in ASVtable_decontam_true_e20.tsv. If you re-export
#      the table with a different column layout, update these.

Rscript owi_renormalize.R \
  -i ASVtable_decontam_true_e20.tsv \
  -o ASVtable_tsc_e20.tsv \
  -c 0.97 -s 2 -e 171
