# 02 — Taxonomic assignment (MARES + MEGAN6)

## Input
- `ASVtable/ASV-frequency-table.tsv`, `ASVseqs/rep-seq-ASV.fasta`
  (from `01_raw_to_asv/`)
- `MARES_NOBAR_16Jan2023` reference database

## Steps
1. **BLAST** each ASV representative sequence against the MARES database
   (`01_blastn_vs_mares.sh`).
2. **MEGAN6** — manual step, not scriptable end-to-end:
   - Download/open [MEGAN6](https://software-ab.cs.uni-tuebingen.de/download/megan6/welcome.html)
   - Import `rep-seq-ASV.fasta` + the BLAST XML output
   - Run the LCA algorithm with the parameters below
   - Export: level = species, format = Txt (csv),
     `readName_to_taxonPathKPCOFGS`
3. Clean the exported taxonomy path (fill "unknown/unclassified/unassigned"
   gaps with the last known rank) — `02_clean_taxonomy_path.sh`.
4. Build the ASV-level phyloseq object from the ASV table + cleaned
   taxonomy + sample metadata + representative sequences —
   `03_build_phyloseq.R`.

## MEGAN6 LCA parameters
| Parameter | Value |
|---|---|
| BLAST e-value (upstream) | 1e-20 |
| min score | 98 |
| max expected | 1×10⁻⁸ |
| min % identity | 80 |
| top % | 10 |
| min support % | 0 (off) |
| min support | 1 |

## Output
- `phyloseq_MARES_full_e20.rds` — ASV-level phyloseq object with taxonomy
  → input for `03_decontam/`

## Reproducibility note
Step 2 (MEGAN6) is a manual GUI step and is the one part of this pipeline
that cannot be fully scripted.
