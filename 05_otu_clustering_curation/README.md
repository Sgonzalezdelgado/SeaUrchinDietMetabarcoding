# 05 — OTU clustering and post-clustering curation

## 5.1 Clustering (vsearch)
ASVs clustered into OTUs at 97% similarity:
```bash
vsearch --cluster_fast ASV_repseqs_tsc_e20.fasta --id 0.97 \
  --centroids OTU97_centroids_e20.fasta --uc OTU97_mapping_e20.uc --threads 35
```
ASV→OTU97 table collapsed in R (see `build_otu97_table.R`) to produce
`phyloseq_OTU97_decontam_tsc_e20.rds` (32,993 OTUs × 170 samples).

## 5.2 Post-clustering curation (LULU)

Self-match list built with vsearch at 84% identity:
```bash
vsearch --usearch_global OTU97_centroids_e20.fasta --db OTU97_centroids_e20.fasta \
  --id 0.84 --self --iddef 1 --maxaccepts 0 --maxrejects 0 \
  --uc match_list_OTU97_e20.uc --threads 35
```

### ⚠️ Known bug in `tobiasgf/lulu`
The official package (commit `50187ba8bb5eb3a36f4509e8dc0a7aeaa63a1e75`,
2023-12-18, `Version 0.1.0`) has an unfixed bug
([issue #8](https://github.com/tobiasgf/lulu/issues/8)): with
`minimum_ratio_type = "min"`, OTUs are only merged when co-occurrence with
their parent equals 1, regardless of `minimum_relative_cooccurence`.

**Fix used:** the community-patched function from
[`jesuszarcero/LULU_corrected`](https://github.com/jesuszarcero/LULU_corrected)
(`lulu2()`, same argument signature as `lulu()`). Only the function
definition is sourced — not the author's example script at the bottom of
the file.

### Parameters
```r
cur <- lulu2(otu_df, match_df,
             minimum_ratio_type = "min",
             minimum_relative_cooccurence = 0.5,
             minimum_match = 0.84)
```
`minimum_match` is set to **0.84** (not the package default of `84`) because
the similarity values in `match_df` are scaled 0–1, not 0–100. Using the
default here would silently reject every match.

### Result
- Corrected (`lulu2()`): 32,993 → **30,386 OTUs** (2,607 merged as artefacts)

## Output
- `phyloseq_OTU97_LULU_e20.rds` — curated OTU table (30,386 OTUs × 170 samples),
  input for step 06


