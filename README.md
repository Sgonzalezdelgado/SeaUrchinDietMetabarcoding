# Sea urchin diet metabarcoding pipeline

Code and analysis pipeline for the diet metabarcoding study of two sea urchin
species (*Arbacia lixula* and *Paracentrotus lividus*) across ambient and
naturally CO2-acidified (vent) localities, in La Palma (Canary Islands) and
Blanes (NW Mediterranean, methodological validation site).

This repository accompanies the manuscript:

> González-Delgado, S., Hernández, J. C., Pérez-Portela, R., Arranz, V. [Year]. Non-invasive DNA metabarcoding reveals distinct diets in two sympatric sea urchins across natural acidified systems. *Journal* (in prep). DOI: TBD

Real locality names and exact coordinates are
reported only in the manuscript's.


## Pipeline overview

| Step | Folder | What it does |
|---|---|---|
| 1 | [`01_raw_to_asv/`](01_raw_to_asv/) | Raw read import, primer trimming (cutadapt), ASV inference (DADA2, QIIME2 2022.11) |
| 2 | [`02_taxonomic_assignment/`](02_taxonomic_assignment/) | BLAST vs. MARES reference database + LCA taxonomic assignment (MEGAN6) |
| 3 | [`03_decontam/`](03_decontam/) | Contaminant removal (R `decontam`, prevalence method) |
| 4 | [`04_tag_switch_correction/`](04_tag_switch_correction/) | Tag-jump/tag-switch correction (`owi_renormalize`) |
| 5 | [`05_otu_clustering_curation/`](05_otu_clustering_curation/) | OTU clustering (vsearch, 97%) and post-clustering curation (LULU) |
| 6 | [`06_taxonomic_curation/`](06_taxonomic_curation/) | Eukaryote restriction, rule-based dietary curation (29 documented steps), species-level collapse |
| 7 | [`07_diversity_analyses/`](07_diversity_analyses/) | Alpha/beta diversity, PERMANOVA, IndVal, interspecific overlap (Objectives 1–3) |


Each folder has its own short `README.md` describing the scripts inside and
the exact parameters used. Pipeline followed for steps 1–2:
[vanearranz/Metabarcoding_CO1_kelpholdfast](https://github.com/vanearranz/Metabarcoding_CO1_kelpholdfast).

## Software and key parameters

| Tool | Use | Notes |
|---|---|---|
| QIIME2 | Import, trimming, ASV pipeline wrapper | version 2022.11 |
| cutadapt (via QIIME2) | Primer trimming | error rate 0.1 |
| DADA2 (via QIIME2) | ASV inference | trunc-len 210/180, chimera method consensus |
| MARES | Reference database for taxonomic assignment | `MARES_NOBAR_16Jan2023` |
| BLAST | Sequence search vs. MARES | `-evalue 1e-20` |
| MEGAN6 | LCA taxonomic assignment | min score 98, max expected 1×10⁻⁸, min identity 80%, top 10% |
| decontam | Contaminant removal | prevalence method, threshold 0.5; version TBD — see `environment/` |
| `owi_renormalize` (Wangensteen, Project Metabarpark) | Tag-switch correction | cumulative frequency cutoff = 0.97 |
| vsearch | OTU clustering | `--id 0.97` |
| `lulu` (tobiasgf/lulu) | Post-clustering curation | **known bug**, see note below |
| R / phyloseq / vegan / lme4 / emmeans / indicspecies / betadisper | Diversity & statistical analyses | see `environment/` |

**Note on `lulu`:** the official `tobiasgf/lulu` package (GitHub, no CRAN
release, stuck at `Version: 0.1.0`) has an unresolved bug
([tobiasgf/lulu#8](https://github.com/tobiasgf/lulu/issues/8)) that causes
OTUs to merge only when co-occurrence equals 1, regardless of the
`minimum_relative_cooccurence` value supplied. This pipeline uses the
community-patched function (`jesuszarcero/LULU_corrected`) instead. Details
and the exact commit used are in `05_otu_clustering_curation/README.md`.

## Data availability

Raw sequencing reads: deposited in SRA under accession `PRJNA1552446`.
This repository contains code only — no raw sequence data.

## Reproducing this pipeline

Run the folders in numeric order (`01` → `07`). Each step's `README.md`
lists its exact inputs and outputs, so you can start from any intermediate
file without re-running everything from raw reads.

## Citation

If you use this code, please cite the manuscript above.
