# Environment

## R statistical environment (recommended)
`environment.yml` in this folder builds a clean conda environment for all
R steps in this repo (02_taxonomic_assignment/03_build_phyloseq.R,
03_decontam, 04_tag_switch_correction, 05/06 curation, 07 diversity
analyses), **separate from the qiime2 environment used in step 01**.

```bash
conda env create -f environment.yml
conda activate r-stats
R.version.string
```

Verified package list is derived directly from the `library()`/setup calls
actually present in the final scripts in this repo — not a generic list.

For the packages installed from GitHub:
```r
packageDescription("lulu")   # RemoteSha field
```
and note the source of `lulu2()` (community patch, not an installable
package): https://github.com/jesuszarcero/LULU_corrected

## Command-line tools
| Tool | Version |
|---|---|
| vsearch | `[run: vsearch --version]` |
| DADA2 | `qiime2-2022.11` |
| MEGAN6 | 6.25.9 |
| QIIME2 | 2022.11 |

## Files to add here
- `environment.yml` — done, see above
