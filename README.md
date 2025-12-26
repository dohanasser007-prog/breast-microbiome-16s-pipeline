# breast-microbiome-16s-pipeline

A reproducible QIIME 2 + R workflow for 16S rRNA amplicon microbiome analysis comparing **Breast_cancer** vs **Healthy_control** samples, created for the Bioinformatics Volunteer Challenge (Task 1–3). 
The repository contains QIIME 2 artifacts/visualizations plus R scripts to generate genus-level figures, diversity outputs, and differential abundance results for reporting. 

---

## Challenge tasks (coverage)

### Task 1 — Data acquisition & scoping
**Goal:** obtain raw FASTQ files and prepare metadata with the key variable `Status` (Breast Cancer vs Healthy Control).   
**In this repo:** input reads are stored in `data/`, and import is driven by manifest files (`manifest*.tsv`) plus `metadata.txt`. 

> Optional (if starting from scratch): FASTQ files can be downloaded from SRA using SRA Toolkit (`prefetch` / `fasterq-dump`) and then organized to match the manifest format. 

### Task 2 — Core bioinformatics analysis
**Goal:** run a standard 16S pipeline including QC, ASV inference, taxonomy, essential figures, diversity testing, and differential abundance.
**In this repo:** QIIME 2 `.qza` artifacts and `.qzv` visualizations capture outputs and provenance for reproducibility.

### Task 3 — Automation & reproducibility
**Goal:** provide scripts and a logical folder structure to reproduce the workflow and generate final outputs. 
**In this repo:** plotting is automated with `scripts/plot_top500_genus.R` (genus composition + log2FC effect size) and `scripts/plot_ancombc_genus.R` (ANCOM‑BC volcano/barplot). 

---

## Project layout

breast_microbiome/
data/ # FASTQ inputs
metadata.txt # SampleID + Status

manifest.tsv # QIIME 2 manifest(s)
manifest_batch1.tsv
manifest_batch2.tsv
manifest_batch2a.tsv
manifest_batch2b.tsv
manifest_single.tsv

refs_silva/ # SILVA reference artifacts (qza)
scripts/
plot_top500_genus.R # FINAL: genus plots from Top500 (QIIME2 exported table + taxonomy)
plot_ancombc_genus.R # FINAL: ANCOM-BC volcano + significant barplot
analysis_microbiome.R # LEGACY: DECIPHER IdTaxa attempt (not used for final figs)

r_input/ # exported inputs for R + helper files
feature-table-top500.tsv
exported-taxonomy-top500/taxonomy.tsv
ancombc_genus_top500_export/ # exported ANCOM-BC slices (csv)
top500_asv_ids.txt
top500_asv_ids_metadata.tsv

core-metrics-no-tree/ # non-phylogenetic diversity outputs
results/
figures/
tables/


---

## Requirements

### QIIME 2 (WSL/Ubuntu)
- QIIME 2 amplicon distribution environment (example used here: `qiime2-amplicon-2024.10`).

### R (Windows or Linux)
- R + packages: `ggplot2` (minimum for the plotting scripts).

---

## QIIME 2 outputs (this run)

### Demultiplexing / QC
- `demux-single-end-amplicon.qza`
- `demux-single-end-amplicon.qzv`
- Batch demux artifacts (examples): `demux-batch1.qza`, `demux-batch2a.qza`, `demux-batch2b.qza`

### Denoising / ASVs (DADA2)
- Merged/all artifacts:
  - `table-all.qza` (+ `table-all.qzv`)
  - `rep-seqs-all.qza` (+ `rep-seqs-all.qzv`)
- Batch-level artifacts:
  - `table-batch1.qza`, `table-batch2a.qza`, `table-batch2b.qza`
  - `rep-seqs-batch1.qza`, `rep-seqs-batch2a.qza`, `rep-seqs-batch2b.qza`
  - `denoising-stats-batch1.qza`, `denoising-stats-batch2a.qza`, `denoising-stats-batch2b.qza`

---

## Diversity analysis (Task 2C)

### Alpha diversity
- Shannon group comparison:
  - `core-metrics-no-tree/shannon-Status-significance.qzv`

### Beta diversity (no-tree)
- Bray–Curtis PCoA:
  - `core-metrics-no-tree/bray_curtis_emperor.qzv`
- Jaccard PCoA:
  - `core-metrics-no-tree/jaccard_emperor.qzv`
- Bray–Curtis group significance:
  - `core-metrics-no-tree/bray-curtis-Status-significance.qzv`

### Weighted UniFrac workaround (filtered table + tree)
Tree building can be memory intensive on large ASV sets, so this run includes a filtered-table workaround (`filt50`) to enable phylogenetic distances under limited RAM. 
Key outputs:
- Filtered artifacts:
  - `table-filt50.qza` (+ `table-filt50.qzv`)
  - `rep-seqs-filt50.qza`
  - `aligned-rep-seqs-filt50.qza`, `masked-aligned-rep-seqs-filt50.qza`
  - `unrooted-tree-filt50.qza`, `rooted-tree-filt50.qza`
- PERMANOVA:
  - `braycurtis-permanova-filt50.qzv`
  - `weighted-unifrac-permanova-filt50.qzv`
- Emperor plots:
  - `braycurtis-emperor-filt50.qzv`
  - `weighted-unifrac-emperor-filt50.qzv`

---

## Taxonomy + genus-level figures (Task 2B)

### Top500 taxonomy strategy (used for genus plots)
To generate genus-level reporting figures reliably on local resources, taxonomy was assigned for the **Top500 ASVs** only and then exported to R for plotting. 
Taxonomy assignment was performed in QIIME 2 using the VSEARCH consensus classifier (`feature-classifier classify-consensus-vsearch`).

Artifacts:
- `rep-seqs-top500.qza`
- `table-top500.qza`
- `taxonomy-top500.qza`
- `taxonomy-top500-hits.qza`

R inputs (already present in this repo):
- `r_input/feature-table-top500.tsv`
- `r_input/exported-taxonomy-top500/taxonomy.tsv`

### Differential abundance (ANCOM‑BC) — FINAL Figure 2
Differential abundance testing was performed at genus level using QIIME 2 `composition ancombc` with `Status` as the formula, producing FDR-corrected q-values (used for significance).   
ANCOM‑BC results are saved as `ancombc_genus_top500.qza` and viewable as `ancombc_genus_top500.qzv`. 

---

## Quickstart

### 1) Activate QIIME 2
conda activate qiime2-amplicon-2024.10
qiime --help


### 2) Export Top500 inputs for R (only if needed)
If `r_input/feature-table-top500.tsv` and `r_input/exported-taxonomy-top500/taxonomy.tsv` are already present, you can skip this. 

qiime tools export
--input-path taxonomy-top500.qza
--output-path r_input/exported-taxonomy-top500

qiime tools export
--input-path table-top500.qza
--output-path r_input/exported-table-top500

biom convert
-i r_input/exported-table-top500/feature-table.biom
-o r_input/feature-table-top500.tsv
--to-tsv


### 3) Generate genus composition + log2FC effect size (R)
Rscript scripts/plot_top500_genus.R

text

### 4) Run ANCOM‑BC (genus level) + export results for R plotting
qiime taxa collapse
--i-table table-top500.qza
--i-taxonomy taxonomy-top500.qza
--p-level 6
--o-collapsed-table genus-top500.qza

qiime feature-table filter-features
--i-table genus-top500.qza
--p-min-frequency 10
--o-filtered-table genus-top500-filt.qza

qiime composition ancombc
--i-table genus-top500-filt.qza
--m-metadata-file metadata.txt
--p-formula "Status"
--o-differentials ancombc_genus_top500.qza

qiime composition tabulate
--i-data ancombc_genus_top500.qza
--o-visualization ancombc_genus_top500.qzv

qiime tools export
--input-path ancombc_genus_top500.qza
--output-path r_input/ancombc_genus_top500_export


### 5) Plot ANCOM‑BC volcano + significant barplot (R)
Rscript scripts/plot_ancombc_genus.R


---

## Outputs

### Final figures (use these in the report)
Saved under `results/figures/`:
- `Figure1_top500_genus_relative_abundance.png`
- `Figure2_ANCOMBC_volcano.png`
- `Figure2_ANCOMBC_top_sig_barplot.png` 

### Legacy figures (kept for reference)
Saved under `results/figures/`:
- `Figure1_genus_relative_abundance.png`
- `Figure2_genus_log2FC_simple.png`
- `Figure2_top500_genus_log2FC_fixedLegend.png`

### Final tables (Top500 genus summaries)
Saved under `results/tables/`:
- `genus_counts_long_top500.tsv`
- `genus_counts_wide_top500.tsv`
- `genus_log2FC_top500.tsv`

### Differential abundance tables (ANCOM‑BC)
Saved under `results/tables/`:
- `ancombc_genus_top500_results.tsv`
- `ancombc_genus_top20_sig.tsv` 

### Other/legacy tables
Saved under `results/tables/`:
- `asv_to_genus.tsv`
- `genus_counts_long.tsv`
- `genus_counts_wide.tsv`
- `genus_log2FC_simple.tsv`

---

## Viewing `.qzv` results
QIIME 2 visualizations can be opened in the browser using QIIME 2 View:
https://view.qiime2.org/ 
If running locally with a display:
qiime tools view core-metrics-no-tree/bray_curtis_emperor.qzv

---

## Notes / limitations
- Genus-level reporting figures are based on **Top500 ASVs** taxonomy assignment due to local resource constraints. 
- Differential abundance (Figure 2) was performed using ANCOM‑BC with FDR-corrected q-values. 
- The older log2FC plot is kept as a descriptive effect-size (legacy), and should not be interpreted as a statistical differential abundance test. 

---

## License
MIT License is recommended for academic code sharing. 
