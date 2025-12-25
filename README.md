# breast-microbiome-16s-pipeline
microbiome analysis (Task 3)
A reproducible QIIME 2 + R workflow for 16S rRNA amplicon analysis of breast microbiome samples (Breast_cancer vs Healthy_control), producing ASVs, diversity metrics (Bray–Curtis, Jaccard, and Weighted UniFrac workaround), and genus-level summary figures/tables. 

## Quickstart (recommended)
From the project root on WSL/Ubuntu:
conda activate qiime2-amplicon-2024.10
chmod +x scripts/run_pipeline.sh
./scripts/run_pipeline.sh 2>&1 | tee logs/pipeline.log
Then on Windows (R installed):
Rscript scripts/analysis_microbiome.R
r_input/feature-table.tsv
r_input/dna-sequences.fasta
metadata.txt
D:/R_microbiome/Silva_SSU_r138_1_training.RData
results/figures
results/tables

Open `.qzv` outputs using QIIME 2 View to export plots/tables for your report:
https://view.qiime2.org 

---

## Project Description
This project implements the Task 3 microbiome pipeline using QIIME 2 (WSL/Ubuntu) for denoising and diversity analysis and an R module (Windows) for taxonomy assignment and publication-style plots. 
QIIME 2 results are stored as `.qza` artifacts and `.qzv` visualizations, designed for reproducibility and easy sharing. 
Weighted UniFrac was computed using a filtered feature table (`table-filt50`) because tree building on the full ASV set was limited by available RAM.

---

## Installation Instructions

### A) QIIME 2 (WSL/Ubuntu)
Install QIIME 2 Amplicon 2024.10 using the official conda-based instructions for Linux/WSL. [web:300][web:301]  
After installation, activate with `conda activate qiime2-amplicon-2024.10` and test with `qiime --help`. 

### B) R (Windows)
Install R/RStudio and required packages for the R module, including `DECIPHER::IdTaxa` for taxonomy assignment.  
Provide the SILVA DECIPHER training set file used by the script (example): `D:/R_microbiome/Silva_SSU_r138_1_training.RData`. 

---

## Getting Started

### Expected folder structure
breast_microbiome/
data/ # FASTQ files
metadata.txt # SampleID + Status
manifest_batch1.tsv
manifest_batch2a.tsv
manifest_batch2b.tsv

scripts/
run_pipeline.sh
analysis_microbiome.R

r_input/ # exported inputs for R
feature-table.biom
feature-table.tsv
dna-sequences.fasta

results/
figures/
tables/

logs/
pipeline.log


### Metadata format
`metadata.txt` must contain two columns:
- `SampleID`
- `Status` with levels: `Breast_cancer` and `Healthy_control`

---

## Usage Guide

### 1) QIIME 2 workflow (WSL)
The `scripts/run_pipeline.sh` script generates the following key outputs (from this run):

**Demultiplexing**
- `demux-single-end-amplicon.qza`
- `demux-single-end-amplicon.qzv`

**Merged DADA2 results**
- `table-all.qza`, `table-all.qzv`
- `rep-seqs-all.qza`, `rep-seqs-all.qzv`
- Batch-level artifacts: `table-batch1.qza`, `table-batch2a.qza`, `table-batch2b.qza`, and their `rep-seqs-*` and `denoising-stats-*` artifacts.

**Non-phylogenetic diversity (core-metrics-no-tree)**
Generated under `core-metrics-no-tree/` (examples from this run):
- `bray_curtis_distance_matrix.qza`, `bray_curtis_pcoa_results.qza`, `bray_curtis_emperor.qzv`
- `bray-curtis-Status-significance.qzv`
- `jaccard_distance_matrix.qza`, `jaccard_pcoa_results.qza`, `jaccard_emperor.qzv`
- `observed_features_vector.qza`, `shannon_vector.qza`, `evenness_vector.qza`
- `shannon-Status-significance.qzv`
- `rarefied_table.qza`

**Weighted UniFrac workaround (filtered table + tree)**
This run includes the tree + phylogenetic beta diversity outputs:
- `table-filt50.qza`, `table-filt50.qzv`, `rep-seqs-filt50.qza`
- `aligned-rep-seqs-filt50.qza`, `masked-aligned-rep-seqs-filt50.qza`
- `unrooted-tree-filt50.qza`, `rooted-tree-filt50.qza`
- `weighted-unifrac-distance-filt50.qza`, `weighted-unifrac-pcoa-filt50.qza`
- `weighted-unifrac-emperor-filt50.qzv`, `weighted-unifrac-permanova-filt50.qzv`
- Companion Bray–Curtis (filt50): `braycurtis-distance-filt50.qza`, `braycurtis-emperor-filt50.qzv`, `braycurtis-permanova-filt50.qzv` 

### 2) R workflow (Windows)
The R script assigns taxonomy using `DECIPHER::IdTaxa` and produces genus-level figures/tables.  
Figure 2 is computed as a simple log2 fold-change from mean relative abundance (effect size), and can be replaced later by ANCOM-BC/DESeq2 if needed. 

### 3) Viewing/exporting results
Use QIIME 2 View to open `.qzv` files and export images/tables for reporting. 

---

## Outputs

### Final figures (R)
Saved under `results/figures/`:
- `Figure1_genus_relative_abundance.png`
- `Figure2_genus_log2FC_simple.png`

### Final tables (R)
Saved under `results/tables/`:
- `genus_counts_long.tsv`
- `genus_log2FC_simple.tsv`
- `alpha_diversity.tsv`
- `alpha_wilcoxon_results.txt`

### Key interactive QIIME 2 visualizations
Open these in QIIME 2 View:
- `core-metrics-no-tree/bray_curtis_emperor.qzv`
- `core-metrics-no-tree/jaccard_emperor.qzv`
- `braycurtis-emperor-filt50.qzv`
- `weighted-unifrac-emperor-filt50.qzv`
- `braycurtis-permanova-filt50.qzv`
- `weighted-unifrac-permanova-filt50.qzv` 

---

## Contributing Guidelines
Contributions are welcome:
- Fork the repo and create a feature branch.
- Open a pull request with a clear description and testing notes.
- Report issues via GitHub Issues and attach relevant logs (e.g., `logs/pipeline.log`).

---

## License Information
MIT License (recommended for academic code sharing).  
Add a `LICENSE` file in the repository root containing the MIT license text.

---

## Acknowledgments and Credits
- QIIME 2 for the microbiome analysis framework and artifact/visualization system. 
- QIIME 2 View for viewing/exporting `.qza/.qzv` results. 
- DECIPHER for taxonomic classification using `IdTaxa`. 

