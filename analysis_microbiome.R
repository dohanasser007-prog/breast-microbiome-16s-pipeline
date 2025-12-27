suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ggplot2)
})

cat(">>> R analysis started\n")

# ---- Inputs (WSL project paths) ----
feature_tsv <- "r_input/feature-table.tsv"
dna_fasta   <- "r_input/dna-sequences.fasta"
metadata_fp <- "metadata.txt"   # use root metadata
silva_rdata <- "D:/R_microbiome/Silva_SSU_r138_1_training.RData"

fig_dir <- "results/figures"
tab_dir <- "results/tables"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(feature_tsv))
stopifnot(file.exists(dna_fasta))
stopifnot(file.exists(metadata_fp))
stopifnot(file.exists(silva_rdata))

# ---- Metadata ----
metadata <- read.table(metadata_fp, header = TRUE, sep = "\t",
                       check.names = FALSE, comment.char = "",
                       stringsAsFactors = FALSE)

if (ncol(metadata) < 2) stop("metadata.txt must have at least 2 columns: SampleID and Status")
metadata <- metadata[, 1:2]
colnames(metadata) <- c("SampleID", "Status")

# ---- Feature table (TSV exported from QIIME2; often includes a comment line) ----
feat_tab_raw <- read.table(feature_tsv, header = TRUE, sep = "\t",
                           check.names = FALSE, comment.char = "",
                           skip = 1, stringsAsFactors = FALSE)

colnames(feat_tab_raw)[1] <- "FeatureID"
rownames(feat_tab_raw) <- feat_tab_raw$FeatureID
feat_tab_raw$FeatureID <- NULL

count_mat <- as.matrix(feat_tab_raw)
mode(count_mat) <- "numeric"
cat("Count matrix dim:", paste(dim(count_mat), collapse = " x "), "\n")

# ---- ASV sequences ----
asv_seqs <- readDNAStringSet(dna_fasta)
cat("ASV sequences:", length(asv_seqs), "\n")

# ---- Ensure matching IDs ----
common_asvs <- intersect(rownames(count_mat), names(asv_seqs))
cat("Common ASVs between table and FASTA:", length(common_asvs), "\n")

if (length(common_asvs) < 1000) {
  warning("Low overlap between feature table and FASTA IDs. Check QIIME2 exports.")
}

count_mat <- count_mat[common_asvs, , drop = FALSE]
asv_seqs  <- asv_seqs[common_asvs]

# ---- Load SILVA training set (DECIPHER) ----
load(silva_rdata)
training_obj <- NULL
if (exists("TrainingSet_16S")) training_obj <- get("TrainingSet_16S")
if (is.null(training_obj) && exists("trainingSet")) training_obj <- get("trainingSet")
if (is.null(training_obj)) stop("No TrainingSet_16S/trainingSet object found in SILVA RData.")

# ---- IdTaxa on abundant ASVs only ----
asv_totals <- rowSums(count_mat)
min_total_reads_for_tax <- 50
keep_ids <- names(asv_totals[asv_totals >= min_total_reads_for_tax])

asv_seqs_keep <- asv_seqs[names(asv_seqs) %in% keep_ids]
cat("ASVs used for taxonomy (>=50 reads):", length(asv_seqs_keep), "\n")

tax_res <- IdTaxa(asv_seqs_keep,
                  trainingSet = training_obj,
                  strand = "both",
                  processors = 2,
                  verbose = TRUE)

tax_genus <- setNames(
  sapply(tax_res, function(x) {
    if (!is.null(x$rank) && "Genus" %in% x$rank) {
      x$taxon[match("Genus", x$rank)]
    } else {
      NA_character_
    }
  }),
  names(asv_seqs_keep)
)

# ---- Map genus to ALL ASVs; missing/NA => Unclassified ----
genus_map <- rep("Unclassified", nrow(count_mat))
names(genus_map) <- rownames(count_mat)
genus_map[names(tax_genus)] <- ifelse(is.na(tax_genus), "Unclassified", tax_genus)
genus_map[is.na(genus_map)] <- "Unclassified"

write.table(
  data.frame(ASV = names(genus_map), Genus = genus_map),
  file = file.path(tab_dir, "asv_to_genus.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- Genus aggregation (fast) ----
genus_levels <- unique(genus_map)
genus_counts_mat <- sapply(genus_levels, function(g) {
  colSums(count_mat[genus_map == g, , drop = FALSE])
})
genus_counts_mat <- t(genus_counts_mat)

genus_counts_df <- data.frame(
  Genus = rownames(genus_counts_mat),
  genus_counts_mat,
  check.names = FALSE,
  row.names = NULL
)

write.table(genus_counts_df,
            file = file.path(tab_dir, "genus_counts_wide.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Long for plotting (safe; avoids reshape() issues) ----
m <- as.matrix(genus_counts_df[, -1, drop = FALSE])
rownames(m) <- genus_counts_df$Genus

tmp <- stack(as.data.frame(t(m), check.names = FALSE))
colnames(tmp) <- c("Count", "Genus")
tmp$SampleID <- rep(rownames(t(m)), times = ncol(t(m)))

genus_long <- tmp[, c("SampleID", "Genus", "Count")]
genus_long$Count <- as.numeric(genus_long$Count)

# Join metadata and keep valid samples
genus_long <- merge(genus_long, metadata, by = "SampleID", all.x = TRUE)
genus_long <- genus_long[!is.na(genus_long$Status), ]

write.table(genus_long,
            file = file.path(tab_dir, "genus_counts_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Figure 1: Relative abundance ----
sample_totals <- tapply(genus_long$Count, genus_long$SampleID, sum)
genus_long$RelAbundance <- genus_long$Count / sample_totals[genus_long$SampleID]

mean_rel <- tapply(genus_long$RelAbundance, genus_long$Genus, mean)
top20 <- names(sort(mean_rel, decreasing = TRUE))[1:min(20, length(mean_rel))]
genus_long$Genus_plot <- ifelse(genus_long$Genus %in% top20, genus_long$Genus, "Other")

p1 <- ggplot(genus_long, aes(x = SampleID, y = RelAbundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Status, scales = "free_x") +
  ylab("Relative abundance") + xlab("Samples") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(fig_dir, "Figure1_genus_relative_abundance.png"),
       plot = p1, width = 10, height = 6, dpi = 300)

# ---- Figure 2: Simple log2FC (effect size) ----
mean_rel_by_status <- aggregate(RelAbundance ~ Status + Genus, data = genus_long, FUN = mean)

cancer  <- subset(mean_rel_by_status, Status == "Breast_cancer",   select = c("Genus", "RelAbundance"))
healthy <- subset(mean_rel_by_status, Status == "Healthy_control", select = c("Genus", "RelAbundance"))
colnames(cancer)[2]  <- "Breast_cancer"
colnames(healthy)[2] <- "Healthy_control"

genus_wide <- merge(cancer, healthy, by = "Genus", all = TRUE)
genus_wide[is.na(genus_wide)] <- 0
genus_wide$log2FC <- log2((genus_wide$Breast_cancer + 1e-6) / (genus_wide$Healthy_control + 1e-6))

write.table(genus_wide,
            file = file.path(tab_dir, "genus_log2FC_simple.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top_da <- head(genus_wide[order(-abs(genus_wide$log2FC)), ], 20)

p2 <- ggplot(top_da, aes(x = reorder(Genus, log2FC), y = log2FC, fill = log2FC > 0)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                    labels = c("FALSE" = "Enriched in Healthy", "TRUE" = "Enriched in Cancer"),
                    name = "") +
  ylab("log2 Fold Change (Cancer / Healthy)") +
  xlab("Genus") +
  theme_bw()

ggsave(file.path(fig_dir, "Figure2_genus_log2FC_simple.png"),
       plot = p2, width = 7, height = 6, dpi = 300)

cat(">>> R analysis finished\n")
