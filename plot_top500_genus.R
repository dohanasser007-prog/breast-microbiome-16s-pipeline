suppressPackageStartupMessages({
  library(ggplot2)
})

cat(">>> plot_top500_genus.R started\n")

# ---- Inputs (from QIIME2 exports) ----
feature_tsv  <- "r_input/feature-table-top500.tsv"
taxonomy_tsv <- "r_input/exported-taxonomy-top500/taxonomy.tsv"
metadata_fp  <- "metadata.txt"

fig_dir <- "results/figures"
tab_dir <- "results/tables"

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(feature_tsv))
stopifnot(file.exists(taxonomy_tsv))
stopifnot(file.exists(metadata_fp))

# ---- Metadata ----
meta <- read.table(metadata_fp, header = TRUE, sep = "\t",
                   check.names = FALSE, comment.char = "",
                   stringsAsFactors = FALSE)
if (ncol(meta) < 2) stop("metadata.txt must have at least 2 columns: SampleID and Status")
meta <- meta[, 1:2]
colnames(meta) <- c("SampleID", "Status")

# ---- Taxonomy (QIIME2 export) ----
tax <- read.table(taxonomy_tsv, header = TRUE, sep = "\t",
                  quote = "", check.names = FALSE, comment.char = "")
stopifnot(all(c("Feature ID", "Taxon") %in% colnames(tax)))

# Parse genus from SILVA-style string containing g__Genus
tax$Genus <- sub(".*;\\s*g__([^;]+).*", "\\1", tax$Taxon)
tax$Genus[tax$Genus == tax$Taxon] <- "Unclassified"
tax$Genus[is.na(tax$Genus) | tax$Genus == "" | tax$Genus == "NA"] <- "Unclassified"

genus_map <- setNames(tax$Genus, tax[["Feature ID"]])

# ---- Feature table TSV (BIOM -> TSV has a comment line) ----
tab <- read.table(feature_tsv, header = TRUE, sep = "\t",
                  check.names = FALSE, comment.char = "",
                  skip = 1, stringsAsFactors = FALSE)

colnames(tab)[1] <- "FeatureID"
rownames(tab) <- tab$FeatureID
tab$FeatureID <- NULL

count_mat <- as.matrix(tab)
mode(count_mat) <- "numeric"

# Keep only ASVs present in taxonomy (Top500 consistency)
common_asvs <- intersect(rownames(count_mat), names(genus_map))
count_mat <- count_mat[common_asvs, , drop = FALSE]
genus_vec <- genus_map[common_asvs]
genus_vec[is.na(genus_vec)] <- "Unclassified"

# ---- Aggregate to genus x sample ----
genus_levels <- unique(genus_vec)
genus_counts <- sapply(genus_levels, function(g) colSums(count_mat[genus_vec == g, , drop = FALSE]))
genus_counts <- t(genus_counts)

genus_counts_df <- data.frame(Genus = rownames(genus_counts), genus_counts, check.names = FALSE)
write.table(genus_counts_df,
            file = file.path(tab_dir, "genus_counts_wide_top500.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Long format ----
m <- as.matrix(genus_counts_df[, -1, drop = FALSE])
rownames(m) <- genus_counts_df$Genus

tmp <- stack(as.data.frame(t(m), check.names = FALSE))
colnames(tmp) <- c("Count", "Genus")
tmp$SampleID <- rep(rownames(t(m)), times = ncol(t(m)))

genus_long <- tmp[, c("SampleID", "Genus", "Count")]
genus_long$Count <- as.numeric(genus_long$Count)

genus_long <- merge(genus_long, meta, by = "SampleID", all.x = TRUE)
genus_long <- genus_long[!is.na(genus_long$Status), ]

write.table(genus_long,
            file = file.path(tab_dir, "genus_counts_long_top500.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Relative abundance ----
sample_totals <- tapply(genus_long$Count, genus_long$SampleID, sum)
genus_long$RelAbundance <- genus_long$Count / sample_totals[genus_long$SampleID]

# ---- Figure 1: Relative abundance (Top20 + Other) ----
mean_rel <- tapply(genus_long$RelAbundance, genus_long$Genus, mean)
top20 <- names(sort(mean_rel, decreasing = TRUE))[1:min(20, length(mean_rel))]
genus_long$Genus_plot <- ifelse(genus_long$Genus %in% top20, genus_long$Genus, "Other")

p1 <- ggplot(genus_long, aes(x = SampleID, y = RelAbundance, fill = Genus_plot)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Status, scales = "free_x") +
  ylab("Relative abundance") + xlab("SampleID") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(file.path(fig_dir, "Figure1_top500_genus_relative_abundance.png"),
       plot = p1, width = 10, height = 6, dpi = 300)

# ---- Figure 2: log2FC (Cancer/Healthy) ----
mean_rel_by_status <- aggregate(RelAbundance ~ Status + Genus, data = genus_long, FUN = mean)

cancer  <- subset(mean_rel_by_status, Status == "Breast_cancer",   select = c("Genus", "RelAbundance"))
healthy <- subset(mean_rel_by_status, Status == "Healthy_control", select = c("Genus", "RelAbundance"))
colnames(cancer)[2]  <- "Breast_cancer"
colnames(healthy)[2] <- "Healthy_control"

wide <- merge(cancer, healthy, by = "Genus", all = TRUE)
wide[is.na(wide)] <- 0
wide$log2FC <- log2((wide$Breast_cancer + 1e-6) / (wide$Healthy_control + 1e-6))
wide$Enrichment <- ifelse(wide$log2FC > 0, "Enriched in Cancer", "Enriched in Healthy")

write.table(wide,
            file = file.path(tab_dir, "genus_log2FC_top500.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

top_da <- head(wide[order(-abs(wide$log2FC)), ], 20)

p2 <- ggplot(top_da, aes(x = reorder(Genus, log2FC), y = log2FC, fill = Enrichment)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(
    name = "",
    values = c("Enriched in Healthy" = "blue", "Enriched in Cancer" = "red")
  ) +
  ylab("log2 Fold Change (Cancer / Healthy)") +
  xlab("Genus")

ggsave(file.path(fig_dir, "Figure2_top500_genus_log2FC_fixedLegend.png"),
       plot = p2, width = 7, height = 6, dpi = 300)

cat(">>> plot_top500_genus.R finished\n")
