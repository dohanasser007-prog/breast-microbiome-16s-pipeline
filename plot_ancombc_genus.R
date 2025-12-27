suppressPackageStartupMessages({
  library(ggplot2)
})

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

ancom_dir <- "r_input/ancombc_genus_top500_export"

lfc_fp <- file.path(ancom_dir, "lfc_slice.csv")
q_fp   <- file.path(ancom_dir, "q_val_slice.csv")
p_fp   <- file.path(ancom_dir, "p_val_slice.csv")
se_fp  <- file.path(ancom_dir, "se_slice.csv")
w_fp   <- file.path(ancom_dir, "w_slice.csv")

stopifnot(file.exists(lfc_fp), file.exists(q_fp), file.exists(p_fp), file.exists(se_fp), file.exists(w_fp))

lfc <- read.csv(lfc_fp, check.names = FALSE)
qv  <- read.csv(q_fp,  check.names = FALSE)
pv  <- read.csv(p_fp,  check.names = FALSE)
se  <- read.csv(se_fp, check.names = FALSE)
w   <- read.csv(w_fp,  check.names = FALSE)

# Each file typically has: id, (Intercept), StatusHealthy_control
needed_cols <- c("id", "StatusHealthy_control")
for (df in list(lfc, qv, pv, se, w)) {
  if (!all(needed_cols %in% colnames(df))) {
    stop("Expected columns not found. Needed: id and StatusHealthy_control")
  }
}

res <- data.frame(
  id = lfc$id,
  lfc = lfc$StatusHealthy_control,
  q_val = qv$StatusHealthy_control,
  p_val = pv$StatusHealthy_control,
  se = se$StatusHealthy_control,
  W = w$StatusHealthy_control,
  stringsAsFactors = FALSE
)

# Try to extract genus name from taxonomy string if present; otherwise keep id
res$Genus <- sub(".*;\\s*g__([^;]+).*", "\\1", res$id)
res$Genus[res$Genus == res$id] <- res$id

res$neglog10q <- -log10(res$q_val + 1e-300)
res$Significant <- is.finite(res$q_val) & (res$q_val < 0.05)

write.table(res, file = "results/tables/ancombc_genus_top500_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Volcano
p_volcano <- ggplot(res, aes(x = lfc, y = neglog10q, color = Significant)) +
  geom_point(alpha = 0.8) +
  theme_bw() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  xlab("ANCOM-BC log fold change (Healthy_control vs Breast_cancer)") +
  ylab("-log10(q-value)") +
  ggtitle("Differential abundance (ANCOM-BC), genus level")

ggsave("results/figures/Figure2_ANCOMBC_volcano.png",
       plot = p_volcano, width = 7, height = 5, dpi = 300)

# Top significant barplot
sig <- res[res$Significant & is.finite(res$lfc), ]
sig <- sig[order(-abs(sig$lfc)), ]
topn <- head(sig, 20)

if (nrow(topn) > 0) {
  topn$Direction <- ifelse(topn$lfc > 0, "Enriched in Healthy", "Enriched in Cancer")
  
  p_bar <- ggplot(topn, aes(x = reorder(Genus, lfc), y = lfc, fill = Direction)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    scale_fill_manual(values = c("Enriched in Healthy" = "blue", "Enriched in Cancer" = "red")) +
    xlab("Genus") +
    ylab("ANCOM-BC log fold change (Healthy_control vs Breast_cancer)") +
    ggtitle("Top significant genera (q < 0.05)")
  
  ggsave("results/figures/Figure2_ANCOMBC_top_sig_barplot.png",
         plot = p_bar, width = 8, height = 6, dpi = 300)
  
  write.table(topn, file = "results/tables/ancombc_genus_top20_sig.tsv",
              sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  message("No significant genera at q < 0.05. Barplot not generated.")
}

cat("Done. Check results/figures and results/tables\n")
