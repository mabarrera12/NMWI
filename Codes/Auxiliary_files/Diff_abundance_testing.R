# Differential abundance testing
library(phyloseq)
library(ANCOMBC)
library(ggplot2)
library(ggrepel)

data <- read.table('raw_counts.tsv', sep = '\t', header = TRUE, check.names = FALSE)
meta <- data[, !grepl("^d__", colnames(data))]
taxa <- data[, grepl("^d__", colnames(data))] # columns are taxa; Data is counts

rownames(taxa) <- rownames(meta) <- data$Sample_id

meta$Category <- factor(meta$Category)
meta$Category <- relevel(meta$Category, ref = "control")


# have to build a phyloseq object
otu <- otu_table(as.matrix(taxa), taxa_are_rows = FALSE)
sam <- sample_data(meta)
ps <- phyloseq(otu, sam)

out <- ancombc2(
  data = ps,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "Category",
  rand_formula = NULL,
  p_adj_method = "holm",
  prv_cut = 0.15,
  lib_cut = 0,
  group = "Category",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE
)

res <- data.frame(taxon = rownames(out$res), out$res)

sig_global <- res[res$q_Categorycase < 0.05, ]
sig_global <- sig_global[order(sig_global$q_Categorycase), ]


res$logFC <- res$lfc_Categorycase
res$negLogQ <- -log10(res$q_Categorycase)
res$negLogQ[is.infinite(res$negLogQ)] <- NA
res$significance <- "Not significant"
res$significance[res$q_Categorycase < 0.05 & res$logFC > 0] <- "Higher in case"
res$significance[res$q_Categorycase < 0.05 & res$logFC < 0] <- "Lower in case"

sif_data <- res[res$q_Categorycase < 0.05 & (res$lfc_Categorycase > 0.5 | res$lfc_Categorycase < -0.5 ), ]
top_taxa <- res[res$q_Categorycase < 0.05, ]
top_taxa <- top_taxa[order(top_taxa$q_Categorycase), ][1:10, ]

ggplot(res, aes(x = logFC, y = negLogQ, color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = top_taxa,
    aes(label = taxon.1),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c(
    "Higher in case" = "red",
    "Lower in case" = "blue",
    "Not significant" = "grey"
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = 'dashed')+
  geom_vline(xintercept = 0.5, linetype = 'dashed')+
  theme_minimal() +
  labs(
    x = "Log Fold Change (case vs control)",
    y = "-log10(q-value)",
    title = "ANCOM-BC2 Volcano Plot"
  )

write.csv(sif_data, "sig_data_genus.csv", row.names = FALSE)


# ------------ Family
t.taxa <- t(taxa)
new_names <- sapply(strsplit(rownames(t.taxa), "\\;"), function(x) {paste(x[-length(x)], collapse = ";")})
tax_grouped <- aggregate(t.taxa, by = list(new_names), FUN = sum)
rownames(tax_grouped) <- tax_grouped$Group.1
tax_grouped <- tax_grouped[,2:ncol(tax_grouped)]
taxa <- t(tax_grouped)

# have to build a phyloseq object
otu <- otu_table(as.matrix(taxa), taxa_are_rows = FALSE)
sam <- sample_data(meta)
ps <- phyloseq(otu, sam)

out <- ancombc2(
  data = ps,
  assay_name = "counts",
  tax_level = NULL,
  fix_formula = "Category",
  rand_formula = NULL,
  p_adj_method = "holm",
  prv_cut = 0.15,
  lib_cut = 0,
  group = "Category",
  struc_zero = TRUE,
  neg_lb = TRUE,
  alpha = 0.05,
  n_cl = 1,
  verbose = TRUE
)

res <- data.frame(taxon = rownames(out$res), out$res)

sig_global <- res[res$q_Categorycase < 0.05, ]
sig_global <- sig_global[order(sig_global$q_Categorycase), ]


res$logFC <- res$lfc_Categorycase
res$negLogQ <- -log10(res$q_Categorycase)
res$negLogQ[is.infinite(res$negLogQ)] <- NA
res$significance <- "Not significant"
res$significance[res$q_Categorycase < 0.05 & res$logFC > 0] <- "Higher in case"
res$significance[res$q_Categorycase < 0.05 & res$logFC < 0] <- "Lower in case"

library(ggplot2)
library(ggrepel)
sif_data <- res[res$q_Categorycase < 0.05 & (res$lfc_Categorycase > 0.5 | res$lfc_Categorycase < -0.5 ), ]
top_taxa <- res[res$q_Categorycase < 0.05, ]
top_taxa <- top_taxa[order(top_taxa$q_Categorycase), ][1:10, ]

ggplot(res, aes(x = logFC, y = negLogQ, color = significance)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_text_repel(
    data = top_taxa,
    aes(label = taxon.1),
    size = 3,
    max.overlaps = 20
  ) +
  scale_color_manual(values = c(
    "Higher in case" = "red",
    "Lower in case" = "blue",
    "Not significant" = "grey"
  )) +geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = -0.5, linetype = 'dashed')+
  geom_vline(xintercept = 0.5, linetype = 'dashed')+
  theme_minimal() +
  labs(
    x = "Log Fold Change (case vs control)",
    y = "-log10(q-value)",
    title = "ANCOM-BC2 Volcano Plot"
  )

write.csv(sif_data, "sig_data_family.csv", row.names = FALSE)
