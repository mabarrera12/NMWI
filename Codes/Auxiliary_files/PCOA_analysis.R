# Load packages -----------------------------------------------------------
library(vegan)
library(ade4)
library(ggplot2)

# Load input data ---------------------------------------------------------
training_set <- read.delim("Data/training_set.tsv", sep = "\t", header = TRUE, check.names = FALSE)

# Select only taxonomy abundance columns (columns starting with "d__")
taxa_mat <- training_set[, grepl("^d__", colnames(training_set))]

# Recode group labels for case/control comparison
# control = TRUE, case = FALSE
group <- c(case = FALSE, control = TRUE)[training_set$Category]

# Build data frame for PERMANOVA
permanova_df <- data.frame(
  Group = group,
  BioProject = training_set$BioProject.Number,
  taxa_mat,
  check.names = FALSE
)

# Identify taxa columns only
taxa_idx <- which(!(colnames(permanova_df) %in% c("Group", "BioProject")))

# Calculate Bray-Curtis distance matrix from relative abundance data
bc_dist <- vegdist(permanova_df[, taxa_idx], method = "bray")

# Run PERMANOVA, stratifying by BioProject
adonis_res <- adonis2(
  bc_dist ~ Group,
  data = permanova_df,
  permutations = 999,
  strata = permanova_df$BioProject
)

print(adonis_res)

# Extract R-squared and p-value for the Group term
adon_df <- as.data.frame(adonis_res)
r2_val <- if ("Group" %in% rownames(adon_df)) adon_df["Group", "R2"] else NA
p_val  <- if ("Group" %in% rownames(adon_df)) adon_df["Group", "Pr(>F)"] else NA

# Run PCoA ---------------------------------------------------------------
pcoa_res <- dudi.pco(bc_dist, scannf = FALSE, nf = 3)

# Keep sample metadata for plotting
pcoa_meta <- data.frame(
  Group_all = as.character(training_set$Category),
  Sample_ID = training_set$Sample_id,
  row.names = rownames(permanova_df)
)

# Percent variance explained by first two axes
eigenvalues <- pcoa_res$eig
variance_explained <- eigenvalues / sum(eigenvalues)
variance1 <- 100 * signif(variance_explained[1], 2)
variance2 <- 100 * signif(variance_explained[2], 2)

# Create plotting data frame
pc_plot_data <- data.frame(
  x = pcoa_res$li$A1,
  y = pcoa_res$li$A2,
  Group = permanova_df$Group,
  Group_all = pcoa_meta$Group_all,
  Sample_ID = pcoa_meta$Sample_ID,
  Phenotype = training_set$Phenotype
)

# Add group centroids
pc_plot_data_centroids <- merge(
  pc_plot_data,
  aggregate(cbind(mean.x = x, mean.y = y) ~ Group, pc_plot_data, mean),
  by = "Group"
)

# Plot PCoA colored by case/control status -------------------------------
fig_case_control <- ggplot(
  pc_plot_data_centroids,
  aes(x, y, color = factor(Group_all))
) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95) +
  labs(
    x = paste0("PCoA1 (", variance1, "%)"),
    y = paste0("PCoA2 (", variance2, "%)")
  ) +
  theme_bw() +
  theme(legend.title = element_blank()) +
  scale_colour_manual(values = c(
    "Healthy" = "steelblue",
    "Nonhealthy" = "orange2",
    "control" = "steelblue",
    "case" = "orange2"
  ))

print(
  fig_case_control +
    annotate(
      "text",
      label = paste0(
        "R-squared = ", signif(r2_val, 4),
        "\n",
        "p-value = ", signif(p_val, 3)
      ),
      x = min(pc_plot_data_centroids$x, na.rm = TRUE),
      y = max(pc_plot_data_centroids$y, na.rm = TRUE),
      hjust = 0,
      vjust = 1,
      color = "black"
    )
)

# Plot PCoA colored by phenotype -----------------------------------------
fig_phenotype <- ggplot(pc_plot_data, aes(x, y, color = factor(Phenotype))) +
  geom_point(size = 2) +
  labs(
    x = paste0("PCoA1 (", variance1, "%)"),
    y = paste0("PCoA2 (", variance2, "%)")
  ) +
  theme_bw() +
  theme(legend.title = element_blank())

print(fig_phenotype)

# Test effect of sampling region and sequencing platform -----------------
permanova_full <- data.frame(
  Group = group,
  BioProject = training_set$BioProject.Number,
  training_set,
  check.names = FALSE
)

# Reuse only taxa columns for Bray-Curtis distance
taxa_idx <- which(colnames(permanova_full) %in% colnames(taxa_mat))
bc_dist <- vegdist(permanova_full[, taxa_idx], method = "bray")

# Test whether 16S region explains microbiome composition
adonis_region <- adonis2(
  bc_dist ~ `16S Region`,
  data = permanova_full,
  permutations = 999
)
print(adonis_region)

# Test whether BioProject and sequencing platform explain variation
adonis_platform <- adonis2(
  bc_dist ~ BioProject + `Sequencing platform`,
  data = permanova_full,
  permutations = 999
)
print(adonis_platform)

# Restrict to healthy/control samples only -------------------------------
healthy_set <- subset(training_set, Category == "control")

# Extract taxa columns from healthy samples
healthy_taxa <- healthy_set[, grepl("^d__", colnames(healthy_set))]

# Build data frame for healthy-only PERMANOVA
permanova_hv <- data.frame(
  BioProject = healthy_set$BioProject.Number,
  HypervariableRegion = healthy_set$`16S Region`,
  healthy_taxa
)

# Identify taxa columns
taxa_idx <- which(!(colnames(permanova_hv) %in% c("BioProject", "HypervariableRegion")))

# Compute Bray-Curtis distances for healthy-only samples
bc_dist_healthy <- vegdist(permanova_hv[, taxa_idx], method = "bray")

# Test BioProject and hypervariable region effects in healthy controls
adonis_hv <- adonis2(
  bc_dist_healthy ~ BioProject + HypervariableRegion,
  data = permanova_hv,
  permutations = 999
)
print(adonis_hv)

# Quick check of healthy samples by hypervariable region
table(permanova_hv$HypervariableRegion)

# Test BioProject and hypervariable region effects in healthy controls
adonis_hv <- adonis2(
  bc_dist_healthy ~  HypervariableRegion,
  data = permanova_hv,
  permutations = 999
)
print(adonis_hv)


disp_hv <- betadisper(
  bc_dist_healthy,
  group = permanova_hv$HypervariableRegion
)

# Permutation test for dispersion differences
perm_disp_hv <- permutest(
  disp_hv,
  permutations = 999
)

print(perm_disp_hv)