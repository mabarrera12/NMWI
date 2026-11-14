# Load packages -----------------------------------------------------------
library(vegan)
library(ade4)
library(ggplot2)
set.seed(1234)

# Load input data ---------------------------------------------------------
training_set <- read.delim("Data/Output/training_set.tsv", sep = "\t", header = TRUE, check.names = FALSE)

# Select only taxonomy abundance columns (columns starting with "d__")
taxa_mat <- training_set[, grepl("^d__", colnames(training_set))]

# Recode group labels for case/control comparison
# control = TRUE, case = FALSE
group <- c(case = FALSE, control = TRUE)[training_set$Category]

# Build data frame for PERMANOVA
permanova_df <- data.frame(
  Group = group,
  BioProject = training_set$`BioProject Number`,
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
  Sample_ID = training_set$`Sample ID`,
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
  Disease = training_set$Disease
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

fig_case_control


# Plot PCoA colored by Disease -----------------------------------------
fig_Disease <- ggplot(pc_plot_data, aes(x, y, color = factor(Disease))) +
  geom_point(size = 2) +
  labs(
    x = paste0("PCoA1 (", variance1, "%)"),
    y = paste0("PCoA2 (", variance2, "%)")
  ) +
  theme_bw() +
  theme(legend.title = element_blank())

print(fig_Disease)
