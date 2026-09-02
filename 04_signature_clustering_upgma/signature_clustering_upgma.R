library(tidyverse)
library(ggdendro)
library(ggplot2)

# CONFIGURATION

input_file <- "ssGSEA_scores_TC_5_public_mouse__with_modules_signatures_scaled.csv"
output_dir <- "signature_clustering_analysis"
dir.create(output_dir, showWarnings = FALSE)

# All 30 module columns (as named in your file)
module_cols <- paste0("Module-Module-", 1:30)

# The 6 signatures to cluster (rename labels as you like for display)
signature_cols <- c("Module-CBC", "Module-EpiHR", "Module-Fetal-Mustata",
                    "Module-OncoFetal", "Module-RSC", "Module-coreHRC")

# Display labels matching your figure style (edit freely)
signature_labels <- c(
  "Module-CBC"           = "White CBC (WNT)",
  "Module-EpiHR"         = "EpiHR (high relapse)",
  "Module-Fetal-Mustata" = "Mustata fetal",
  "Module-OncoFetal"     = "Han oncofetal",
  "Module-RSC"           = "White RSC",
  "Module-coreHRC"       = "coreHRC (EMP1)"
)

cor_method   <- "pearson"   # correlation method for building the 30-program profile
clust_method <- "average"    # "average" = UPGMA

# Distance threshold lines to draw on the dendrogram.
# Set to NULL to skip. Set to "auto" to pick the 2 largest merge-height gaps.
# Or supply manual numeric values once you've reviewed the merge table, e.g.
#   distance_cutoffs <- c(0.14, 0.38)
distance_cutoffs <- "auto"

# TIFF output settings
tiff_dpi         <- 600
tiff_compression <- "lzw"
tiff_type        <- "cairo"

open_tiff <- function(filepath, width, height, dpi, compression, type) {
  tryCatch({
    tiff(filepath, width = width, height = height, units = "in",
         res = dpi, compression = compression, type = type)
  }, error = function(e) {
    message("  [cairo unavailable, falling back to default type] ", basename(filepath))
    tiff(filepath, width = width, height = height, units = "in",
         res = dpi, compression = compression)
  })
}

# STEP 1: Load Data

cat(paste(rep("=", 80), collapse = ""), "\n")
print("STEP 1: Loading ssGSEA score matrix...")
cat(paste(rep("=", 80), collapse = ""), "\n")

raw <- read.csv(input_file, row.names = 1, check.names = FALSE)
cat(paste("✓ Matrix dimensions:", nrow(raw), "spots ×", ncol(raw), "features\n"))

stopifnot(all(signature_cols %in% colnames(raw)))
stopifnot(all(module_cols %in% colnames(raw)))

# STEP 2: Build the 30-program correlation profile for each signature

cat("\nSTEP 2: Computing each signature's correlation profile against 30 modules...\n")

mat_clean <- raw %>% drop_na(all_of(c(signature_cols, module_cols)))
cat(paste("✓ Spots retained after NA removal:", nrow(mat_clean), "/", nrow(raw), "\n"))

# For each signature, compute its correlation (rho) against each of the 30 modules
# Result: a (6 signatures) x (30 modules) matrix of rho values
profile_matrix <- sapply(signature_cols, function(sig) {
  sapply(module_cols, function(mod) {
    cor(mat_clean[[sig]], mat_clean[[mod]], method = cor_method, use = "complete.obs")
  })
})

profile_matrix <- t(profile_matrix)  # rows = signatures, cols = modules
rownames(profile_matrix) <- signature_labels[signature_cols]

cat(paste("✓ Correlation profile matrix:", nrow(profile_matrix), "signatures ×",
          ncol(profile_matrix), "modules\n"))

write.csv(profile_matrix, file.path(output_dir, "signature_30module_correlation_profiles.csv"))

# STEP 3: Distance between signatures' correlation profiles (1 - Pearson r)

cat("\nSTEP 3: Computing inter-signature distance (1 - correlation) and UPGMA clustering...\n")

# Distance = 1 - Pearson correlation BETWEEN each pair of signatures' 30-module
# profiles. This compares profile SHAPE (do two signatures relate to the same
# modules in the same pattern), ignoring absolute differences in magnitude.
# Identical profile shape -> distance 0. Opposite profile shape -> distance 2.
profile_cor <- cor(t(profile_matrix), method = "pearson")
dist_matrix <- as.dist(1 - profile_cor)

cat("\nPairwise correlation between signature profiles:\n")
print(round(profile_cor, 4))

cat("\nPairwise distance (1 - correlation):\n")
dist_mat_full <- round(as.matrix(dist_matrix), 4)
print(dist_mat_full)

# Save full pairwise distance matrix as CSV
write.csv(dist_mat_full, file.path(output_dir, "signature_distance_matrix_correlation.csv"))
cat("\n✓ Full pairwise distance matrix saved: signature_distance_matrix_correlation.csv\n")

# UPGMA = hierarchical clustering with "average" linkage
hc <- hclust(dist_matrix, method = clust_method)

# STEP 4: Decode and save the exact merge order / merge heights

cat("\nSTEP 4: Decoding exact merge order and heights...\n")

# hc$merge: each row = one merge step. Negative numbers = original signatures
# (by index), positive numbers = the cluster formed at that earlier step.
# hc$height: the distance at which that merge happened.
n <- length(hc$labels)
cluster_members <- vector("list", n - 1)

describe_member <- function(idx) {
  if (idx < 0) {
    hc$labels[-idx]
  } else {
    paste0("{", paste(cluster_members[[idx]], collapse = ", "), "}")
  }
}

merge_table <- tibble(
  step = seq_len(nrow(hc$merge)),
  height = round(hc$height, 4),
  member_1 = sapply(hc$merge[, 1], describe_member),
  member_2 = sapply(hc$merge[, 2], describe_member)
) %>%
  mutate(merged_into = paste0("{", member_1, " + ", member_2, "}"))

for (i in seq_len(nrow(hc$merge))) {
  m1 <- if (hc$merge[i, 1] < 0) hc$labels[-hc$merge[i, 1]] else cluster_members[[hc$merge[i, 1]]]
  m2 <- if (hc$merge[i, 2] < 0) hc$labels[-hc$merge[i, 2]] else cluster_members[[hc$merge[i, 2]]]
  cluster_members[[i]] <- c(m1, m2)
}

cat("\nFull merge sequence (in the order UPGMA joined signatures):\n")
print(merge_table %>% select(step, height, merged_into), width = Inf)

write.csv(merge_table %>% select(step, height, merged_into),
          file.path(output_dir, "signature_merge_order.csv"), row.names = FALSE)
cat("\n✓ Exact merge order/heights saved: signature_merge_order.csv\n")
cat("  (Use this table — not visual inspection of the dendrogram — to state\n")
cat("   exactly which signatures join at which distance.)\n")

# ---- Resolve distance_cutoffs: auto-detect from tree if requested ----
if (!is.null(distance_cutoffs) && identical(distance_cutoffs, "auto")) {

  merge_heights <- sort(hc$height)
  gaps <- diff(merge_heights)

  n_cuts <- min(2, length(gaps))
  top_gap_idx <- order(gaps, decreasing = TRUE)[1:n_cuts]
  distance_cutoffs <- sort(merge_heights[top_gap_idx] + gaps[top_gap_idx] / 2)

  cat("\n✓ Auto-selected cutoffs based on largest gaps in merge heights:\n")
  print(round(distance_cutoffs, 4))
}

# STEP 5: Plot dendrogram

cat("\nSTEP 5: Plotting dendrogram...\n")

dend_data <- dendro_data(as.dendrogram(hc))

y_range <- max(dend_data$segments$y, dend_data$segments$yend)
label_offset <- y_range * 0.035

p <- ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               color = "black", linewidth = 0.5) +
  geom_point(data = dend_data$labels,
             aes(x = x, y = 0),
             size = 5, color = "#8e6bb5") +
  geom_text(data = dend_data$labels,
            aes(x = x, y = -y_range * 0.05, label = label),
            angle = 30, hjust = 1, vjust = 1, size = 3.2,
            lineheight = 0.85) +
  labs(title = "Unsupervised clustering of six signatures",
       subtitle = "UPGMA on 30-program correlation profiles  |  distance = 1 - correlation",
       y = "Distance", x = NULL) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.18, 0.08))) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.margin = margin(t = 10, r = 20, b = 70, l = 10))

if (!is.null(distance_cutoffs)) {
  cutoff_colors <- c("#9b59b6", "#3498db")
  for (i in seq_along(distance_cutoffs)) {
    cutoff_label <- paste0("d = ", round(distance_cutoffs[i], 2))
    p <- p + geom_hline(yintercept = distance_cutoffs[i],
                        linetype = "dashed",
                        color = cutoff_colors[((i - 1) %% length(cutoff_colors)) + 1],
                        linewidth = 0.6) +
      annotate("label", x = -Inf, y = distance_cutoffs[i] + label_offset,
               label = cutoff_label,
               color = cutoff_colors[((i - 1) %% length(cutoff_colors)) + 1],
               hjust = 0, size = 4, fontface = "bold",
               label.size = 0, fill = "white")
  }
}

dendrogram_file <- file.path(output_dir, "signature_UPGMA_dendrogram_correlation.tiff")
open_tiff(dendrogram_file, width = 9, height = 7,
         dpi = tiff_dpi, compression = tiff_compression, type = tiff_type)
print(p)
dev.off()
cat(paste("✓ Dendrogram saved:", basename(dendrogram_file), "\n"))

pdf(file.path(output_dir, "signature_UPGMA_dendrogram_correlation.pdf"), width = 9, height = 7)
print(p)
dev.off()

# STEP 6: Cluster membership at each cutoff

cat("\nSTEP 6: Cluster membership at specified distance cutoffs...\n")

if (!is.null(distance_cutoffs)) {
  for (cutoff in distance_cutoffs) {
    clusters <- cutree(hc, h = cutoff)
    cat(paste0("\nAt d = ", round(cutoff, 2), ":\n"))
    for (cl in sort(unique(clusters))) {
      members <- names(clusters)[clusters == cl]
      cat(paste0("  Cluster ", cl, ": ", paste(members, collapse = ", "), "\n"))
    }
  }
}

# DONE

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
print("✓ ANALYSIS COMPLETE!")
cat(paste(rep("=", 80), collapse = ""), "\n")
print(paste("Output directory:", output_dir))
print("Files generated:")
print("  1. signature_30module_correlation_profiles.csv  (6 signatures x 30 modules, rho)")
print("  2. signature_distance_matrix_correlation.csv    (6x6 pairwise distance matrix)")
print("  3. signature_merge_order.csv                    (exact UPGMA merge steps/heights)")
print("  4. signature_UPGMA_dendrogram_correlation.tiff   (600 DPI)")
print("  5. signature_UPGMA_dendrogram_correlation.pdf")
cat(paste(rep("=", 80), collapse = ""), "\n")
