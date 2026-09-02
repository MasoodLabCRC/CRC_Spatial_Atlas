## =============================================================
## 02_comet_visium_scatter.R
## Cross-platform validation of spatial ecotypes -- Extended Data Fig. 1e
## Visium SE vs COMET pSE scatter: matched CC pairs x cell types (10 x 6 = 60),
## CC matching via the optimal (Hungarian) assignment from
## 01_comet_visium_alignment.R. Both platforms renormalized to common
## cell types only, zero-signal COMET pseudospots excluded before renorm.
##
## Inputs (place in ./data/):
##   k10/ISCHIA_k10_only_spatial_Assay.rds  - Visium ecotype assignments
##   COMET_pseudospots_seurat.rds           - COMET pseudospot ecotype assignments
##   RCTD_norm_weights.csv                  - Visium RCTD cell-type proportions
## Also requires results/CC_best_match_assignment.csv from
## 01_comet_visium_alignment.R (run that script first).
## =============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

dir.create("results", showWarnings = FALSE)

## ---------------------------------------------------------------
## 1. LOAD OBJECTS
## ---------------------------------------------------------------

visium <- readRDS("data/k10/ISCHIA_k10_only_spatial_Assay.rds")
seu    <- readRDS("data/COMET_pseudospots_seurat.rds")

cat("Visium CC distribution:\n"); print(table(visium$CompositionCluster_CC))
cat("\nCOMET CC distribution:\n");  print(table(seu$CompositionCluster_CC))

## ---------------------------------------------------------------
## 2. LOAD PRE-COMPUTED OPTIMAL CC MATCH
## ---------------------------------------------------------------

optimal_match <- read.csv("results/CC_best_match_assignment.csv") %>%
  select(Visium_CC = Visium_CC, COMET_CC = Best_COMET_match, Pearson_r = Pearson_r)

cat("\nOptimal CC matches loaded:\n")
print(optimal_match)

## ---------------------------------------------------------------
## 3. LOAD RCTD PROPORTION MATRICES
## ---------------------------------------------------------------

a <- read.csv("data/RCTD_norm_weights.csv", row.names = 1)

deconv.mat <- as.data.frame(t(as.matrix(
  GetAssayData(seu, assay = "predictions", slot = "data")
)))
stopifnot(all(rownames(deconv.mat) == colnames(seu)))

## ---------------------------------------------------------------
## 4. HARMONIZE SAMPLE NAMES & SUBSET VISIUM TO MATCHED SAMPLES
## ---------------------------------------------------------------

visium$orig.ident_clean <- gsub("_", "-", visium$orig.ident)
comet_samples   <- unique(seu$sample)
matched_samples <- intersect(unique(visium$orig.ident_clean), comet_samples)
cat("\nMatched samples:\n"); print(matched_samples)

visium_sub    <- subset(visium, subset = orig.ident_clean %in% matched_samples)
a_sub         <- a[colnames(visium_sub), ]
visium_cc_sub <- visium_sub$CompositionCluster_CC[rownames(a_sub)]
comet_cc      <- seu$CompositionCluster_CC[rownames(deconv.mat)]

stopifnot(!any(is.na(visium_cc_sub)))
stopifnot(!any(is.na(comet_cc)))

## ---------------------------------------------------------------
## 5. HARMONIZE COLUMN NAMES & FIND COMMON CELL TYPES
## ---------------------------------------------------------------

colnames(a_sub)      <- gsub("-", "_", colnames(a_sub))
colnames(deconv.mat) <- gsub("-", "_", colnames(deconv.mat))

common_cols <- intersect(colnames(a_sub), colnames(deconv.mat))
cat("\nCommon cell-type columns:\n"); print(common_cols)

visium_comp <- a_sub[, common_cols]
comet_comp  <- deconv.mat[, common_cols]

## ---------------------------------------------------------------
## 6. RENORMALIZE BOTH TO COMMON COLUMNS (sum to 1)
## ---------------------------------------------------------------

visium_comp_norm <- visium_comp / rowSums(visium_comp)

keep              <- rowSums(comet_comp) > 0
comet_comp_norm   <- comet_comp[keep, ] / rowSums(comet_comp[keep, ])
comet_cc_norm     <- comet_cc[keep]

cat("\nRow sums after renorm — Visium:",
    round(mean(rowSums(visium_comp_norm)), 3),
    " COMET:", round(mean(rowSums(comet_comp_norm)), 3), "\n")

## ---------------------------------------------------------------
## 7. COMPUTE PER-CC CENTROIDS
## ---------------------------------------------------------------

visium_centroid <- aggregate(visium_comp_norm,
                             by = list(CC = visium_cc_sub), FUN = mean)
comet_centroid  <- aggregate(comet_comp_norm,
                             by = list(CC = comet_cc_norm), FUN = mean)

rownames(visium_centroid) <- visium_centroid$CC; visium_centroid$CC <- NULL
rownames(comet_centroid)  <- comet_centroid$CC;  comet_centroid$CC  <- NULL

visium_centroid <- visium_centroid[order(rownames(visium_centroid)), ]
comet_centroid  <- comet_centroid[order(rownames(comet_centroid)),  ]

cat("\nVisium centroids:\n"); print(round(visium_centroid, 3))
cat("\nCOMET centroids:\n");  print(round(comet_centroid,  3))

## ---------------------------------------------------------------
## 8. FLATTEN CENTROIDS TO LONG FORMAT
## ---------------------------------------------------------------

visium_long <- visium_centroid %>%
  tibble::rownames_to_column("Visium_CC") %>%
  pivot_longer(-Visium_CC, names_to = "CellType", values_to = "Visium_prop")

comet_long <- comet_centroid %>%
  tibble::rownames_to_column("COMET_CC") %>%
  pivot_longer(-COMET_CC, names_to = "CellType", values_to = "COMET_prop")

## ---------------------------------------------------------------
## 9. BUILD SCATTER DATA USING OPTIMAL MATCHED CC PAIRS
##    Visium CC1 -> COMET CC6 (optimal), not CC1 -> CC1
## ---------------------------------------------------------------

scatter_df <- optimal_match %>%
  left_join(visium_long, by = "Visium_CC") %>%
  left_join(comet_long,  by = c("COMET_CC", "CellType")) %>%
  mutate(Pair = paste0(Visium_CC, " \u2194 ", COMET_CC))

cat("\nScatter points:", nrow(scatter_df),
    "(expected:", nrow(optimal_match) * length(common_cols), ")\n")
write.csv(scatter_df, "results/SE_pSE_scatter_data.csv", row.names = FALSE)

## ---------------------------------------------------------------
## 10. PEARSON CORRELATION ACROSS ALL MATCHED POINTS
## ---------------------------------------------------------------

cor_test <- cor.test(scatter_df$Visium_prop, scatter_df$COMET_prop,
                     method = "pearson")
r_val   <- round(cor_test$estimate, 2)
p_val   <- cor_test$p.value
p_label <- ifelse(p_val < 1e-6,
                  "P < 1e\u22126",
                  paste0("P = ", formatC(p_val, format = "e", digits = 1)))

cat("\nOverall Pearson R:", r_val, "  P:", p_val, "\n")

## ---------------------------------------------------------------
## 11. SCATTER PLOT — single purple color
## ---------------------------------------------------------------

p <- ggplot(scatter_df, aes(x = Visium_prop, y = COMET_prop)) +
  geom_point(color = "#7B2D8B", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linewidth = 0.8) +
  annotate("text",
           x     = max(scatter_df$Visium_prop, na.rm = TRUE) * 0.98,
           y     = min(scatter_df$COMET_prop,  na.rm = TRUE) +
                   diff(range(scatter_df$COMET_prop, na.rm = TRUE)) * 0.08,
           label = paste0("R = ", r_val, "\n", p_label),
           hjust = 1, size = 4,
           color = "#7B2D8B", fontface = "italic") +
  scale_x_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.01, 0.05))) +
  scale_y_continuous(limits = c(0, NA),
                     expand = expansion(mult = c(0.01, 0.05))) +
  labs(
    title    = paste0("COMET vs Visium ecotype composition\n",
                      "(n = ", length(matched_samples), " samples)"),
    subtitle = paste0(nrow(scatter_df), " points  |  ",
                      length(common_cols), " cell types  \u00d7  ",
                      nrow(optimal_match), " optimal CC pairs"),
    x        = "Visium (RNA)",
    y        = "COMET (protein)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, size = 11),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    axis.title    = element_text(face = "bold"),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave("results/SE_pSE_scatter.pdf", plot = p, width = 5, height = 5)
ggsave("results/SE_pSE_scatter.png", plot = p, width = 5, height = 5, dpi = 300)
print(p)

cat("\nOutputs saved:\n")
cat(" - results/SE_pSE_scatter.pdf / .png\n")
cat(" - results/SE_pSE_scatter_data.csv\n")
