## =============================================================
## 01_comet_visium_alignment.R
## Cross-platform validation of spatial ecotypes (Fig. 1e, Extended Data Fig. 1e)
## Stage 1: Visium SE vs COMET pSE -- centroids, similarity matrices, and
## optimal one-to-one ecotype matching (Hungarian algorithm).
##
## Output CC_best_match_assignment.csv feeds into
## 02_comet_visium_scatter.R (Extended Data Fig. 1e).
##
## Inputs (place in ./data/):
##   k10/ISCHIA_k10_only_spatial_Assay.rds  - Visium ecotype assignments
##   COMET_pseudospots_seurat.rds           - COMET pseudospot ecotype assignments
##   RCTD_norm_weights.csv                  - Visium RCTD cell-type proportions
## =============================================================

library(Seurat)
library(dplyr)
library(clue)

dir.create("results", showWarnings = FALSE)

## ---------------------------------------------------------------
## 1. LOAD OBJECTS (after ISCHIA has been run on both)
## ---------------------------------------------------------------

visium <- readRDS("data/CRC_spatial_atlas.rds")
seu    <- readRDS("data/COMET_pseudospots_seurat.rds")

cat("Visium CC distribution:\n"); print(table(visium$CompositionCluster_CC))
cat("\nCOMET CC distribution:\n");  print(table(seu$CompositionCluster_CC))

## ---------------------------------------------------------------
## 2. LOAD RCTD PROPORTION MATRICES
## ---------------------------------------------------------------

a <- read.csv("data/RCTD_norm_weights.csv", row.names = 1)

deconv.mat <- as.data.frame(t(as.matrix(
  GetAssayData(seu, assay = "predictions", slot = "data")
)))
stopifnot(all(rownames(deconv.mat) == colnames(seu)))
cat("\nCOMET deconv.mat dimensions:", dim(deconv.mat), "\n")

## ---------------------------------------------------------------
## 3. HARMONIZE SAMPLE NAMES & SUBSET VISIUM TO MATCHED SAMPLES
## ---------------------------------------------------------------

visium$orig.ident_clean <- gsub("_", "-", visium$orig.ident)
stopifnot(length(unique(visium$orig.ident)) == length(unique(visium$orig.ident_clean)))

comet_samples   <- unique(seu$sample)
matched_samples <- intersect(unique(visium$orig.ident_clean), comet_samples)
cat("\nMatched samples:\n"); print(matched_samples)
stopifnot(length(matched_samples) > 0)

visium_sub <- subset(visium, subset = orig.ident_clean %in% matched_samples)
cat("\nVisium subset CC distribution:\n"); print(table(visium_sub$CompositionCluster_CC))

a_sub <- a[colnames(visium_sub), ]
visium_cc_sub <- visium_sub$CompositionCluster_CC[rownames(a_sub)]
comet_cc      <- seu$CompositionCluster_CC[rownames(deconv.mat)]

stopifnot(!any(is.na(visium_cc_sub)))
stopifnot(!any(is.na(comet_cc)))

## ---------------------------------------------------------------
## 4. HARMONIZE COLUMN NAMES & FIND COMMON CELL TYPES
## ---------------------------------------------------------------

colnames(a_sub)      <- gsub("-", "_", colnames(a_sub))
colnames(deconv.mat) <- gsub("-", "_", colnames(deconv.mat))

common_cols <- intersect(colnames(a_sub), colnames(deconv.mat))
cat("\nCommon cell-type columns:\n"); print(common_cols)

visium_comp <- a_sub[, common_cols]
comet_comp  <- deconv.mat[, common_cols]

## ---------------------------------------------------------------
## 5. RENORMALIZE BOTH PLATFORMS TO COMMON CELL TYPES (sum to 1)
##    COMET pseudospots with zero total signal across common_cols
##    are excluded before renormalizing (per Methods).
## ---------------------------------------------------------------

visium_comp_norm <- visium_comp / rowSums(visium_comp)

keep            <- rowSums(comet_comp) > 0
comet_comp_norm <- comet_comp[keep, ] / rowSums(comet_comp[keep, ])
comet_cc_norm   <- comet_cc[keep]

cat("\nRow sums after renorm -- Visium:",
    round(mean(rowSums(visium_comp_norm)), 3),
    " COMET:", round(mean(rowSums(comet_comp_norm)), 3), "\n")
cat("COMET pseudospots excluded (zero signal):", sum(!keep), "of", length(keep), "\n")

## ---------------------------------------------------------------
## 6. COMPUTE PER-CC CENTROIDS
## ---------------------------------------------------------------

visium_centroid <- aggregate(visium_comp_norm, by = list(CC = visium_cc_sub), FUN = mean)
comet_centroid  <- aggregate(comet_comp_norm,  by = list(CC = comet_cc_norm), FUN = mean)

rownames(visium_centroid) <- visium_centroid$CC; visium_centroid$CC <- NULL
rownames(comet_centroid)  <- comet_centroid$CC;  comet_centroid$CC  <- NULL

visium_centroid <- visium_centroid[order(rownames(visium_centroid)), ]
comet_centroid  <- comet_centroid[order(rownames(comet_centroid)),  ]

cat("\nVisium centroids:\n"); print(round(visium_centroid, 3))
cat("\nCOMET centroids:\n");  print(round(comet_centroid, 3))

saveRDS(visium_centroid, "results/visium_centroid.rds")
saveRDS(comet_centroid,  "results/comet_centroid.rds")

## ---------------------------------------------------------------
## 7. SIMILARITY MATRIX (Pearson)
## ---------------------------------------------------------------

sim_mat_pearson <- cor(t(visium_centroid), t(comet_centroid), method = "pearson")
cat("\nPearson correlation matrix:\n"); print(round(sim_mat_pearson, 2))

## ---------------------------------------------------------------
## 8. OPTIMAL 1:1 ASSIGNMENT (Hungarian algorithm)
## ---------------------------------------------------------------

sim_mat_shifted <- sim_mat_pearson + 1  # shift to non-negative for solve_LSAP
assignment <- solve_LSAP(sim_mat_shifted, maximum = TRUE)

best_match_df <- data.frame(
  Visium_CC        = rownames(sim_mat_pearson),
  Best_COMET_match = colnames(sim_mat_pearson)[assignment],
  Pearson_r        = sim_mat_pearson[cbind(seq_along(assignment), assignment)]
)
cat("\nOptimal 1:1 CC matches:\n"); print(best_match_df)
write.csv(best_match_df, "results/CC_best_match_assignment.csv", row.names = FALSE)

cat("\nDone. Output: results/CC_best_match_assignment.csv\n")
cat("Feed this into 02_comet_visium_scatter.R for Extended Data Fig. 1e.\n")
