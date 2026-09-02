# core_to_edge_region_gradient.R
#
# Core-to-edge module gradient analysis (Fig. 2b, 2d).
# Spots are binned into four tissue regions by tumor-epithelium fraction
# (Core/Inner margin/Outer margin/Edge). For each of the 30 spatial modules,
# ssGSEA scores are tested against ordinal region rank by Spearman
# correlation (BH-corrected across 30 modules).
#
# Panel A: mean module activity per region (heatmap) + per-module Spearman
#          rho, coloured by peak region.
# Panel B: per-tumor core-to-edge gradient, reproducibility across tumors.
#
# Inputs (place in ./data/):
#   Tumor_region_meta.csv   - per-spot metadata: barcode, orig.ident, Tumor_region
#   ssgsea_scores_scaled.csv - per-spot ssGSEA scores, Module_1..Module_30

library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(tibble)
library(gridExtra)
library(grid)
library(conflicted)

conflict_prefer("select",    "dplyr")
conflict_prefer("filter",    "dplyr")
conflict_prefer("mutate",    "dplyr")
conflict_prefer("summarise", "dplyr")

# ---- output folder ----
out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ---- color scheme ----
heatmap_colors <- colorRampPalette(c(
  "#2166ac",   # blue = low/negative
  "white",     # white = zero
  "#b2182b"    # red = high/positive
))(100)

# ---- module annotations ----
module_labels_raw <- c(
  M1="Absorptive progenitor epithelial",       M2="EMP1+ injury-associated epithelial",
  M3="HER2/fetal-like epithelial stress",      M4="Squamous-like epithelial differentiation",
  M5="IFN-y/checkpoint immune",                M6="CAF/stromal ECM-remodelling",
  M7="SPP1+ TAM/myeloid",                      M8="PROX1+ neuroendocrine-like",
  M9="RTK/MAPK-feedback epithelial adaptation",M10="Quiescent/dormant epithelial",
  M11="ISC-like stem-cell",                    M12="Secretory intestine/TFF",
  M13="Inflammatory epithelial stress",        M14="mTORC1/lipogenic metabolic",
  M15="Hypoxia/endoderm-like metabolic",       M16="Goblet/secretory epithelial",
  M17="Invasive epithelial/partial EMT-like",  M18="Notch/TGF-b inflammatory signalling",
  M19="Ribosomal/translation proliferative",   M20="EGFR-MYC proliferative epithelial",
  M21="Active CBC/transit-amplifying",         M22="Differentiated neuroendocrine-like",
  M23="IFN-response inflammatory",             M24="Hypoxia/glycolytic stress",
  M25="OXPHOS-associated metabolic",           M26="Absorptive intestine/lipid-catabolic",
  M27="Paneth/intestinal secretory",           M28="L1CAM+ oncofetal/injury-regenerative",
  M29="Cytoskeletal/junctional adaptation",    M30="Smooth muscle/pericyte vascular"
)
names(module_labels_raw) <- gsub("^M", "Module_", names(module_labels_raw))
module_labels <- module_labels_raw

# ---- load and merge ----
meta    <- read.csv("data/Tumor_region_meta.csv")
modules <- read.csv("data/ssgsea_scores_scaled.csv", row.names = 1)

colnames(modules) <- gsub("^Module\\.", "Module_", colnames(modules))
colnames(modules) <- gsub("^Module_Module\\.", "Module_", colnames(modules))
colnames(meta)[1] <- "barcode"
modules <- rownames_to_column(modules, var = "barcode")

module_cols  <- paste0("Module_", 1:30)
region_order <- c("Core", "Inner margin", "Outer margin", "Edge")

merged <- merge(
  meta[, c("barcode", "orig.ident", "Tumor_region")],
  modules[, c("barcode", module_cols)],
  by = "barcode"
)

merged$Tumor_region <- factor(merged$Tumor_region, levels = region_order)
merged              <- merged %>% filter(!is.na(Tumor_region))
merged$region_rank  <- c("Core" = 1, "Inner margin" = 2, "Outer margin" = 3, "Edge" = 4)[
  as.character(merged$Tumor_region)]

cat("Merged spots:", nrow(merged), "\n")
cat("Tumors:", length(unique(merged$orig.ident)), "\n")
print(table(merged$Tumor_region))

# ---- Spearman correlation vs. region rank (per module) ----
spearman_results <- lapply(module_cols, function(mod) {
  test <- cor.test(merged[[mod]], merged$region_rank, method = "spearman", exact = FALSE)
  data.frame(module = mod, rho = test$estimate, pvalue = test$p.value)
})

spearman_df             <- bind_rows(spearman_results)
spearman_df$q_value     <- p.adjust(spearman_df$pvalue, method = "BH")
spearman_df$significant <- spearman_df$q_value < 0.05

# Data-driven module order: most core-enriched -> most edge-enriched
module_rho_order <- spearman_df %>% arrange(rho) %>% pull(module)

cat("\nTop 5 core-enriched modules:\n")
print(module_labels[head(module_rho_order, 5)])
cat("\nTop 5 edge-enriched modules:\n")
print(module_labels[tail(module_rho_order, 5)])

spearman_table <- spearman_df %>%
  mutate(Module_Name = module_labels[module]) %>%
  select(Module = module, Module_Name, rho, pvalue, q_value, significant) %>%
  arrange(rho)

write.csv(spearman_table, file.path(out_dir, "Spearman_correlation_table.csv"), row.names = FALSE)
cat("\nSaved Spearman correlation table to:", file.path(out_dir, "Spearman_correlation_table.csv"), "\n")

# ---- Panel A left: mean z-score heatmap ----
mean_z <- merged %>%
  group_by(Tumor_region) %>%
  summarise(across(all_of(module_cols), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("Tumor_region") %>%
  t() %>%
  as.data.frame()

mean_z_ordered <- mean_z[module_rho_order, ]
stopifnot(rownames(mean_z_ordered) == module_rho_order)
rownames(mean_z_ordered) <- module_labels[rownames(mean_z_ordered)]

panel_A_left <- pheatmap(
  mean_z_ordered[, region_order],
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = heatmap_colors, breaks = seq(-0.5, 0.5, length.out = 101),
  border_color = "white",
  main = "Module activity across the core-to-edge axis\n(mean z-score per bin)",
  angle_col = 45, fontsize_row = 8, fontsize_col = 10,
  legend_breaks = c(-0.4, -0.2, 0, 0.2, 0.4),
  labels_col = c("Core", "Inner", "Outer", "Edge"),
  cellwidth = 40, cellheight = 12, silent = TRUE
)

# ---- Panel A right: Spearman bar chart ----
peak_region <- merged %>%
  group_by(Tumor_region) %>%
  summarise(across(all_of(module_cols), \(x) mean(x, na.rm = TRUE))) %>%
  pivot_longer(-Tumor_region, names_to = "module", values_to = "mean_score") %>%
  group_by(module) %>%
  slice_max(mean_score, n = 1) %>%
  dplyr::select(module, peak_region = Tumor_region)

spearman_df <- left_join(spearman_df, peak_region, by = "module")

region_colors <- c(
  "Core"         = "#b2182b",
  "Inner margin" = "#d6604d",
  "Outer margin" = "#f4a582",
  "Edge"         = "#fddbc7"
)

spearman_df$label <- factor(module_labels[spearman_df$module], levels = rownames(mean_z_ordered))
max_rho <- ceiling(max(abs(spearman_df$rho), na.rm = TRUE) * 10) / 10

panel_A_right <- ggplot(spearman_df, aes(x = rho, y = label, fill = peak_region)) +
  geom_col(width = 0.7) +
  geom_text(
    data = spearman_df %>% filter(significant),
    aes(x = ifelse(rho > 0, rho + 0.01, rho - 0.01), label = "*"),
    hjust = ifelse((spearman_df %>% filter(significant))$rho > 0, 0, 1),
    size = 4, color = "black"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = region_colors, na.value = "grey70", name = "Peak region") +
  scale_x_continuous(
    limits = c(-max_rho - 0.05, max_rho + 0.05),
    breaks = round(seq(-max_rho, max_rho, length.out = 5), 1)
  ) +
  labs(
    title = "Core-to-edge coordinate per module\n(Spearman rho; * q<0.05)",
    x = "<- core          edge-trend (rho)          interface ->", y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 7), legend.position = "bottom",
        plot.title = element_text(size = 9, hjust = 0.5))

# ---- Panel B: per-tumor gradient heatmap ----
valid_tumors <- merged %>%
  group_by(orig.ident) %>%
  summarise(n_regions = n_distinct(Tumor_region)) %>%
  filter(n_regions == 4) %>%
  pull(orig.ident)

cat("\nTumors with all 4 regions:", length(valid_tumors), "\n")

tumor_gradients <- merged %>%
  filter(orig.ident %in% valid_tumors) %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(module_cols), \(x) cor(x, region_rank, method = "spearman", use = "complete.obs"))) %>%
  column_to_rownames("orig.ident") %>%
  t() %>%
  as.data.frame()

tumor_gradients_ordered <- tumor_gradients[module_rho_order, ]
stopifnot(rownames(tumor_gradients_ordered) == module_rho_order)
rownames(tumor_gradients_ordered) <- module_labels[rownames(tumor_gradients_ordered)]

col_annotation <- data.frame(
  Position = sapply(colnames(tumor_gradients_ordered), function(s) {
    as.character(merged %>% filter(orig.ident == s) %>%
                   count(Tumor_region) %>% slice_max(n, n = 1) %>% pull(Tumor_region))
  })
)
rownames(col_annotation) <- colnames(tumor_gradients_ordered)
ann_colors <- list(Position = region_colors)

panel_B <- pheatmap(
  tumor_gradients_ordered,
  cluster_rows = FALSE, cluster_cols = TRUE,
  color = heatmap_colors, breaks = seq(-0.45, 0.45, length.out = 101),
  border_color = "white",
  annotation_col = col_annotation, annotation_colors = ann_colors,
  main = "Core-to-edge module gradient is reproducible across tumors\n(Each column = one tumor)",
  angle_col = 90, fontsize_row = 7, fontsize_col = 6,
  cellwidth = 12, cellheight = 12, silent = TRUE
)

# ---- save outputs ----
png(file.path(out_dir, "Panel_A_left.png"),  width = 5500, height = 3200, res = 300)
grid::grid.draw(panel_A_left$gtable)
dev.off()

png(file.path(out_dir, "Panel_A_right.png"), width = 4000, height = 3200, res = 300)
print(panel_A_right)
dev.off()

png(file.path(out_dir, "Panel_B.png"), width = 5500, height = 3200, res = 300)
grid::grid.draw(panel_B$gtable)
dev.off()

cat("\nAll panels + table saved to:", out_dir, "\n")
list.files(out_dir)
