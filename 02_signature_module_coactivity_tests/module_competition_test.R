#!/usr/bin/env Rscript
# ============================================================================
# visium_module_competition_test.R 
#
# ============================================================================

suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(data.table) })
set.seed(1)

visium_obj <- readRDS("CRC_spatial_atlas.rds")

ssgsea <- readRDS("ssGSEA_scores_TC_5_public_mouse__with_modules_signatures.rds")
ssgsea <- as.data.frame(as.matrix(t(ssgsea))) 
# 0. PARAMETERS

TARGET     <- "M2"
CONVERGING <- c("coreHRC", "EpiHR", "Fetal_Mustata", "OncoFetal", "RSC")
CONTROL    <- c("CBC")
STEM       <- c("M11", "M21")

OUT_DIR <- "VisiumDiscovery"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# 1. LOAD Visium OBJECT + ssGSEA SCORES

stopifnot(exists("visium_obj"), exists("ssgsea"))

module_cols <- paste0("Module_", 1:30)
if (!all(module_cols %in% colnames(ssgsea))) module_cols <- paste0("Module-", 1:30)
sig_cols_raw <- c(CONVERGING, CONTROL)
stopifnot(all(module_cols %in% colnames(ssgsea)))
stopifnot(all(sig_cols_raw %in% colnames(ssgsea)))

md <- visium_obj@meta.data
common <- intersect(rownames(md), rownames(ssgsea))
cat(sprintf("Matched spots: %d of %d (object), %d (ssGSEA)\n", length(common), nrow(md), nrow(ssgsea)))

md  <- md[common, ]
MOD <- as.matrix(ssgsea[common, module_cols]); colnames(MOD) <- paste0("M", 1:30)
SIG <- as.matrix(ssgsea[common, sig_cols_raw])

# 2. DEPTH + BATCH RESIDUALIZATION -> COA (signature x module correlation matrix)

nC_col <- grep("^nCount_",  colnames(md), value=TRUE)[1]
nF_col <- grep("^nFeature_",colnames(md), value=TRUE)[1]
if (is.na(nC_col) || is.na(nF_col)) stop("Could not find nCount_/nFeature_ columns.")
nC <- md[[nC_col]]; nF <- md[[nF_col]]

force_col <- Sys.getenv("BATCH_COL","")
if (nzchar(force_col) && force_col %in% names(md)) { pick_samp <- force_col } else {
  pref <- c("patient","donor","subject","sample","orig\\.ident","slide","section","slice","tissue","run")
  cand <- unique(unlist(lapply(pref, function(p) grep(p, names(md), value=TRUE, ignore.case=TRUE))))
  cand <- c(cand, setdiff(grep("fov|image", names(md), value=TRUE, ignore.case=TRUE), cand))
  pick_samp <- NULL
  for (cc in cand){ nl<-nlevels(factor(md[[cc]])); if(nl>=2 && nl<=100){pick_samp<-cc; break} }
}
samp <- if (!is.null(pick_samp)) factor(md[[pick_samp]]) else factor(rep(1,nrow(md)))
cat(sprintf("Batch column: %s (%d levels)\n", ifelse(is.null(pick_samp),"<none>",pick_samp), nlevels(samp)))

scale01 <- function(v){ v<-as.numeric(v); s<-sd(v); if(!is.finite(s)||s==0) v-mean(v) else (v-mean(v))/s }
depth_cov <- cbind(scale01(rank(log1p(as.numeric(nC)))), scale01(rank(log1p(as.numeric(nF)))))
Z    <- if (nlevels(samp)>=2) Matrix::sparse.model.matrix(~ samp + depth_cov) else Matrix::sparse.model.matrix(~ depth_cov)
XtX  <- as.matrix(Matrix::crossprod(Z)); ridge0 <- 1e-6*mean(diag(XtX))
solve_reg <- function(A,lam){ tryCatch(solve(A+diag(lam,ncol(A))), error=function(e) NULL) }
XtXi <- solve_reg(XtX,ridge0)
if(is.null(XtXi)) XtXi<-solve_reg(XtX,1e-3*mean(diag(XtX)))
if(is.null(XtXi)) XtXi<-solve_reg(XtX,1e-1*mean(diag(XtX)))
if(is.null(XtXi)){ ev<-eigen(XtX+diag(ridge0,ncol(XtX)),symmetric=TRUE)
  di<-ifelse(ev$values>max(ev$values)*1e-8,1/ev$values,0); XtXi<-ev$vectors%*%(di*t(ev$vectors)) }
residZ <- function(v){ r<-rank(v); b<-XtXi%*%as.numeric(Matrix::crossprod(Z,r)); as.numeric(r-as.numeric(Z%*%b)) }

cat("Residualizing module + signature scores...\n")
MODp <- apply(MOD, 2, residZ)
SIGp <- apply(SIG, 2, residZ)

prho <- function(a,b){ s<-sqrt(sum(a^2)*sum(b^2)); if(s==0) NA_real_ else sum(a*b)/s }
sig_lab <- colnames(SIGp); mod_lab <- colnames(MODp)
COA <- outer(sig_lab, mod_lab, Vectorize(function(s,m) prho(SIGp[,s], MODp[,m])))
dimnames(COA) <- list(sig_lab, mod_lab)

write.csv(as.data.frame(COA), file.path(OUT_DIR, "visium_coactivity_matrix.csv"))

# 3. MODULE-COMPETITION TEST 

NMOD <- ncol(COA)

summ <- function(s){
  r  <- COA[s, ]
  ord<- order(r, decreasing=TRUE)
  am <- colnames(COA)[ord[1]]
  m2rank <- match(TARGET, colnames(COA)[ord])
  runner <- colnames(COA)[ord[ if(am==TARGET) 2 else 1]]
  data.table::data.table(
    signature = s, argmax = am, M2_rank = m2rank,
    rho_M2 = round(r[[TARGET]],3),
    rho_top = round(r[[ord[1]]],3),
    gap_M2_vs_runnerup = round(r[[TARGET]] - r[[ setdiff(ord[1:2], match(TARGET,colnames(COA)))[1] ]],3))
}
tab <- data.table::rbindlist(lapply(rownames(COA), summ))
cat("\n== per-signature ranking ==\n"); print(tab, row.names=FALSE)
write.csv(tab, file.path(OUT_DIR, "visium_module_competition_ranking.csv"), row.names = FALSE)

n  <- length(CONVERGING)
k  <- sum(tab[signature %in% CONVERGING, argmax] == TARGET)
p0 <- 1/NMOD

# Test A -- pre-specified target
pA <- pbinom(k-1, n, p0, lower.tail=FALSE)
cat(sprintf("\n[A] argmax on pre-specified %s: %d/%d signatures | Binom(n,1/%d) P = %.2e\n",
            TARGET, k, n, NMOD, pA))

# Test B -- module-agnostic, conservative union bound
pB <- min(1, NMOD * pbinom(k-1, n, p0, lower.tail=FALSE))
cat(sprintf("[B] any-module concordance (conservative upper bound) P <= %.2e\n", pB))

# Test C -- Fisher-combined over full M2 ranks
pC_each <- tab[signature %in% CONVERGING, M2_rank] / NMOD
X  <- -2*sum(log(pC_each)); pC <- pchisq(X, df=2*n, lower.tail=FALSE)
cat(sprintf("[C] Fisher over M2-ranks (%s): X2(%d)=%.1f  P = %.2e\n",
            paste(round(pC_each*NMOD), collapse=","), 2*n, X, pC))

cat("\n== control ==\n"); print(tab[signature %in% CONTROL], row.names=FALSE)

summary_stats <- data.frame(
  test = c("A_prespecified_argmax", "B_any_module_upper_bound", "C_fisher_combined_rank"),
  statistic = c(k, NA, round(X,2)),
  p_value = c(pA, pB, pC)
)
write.csv(summary_stats, file.path(OUT_DIR, "visium_module_competition_pvalues.csv"), row.names = FALSE)

# 4. FIGURES

library(ggplot2)
blue_red <- c("#2166ac", "#b2182b")

tab_plot <- tab
tab_plot$group <- ifelse(tab_plot$signature %in% CONTROL, "CBC (control)", "Converging")

# ---- FIG A: M2 rank per signature (lollipop, matches earlier argmax figure style) ----
tab_plot$signature <- factor(tab_plot$signature, levels = tab_plot$signature[order(-tab_plot$M2_rank)])

pA_fig <- ggplot(tab_plot, aes(x = M2_rank, y = signature, color = group)) +
  geom_segment(aes(x = 0, xend = M2_rank, y = signature, yend = signature), linewidth = 0.8) +
  geom_point(size = 4) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = paste0("argmax: ", argmax)), hjust = -0.15, size = 3.2, color = "black") +
  scale_color_manual(values = c("Converging" = blue_red[2], "CBC (control)" = blue_red[1]), name = NULL) +
  scale_x_reverse(limits = c(max(tab_plot$M2_rank) + 5, 0.5)) +
  labs(title = "M2's rank among 30 modules per signature",
       subtitle = sprintf("Test A: %d/%d converging hit M2 exactly (argmax), P = %.2e", k, n, pA),
       x = "M2 rank (1 = best)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "FIG_A_M2_rank_per_signature.pdf"), pA_fig, width = 8, height = 5, device = cairo_pdf)

# ---- FIG B: rho_M2 vs rho_top (gap visualization) ----
tab_plot2 <- tab
tab_plot2$group <- ifelse(tab_plot2$signature %in% CONTROL, "CBC (control)", "Converging")
tab_plot2$signature <- factor(tab_plot2$signature, levels = tab_plot2$signature[order(tab_plot2$rho_M2)])

pB_fig <- ggplot(tab_plot2) +
  geom_segment(aes(x = rho_M2, xend = rho_top, y = signature, yend = signature), color = "grey60", linewidth = 0.6) +
  geom_point(aes(x = rho_top, y = signature), color = "grey40", size = 3, shape = 1) +
  geom_point(aes(x = rho_M2, y = signature, color = group), size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Converging" = blue_red[2], "CBC (control)" = blue_red[1]), name = "rho_M2") +
  labs(title = "rho_M2 (filled) vs top-ranked module's rho (open circle)",
       subtitle = "Small gap = M2 IS the top module; large gap = another module wins",
       x = "Spearman rho", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(file.path(OUT_DIR, "FIG_B_rho_M2_vs_top.pdf"), pB_fig, width = 8, height = 5, device = cairo_pdf)

# ---- FIG C: p-value summary bar (log scale) ----
summary_stats$test_label <- c("A: pre-specified argmax", "B: any-module (conservative)", "C: Fisher-combined rank")
summary_stats$test_label <- factor(summary_stats$test_label, levels = rev(summary_stats$test_label))

pC_fig <- ggplot(summary_stats, aes(x = -log10(p_value), y = test_label)) +
  geom_col(fill = blue_red[2], width = 0.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dotted", color = "grey40") +
  geom_text(aes(label = sprintf("P = %.2e", p_value)), hjust = -0.1, size = 3.5) +
  labs(title = "Module-competition test: significance across three formulations",
       subtitle = "Dotted line = P = 0.05 threshold",
       x = "-log10(P)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) +
  xlim(0, max(-log10(summary_stats$p_value)) * 1.3)

ggsave(file.path(OUT_DIR, "FIG_C_pvalue_summary.pdf"), pC_fig, width = 8, height = 4, device = cairo_pdf)

# ---- FIG D: full co-activity heatmap, M2 column boxed ----
coa_long <- as.data.frame(COA) %>%
  tibble::rownames_to_column("signature") %>%
  tidyr::pivot_longer(-signature, names_to = "module", values_to = "rho")
coa_long$module <- factor(coa_long$module, levels = paste0("M", 1:30))
coa_long$signature <- factor(coa_long$signature, levels = rev(c(CONVERGING, CONTROL)))

pD_fig <- ggplot(coa_long, aes(x = module, y = signature, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_tile(data = subset(coa_long, module == TARGET), fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_gradient2(low = blue_red[1], mid = "white", high = blue_red[2], midpoint = 0, name = "rho") +
  labs(title = "Signature x module co-activity (depth+batch residualized)",
       subtitle = paste0("Black outline = ", TARGET),
       x = "Program", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        panel.grid = element_blank())

ggsave(file.path(OUT_DIR, "FIG_D_coactivity_heatmap.pdf"), pD_fig, width = 11, height = 5, device = cairo_pdf)

cat("\nDone. Outputs in:", normalizePath(OUT_DIR), "\n")
cat("  visium_coactivity_matrix.csv\n  visium_module_competition_ranking.csv\n  visium_module_competition_pvalues.csv\n")
cat("  FIG_A_M2_rank_per_signature.pdf\n  FIG_B_rho_M2_vs_top.pdf\n  FIG_C_pvalue_summary.pdf\n  FIG_D_coactivity_heatmap.pdf\n")