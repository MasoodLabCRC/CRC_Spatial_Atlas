#!/usr/bin/env Rscript
# survival_cox_analysis.R
#
# Prognostic association of M2 activity (Fig. 6b-d).
# Discovery meta-cohort and GSE17536 are restricted to Stage II/III
# (Stage I dropped: 0-1 recurrence events, unstable reference category;
# Stage IV dropped: sparse/absent across Discovery sub-cohorts).
# TCGA-COADREAD is retained at its full Stage I-IV range, since its
# per-stage hazard ratio confidence intervals are stable across all four
# levels (unlike Discovery's).
#
# Discovery: age+stage+CMS+proliferation, strata(cohort), stage restricted to II/III
# GSE17536:  age+stage+proliferation, stage from 'stageall' column, restricted to II/III
#            (CMS not available for this cohort)
# TCGA:      age+stage+CMS+proliferation, full Stage I-IV
#
# Inputs (place in ./data/):
#   program_genes.csv        - gene membership per program (incl. M2)
#   GSE17536_eset.rds        - GSE17536 ExpressionSet
#
# TCGA-COADREAD expression/survival data are downloaded automatically via
# UCSCXenaTools on first run and cached under results/xena_gdc/.
#
# Outputs -> results/

suppressPackageStartupMessages({
  library(survival); library(Biobase); library(ExperimentHub)
  library(org.Hs.eg.db); library(ggplot2); library(dplyr)
})

PROG          <- "data/program_genes.csv"
GSE17536_RDS  <- "data/GSE17536_eset.rds"

OUTDIR <- "results"; dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
XDIR   <- file.path(OUTDIR, "xena_gdc"); dir.create(XDIR, showWarnings = FALSE, recursive = TRUE)

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
chr <- function(x) trimws(as.character(x))

# ---- gene sets --------------------------------------------------------------
pg <- read.csv(PROG, stringsAsFactors = FALSE); names(pg) <- tolower(names(pg))
m2genes <- unique(toupper(trimws(pg$gene[trimws(pg$program) == "M2"])))

g1s_genes <- c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
  "CDCA7","DTL","PRIM1","UHRF1","MLF1IP","HELLS","RFC2","RPA2","NASP",
  "RAD51AP1","GMNN","WDR76","SLBP","CCNE2","UBR7","POLD3","MSH2","ATAD2",
  "RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN","DSCC1","BLM","CASP8AP2",
  "USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8")
g2m_genes <- c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
  "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","FAM64A","SMC4","CCNB2",
  "CKAP2L","CKAP2","AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1",
  "KIF20B","HJURP","CDCA3","HN1","CDC20","TTK","CDC25C","KIF2C","RANGAP1",
  "NCAPD2","DLGAP5","CDCA2","CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1",
  "ANLN","LBR","CKAP5","CENPE","CTCF","NEK2","G2E3","GAS2L3","CBX5","CENPA")
prolif_genes <- unique(c(g1s_genes, g2m_genes))

score_signature <- function(E, genes) {
  g <- intersect(genes, rownames(E))
  if (length(g) < 3) return(setNames(rep(NA_real_, ncol(E)), colnames(E)))
  Z <- (E[g,] - rowMeans(E[g,])) / apply(E[g,], 1, sd)
  setNames(colMeans(Z), colnames(E))
}

# ---- helpers -----------------------------------------------------------------
stratify_tertile <- function(score) {
  q <- quantile(score, c(0.40, 0.60), na.rm = TRUE)
  cut(score, c(-Inf, q, Inf), labels = c("M2-low", "M2-mid", "M2-high"))
}

clean_stage <- function(x) {
  s <- toupper(chr(x)); s <- gsub("STAGE\\s*", "", s); s <- gsub("[ABC]$", "", s)
  s[s %in% c("", "NA", "N/A", "X", "UNK", "UNKNOWN", "[DISCREPANCY]", "NULL", "0")] <- NA
  factor(s)
}

stage_num_to_roman <- function(x) {
  x <- chr(x); map <- c("1" = "I", "2" = "II", "3" = "III", "4" = "IV")
  factor(unname(map[x]))
}

code_event <- function(v) {
  v <- tolower(chr(v)); e <- rep(NA_real_, length(v))
  e[grepl("recur|relaps|dead|decease|event|progress|yes", v)] <- 1
  e[grepl("^no|norecur|free|censor|living|alive|^0$", v)] <- 0
  e
}

run_cox <- function(d, formula_str, label_row = "grpM2-high") {
  cx <- summary(coxph(as.formula(formula_str), data = d))
  gi <- grep(paste0("^", label_row), rownames(cx$conf.int))
  if (length(gi) == 0) return(NULL)
  list(HR = cx$conf.int[gi, "exp(coef)"], lo = cx$conf.int[gi, "lower .95"],
       hi = cx$conf.int[gi, "upper .95"], p = cx$coefficients[gi, "Pr(>|z|)"])
}

extract_all_terms <- function(fit, cohort_name) {
  cs <- summary(fit)$conf.int; co <- summary(fit)$coefficients
  data.frame(cohort = cohort_name, term = rownames(cs),
             HR = cs[, "exp(coef)"], lo = cs[, "lower .95"], hi = cs[, "upper .95"],
             p = co[, "Pr(>|z|)"], row.names = NULL, stringsAsFactors = FALSE)
}

pretty_label <- function(x) {
  x <- gsub("^cms", "CMS_", x); x <- gsub("^stage", "Stage_", x); x <- gsub("^grp", "", x)
  x <- gsub("^proliferation$", "Proliferation", x); x <- gsub("^age$", "Age", x); x
}

plot_forest_covariates <- function(fit, keep, outfile, title = "Hazard Ratio", xlim = NULL) {
  cs <- summary(fit)$conf.int; co <- summary(fit)$coefficients
  terms <- rownames(cs); keep_pattern <- paste0("^(", paste(keep, collapse = "|"), ")")
  sel <- grepl(keep_pattern, terms)
  df <- data.frame(term = terms[sel], HR = cs[sel, "exp(coef)"], lo = cs[sel, "lower .95"],
                    hi = cs[sel, "upper .95"], p = co[sel, "Pr(>|z|)"], stringsAsFactors = FALSE)
  df$label <- factor(pretty_label(df$term), levels = rev(pretty_label(df$term)))
  p <- ggplot(df, aes(x = HR, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.7, color = "black") +
    geom_point(size = 3, color = "black") + labs(x = title, y = NULL) + theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank())
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  ggsave(outfile, p, width = 6, height = 0.6 * nrow(df) + 1.5, limitsize = FALSE)
  cat("Saved:", normalizePath(outfile), "\n")
}

main_results <- list(); all_terms <- list()

# DISCOVERY -- restricted to Stage II/III

cat("\n============ DISCOVERY COHORT (Stage II/III only) ============\n")
es <- query(ExperimentHub(), "mcsurvdata")[["EH1498"]]
E  <- exprs(es); ph <- pData(es)
sym <- mapIds(org.Hs.eg.db, keys = rownames(E), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
k <- !is.na(sym) & !duplicated(sym); E <- E[k,]; rownames(E) <- toupper(sym[k])
mss_mask <- toupper(chr(ph$msi)) == "MSS"; E <- E[, mss_mask]; ph <- ph[mss_mask, ]

M2_disc <- score_signature(E, m2genes)
Prolif_disc <- score_signature(E, prolif_genes)

d <- data.frame(
  sample_id = colnames(E), time = num(ph$tev), status = num(ph$evn),
  grp  = stratify_tertile(M2_disc), age  = num(ph$age), stage = clean_stage(ph$stage),
  cms = factor(chr(ph$cms)), cohort = factor(chr(ph$dataset)),
  proliferation = as.numeric(scale(Prolif_disc)), stringsAsFactors = FALSE
)
d <- d[is.finite(d$time) & d$time > 0 & !is.na(d$status), ]

d2 <- d[d$grp %in% c("M2-low","M2-high"), ]; d2$grp <- droplevels(d2$grp)
n_before <- nrow(d2)
d2 <- d2[d2$stage %in% c("II","III"), ]; d2$stage <- droplevels(d2$stage)
cat(sprintf("[Discovery] Stage II/III restriction: %d -> %d patients (%.1f%% kept)\n",
            n_before, nrow(d2), 100 * nrow(d2) / n_before))

sf  <- survfit(Surv(time, status) ~ grp, data = d[d$stage %in% c("II","III") & !is.na(d$stage), ])
lr  <- survdiff(Surv(time, status) ~ grp, data = d[d$stage %in% c("II","III") & !is.na(d$stage), ])
lrp <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
pdf(file.path(OUTDIR, "M2_kaplanmeier_discovery.pdf"), width = 6.2, height = 5.2)
plot(sf, col = c("#2a78d6","#B4B2A9","#e34948"), lwd = 2, xlab = "years",
     ylab = "recurrence-free survival",
     main = sprintf("Discovery (Stage II/III only)   log-rank P = %.2g", lrp))
legend("bottomleft", bty = "n", lwd = 2, col = c("#e34948","#B4B2A9","#2a78d6"),
       legend = c("M2-high","M2-mid","M2-low"))
dev.off()

fit_disc <- coxph(Surv(time,status) ~ grp + age + stage + cms + proliferation + strata(cohort), data = d2)
cox_disc <- run_cox(d2, "Surv(time,status) ~ grp + age + stage + cms + proliferation + strata(cohort)")
cat(sprintf("[Cox] Discovery (Stage II/III), M2-high vs low: HR = %.2f (%.2f-%.2f), P = %.2g\n",
            cox_disc$HR, cox_disc$lo, cox_disc$hi, cox_disc$p))

main_results$discovery <- data.frame(
  cohort = sprintf("Discovery (n=%d, Stage II/III only)\nage+stage+CMS+proliferation", nrow(d2)),
  HR = cox_disc$HR, lo = cox_disc$lo, hi = cox_disc$hi, p = cox_disc$p)
all_terms$discovery <- extract_all_terms(fit_disc, "Discovery (Stage II/III)")
plot_forest_covariates(fit_disc, keep = c("grp","age","stage","cms","proliferation"),
  outfile = file.path(OUTDIR, "M2_forest_covariates_discovery.pdf"),
  title = "Discovery (Stage II/III only) - Hazard Ratio")

# GSE17536 -- restricted to Stage II/III (using 'stageall' column)

if (file.exists(GSE17536_RDS)) {
  cat("\n============ GSE17536 (Stage II/III only) ============\n")
  g <- readRDS(GSE17536_RDS)
  Eg <- exprs(g); rownames(Eg) <- toupper(rownames(Eg)); phg <- pData(g)

  M2g <- score_signature(Eg, m2genes)
  Prolif_g <- score_signature(Eg, prolif_genes)
  rt <- intersect(c("days_to_recurrence","days_to_tumor_recurrence"), colnames(phg))[1]
  ag <- intersect(c("age_at_initial_pathologic_diagnosis","age"), colnames(phg))[1]
  st <- intersect(c("stageall","summarystage","tumorstage","stage"), colnames(phg))[1]

  stage_clean_g <- if (!is.na(st) && st == "stageall") stage_num_to_roman(phg[[st]]) else clean_stage(phg[[st]])

  dg <- data.frame(sample_id = colnames(Eg), time = num(phg[[rt]]),
                    status = code_event(phg$recurrence_status),
                    grp = stratify_tertile(M2g), age = num(phg[[ag]]),
                    stage = stage_clean_g,
                    proliferation = as.numeric(scale(Prolif_g)))
  dg <- dg[is.finite(dg$time) & dg$time > 0 & !is.na(dg$status), ]

  dg2 <- dg[dg$grp %in% c("M2-low","M2-high"), ]; dg2$grp <- droplevels(dg2$grp)
  n_before_g <- nrow(dg2)
  dg2 <- dg2[dg2$stage %in% c("II","III"), ]; dg2$stage <- droplevels(dg2$stage)
  cat(sprintf("[GSE17536] Stage II/III restriction: %d -> %d patients (%.1f%% kept)\n",
              n_before_g, nrow(dg2), 100 * nrow(dg2) / n_before_g))

  dg2_full <- dg2[complete.cases(dg2[, c("time","status","grp","age","stage","proliferation")]), ]
  fit_g <- coxph(Surv(time,status) ~ grp + age + stage + proliferation, data = dg2_full)
  cox_g <- run_cox(dg2_full, "Surv(time,status) ~ grp + age + stage + proliferation")

  pdf(file.path(OUTDIR, "M2_kaplanmeier_GSE17536.pdf"), width = 6.2, height = 5.2)
  plot(survfit(Surv(time,status) ~ grp, data = dg[dg$stage %in% c("II","III") & !is.na(dg$stage), ]),
       col = c("#2a78d6","#B4B2A9","#e34948"), lwd = 2, xlab = "months",
       ylab = "recurrence-free survival",
       main = sprintf("GSE17536 (Stage II/III only)   n=%d", nrow(dg2_full)))
  legend("bottomleft", bty = "n", lwd = 2, col = c("#e34948","#B4B2A9","#2a78d6"),
         legend = c("M2-high","M2-mid","M2-low"))
  dev.off()

  cat(sprintf("[GSE17536] age+stage+proliferation (Stage II/III): HR %.2f (%.2f-%.2f), P = %.2g\n",
              cox_g$HR, cox_g$lo, cox_g$hi, cox_g$p))

  main_results$gse17536 <- data.frame(
    cohort = sprintf("GSE17536 (n=%d, Stage II/III only)\nage+stage+proliferation", nrow(dg2_full)),
    HR = cox_g$HR, lo = cox_g$lo, hi = cox_g$hi, p = cox_g$p)
  all_terms$gse17536 <- extract_all_terms(fit_g, "GSE17536 (Stage II/III)")
  plot_forest_covariates(fit_g, keep = c("grp","age","stage","proliferation"),
    outfile = file.path(OUTDIR, "M2_forest_covariates_GSE17536.pdf"),
    title = "GSE17536 (Stage II/III only) - Hazard Ratio")
} else cat(sprintf("\n[GSE17536] %s not found - skipping\n", GSE17536_RDS))

# TCGA-COADREAD -- UNCHANGED, full Stage I-IV (stage coefficients are stable)

if (requireNamespace("UCSCXenaTools", quietly = TRUE)) {
  cat("\n============ TCGA-COADREAD (full Stage I-IV, unchanged) ============\n")
  suppressPackageStartupMessages(library(UCSCXenaTools))

  one <- function(pat, hub = "gdcHub") {
    hit <- XenaData[XenaData$XenaHostNames == hub & grepl(pat, XenaData$XenaDatasets), ]
    ds <- hit$XenaDatasets[1]
    dl <- XenaDownload(XenaQuery(XenaGenerate(subset = XenaDatasets == ds)), destdir = XDIR, force = FALSE)
    dl$destfiles[file.exists(dl$destfiles)][1]
  }
  rd <- function(p) { con <- if (grepl("\\.gz$", p)) gzfile(p) else p
    read.delim(con, row.names = 1, check.names = FALSE) }

  Ec <- as.matrix(rd(one("TCGA-COAD\\.star_fpkm\\.tsv$")))
  Er <- as.matrix(rd(one("TCGA-READ\\.star_fpkm\\.tsv$")))
  cg <- intersect(rownames(Ec), rownames(Er)); Et <- cbind(Ec[cg,], Er[cg,])
  cdr <- rd(one("Survival_SupplementalTable", "pancanAtlasHub"))
  ens <- sub("\\..*", "", rownames(Et))
  sym_t <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  keep <- !is.na(sym_t) & !duplicated(sym_t); Et <- Et[keep,]; rownames(Et) <- toupper(sym_t[keep])
  colnames(Et) <- substr(colnames(Et), 1, 15); rownames(cdr) <- substr(rownames(cdr), 1, 15)
  Et <- Et[, grepl("-01$", colnames(Et)) & !duplicated(colnames(Et)), drop = FALSE]
  cm <- intersect(colnames(Et), rownames(cdr)); Et <- Et[, cm]; cdr <- cdr[cm, ]

  M2t <- score_signature(Et, m2genes); Prolif_t <- score_signature(Et, prolif_genes)
  agec <- intersect(c("age_at_initial_pathologic_diagnosis","age"), colnames(cdr))[1]
  stc  <- intersect(c("ajcc_pathologic_tumor_stage","clinical_stage","tumor_stage"), colnames(cdr))[1]

  cms_col_t <- NA
  if (file.exists(CMS_LABELS_FILE)) {
    cms_all <- read.delim(CMS_LABELS_FILE, stringsAsFactors = FALSE)
    cms_tcga <- cms_all[grepl("tcga", cms_all$dataset, ignore.case = TRUE), ]
    cms_tcga$patient_barcode <- toupper(chr(cms_tcga$sample))
    lab <- chr(cms_tcga$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples)
    lab[lab %in% c("UNK", "NOLBL", "", "NA")] <- NA
    cms_lookup <- setNames(lab, cms_tcga$patient_barcode)
    cms_col_t <- "cms_synapse"
  }
  patient_bc_t <- setNames(substr(colnames(Et), 1, 12), colnames(Et))
  cms_vec_t <- if (!is.na(cms_col_t)) factor(unname(cms_lookup[patient_bc_t])) else NA

  dt <- data.frame(sample_id = colnames(Et), time = num(cdr$PFI.time)/30.44, status = num(cdr$PFI),
                    grp = stratify_tertile(M2t), age = num(cdr[[agec]]), stage = clean_stage(cdr[[stc]]),
                    proliferation = as.numeric(scale(Prolif_t)), cms = cms_vec_t)
  dt <- dt[is.finite(dt$time) & dt$time > 0 & !is.na(dt$status), ]

  pdf(file.path(OUTDIR, "M2_kaplanmeier_TCGA.pdf"), width = 6.2, height = 5.2)
  plot(survfit(Surv(time,status) ~ grp, data = dt), col = c("#2a78d6","#B4B2A9","#e34948"), lwd = 2,
       xlab = "months", ylab = "progression-free interval",
       main = sprintf("TCGA-COADREAD (full Stage I-IV, n=%d)", nrow(dt)))
  legend("bottomleft", bty = "n", lwd = 2, col = c("#e34948","#B4B2A9","#2a78d6"),
         legend = c("M2-high","M2-mid","M2-low"))
  dev.off()

  dt2 <- dt[dt$grp %in% c("M2-low","M2-high"), ]; dt2$grp <- droplevels(dt2$grp)
  if (!is.na(cms_col_t)) {
    dt2_full <- dt2[complete.cases(dt2[, c("time","status","grp","age","stage","cms","proliferation")]), ]
    fit_t <- coxph(Surv(time,status) ~ grp + age + stage + cms + proliferation, data = dt2_full)
    keep_t <- c("grp","age","stage","cms","proliferation"); adj_label_t <- "age+stage+CMS+proliferation"
  } else {
    dt2_full <- dt2[complete.cases(dt2[, c("time","status","grp","age","stage","proliferation")]), ]
    fit_t <- coxph(Surv(time,status) ~ grp + age + stage + proliferation, data = dt2_full)
    keep_t <- c("grp","age","stage","proliferation"); adj_label_t <- "age+stage+proliferation (no CMS)"
  }
  cox_t <- run_cox(dt2_full, paste("Surv(time,status) ~", as.character(formula(fit_t))[3]))
  cat(sprintf("[TCGA-GDC] %s: HR %.2f (%.2f-%.2f), P = %.2g\n", adj_label_t, cox_t$HR, cox_t$lo, cox_t$hi, cox_t$p))

  main_results$tcga <- data.frame(
    cohort = sprintf("TCGA-COADREAD (n=%d, full Stage I-IV)\n%s", nrow(dt2_full), adj_label_t),
    HR = cox_t$HR, lo = cox_t$lo, hi = cox_t$hi, p = cox_t$p)
  all_terms$tcga <- extract_all_terms(fit_t, "TCGA-COADREAD (full I-IV)")
  plot_forest_covariates(fit_t, keep = keep_t,
    outfile = file.path(OUTDIR, "M2_forest_covariates_TCGA.pdf"),
    title = "TCGA-COADREAD (full Stage I-IV) - Hazard Ratio")
} else cat("\n[TCGA] install UCSCXenaTools for TCGA validation\n")

# CROSS-COHORT FOREST PLOT + CSVs

all_results <- do.call(rbind, main_results)
all_results$cohort <- factor(all_results$cohort, levels = rev(all_results$cohort))
write.csv(all_results, file.path(OUTDIR, "M2_forest_stageII_III.csv"), row.names = FALSE)

p_forest <- ggplot(all_results, aes(x = HR, y = cohort)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.15, linewidth = 0.7, color = "#2a78d6") +
  geom_point(size = 3, color = "#e34948") +
  geom_text(aes(label = sprintf("HR=%.2f (%.2f-%.2f), P=%.2g", HR, lo, hi, p)),
            hjust = -0.1, size = 3.2, nudge_y = 0.15) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale), M2-high vs M2-low", y = NULL,
       title = "M2 prognostic value -- Discovery/GSE17536 restricted to Stage II/III, TCGA full I-IV") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))
ggsave(file.path(OUTDIR, "M2_forest_plot_stageII_III.pdf"), p_forest, width = 10, height = 4.5)

all_terms_df <- do.call(rbind, all_terms)
write.csv(all_terms_df, file.path(OUTDIR, "M2_all_covariates_stageII_III.csv"), row.names = FALSE)

cat("\n============ DONE (Stage II/III restriction variant) ============\n")
cat("Outputs in:", normalizePath(OUTDIR), "\n")
