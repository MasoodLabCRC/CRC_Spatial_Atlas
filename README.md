# CRC Spatial Atlas

Analysis code accompanying **"Tissue position resolves colorectal cancer plasticity states to three compartments."**

This repository contains the analysis scripts only. Data are not hosted here — see **Data availability** below.

## Repository structure

| Folder | Script | Purpose |
|---|---|---|
| `01_core_to_edge_region_gradient/` | `core_to_edge_region_gradient.R` | Core-to-edge module gradient analysis (Fig. 2b, 2d) |
| `02_signature_module_coactivity_tests/` | `module_competition_test.R` | Signature x module co-activity / module-competition test (Fig. 3c–e, Ext. Data Fig. 3) |
| `03_cross_tumor_reproducibility/` | `cross_tumor_reproducibility.py` | Cross-tumor reproducibility of module region profiles (Fig. 2e, Supp. Table 9) |
| `04_signature_clustering_upgma/` | `signature_clustering_upgma.R` | UPGMA clustering of plasticity signatures (Fig. 3f) |
| `05_comet_visium_alignment/` | `01_comet_visium_alignment.R`, `02_comet_visium_scatter.R` | COMET–Visium ecotype cross-platform validation (Fig. 1e, Ext. Data Fig. 1e) |
| `06_spatial_context_map/` | `spatial_context_map.py` | Spatial context / ecotype adjacency map (Fig. 1f) |
| `07_survival_cox_analysis/` | `survival_cox_analysis.R` | Prognostic association of M2 activity across survival cohorts (Fig. 6b–d) |

## Data availability

Data used by these scripts are deposited externally, not in this repository:

- Spatial transcriptomic and COMET imaging data: Zenodo (DOI to be added on publication)
- Raw sequencing data: GEO (accession to be added on publication)
- Mouse scRNA-seq data: GEO accession GSE290028
- TCGA-COADREAD expression/survival data: downloaded automatically at runtime via `UCSCXenaTools` (see `07_survival_cox_analysis/survival_cox_analysis.R`)
- GSE17536: publicly available via GEO; loaded as a Bioconductor `ExpressionSet`
- Bulk Caris Life Sciences data: available from the corresponding author upon reasonable request under a data-use agreement

Each script's header comment lists its expected input files and where to place them locally (typically a `data/` subfolder, which is intentionally excluded from this repository — see `.gitignore`).

## Requirements

R (Seurat, dplyr, tidyr, ggplot2, survival, data.table, clue, Matrix, Biobase, ExperimentHub, org.Hs.eg.db, UCSCXenaTools) and Python 3.8 (pandas, numpy, scipy, statsmodels, matplotlib). See individual script headers for the exact packages each one needs.
