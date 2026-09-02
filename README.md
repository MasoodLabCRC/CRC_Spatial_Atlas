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

- Spatial transcriptomic and COMET imaging data: Zenodo (https://zenodo.org/records/22070345?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6IjFlY2ZjMDc4LTUwMjUtNGNlMi05MDFlLTMwNjhmY2EzNzBmZiIsImRhdGEiOnt9LCJyYW5kb20iOiI0ZGRlYWI2MWNmYzJjNmU4NzFlMzdhNTI5YjdlMWFjZSJ9.wmrcZXu8IKzf_So7ixvnEXHtJ8_fHE-6eV3otRTWb0GnX1XiawaRytXHKTuuCfvA_kYdD734NSxowq2svyTSHA)
- Raw sequencing data: GEO (accession to be added on publication)
- Mouse scRNA-seq data: GEO accession GSE290028
- TCGA-COADREAD expression/survival data: downloaded automatically at runtime via `UCSCXenaTools` (see `07_survival_cox_analysis/survival_cox_analysis.R`)
- GSE17536: publicly available via GEO; loaded as a Bioconductor `ExpressionSet`
- Bulk Caris Life Sciences data: available from the corresponding author upon reasonable request under a data-use agreement

Each script's header comment lists its expected input files and where to place them locally (typically a `data/` subfolder, which is intentionally excluded from this repository — see `.gitignore`).

## 1. System requirements

**Software dependencies and versions actually required to run this repository's code:**
- R 4.3.3
- Python 3.8.15
- Seurat v4.4.0
- ggplot2 (version not separately pinned in the manuscript)
- R `survival` package (version not separately pinned in the manuscript)
- Additional packages used by individual scripts: dplyr, tidyr, pheatmap, tibble, gridExtra, clue, data.table, Matrix, ggdendro, tidyverse, cowplot, patchwork, ggpubr, ggthemes, ggforce, pals, scales, readr, conflicted, Biobase, ExperimentHub, UCSCXenaTools, org.Hs.eg.db (R); pandas, numpy, scipy, statsmodels, matplotlib, scanpy, networkx, python-louvain (Python). See `environment.yml` for the full pinned list.


**Operating systems tested:** Linux (Ubuntu 22.04, via institutional HPC cluster with conda-managed R/Python environments). Not tested on Windows or macOS, but no OS-specific code is used.

**Non-standard hardware:** None required. All scripts run on a standard desktop/laptop; no GPU needed.

## 2. Installation guide

```bash
git clone [REPO_URL] CRC_Spatial_Atlas
cd CRC_Spatial_Atlas
conda env create -f environment.yml
conda activate crc-spatial-atlas
```

**Typical install time:** ~10-15 minutes on a normal desktop computer (dominated by R package compilation).

## 3. Demo

Each folder's `data/` subfolder contains real (non-simulated) input data — subsets of the study's spatial transcriptomic and single-cell data sufficient to run that component's analysis end-to-end without needing the full 24-section, 264,255-spot atlas.

**To run any component**, e.g. folder 01:
```bash
cd 01_core_to_edge_region_gradient
Rscript core_to_edge_region_gradient.R
```

**Expected output:** each script prints progress/summary statistics to the console and writes figures (PDF/PNG) and result tables (CSV) to a `results/` subfolder. See each folder's script header comment for its specific outputs.


## 4. Instructions for use

To run any script on your own data, replace the files in that folder's `data/` subfolder with your own, matching the same column names/structure documented in that script's header comment (input file names and expected columns are listed at the top of each `.R`/`.py` file).

### Reproducing manuscript results

Each script's header comment states which figure/table it reproduces. Running a script on the data provided in `data/` reproduces that figure/panel as shown in the manuscript.

## License

This project is licensed under the MIT License — see `LICENSE` for details.

## Contact

Corresponding author: Ashiq Masood (asmasood@iu.edu), Indiana University School of Medicine.
