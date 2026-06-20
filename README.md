# SVG-Rotation-Bench

**Benchmarking Rotation Invariance and Performance of Kernel-Based Methods for Spatially Variable Genes (SVGs) Detection**

## Overview

Advances in spatially resolved transcriptomics (SRT) require robust computational tools to identify Spatially Variable Genes (SVGs). However, many kernel-based methods exhibit a critical, yet overlooked, technical vulnerability: **rotation variance**. In laboratory practice, tissue sections are randomly positioned on slides, altering the absolute spatial coordinates. Theoretically, the biological signal remains unchanged, but algorithms relying on axis-dependent transformations may produce highly discordant results when the coordinate system is rotated.

This repository provides a systematic benchmarking framework designed to evaluate how the underlying statistical formulation of kernel-based methods affects rotation invariance, classification accuracy, ranking power, and computational scalability.

## Scope and Included Methods

We benchmark 5 representative kernel-based methods across two main statistical schools:

**1. Dependence Tests** (Evaluating global independence between expression and location):
- **SPARK-X:** Non-parametric covariance test.
- **Moran's I:** Spatial autocorrelation statistic (via squidpy).
- **SMASH:** Generalized non-parametric method incorporating distance-based kernels.

**2. Random-Effect Regression** (Modeling spatial variation as a Gaussian Process):
- **SpatialDE:** Linear mixed model using standard Gaussian Processes (GP).
- **nnSVG:** Nearest-Neighbor Gaussian Process (NNGP) for linear scalability.
- **BOOST-GP:** Bayesian hierarchical model with GP-based overdispersion.

## Pipeline

The benchmark is organized into four sequential stages:

```
src/
├── 01_simulation/
│   ├── scripts/               # R simulation scripts
│   └── outputs/scDesign3/
│       ├── data/              # counts.csv, location.csv, sim_sce.rds
│       └── figures/           # Diagnostic plots (.png, .pdf)
├── 02_rotation/
│   ├── scripts/               # R rotation script + Python .h5ad converter
│   └── outputs/
│       ├── locations/         # rotated_locations_{angle}.csv
│       └── anndata/
│           ├── data/          # scdesign3_angle{angle}.h5ad
│           └── figures/       # Verification plots (.png)
├── 03_benchmark/              # Run SVG detection methods across all rotations
│   ├── scripts/               # Per-method benchmark scripts
│   ├── tools/                 # External tools (BOOST-GP, SMASH)
│   └── outputs/{method}/      # Raw results per method/angle
└── 04_metrics/                # Collect and compare results across methods/rotations
    ├── scripts/
    └── outputs/
        ├── comparison/        # Cross-method figures and CSVs
        └── {method}/          # Per-method figures and CSVs
```

## Setup

### 1. Python environment (conda)

The Python-based methods and data conversion scripts require the following environment:

```bash
conda env create -f environment.yml
conda activate svg-rotation-bench
```

This installs: `numpy`, `scipy`, `pandas`, `matplotlib`, `h5py`, `scikit-learn`, `scikit-image`, `geopandas`, `scanpy`, `squidpy`, `anndata`, `statsmodels`, `tqdm`, plus the PyPI packages `SpatialDE` and `NaiveDE`.

### 2. R dependencies

The R-based methods and metrics require packages from both **CRAN** and **Bioconductor**.

#### CRAN packages

```r
install.packages(c(
  "ggplot2", "reshape2", "RColorBrewer", "patchwork",
  "cowplot", "dplyr", "scales", "MASS", "ggVennDiagram"
))
```

#### Bioconductor packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SingleCellExperiment", "SpatialExperiment", "scran",
  "nnSVG", "anndata"
))
```

#### Seurat (CRAN / remotes)

```r
install.packages("Seurat")
```

#### SPARK-X (GitHub)

SPARK-X is not available on CRAN or Bioconductor. Install from the authors' GitHub repository:

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("xzhoulab/SPARK")
```

> **Note:** `SPARK::sparkx()` is the function used in the benchmark. The package may have additional system dependencies (e.g., a working C++ compiler for Rcpp).
> 
#### scDesign3 (GitHub)

scDesign3 is not available on CRAN or Bioconductor. Install from the authors' GitHub repository:

```r
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("SONGDONGYUAN1994/scDesign3")
```

### 3. External tools (manual download)

Two methods require manual cloning/downloading into `src/03_benchmark/tools/`:

#### SMASH

Clone the SMASH repository into the tools directory:

```bash
cd src/03_benchmark/tools/
git clone https://github.com/your-org/SMASH.git  # Replace with actual URL
```

SMASH is imported directly from this local path by `run_smash.py`. It requires `numpy`, `scipy`, `pandas`, `scikit-learn`, and `tqdm` (all provided by the conda environment).

> **Compatibility note:** SMASH uses the deprecated BLAS function `sgemm` which fails with float64 arrays on modern scipy. The benchmark script (`run_smash.py`) includes a monkey-patch to fix this issue automatically.

#### BOOST-GP

Clone the BOOST-GP repository into the tools directory:

```bash
cd src/03_benchmark/tools/
git clone https://github.com/your-org/BOOST-GP.git  # Replace with actual URL
```

The benchmark script (`run_boostgp.R`) sources `R/boost.gp.R` from this directory. BOOST-GP uses MCMC sampling and is **computationally expensive** (~hours per dataset). The benchmark uses `iter=100, burn=50` for feasibility.

### 4. SpatialDE compatibility patch

SpatialDE is unmaintained and incompatible with `scipy >= 1.12` (removed `scipy.misc.derivative` and `scipy.arange`). The benchmark script (`run_spatialde.py`) includes runtime compatibility patches to address these issues. No manual intervention is required.

### 5. Full dependency summary

| Package | Source | Used by | Notes |
|---------|--------|---------|-------|
| Python `numpy` | conda | All Python scripts | |
| Python `scipy` | conda | SpatialDE, SMASH | |
| Python `pandas` | conda | All Python scripts | |
| Python `matplotlib` | conda | SMASH | |
| Python `h5py` | conda | anndata | |
| Python `scikit-learn` | conda | SMASH | |
| Python `scanpy` | conda | Moran's I, SpatialDE, SMASH | |
| Python `squidpy` | conda | Moran's I | |
| Python `anndata` | conda | All Python scripts | |
| Python `statsmodels` | conda | SMASH (FDR correction) | |
| Python `tqdm` | conda | SMASH | |
| Python `SpatialDE` | pip | SpatialDE | Unmaintained; patched at runtime |
| Python `NaiveDE` | pip | SpatialDE | Companion to SpatialDE |
| R `Seurat` | CRAN | Simulation | |
| R `SingleCellExperiment` | Bioc | Simulation, nnSVG, BOOST-GP | |
| R `SpatialExperiment` | Bioc | nnSVG | |
| R `scran` | Bioc | nnSVG | |
| R `nnSVG` | Bioc | nnSVG | |
| R `scDesign3` | Bioc | Simulation | |
| R `anndata` (R pkg) | Bioc | nnSVG, BOOST-GP | R interface to .h5ad files |
| R `SPARK` | GitHub | SPARK-X | `devtools::install_github("xzhoulab/SPARK")` |
| R `ggplot2` | CRAN | Metrics, simulation | |
| R `reshape2` | CRAN | Metrics | |
| R `RColorBrewer` | CRAN | Metrics | |
| R `patchwork` | CRAN | Metrics | |
| R `cowplot` | CRAN | Simulation | |
| R `dplyr` | CRAN | Simulation | |
| R `scales` | CRAN | Simulation | |
| R `MASS` | CRAN | SPARK-X | |
| R `ggVennDiagram` | CRAN | Metrics (optional) | For Venn diagrams; skip if unavailable |

## Stage 1: Simulation (`src/01_simulation/`)

Realistic synthetic data is generated using `scDesign3`, fitted on the **10x Visium mouse brain** reference dataset (`data/VISIUM_sce.rds`).

- **Marginal Modeling:** Learns the true spatial mean $\mu_s(s)$ using a 2D Gaussian Process.
- **Joint Modeling:** Preserves gene-gene correlation using a Gaussian Copula.
- **Ground Truth Generation:** True SVGs are defined by mixing the spatial signal with a randomized null signal $\mu_{ns}(s)$ using a mixing parameter $\alpha$:  
  $\mu(s) = \alpha \cdot \mu_s(s) + (1 - \alpha) \cdot \mu_{ns}(s)$

**Scripts:**
- `scripts/run_scdesign3_fit.R` — Runs scDesign3 fitting (GP marginal models + Gaussian copula) and saves fitted parameters to `outputs/scDesign3/data/scdesign3_fitted.rds`. Runs once; also generates diagnostic plots.
- `scripts/run_scdesign3_sweep.R` — Loads the fitted parameters and generates the alpha-sweep count matrix at 21 mixing levels (0.00–1.00 in 0.05 steps). Outputs `counts.csv` (1050 features: 50 genes × 21 alpha levels) and `location.csv` (2696 spots) to `outputs/scDesign3/data/`.

## Stage 2: Rotation (`src/02_rotation/`)

The generated count matrix is frozen. To test geometrical robustness, the original spatial coordinates matrix $S \in \mathbb{R}^{n \times 2}$ is multiplied by an orthogonal rotation matrix $R$ at specific angles $\theta \in \{0^\circ, 30^\circ, 45^\circ, 60^\circ\}$:

$$S^* = S R^T$$

Where $R$ is the standard 2D rotation matrix:

$$R = \begin{bmatrix}
\cos(\frac{\theta}{180}\pi) & -\sin(\frac{\theta}{180}\pi) \\
\sin(\frac{\theta}{180}\pi) &  \cos(\frac{\theta}{180}\pi)
\end{bmatrix}$$

Rotation is centered on the tissue's center of mass to preserve the overall spatial distribution.

**Scripts:**
- `scripts/generate_rotated_locations.R` — Reads the original `location.csv` and outputs per-angle rotated coordinate files to `outputs/locations/`.
- `scripts/convert_to_anndata.py` — Builds separate AnnData `.h5ad` files for each rotation angle, combining the fixed `counts.csv` with the rotated spatial coordinates. `.h5ad` files saved to `outputs/anndata/data/`, verification plots to `outputs/anndata/figures/`.

## Stage 3: Benchmark (`src/03_benchmark/`)

Each SVG detection method is implemented as a standalone script under `scripts/`. Each script reads the `.h5ad` or CSV files from Stage 2 and runs the method across all four rotation angles. Raw outputs (p-values, test statistics per gene) are saved to `outputs/{method}/`.

**Scripts:**
- `scripts/run_sparkx.R` — SPARK-X non-parametric covariance test per angle.
- `scripts/run_nnsvg.R` — nnSVG nearest-neighbor Gaussian process per angle.
- `scripts/run_spatialde.py` — SpatialDE Gaussian process regression per angle.
- `scripts/run_moransi.py` — Moran's I spatial autocorrelation per angle.
- `scripts/run_smash.py` — SMASH non-parametric kernel test per angle.
- `scripts/run_boostgp.R` — BOOST-GP Bayesian GP model per angle (very slow; not run by default).

## Stage 4: Metrics (`src/04_metrics/`)

Raw benchmark outputs from Stage 3 are evaluated using classification and ranking metrics across rotation angles.

**Scripts:**
- `scripts/compute_all_metrics.R` — Reads `.rds` results per method/angle, computes:
  - **Threshold-independent:** AUC-ROC, auPRC, Kendall's $\tau_\alpha$ (alpha vs score)
  - **Rotation invariance:** Kendall's $\tau_{rotation}$ (cross-angle score consistency)
  - **Classification metrics (adj p < 0.05):** Sensitivity, Specificity, Precision, FDR, F1, MCC, Balanced Accuracy
  - **Set overlap:** Jaccard index, consistency rate, Venn diagrams
  - **Classification stability:** Spearman $\rho$ between angle and classification metric

Outputs CSVs and publication-ready figures (PNG 300dpi + PDF) to `outputs/comparison/` and `outputs/{method}/`.
