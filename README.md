# SVG-Rotation-Bench

**Benchmarking Rotation Invariance and Performance of Kernel-Based Methods for Spatially Variable Genes (SVGs) Detection**

## Overview

Advances in spatially resolved transcriptomics (SRT) require robust computational tools to identify Spatially Variable Genes (SVGs). However, many kernel-based methods exhibit a critical, yet overlooked, technical vulnerability: **rotation variance**. In laboratory practice, tissue sections are randomly positioned on slides, altering the absolute spatial coordinates. Theoretically, the biological signal remains unchanged, but algorithms relying on axis-dependent transformations may produce highly discordant results when the coordinate system is rotated.

This repository provides a systematic benchmarking framework designed to evaluate how the underlying statistical formulation of kernel-based methods affects rotation invariance, classification accuracy, ranking power, and computational scalability.

## Scope and Included Methods

We benchmark 4 representative kernel-based methods across two main statistical schools:

**1. Dependence Tests** (Evaluating global independence between expression and location):
- **SPARK-X:** Non-parametric covariance test.
- **SMASH:** Generalized non-parametric method incorporating distance-based kernels.

**2. Random-Effect Regression** (Modeling spatial variation as a Gaussian Process):
- **SpatialDE:** Linear mixed model using standard Gaussian Processes (GP).
- **nnSVG:** Nearest-Neighbor Gaussian Process (NNGP) for linear scalability.

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
├── 03_benchmark/         # Run SVG detection methods across all rotations
└── 04_metrics/           # Collect and compare results across methods/rotations
```

## Setup

### Python environment (conda)

```bash
conda env create -f environment.yml
conda activate svg-rotation-bench
```

### R dependencies

The simulation scripts require the following R packages:

```
Seurat, SingleCellExperiment, scDesign3, scales, ggplot2, cowplot, dplyr, glue, SPARK
```

## Stage 1: Simulation (`src/01_simulation/`)

Realistic synthetic data is generated using `scDesign3`, fitted on the **10x Visium mouse brain** reference dataset (`data/VISIUM_sce.rds`).

- **Marginal Modeling:** Learns the true spatial mean $\mu_s(s)$ using a 2D Gaussian Process.
- **Joint Modeling:** Preserves gene-gene correlation using a Gaussian Copula.
- **Ground Truth Generation:** True SVGs are defined by mixing the spatial signal with a randomized null signal $\mu_{ns}(s)$ using a mixing parameter $\alpha$:  
  $\mu(s) = \alpha \cdot \mu_s(s) + (1 - \alpha) \cdot \mu_{ns}(s)$

**Scripts:**
- `scripts/run_scdesign3.R` — Runs scDesign3 fitting and simulation. Outputs `counts.csv` (300 features: 50 genes × 6 alpha levels) and `location.csv` (2696 spots with spatial coordinates) to `outputs/scDesign3/data/`. Diagnostic plots saved to `outputs/scDesign3/figures/`.

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

## Stage 4: Metrics (`src/04_metrics/`)

Performance is evaluated using the following metrics:
- **Classification Accuracy:** Area under the Precision-Recall curve (**auPRC**) for classifying true SVGs vs. non-SVGs.
- **Ranking Accuracy:** **Kendall's rank correlation coefficient ($\tau$)** to evaluate the biological prioritization of genes.
- **Computational Scalability:** Peak RAM usage (GB) and total running time.
