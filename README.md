# SVG-Rotation-Bench 🧬🔄

**Benchmarking Rotation Invariance and Performance of Kernel-Based Methods for Overall Spatially Variable Genes (SVGs) Detection**

## 📌 Overview
Advances in spatially resolved transcriptomics (SRT) require robust computational tools to identify Spatially Variable Genes (SVGs). However, many kernel-based methods exhibit a critical, yet overlooked, technical vulnerability: **rotation variance**. In laboratory practice, tissue sections are randomly positioned on slides, altering the absolute spatial coordinates. Theoretically, the biological signal remains unchanged, but algorithms relying on axis-dependent transformations may produce highly discordant results when the coordinate system is rotated.

This repository provides a systematic benchmarking framework designed to evaluate how the underlying statistical formulation of kernel-based methods—specifically **Dependence Tests vs. Random-Effect Regression**—affects rotation invariance, classification accuracy, ranking power, and computational scalability.

## 🔬 Scope and Included Methods
We benchmark 4 representative kernel-based methods across two main statistical schools:

**1. Dependence Tests** (Evaluating global independence between expression and location):
*   **SPARK-X:** Non-parametric covariance test.
*   **SMASH:** Generalized non-parametric method incorporating distance-based kernels.

**2. Random-Effect Regression** (Modeling spatial variation as a Gaussian Process):
*   **SpatialDE:** Linear mixed model using standard Gaussian Processes (GP).
*   **nnSVG:** Nearest-Neighbor Gaussian Process (NNGP) for linear scalability.

## ⚙️ Methodology

### 1. Realistic Data Simulation (`scDesign3`)
To avoid artificial patterns, we use `scDesign3` to simulate realistic benchmarking data based on the **10x Visium Human Dorsolateral Prefrontal Cortex (DLPFC)** dataset.
*   **Marginal Modeling:** Learns the true spatial mean $\mu_s(s)$ using a 2D Gaussian Process.
*   **Joint Modeling:** Preserves gene-gene correlation using a Gaussian Copula.
*   **Ground Truth Generation:** We define true SVGs and non-SVGs by mixing the spatial signal with a randomized null signal $\mu_{ns}(s)$ using a mixing parameter $\alpha$:  
    $\mu(s) = \alpha \cdot \mu_s(s) + (1 - \alpha) \cdot \mu_{ns}(s)$

### 2. Rotation Invariance Testing
The generated count matrices are frozen. To test geometrical robustness, the original spatial coordinates matrix $S \in \mathbb{R}^{n \times 2}$ is multiplied by an orthogonal rotation matrix $R$ at specific angles $\theta \in \{0^\circ, 30^\circ, 45^\circ, 60^\circ, 90^\circ\}$:
$$S^* = S R^T$$
The algorithms are evaluated on their ability to retain consistent results using the rotated coordinates $S^*$.

### 3. Evaluation Metrics
Performance is evaluated using the following metrics:
*   **Classification Accuracy:** Area under the Precision-Recall curve (**auPRC**) for classifying true SVGs vs. non-SVGs.
*   **Ranking Accuracy:** **Kendall's rank correlation coefficient ($\tau$)** to evaluate the biological prioritization of genes.
*   **Computational Scalability:** Peak RAM usage (GB) and total running time.
