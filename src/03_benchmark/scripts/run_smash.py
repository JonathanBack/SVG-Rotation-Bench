# ==============================================================================
# run_smash.py
# Benchmarks the SMASH method for SVG detection across rotated datasets.
# For each angle: loads AnnData, extracts counts and coordinates, runs
# SMASH (mean + covariance tests), applies Benjamini-Hochberg FDR correction,
# and saves results via RDS wrapper. Includes a BLAS compatibility fix.
# Output: scdesign3_angle{angle}_results.rds and runtime CSV.
# ==============================================================================

import os
import sys
import time
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
from statsmodels.stats.multitest import multipletests

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(BENCHMARK_DIR))

# Add SMASH tool to Python path and import
sys.path.insert(0, os.path.join(BENCHMARK_DIR, "tools", "SMASH"))
import SMASH
import Testing_functions as _tf

# --- BLAS compatibility fix: monkey-patch sgemm to handle transposition flags ---
# The distributed SMASH code uses scipy.linalg.blas.sgemm which may not support
# the trans_a/trans_b parameters; this wrapper provides that interface.
def _sgemm_fixed(alpha=1, a=None, b=None, trans_a=0, trans_b=0):
    if trans_a:
        a = a.T
    if trans_b:
        b = b.T
    return alpha * (a @ b)

_tf.sgemm = _sgemm_fixed

H5AD_TEMPLATE = os.path.join(
    PROJECT_ROOT, "src", "02_rotation", "outputs", "anndata",
    "data", "scdesign3_angle{angle}.h5ad"
)
OUTPUT_DIR = os.path.join(BENCHMARK_DIR, "outputs", "smash")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ANGLES = [0, 30, 45, 60]


# Helper: save a DataFrame as an R .rds file via CSV round-trip with Rscript
def save_rds(df, rds_path):
    pvals_csv = rds_path.replace(".rds", "_tmp.csv")
    df.to_csv(pvals_csv, index=False)
    r_code = (
        'df <- read.csv("' + pvals_csv + '");'
        'rownames(df) <- df$feature;'
        'df$feature <- NULL;'
        'result <- list(res_mtest = df);'
        'saveRDS(result, "' + rds_path + '")'
    )
    subprocess.run(["Rscript", "-e", r_code], check=True)
    os.remove(pvals_csv)


for angle in ANGLES:
    rds_file = os.path.join(OUTPUT_DIR, f"scdesign3_angle{angle}_results.rds")
    runtime_file = os.path.join(OUTPUT_DIR, f"scdesign3_angle{angle}_runtime.csv")
    results_csv = os.path.join(OUTPUT_DIR, f"scdesign3_angle{angle}_results.csv")

    if os.path.exists(rds_file) and os.path.exists(runtime_file):
        print(f"Skipping angle {angle} -- outputs already exist")
        continue

    print(f"Running SMASH for angle = {angle}")
    sys.stdout.flush()

    # --- Load AnnData with rotated spatial coordinates ---
    h5ad_path = H5AD_TEMPLATE.format(angle=angle)
    adata = sc.read_h5ad(h5ad_path)

    # Extract raw counts (cells x genes) from the stored layer
    counts_mat = adata.layers["counts"]
    if hasattr(counts_mat, "toarray"):
        counts_mat = counts_mat.toarray()
    Y = pd.DataFrame(
        counts_mat,
        index=adata.obs_names,
        columns=adata.var_names,
    )

    # Extract 2D spatial coordinates
    Cords = pd.DataFrame(
        adata.obsm["spatial"],
        index=adata.obs_names,
        columns=["X", "Y"],
    )

    # --- Run SMASH: mean-only=False enables both mean and covariance tests ---
    t_start = time.time()
    result = SMASH.SMASH(
        Y, Cords, len_l=10,
        mean_only="False", kernel_covariance="All", forcePD="False"
    )
    elapsed = time.time() - t_start

    smash_df = result["SMASH"]
    smash_df.to_csv(results_csv, index=False)

    # Apply Benjamini-Hochberg FDR correction to raw p-values
    _, adj_pvals, _, _ = multipletests(smash_df["p-val"], method="fdr_bh")

    # Rename columns for downstream metric script compatibility
    out_df = smash_df.rename(columns={"Gene": "feature", "p-val": "pval"}).copy()
    out_df["adjusted_pval"] = adj_pvals

    save_rds(out_df, rds_file)

    pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
        runtime_file, index=False
    )

    print(f"Saved {rds_file} ({elapsed:.1f}s)")

print("SMASH benchmark complete.")
