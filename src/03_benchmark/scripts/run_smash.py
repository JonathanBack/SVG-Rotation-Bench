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

sys.path.insert(0, os.path.join(BENCHMARK_DIR, "tools", "SMASH"))
import SMASH
import Testing_functions as _tf

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


def save_rds(features, adjusted_pvals, rds_path):
    pvals_csv = rds_path.replace(".rds", "_tmp_pvals.csv")
    df = pd.DataFrame({"feature": features, "adjusted_pval": adjusted_pvals})
    df.to_csv(pvals_csv, index=False)
    r_code = (
        'df <- read.csv("' + pvals_csv + '");'
        'result <- list(res_mtest = data.frame('
        'adjustedPval = df$adjusted_pval, row.names = df$feature));'
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

    h5ad_path = H5AD_TEMPLATE.format(angle=angle)
    adata = sc.read_h5ad(h5ad_path)

    counts_mat = adata.layers["counts"]
    if hasattr(counts_mat, "toarray"):
        counts_mat = counts_mat.toarray()
    Y = pd.DataFrame(
        counts_mat,
        index=adata.obs_names,
        columns=adata.var_names,
    )

    Cords = pd.DataFrame(
        adata.obsm["spatial"],
        index=adata.obs_names,
        columns=["X", "Y"],
    )

    t_start = time.time()
    result = SMASH.SMASH(
        Y, Cords, len_l=10,
        mean_only="False", kernel_covariance="All", forcePD="False"
    )
    elapsed = time.time() - t_start

    smash_df = result["SMASH"]
    smash_df.to_csv(results_csv, index=False)

    features = list(smash_df["Gene"])
    _, adj_pvals, _, _ = multipletests(smash_df["p-val"], method="fdr_bh")

    save_rds(features, list(adj_pvals), rds_file)

    pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
        runtime_file, index=False
    )

    print(f"Saved {rds_file} ({elapsed:.1f}s)")

print("SMASH benchmark complete.")
