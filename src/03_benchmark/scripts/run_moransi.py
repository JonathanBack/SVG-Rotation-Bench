# ==============================================================================
# run_moransi.py
# Benchmarks Moran's I (via Squidpy) for SVG detection across rotated datasets.
# For each angle: loads AnnData, builds spatial neighbor graph (Delaunay),
# runs Moran's I with 100 permutations, and saves results via RDS wrapper.
# Output: scdesign3_angle{angle}_results.rds and runtime CSV.
# ==============================================================================

import os
import sys
import time
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(BENCHMARK_DIR))

H5AD_TEMPLATE = os.path.join(
    PROJECT_ROOT, "src", "02_rotation", "outputs", "anndata",
    "data", "scdesign3_angle{angle}.h5ad"
)
OUTPUT_DIR = os.path.join(BENCHMARK_DIR, "outputs", "moransi")
os.makedirs(OUTPUT_DIR, exist_ok=True)

ANGLES = [0, 30, 45, 60]


# Helper: save a DataFrame as an R .rds file via CSV round-trip with Rscript.
# This keeps the output format consistent with R-based benchmark runners.
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

    print(f"Running Moran's I for angle = {angle}")
    sys.stdout.flush()

    # --- Load AnnData with rotated spatial coordinates ---
    h5ad_path = H5AD_TEMPLATE.format(angle=angle)
    adata = sc.read_h5ad(h5ad_path)

    # Build spatial neighbor graph using Delaunay triangulation
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # --- Run Moran's I with 100 permutations, parallelized over 10 cores ---
    t_start = time.time()
    sq.gr.spatial_autocorr(
        adata, mode="moran", n_perms=100, n_jobs=10, genes=adata.var_names
    )
    elapsed = time.time() - t_start

    # Extract results and attach gene-level metadata from AnnData
    df_res = adata.uns["moranI"]
    df_res = df_res.loc[adata.var_names]
    df_res[["gene", "spatial_var"]] = adata.var[["gene", "spatial_var"]]
    df_res.to_csv(results_csv)

    save_rds(df_res.reset_index().rename(columns={df_res.index.name or "": "feature"}), rds_file)

    pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
        runtime_file, index=False
    )

    print(f"Saved {rds_file} ({elapsed:.1f}s)")

print("Moran's I benchmark complete.")
