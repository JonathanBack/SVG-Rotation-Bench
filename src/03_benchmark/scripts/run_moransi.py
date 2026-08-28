# ==============================================================================
# run_moransi.py
# Benchmarks Moran's I (via Squidpy) for SVG detection across rotated datasets.
# For each (slice, mode, angle): loads AnnData, builds spatial neighbor graph
# (Delaunay), runs Moran's I with 100 permutations, and saves results via RDS
# wrapper.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/02_rotation/outputs/{slice}/{mode}/anndata/data/scdesign3_angle{angle}.h5ad
# Output: src/03_benchmark/outputs/{slice}/{mode}/moransi/scdesign3_angle{angle}_results.rds
#         src/03_benchmark/outputs/{slice}/{mode}/moransi/scdesign3_angle{angle}_runtime.csv
# ==============================================================================

import os
import sys
import time
import subprocess
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq

# --- Select which slices and modes to run ---
SLICES = ["anterior1", "anterior2", "posterior1", "posterior2"]
MODES = ["simulated", "whole"]

ANGLES = [0, 30, 45, 60]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(BENCHMARK_DIR))

OUTPUTS_ROOT = os.path.join(BENCHMARK_DIR, "outputs")


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


if __name__ == "__main__":
    for slice_name in SLICES:
        for mode in MODES:

            h5ad_template = os.path.join(
                PROJECT_ROOT, "src", "02_rotation", "outputs", slice_name, mode,
                "anndata", "data", "scdesign3_angle{angle}.h5ad"
            )
            output_dir = os.path.join(OUTPUTS_ROOT, slice_name, mode, "moransi")
            os.makedirs(output_dir, exist_ok=True)

            if not os.path.exists(os.path.dirname(h5ad_template.format(angle=0))):
                print(f"Skipping {slice_name}/{mode} -- AnnData dir not found")
                continue

            print(f"\n=== Moran's I: {slice_name} / {mode} ===")

            for angle in ANGLES:
                rds_file = os.path.join(
                    output_dir, f"scdesign3_angle{angle}_results.rds"
                )
                runtime_file = os.path.join(
                    output_dir, f"scdesign3_angle{angle}_runtime.csv"
                )
                results_csv = os.path.join(
                    output_dir, f"scdesign3_angle{angle}_results.csv"
                )

                if os.path.exists(rds_file) and os.path.exists(runtime_file):
                    print(f"  Skipping angle {angle} -- outputs already exist")
                    continue

                h5ad_path = h5ad_template.format(angle=angle)
                if not os.path.exists(h5ad_path):
                    print(f"  Missing AnnData for angle {angle}, skipping")
                    continue

                print(f"  Running Moran's I for angle = {angle}")
                sys.stdout.flush()

                # --- Load AnnData with rotated spatial coordinates ---
                adata = sc.read_h5ad(h5ad_path)

                # --- Filter zero-count genes (harmonize gene universe across methods) ---
                counts_mat = adata.layers["counts"]
                if hasattr(counts_mat, "toarray"):
                    counts_mat = counts_mat.toarray()
                adata = adata[:, counts_mat.sum(axis=0) > 0].copy()

                # Build spatial neighbor graph using Delaunay triangulation
                sq.gr.spatial_neighbors_delaunay(adata)

                # --- Run Moran's I with 100 permutations, parallelized over 10 cores
                t_start = time.time()
                sq.gr.spatial_autocorr(
                    adata, mode="moran", n_perms=100, n_jobs=10,
                    genes=adata.var_names
                )
                elapsed = time.time() - t_start

                # Extract results and attach gene-level metadata from AnnData
                df_res = adata.uns["moranI"]
                df_res = df_res.loc[adata.var_names]
                # gene / spatial_var columns may be absent for whole mode; guard
                for col in ("gene", "spatial_var"):
                    if col in adata.var.columns:
                        df_res[col] = adata.var[col]
                df_res.to_csv(results_csv)

                save_rds(
                    df_res.reset_index().rename(
                        columns={df_res.index.name or "": "feature"}
                    ),
                    rds_file,
                )

                pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
                    runtime_file, index=False
                )

                print(f"    Saved {rds_file} ({elapsed:.1f}s)")

    print("\nMoran's I benchmark complete.")