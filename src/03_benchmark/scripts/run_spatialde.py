# ==============================================================================
# run_spatialde.py
# Benchmarks the SpatialDE method for SVG detection across rotated datasets.
# For each (slice, mode, angle): loads AnnData, applies NaiveDE normalization
# (stabilize + regress out total counts), runs SpatialDE, and saves results via
# RDS wrapper. Includes SciPy/NumPy compatibility shims for older SpatialDE.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/02_rotation/outputs/{slice}/{mode}/anndata/data/scdesign3_angle{angle}.h5ad
# Output: src/03_benchmark/outputs/{slice}/{mode}/spatialde/scdesign3_angle{angle}_results.rds
#         src/03_benchmark/outputs/{slice}/{mode}/spatialde/scdesign3_angle{angle}_runtime.csv
# ==============================================================================

import os
import sys
import time
import subprocess
import numpy as np
import pandas as pd
import scipy.misc
import scipy.optimize
import scanpy as sc
import NaiveDE

import scipy
import numpy as np

# --- SciPy/NumPy compatibility shims for SpatialDE ---
# SpatialDE depends on scipy.misc.derivative (removed in scipy 1.10+) and
# expects certain np namespace objects to be accessible via scipy.*.
# These patches restore the required API surface.

# Copy all public NumPy symbols into scipy namespace
for name in dir(np):
    if not name.startswith("_") and not hasattr(scipy, name):
        try:
            setattr(scipy, name, getattr(np, name))
        except (TypeError, AttributeError):
            pass

# Replace deprecated scipy.misc.derivative with approx_fprime
scipy.misc.derivative = (
    lambda func, x0, dx=1.0, n=1, args=(), order=3:
    scipy.optimize.approx_fprime(x0, func, abs(dx), *args)[0]
)

import SpatialDE

# --- Pandas compatibility shim for SpatialDE qvalue ---
# SpatialDE's qvalue() calls pv.ravel(), but modern pandas Series no longer
# has this method. Convert to numpy array before passing to the original.
_original_qvalue = SpatialDE.util.qvalue
SpatialDE.util.qvalue = lambda pv: _original_qvalue(np.asarray(pv))

# --- Select which slices and modes to run ---
SLICES = ["anterior1", "anterior2", "posterior1", "posterior2"]
MODES = ["simulated", "whole"]

ANGLES = [0, 30, 45, 60]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(BENCHMARK_DIR))

OUTPUTS_ROOT = os.path.join(BENCHMARK_DIR, "outputs")


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


for slice_name in SLICES:
    for mode in MODES:

        h5ad_template = os.path.join(
            PROJECT_ROOT, "src", "02_rotation", "outputs", slice_name, mode,
            "anndata", "data", "scdesign3_angle{angle}.h5ad"
        )
        output_dir = os.path.join(OUTPUTS_ROOT, slice_name, mode, "spatialde")
        os.makedirs(output_dir, exist_ok=True)

        if not os.path.exists(os.path.dirname(h5ad_template.format(angle=0))):
            print(f"Skipping {slice_name}/{mode} -- AnnData dir not found")
            continue

        print(f"\n=== SpatialDE: {slice_name} / {mode} ===")

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

            print(f"  Running SpatialDE for angle = {angle}")
            sys.stdout.flush()

            # --- Load AnnData with rotated spatial coordinates ---
            adata = sc.read_h5ad(h5ad_path)
            sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[10])

            # Extract count matrix and total counts per cell
            counts = sc.get.obs_df(
                adata, keys=list(adata.var_names), use_raw=False, layer="counts"
            )
            total_counts = sc.get.obs_df(adata, keys=["total_counts"])

            # --- NaiveDE normalization pipeline ---
            # 1. Variance-stabilizing transformation
            # 2. Regress out log total_counts to remove library size effects
            t_start = time.time()
            norm_expr = NaiveDE.stabilize(counts.T).T
            resid_expr = NaiveDE.regress_out(
                total_counts, norm_expr.T, "np.log(total_counts)"
            ).T

            # --- Run SpatialDE on the residual expression ---
            df_res = SpatialDE.run(adata.obsm["spatial"], resid_expr)
            elapsed = time.time() - t_start

            # Index results by gene name and attach metadata from AnnData
            df_res.set_index("g", inplace=True)
            df_res = df_res.loc[adata.var_names]
            # gene / spatial_var columns may be absent for whole mode; guard
            for col in ("gene", "spatial_var"):
                if col in adata.var.columns:
                    df_res[col] = adata.var[col]
            df_res.to_csv(results_csv)

            # Deduplicate before saving RDS (SpatialDE may produce duplicate rows)
            dedup = df_res.copy()
            dedup.insert(0, "feature", dedup.index)
            dedup = dedup.drop_duplicates(subset="feature", keep="first")

            save_rds(dedup, rds_file)

            pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
                runtime_file, index=False
            )

            print(f"    Saved {rds_file} ({elapsed:.1f}s)")

print("\nSpatialDE benchmark complete.")