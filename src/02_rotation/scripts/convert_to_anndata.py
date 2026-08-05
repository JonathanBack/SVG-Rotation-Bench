# ==============================================================================
# convert_to_anndata.py
# Builds AnnData (.h5ad) objects for each (slice, mode, angle) by combining the
# simulated/whole counts with the rotated spatial locations. Performs basic QC,
# normalization, and stores raw counts in a layer.
#
# Edit SLICES and MODES below to control scope.
# Input:  src/01_simulation/outputs/scDesign3/{slice}/{mode}/data/counts.csv
#         src/02_rotation/outputs/{slice}/{mode}/locations/rotated_locations_{angle}.csv
# Output: src/02_rotation/outputs/{slice}/{mode}/anndata/data/scdesign3_angle{angle}.h5ad
# ==============================================================================

import os
import matplotlib
matplotlib.use("Agg")  # non-interactive backend for headless execution

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scipy as sp

# --- Select which slices and modes to run ---
SLICES = ["anterior1", "anterior2", "posterior1", "posterior2"]
MODES = ["simulated", "whole"]

ANGLES = [0, 30, 45, 60]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROTATION_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(ROTATION_DIR))

COUNTS_TEMPLATE = os.path.join(
    PROJECT_ROOT, "src", "01_simulation", "outputs", "scDesign3",
    "{slice}", "{mode}", "data", "counts.csv"
)


def parse_var(feature_names, mode):
    """Build var metadata depending on mode.

    simulated: feature names follow "gene_alpha" format; split into gene and
               numeric alpha suffix.
    whole:     feature names are plain gene symbols; record gene name and NA
               spatial_var.
    """
    df_var = pd.DataFrame(data={"feature_name": feature_names})

    if mode == "simulated":
        split = df_var["feature_name"].str.rsplit("_", n=1, expand=True)
        # Guard: only treat suffix as alpha if it is numeric
        is_numeric = split[1].apply(lambda x: x is not None and _is_numeric(x))
        df_var["gene"] = split[0]
        df_var.loc[is_numeric, "spatial_var"] = split.loc[is_numeric, 1]
        df_var.loc[~is_numeric, "spatial_var"] = np.nan
        df_var.loc[~is_numeric, "gene"] = df_var.loc[~is_numeric, "feature_name"]
    else:
        df_var["gene"] = df_var["feature_name"]
        df_var["spatial_var"] = np.nan

    df_var = df_var.set_index("feature_name", drop=True)
    return df_var


def _is_numeric(s):
    try:
        float(s)
        return True
    except (TypeError, ValueError):
        return False


for slice_name in SLICES:
    for mode in MODES:
        counts_path = COUNTS_TEMPLATE.format(slice=slice_name, mode=mode)
        locations_dir = os.path.join(
            ROTATION_DIR, "outputs", slice_name, mode, "locations"
        )
        anndata_data_dir = os.path.join(
            ROTATION_DIR, "outputs", slice_name, mode, "anndata", "data"
        )
        os.makedirs(anndata_data_dir, exist_ok=True)

        if not os.path.exists(counts_path):
            print(f"Skipping {slice_name}/{mode} -- counts.csv not found")
            continue
        if not os.path.exists(locations_dir):
            print(f"Skipping {slice_name}/{mode} -- locations dir not found")
            continue

        print(f"\n=== Processing slice={slice_name} mode={mode} ===")

        # Count matrix is genes x cells (transpose to cells x genes)
        df_count = pd.read_csv(counts_path, index_col=0).transpose()

        # Pre-build var metadata from feature names (shared across angles)
        df_var = parse_var(df_count.columns, mode)

        for angle in ANGLES:
            h5ad_path = os.path.join(
                anndata_data_dir, f"scdesign3_angle{angle}.h5ad"
            )
            loc_path = os.path.join(
                locations_dir, f"rotated_locations_{angle}.csv"
            )

            if not os.path.exists(loc_path):
                print(f"  [angle={angle}] missing locations, skipping")
                continue

            # --- Skip rebuild if AnnData already exists ---
            if os.path.exists(h5ad_path):
                print(f"  [angle={angle}] {h5ad_path} already exists, skipping")
                continue

            print(f"  [angle={angle}] Building AnnData from {loc_path}")

            df_loc = pd.read_csv(loc_path, index_col=0)

            # Construct AnnData with sparse counts, spatial coords in .obsm
            counts = sp.sparse.csr_matrix(df_count.values)

            # Build obs metadata; preserve cell_type if present, else constant
            obs_cols = []
            if "cell_type" in df_loc.columns:
                obs_cols.append("cell_type")
            else:
                df_loc = df_loc.copy()
                df_loc["cell_type"] = "cell_type_1"
                obs_cols.append("cell_type")
            obs_cols += ["spatial1", "spatial2"]

            adata = ad.AnnData(
                counts,
                obs=df_loc[obs_cols],
                obsm={"spatial": df_loc[["spatial1", "spatial2"]].values},
                var=df_var.copy(),
                dtype=np.float32,
            )

            # QC, store raw counts, normalize for any downstream visualization
            sc.pp.calculate_qc_metrics(adata, percent_top=[10])
            adata.layers["counts"] = adata.X.copy()
            adata.uns["spatial"] = {"tissue": {}}

            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)

            adata.write_h5ad(h5ad_path)
            print(f"  [angle={angle}] Saved {h5ad_path}")

print("\nAll rotation .h5ad files built.")