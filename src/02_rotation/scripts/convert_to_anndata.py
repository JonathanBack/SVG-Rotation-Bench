import os
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import anndata as ad
import scipy as sp

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROTATION_DIR = os.path.dirname(SCRIPT_DIR)

LOCATIONS_DIR = os.path.join(ROTATION_DIR, "outputs", "locations")
ANNDATA_DATA_DIR = os.path.join(ROTATION_DIR, "outputs", "anndata", "data")
ANNDATA_FIGURES_DIR = os.path.join(ROTATION_DIR, "outputs", "anndata", "figures")

COUNTS_PATH = os.path.join(
    ROTATION_DIR, "..", "01_simulation", "outputs", "scDesign3", "data", "counts.csv"
)

ANGLES = [0, 30, 45, 60]

Ttr_features = [f"Ttr_{a}" for a in ["1", "0.8", "0.6", "0.4", "0.2", "0"]]
S100a5_features = [
    f"S100a5_{a}" for a in ["1", "0.8", "0.6", "0.4", "0.2", "0"]
]

for angle in ANGLES:
    h5ad_path = os.path.join(ANNDATA_DATA_DIR, f"scdesign3_angle{angle}.h5ad")
    loc_path = os.path.join(LOCATIONS_DIR, f"rotated_locations_{angle}.csv")

    if os.path.exists(h5ad_path):
        print(f"[angle={angle}] {h5ad_path} already exists, loading")
        adata = sc.read_h5ad(h5ad_path)
    else:
        print(f"[angle={angle}] Building AnnData from {loc_path}")

        df_loc = pd.read_csv(loc_path, index_col=0)
        df_count = pd.read_csv(COUNTS_PATH, index_col=0).transpose()

        df_var = pd.DataFrame(data={"feature_name": df_count.columns})
        df_var[["gene", "spatial_var"]] = df_var["feature_name"].str.rsplit(
            "_", n=1, expand=True
        )
        df_var = df_var.set_index("feature_name", drop=True)

        counts = sp.sparse.csr_matrix(df_count.values)
        adata = ad.AnnData(
            counts,
            obs=df_loc[["cell_type", "spatial1", "spatial2"]],
            obsm={"spatial": df_loc[["spatial1", "spatial2"]].values},
            var=df_var,
            dtype=np.float32,
        )

        sc.pp.calculate_qc_metrics(adata, percent_top=[10])
        adata.layers["counts"] = adata.X.copy()
        adata.uns["spatial"] = {"tissue": {}}

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

        adata.write_h5ad(h5ad_path)
        print(f"[angle={angle}] Saved {h5ad_path}")

    print(f"[angle={angle}] {adata}")

    sq.pl.spatial_scatter(
        adata, color=Ttr_features, library_id="tissue", ncols=6, shape=None, size=30
    )
    import matplotlib.pyplot as plt
    plt.savefig(os.path.join(ANNDATA_FIGURES_DIR, f"viz_angle{angle}_Ttr_spatial.png"), dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.violin(adata, keys=Ttr_features)
    plt.savefig(os.path.join(ANNDATA_FIGURES_DIR, f"viz_angle{angle}_Ttr_violin.png"), dpi=150, bbox_inches="tight")
    plt.close()

    sq.pl.spatial_scatter(
        adata, color=S100a5_features, library_id="tissue", ncols=6, shape=None, size=30
    )
    plt.savefig(os.path.join(ANNDATA_FIGURES_DIR, f"viz_angle{angle}_S100a5_spatial.png"), dpi=150, bbox_inches="tight")
    plt.close()

    sc.pl.violin(adata, keys=S100a5_features)
    plt.savefig(os.path.join(ANNDATA_FIGURES_DIR, f"viz_angle{angle}_S100a5_violin.png"), dpi=150, bbox_inches="tight")
    plt.close()

print("All rotation .h5ad files and visualizations saved.")
