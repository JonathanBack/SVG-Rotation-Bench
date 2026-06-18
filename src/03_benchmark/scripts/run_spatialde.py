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

for name in dir(np):
    if not name.startswith("_") and not hasattr(scipy, name):
        try:
            setattr(scipy, name, getattr(np, name))
        except (TypeError, AttributeError):
            pass

scipy.misc.derivative = (
    lambda func, x0, dx=1.0, n=1, args=(), order=3:
    scipy.optimize.approx_fprime(x0, func, abs(dx), *args)[0]
)

import SpatialDE

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BENCHMARK_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(BENCHMARK_DIR))

H5AD_TEMPLATE = os.path.join(
    PROJECT_ROOT, "src", "02_rotation", "outputs", "anndata",
    "data", "scdesign3_angle{angle}.h5ad"
)
OUTPUT_DIR = os.path.join(BENCHMARK_DIR, "outputs", "spatialde")
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

    print(f"Running SpatialDE for angle = {angle}")
    sys.stdout.flush()

    h5ad_path = H5AD_TEMPLATE.format(angle=angle)
    adata = sc.read_h5ad(h5ad_path)
    sc.pp.calculate_qc_metrics(adata, inplace=True, percent_top=[10])

    counts = sc.get.obs_df(
        adata, keys=list(adata.var_names), use_raw=False, layer="counts"
    )
    total_counts = sc.get.obs_df(adata, keys=["total_counts"])

    t_start = time.time()
    norm_expr = NaiveDE.stabilize(counts.T).T
    resid_expr = NaiveDE.regress_out(
        total_counts, norm_expr.T, "np.log(total_counts)"
    ).T
    df_res = SpatialDE.run(adata.obsm["spatial"], resid_expr)
    elapsed = time.time() - t_start

    df_res.set_index("g", inplace=True)
    df_res = df_res.loc[adata.var_names]
    df_res[["gene", "spatial_var"]] = adata.var[["gene", "spatial_var"]]
    df_res.to_csv(results_csv)

    features = list(df_res.index)
    adjusted_pvals = list(df_res["qval"])

    dedup = pd.DataFrame({"feature": features, "adjusted_pval": adjusted_pvals})
    dedup = dedup.drop_duplicates(subset="feature", keep="first")

    save_rds(list(dedup["feature"]), list(dedup["adjusted_pval"]), rds_file)

    pd.DataFrame({"angle": [angle], "elapsed_sec": [elapsed]}).to_csv(
        runtime_file, index=False
    )

    print(f"Saved {rds_file} ({elapsed:.1f}s)")

print("SpatialDE benchmark complete.")
