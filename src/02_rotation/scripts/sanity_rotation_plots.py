# ==============================================================================
# sanity_rotation_plots.py
# Standalone visualization of rotated tissue slices. For each (slice, mode),
# plots the tissue at the 4 rotation angles (0, 30, 45, 60) side-by-side,
# colored by total counts, so rotation consistency can be visually verified.
#
# Optional: set ALL_SLICES_FIGURE = True to also generate one combined figure
# with all slices x all angles.
#
# Edit SLICES and MODES below to control scope.
# Output: src/02_rotation/outputs/{slice}/{mode}/anndata/figures/sanity_rotation.png
# ==============================================================================

import os
import matplotlib
matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

# --- Select which slices and modes to run ---
SLICES = ["anterior1", "anterior2", "posterior1", "posterior2"]
MODES = ["simulated", "whole"]

ANGLES = [0, 30, 45, 60]

# Toggle to also produce an all-slices x all-angles combined figure
ALL_SLICES_FIGURE = True

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROTATION_DIR = os.path.dirname(SCRIPT_DIR)
PROJECT_ROOT = os.path.dirname(os.path.dirname(ROTATION_DIR))

ANNDATA_TEMPLATE = os.path.join(
    ROTATION_DIR, "outputs", "{slice}", "{mode}", "anndata",
    "data", "scdesign3_angle{angle}.h5ad"
)
FIGURES_DIR_TEMPLATE = os.path.join(
    ROTATION_DIR, "outputs", "{slice}", "{mode}", "anndata", "figures"
)


def plot_angle_grid(adatas_by_angle, slice_name, mode, out_dir):
    """Plot one (slice, mode) as a 1x4 grid of angles, colored by total_counts."""
    fig, axes = plt.subplots(1, len(ANGLES), figsize=(5 * len(ANGLES), 5))
    if len(ANGLES) == 1:
        axes = [axes]

    for ax, angle in zip(axes, ANGLES):
        if angle not in adatas_by_angle:
            ax.set_title(f"{angle}° (missing)")
            ax.axis("off")
            continue
        adata = adatas_by_angle[angle]
        coords = adata.obsm["spatial"]
        # Use total_counts if available (added by convert_to_anndata.py)
        if "total_counts" in adata.obs.columns:
            color = adata.obs["total_counts"].values
        else:
            color = np.log1p(np.asarray(adata.X.sum(axis=1)).flatten())
        sc_plot = ax.scatter(
            coords[:, 0], coords[:, 1], c=color, s=4,
            cmap="magma", edgecolors="none"
        )
        ax.set_title(f"{angle}°")
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(
        f"Rotation sanity: {slice_name} / {mode}", fontsize=14, y=1.02
    )
    fig.tight_layout()
    out_path_png = os.path.join(out_dir, "sanity_rotation.png")
    fig.savefig(out_path_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out_path_png}")


for mode in MODES:

    all_slices_adatas = {}

    for slice_name in SLICES:
        out_dir = FIGURES_DIR_TEMPLATE.format(slice=slice_name, mode=mode)
        os.makedirs(out_dir, exist_ok=True)

        # Load all angles for this (slice, mode)
        adatas_by_angle = {}
        for angle in ANGLES:
            h5ad_path = ANNDATA_TEMPLATE.format(
                slice=slice_name, mode=mode, angle=angle
            )
            if not os.path.exists(h5ad_path):
                print(f"Missing {h5ad_path} -- skipping angle {angle}")
                continue
            adatas_by_angle[angle] = sc.read_h5ad(h5ad_path)

        if not adatas_by_angle:
            print(f"No AnnData for {slice_name}/{mode} -- skipping")
            continue

        print(f"\n=== Plotting {slice_name} / {mode} ===")
        plot_angle_grid(adatas_by_angle, slice_name, mode, out_dir)
        all_slices_adatas[slice_name] = adatas_by_angle

    # Optional combined figure: slices (rows) x angles (cols)
    if ALL_SLICES_FIGURE and all_slices_adatas:
        n_slices = len(all_slices_adatas)
        fig, axes = plt.subplots(
            n_slices, len(ANGLES),
            figsize=(5 * len(ANGLES), 5 * n_slices),
            squeeze=False
        )
        for i, (slice_name, adatas_by_angle) in enumerate(
            all_slices_adatas.items()
        ):
            for j, angle in enumerate(ANGLES):
                ax = axes[i][j]
                if angle not in adatas_by_angle:
                    ax.set_title(f"{slice_name} {angle}° (missing)")
                    ax.axis("off")
                    continue
                adata = adatas_by_angle[angle]
                coords = adata.obsm["spatial"]
                if "total_counts" in adata.obs.columns:
                    color = adata.obs["total_counts"].values
                else:
                    color = np.log1p(np.asarray(adata.X.sum(axis=1)).flatten())
                ax.scatter(
                    coords[:, 0], coords[:, 1], c=color, s=4,
                    cmap="magma", edgecolors="none"
                )
                ax.set_title(f"{slice_name} {angle}°")
                ax.set_aspect("equal")
                ax.set_xticks([])
                ax.set_yticks([])

        fig.suptitle(f"All slices x angles: {mode}", fontsize=16, y=1.00)
        fig.tight_layout()
        combined_dir = os.path.join(ROTATION_DIR, "outputs", "all_slices", mode)
        os.makedirs(combined_dir, exist_ok=True)
        combined_path = os.path.join(combined_dir, "sanity_rotation_all_slices.png")
        fig.savefig(combined_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\nSaved combined figure: {combined_path}")

print("\nAll rotation sanity plots complete.")
