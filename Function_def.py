######################################################## Import libraries ########################################################
import os
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
######################################################## Functions ########################################################
# Custom color palette for Ontogeny
my_palette = [
    "#66C2A5", "#E78AC3", "#ffc773",
    "#1B9E77", "#7570B3", "#FEE090",
    "#A6D854", "#8DA0CB", "#FC8D62",
    "#8FBC94", "#b0a4e3", "#ffa631",
    "#0aa344", "#e4c6d0", "#ffa400",
    "#519a73", "#4b5cc4", "#eedeb0",
    "#549688", "#ffb3a7", "#b35c44",
    "#7fecad", "#a1afc9", "#a78e44",
    "#519a73", "#2e4e7e", "#955539"
]

# Export CSV file for Loupe Browser
# Convert barcode format from AAACCCAAGTATCTGC-1_1 to AAACCCAAGTATCTGC-1
# or from TTTGTTGTCTGGTGCG-1_8 to TTTGTTGTCTGGTGCG-8
def convert_barcode_format(barcode):
    """
    Convert barcode format for Loupe Browser.
    Examples:
    - AAACCCAAGTATCTGC-1_1 -> AAACCCAAGTATCTGC-1
    - TTTGTTGTCTGGTGCG-1_8 -> TTTGTTGTCTGGTGCG-8
    """
    if '_' in barcode:
        # Split by underscore
        parts = barcode.split('_')
        if len(parts) == 2:
            # Get the part before underscore and the part after
            prefix = parts[0]  # e.g., AAACCCAAGTATCTGC-1
            suffix = parts[1]  # e.g., 1 or 8
            # Replace the number after the last dash with the suffix
            if '-' in prefix:
                prefix_parts = prefix.rsplit('-', 1)
                new_barcode = f"{prefix_parts[0]}-{suffix}"
                return new_barcode
    return barcode


def plot_umap_by_stage(
    adata,
    output_dir,
    stage_group_col="stage_grouped",
    identity_col="identity",
    stage_order=None,
    identity_order=None,
    color_palette=None,
    figure_prefix="umap_faceted_by_stage_grouped_identity",
    fig_width_per_panel=6,
    fig_height=6,
    single_plot_width=10,
    point_size=16,
    legend_marker_size=4,
):
    """
    Create UMAP plots: (1) one faceted figure (all time points in one row),
    (2) separate figures per time point. Saves PDF/PNG for faceted and PDF per stage.
    Returns (faceted_fig, list_of_individual_figs). Caller should close figures when done.
    """
    umap_coords = adata.obsm["X_umap"]
    available_groups = adata.obs[stage_group_col].unique()
    stage_order = stage_order or []
    stage_groups = [g for g in stage_order if g in available_groups]
    if not stage_groups:
        stage_groups = sorted(available_groups)
    n_groups = len(stage_groups)

    if identity_order is not None:
        available_identities = adata.obs[identity_col].unique()
        identity_values = [x for x in identity_order if x in available_identities]
    else:
        identity_counts = adata.obs[identity_col].value_counts()
        identity_values = identity_counts.index.tolist()
    identity_counts = adata.obs[identity_col].value_counts()

    if color_palette is None:
        tableau10 = plt.cm.tab10
        identity_colors = {identity: tableau10(i % 10) for i, identity in enumerate(identity_values)}
    else:
        identity_colors = {
            identity: color_palette[i % len(color_palette)]
            for i, identity in enumerate(identity_values)
        }

    os.makedirs(output_dir, exist_ok=True)
    ncols = n_groups
    nrows = 1

    # ----- Faceted figure -----
    fig_facet, axes = plt.subplots(nrows, ncols, figsize=(fig_width_per_panel * ncols, fig_height))
    if n_groups == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    handles_dict = {}
    for idx, stage in enumerate(stage_groups):
        ax = axes[idx]
        mask = adata.obs[stage_group_col] == stage
        for identity in identity_values:
            identity_mask = (adata.obs[identity_col] == identity) & mask
            if identity_mask.sum() > 0:
                scatter = ax.scatter(
                    umap_coords[identity_mask, 0],
                    umap_coords[identity_mask, 1],
                    c=[identity_colors[identity]],
                    label=identity if idx == 0 else "",
                    s=point_size,
                    alpha=1,
                    rasterized=True,
                )
                if idx == 0 and identity not in handles_dict:
                    handles_dict[identity] = scatter
        ax.set_title(f"{stage}", fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=10)
        ax.set_ylabel("UMAP 2", fontsize=10)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    labels_list = [f"{identity} (n={identity_counts.get(identity, 0)})" for identity in identity_values]
    use_line2d = color_palette is not None
    if use_line2d:
        handles_list = [
            Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=identity_colors[identity],
                markersize=legend_marker_size,
                markeredgewidth=0,
                label=identity,
            )
            for identity in identity_values
        ]
    else:
        handles_list = [handles_dict[identity] for identity in identity_values]

    axes[-1].legend(
        handles_list, labels_list, loc="center left", bbox_to_anchor=(1, 0.5),
        fontsize=10, frameon=False, markerscale=2, ncol=1,
    )
    plt.tight_layout(rect=[0, 0, 0.95, 1])

    facet_pdf = os.path.join(output_dir, f"{figure_prefix}.pdf")
    facet_png = os.path.join(output_dir, f"{figure_prefix}.png")
    fig_facet.savefig(facet_pdf, dpi=300, bbox_inches="tight", format="pdf")
    fig_facet.savefig(facet_png, dpi=300, bbox_inches="tight", format="png")

    # ----- Separate figures by time point -----
    individual_figs = []
    for stage in stage_groups:
        fig_single, ax_single = plt.subplots(1, 1, figsize=(single_plot_width, fig_height))
        mask = adata.obs[stage_group_col] == stage
        stage_identity_counts = adata.obs[mask][identity_col].value_counts()
        for identity in identity_values:
            identity_mask = (adata.obs[identity_col] == identity) & mask
            if identity_mask.sum() > 0:
                ax_single.scatter(
                    umap_coords[identity_mask, 0],
                    umap_coords[identity_mask, 1],
                    c=[identity_colors[identity]],
                    label=identity,
                    s=point_size,
                    alpha=1,
                    rasterized=True,
                )
        ax_single.set_title(f"{stage}", fontsize=14, fontweight="bold")
        ax_single.set_xlabel("UMAP 1", fontsize=12)
        ax_single.set_ylabel("UMAP 2", fontsize=12)
        ax_single.grid(False)
        ax_single.spines["top"].set_visible(False)
        ax_single.spines["right"].set_visible(False)
        # Legend: only include identities with n > 0 in this stage (same colors as faceted figure)
        identities_present = [identity for identity in identity_values if stage_identity_counts.get(identity, 0) > 0]
        stage_handles = [
            Line2D(
                [0], [0], marker="o", color="w",
                markerfacecolor=identity_colors[identity],
                markersize=legend_marker_size,
                markeredgewidth=0,
                label=identity,
            )
            for identity in identities_present
        ]
        stage_labels = [f"{identity} (n={stage_identity_counts[identity]})" for identity in identities_present]
        ax_single.legend(
            stage_handles, stage_labels, loc="center left", bbox_to_anchor=(1, 0.5),
            fontsize=10, frameon=False, markerscale=2,
        )
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        stage_filename = stage.replace(" ", "_").replace("/", "_")
        stage_pdf = os.path.join(output_dir, f"umap_{stage_filename}_{identity_col}.pdf")
        fig_single.savefig(stage_pdf, dpi=300, bbox_inches="tight", format="pdf")
        individual_figs.append(fig_single)

    return fig_facet, individual_figs



def plot_score(
    adata,
    color_values,
    output_dir,
    color_label,
    stage_order=None,
    figure_prefix="umap_faceted_by_stage_grouped_score",
    fig_width_per_panel=6,
    fig_height=6,
    single_plot_width=10,
    cmap="viridis",
    point_size=16,
):
    """
    Create UMAP plots colored by continuous values (e.g. score or expression), faceted by stage_grouped.
    Same format as plot_umap_by_stage: (1) one faceted figure, (2) separate figures per time point.
    Uses viridis colormap. Returns (faceted_fig, list_of_individual_figs). Caller should close figures when done.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with UMAP coordinates in obsm['X_umap'] and 'stage_grouped' in obs
    color_values : array-like
        Array of values to color cells by (length should match adata.n_obs)
    output_dir : str
        Output directory for PDF/PNG files
    color_label : str
        Label for the colorbar (and used for per-stage filenames)
    stage_order : list, optional
        Order of stage_grouped categories. If None, uses all available groups sorted.
    figure_prefix : str
        Basename for faceted files (default: umap_faceted_by_stage_grouped_score)
    fig_width_per_panel : float
        Width per panel in inches for faceted figure
    fig_height : float
        Figure height in inches
    single_plot_width : float
        Width for each individual stage figure
    cmap : str
        Colormap for continuous values (default: viridis)
    point_size : float
        Size of scatter points
    """
    umap_coords = adata.obsm["X_umap"]
    stage_col = "stage_grouped"

    if stage_order is None:
        stage_order = sorted(adata.obs[stage_col].unique())
    available_groups = adata.obs[stage_col].unique()
    stage_groups = [g for g in stage_order if g in available_groups]
    if not stage_groups:
        stage_groups = sorted(available_groups)
    n_groups = len(stage_groups)

    if n_groups == 0:
        print(f"Warning: No matching stage_grouped groups found. Available: {available_groups}")
        return None, []

    os.makedirs(output_dir, exist_ok=True)
    ncols = n_groups
    nrows = 1

    # ----- Faceted figure -----
    fig_facet, axes = plt.subplots(nrows, ncols, figsize=(fig_width_per_panel * ncols, fig_height))
    if n_groups == 1:
        axes = [axes]
    else:
        axes = axes.flatten()

    for idx, stage in enumerate(stage_groups):
        ax = axes[idx]
        mask = adata.obs[stage_col] == stage
        scatter = ax.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=color_values[mask],
            s=point_size,
            alpha=1,
            rasterized=True,
            cmap=cmap,
        )
        ax.set_title(f"{stage}", fontsize=12, fontweight="bold")
        ax.set_xlabel("UMAP 1", fontsize=10)
        ax.set_ylabel("UMAP 2", fontsize=10)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(color_label, fontsize=9)

    plt.tight_layout(rect=[0, 0, 1, 1])
    facet_pdf = os.path.join(output_dir, f"{figure_prefix}.pdf")
    facet_png = os.path.join(output_dir, f"{figure_prefix}.png")
    fig_facet.savefig(facet_pdf, dpi=300, bbox_inches="tight", format="pdf")
    fig_facet.savefig(facet_png, dpi=300, bbox_inches="tight", format="png")

    # ----- Separate figures by time point -----
    color_label_safe = color_label.replace(" ", "_")
    individual_figs = []
    for stage in stage_groups:
        fig_single, ax_single = plt.subplots(1, 1, figsize=(single_plot_width, fig_height))
        mask = adata.obs[stage_col] == stage
        scatter = ax_single.scatter(
            umap_coords[mask, 0],
            umap_coords[mask, 1],
            c=color_values[mask],
            s=point_size,
            alpha=1,
            rasterized=True,
            cmap=cmap,
        )
        ax_single.set_title(f"{stage}", fontsize=14, fontweight="bold")
        ax_single.set_xlabel("UMAP 1", fontsize=12)
        ax_single.set_ylabel("UMAP 2", fontsize=12)
        ax_single.grid(False)
        ax_single.spines["top"].set_visible(False)
        ax_single.spines["right"].set_visible(False)
        cbar = plt.colorbar(scatter, ax=ax_single)
        cbar.set_label(color_label, fontsize=10)
        plt.tight_layout(rect=[0, 0, 1, 1])
        stage_filename = stage.replace(" ", "_").replace("/", "_")
        stage_pdf = os.path.join(output_dir, f"umap_{stage_filename}_{color_label_safe}.pdf")
        fig_single.savefig(stage_pdf, dpi=300, bbox_inches="tight", format="pdf")
        individual_figs.append(fig_single)

    # ----- Export two CSVs for Loupe Browser: UMAP coordinates + score/expression -----
    cell_names = adata.obs["CellName"].values if "CellName" in adata.obs.columns else adata.obs_names
    barcodes = [convert_barcode_format(b) for b in cell_names]

    umap_csv = os.path.join(output_dir, "loupe_browser_umap.csv")
    umap_df = pd.DataFrame({
        "Barcode": barcodes,
        "UMAP_1": umap_coords[:, 0],
        "UMAP_2": umap_coords[:, 1],
    })
    umap_df.to_csv(umap_csv, index=False)
    print(f"Saved Loupe Browser UMAP CSV to: {umap_csv}")

    score_csv = os.path.join(output_dir, f"loupe_browser_{color_label_safe}.csv")
    color_values_flat = np.asarray(color_values).ravel()
    score_df = pd.DataFrame({
        "Barcode": barcodes,
        color_label: color_values_flat,
    })
    score_df.to_csv(score_csv, index=False)
    print(f"Saved Loupe Browser score/expression CSV to: {score_csv}")

    print(f"Saved faceted plot to: {facet_pdf}")
    return fig_facet, individual_figs


def convert_barcode_format(barcode):
    """
    Convert barcode format for Loupe Browser.
    Examples:
    - AAACCCAAGTATCTGC-1_1 -> AAACCCAAGTATCTGC-1
    - TTTGTTGTCTGGTGCG-1_8 -> TTTGTTGTCTGGTGCG-8
    """
    if "_" in str(barcode):
        parts = str(barcode).split("_", 1)
        if len(parts) == 2:
            prefix, suffix = parts[0], parts[1]
            if "-" in prefix:
                return f"{prefix.rsplit('-', 1)[0]}-{suffix}"
    return str(barcode)


def plot_cell_proportion(
    adata,
    output_dir,
    stage_order,
    stage_colors,
    output_filename="representations_parasite.pdf",
    identity_col="New_name",
    stage_col="stage_grouped",
    figsize=(12, 5),
    dpi=300,
    export_csv=True,
    csv_filename=None,
):
    """
    Stacked horizontal bar plot: % of Cells and Number of Cells per identity (rows) and stage (colors).
    Saves PDF to output_dir. Optionally exports a CSV with n_cells and pct_cells per (identity, stage).
    Rows sorted by total cell count (highest at top).
    """
    df = adata.obs[[identity_col, stage_col]].copy()
    df = df[df[identity_col].notna()]
    counts_df = df.groupby([identity_col, stage_col]).size().reset_index(name="n_cells")
    new_name_totals = df.groupby(identity_col).size().reset_index(name="total_cells_per_new_name")
    counts_df = counts_df.merge(new_name_totals, on=identity_col)
    counts_df["pct_cells"] = (counts_df["n_cells"] / counts_df["total_cells_per_new_name"]) * 100
    representation_long = pd.melt(
        counts_df,
        id_vars=[identity_col, stage_col],
        value_vars=["pct_cells", "n_cells"],
        var_name="measure",
        value_name="value",
    )
    representation_long["measure"] = representation_long["measure"].replace(
        {"pct_cells": "% of Cells", "n_cells": "Number of Cells"}
    )
    representation_long[stage_col] = pd.Categorical(
        representation_long[stage_col], categories=stage_order, ordered=True
    )
    n_cells_data = representation_long[representation_long["measure"] == "Number of Cells"]
    total_cells_per_identity = n_cells_data.groupby(identity_col)["value"].sum().sort_values(ascending=False)
    identity_order = total_cells_per_identity.index.tolist()
    representation_long[identity_col] = pd.Categorical(
        representation_long[identity_col], categories=identity_order, ordered=True
    )
    sns.set_style("white")
    fig, axes = plt.subplots(1, 2, figsize=figsize, sharey=True)
    measures = ["% of Cells", "Number of Cells"]
    for idx, measure in enumerate(measures):
        ax = axes[idx]
        measure_data = representation_long[representation_long["measure"] == measure]
        pivot_data = measure_data.pivot_table(
            index=identity_col, columns=stage_col, values="value", fill_value=0
        )
        pivot_data = pivot_data[[s for s in stage_order if s in pivot_data.columns]]
        pivot_data = pivot_data.reindex(identity_order[::-1])
        left = np.zeros(len(pivot_data))
        for stage in stage_order:
            if stage in pivot_data.columns:
                ax.barh(
                    range(len(pivot_data)),
                    pivot_data[stage],
                    left=left,
                    label=stage,
                    color=stage_colors.get(stage, "#808080"),
                    height=1.0,
                    edgecolor="none",
                )
                left += pivot_data[stage]
        ax.set_yticks(range(len(pivot_data)))
        ax.set_yticklabels(pivot_data.index)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("black")
        ax.spines["bottom"].set_color("black")
        ax.grid(False)
        ax.set_title(measure, size=15, pad=10)
        ax.set_xlim(left=0)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=stage_colors[stage], label=stage)
        for stage in stage_order if stage in stage_colors
    ]
    fig.legend(
        handles=handles,
        labels=[h.get_label() for h in handles],
        title="Infection Stage",
        loc="lower center",
        ncol=len(stage_order),
        bbox_to_anchor=(0.5, -0.05),
        frameon=False,
        fontsize=12,
        title_fontsize=12,
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1])
    os.makedirs(output_dir, exist_ok=True)
    outpath = os.path.join(output_dir, output_filename)
    fig.savefig(outpath, dpi=dpi, bbox_inches="tight", format="pdf")
    plt.close(fig)

    csv_path = None
    if export_csv:
        if csv_filename is None:
            csv_filename = os.path.splitext(output_filename)[0] + "_proportion.csv"
        csv_path = os.path.join(output_dir, csv_filename)
        export_df = counts_df[[identity_col, stage_col, "n_cells", "pct_cells"]].copy()
        export_df = export_df.sort_values([identity_col, stage_col])
        export_df.to_csv(csv_path, index=False)
    return outpath, csv_path

