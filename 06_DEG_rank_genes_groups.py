# conda activate IntestinalFibroblasts_scRNAseq
######################################################## Import libraries ########################################################
import argparse
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

try:
    _script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _script_dir = os.path.join(os.getcwd(), "Script")

sys.path.append(_script_dir)
from Function_def import my_palette

######################################################## Parse command line arguments ########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="DEG analysis: rank_genes_groups by a grouping column, save Excel and dotplot."
    )
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        default="Data/Parasite_Clustered.h5ad",
        help="Path to input h5ad file (default: Data/Parasite_Clustered.h5ad)",
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="Results/06_DEG_rank_genes_groups",
        help="Output directory for Excel and plots (default: Results/06_DEG_rank_genes_groups)",
    )
    parser.add_argument(
        "--groupby",
        type=str,
        default="stage_grouped",
        help="Column in adata.obs to group by for rank_genes_groups (default: stage_grouped)",
    )
    parser.add_argument(
        "--method",
        type=str,
        default="wilcoxon",
        choices=["wilcoxon", "logreg", "t-test", "t-test_overestim_var"],
        help="Method for rank_genes_groups (default: wilcoxon)",
    )
    parser.add_argument(
        "--n_genes",
        type=int,
        default=10,
        help="Number of top genes per group to show in plot (default: 10)",
    )
    return parser.parse_args()


######################################################## Parse command line and load data ########################################################
args = parse_args()

input_path = args.input.strip() if args.input and str(args.input).strip() else None
output_dir = args.output_dir.strip() if args.output_dir and str(args.output_dir).strip() else "Results/06_DEG_rank_genes_groups"
os.makedirs(output_dir, exist_ok=True)
print(f"Output directory: {output_dir}")

if not input_path or not os.path.isfile(input_path):
    raise FileNotFoundError(f"Input h5ad not found: {input_path}")

adata = sc.read_h5ad(input_path)
print(f"Loaded AnnData: {adata.shape[0]} cells, {adata.shape[1]} genes")
if args.groupby not in adata.obs.columns:
    raise ValueError(f"Column '{args.groupby}' not found in adata.obs. Available: {list(adata.obs.columns)}")
print(f"Grouping by: {args.groupby}")
print(adata.obs[args.groupby].value_counts())

# Reorder groupby column by known category sets: ontogeny or infection
actual = set(adata.obs[args.groupby].astype(str).dropna().unique())
expected_ontogeny = {'after_weaning', 'embryonic', 'shortly_after_birth', 'uninfected_adult'}
order_ontogeny = ['embryonic', 'shortly_after_birth', 'after_weaning', 'uninfected_adult']
expected_infection = {'uninfected_adult', 'early', 'peak', 'repair'}
order_infection = ['uninfected_adult', 'early', 'peak', 'repair']

if actual == expected_ontogeny:
    cat_order = order_ontogeny
elif actual == expected_infection:
    cat_order = order_infection
else:
    cat_order = None

if cat_order is not None:
    adata.obs[args.groupby] = pd.Categorical(
        adata.obs[args.groupby].astype(str),
        categories=cat_order,
        ordered=True,
    )
    print(f"Reordered {args.groupby} to: {adata.obs[args.groupby].cat.categories}")
else:
    print(f"{args.groupby} values are {sorted(actual)}; no predefined order, skipping reorder.")

# Remove groups with too few cells (rank_genes_groups needs >= 2 samples per group)
min_cells_per_group = 2
counts = adata.obs[args.groupby].value_counts()
valid_groups = counts[counts >= min_cells_per_group].index.tolist()
dropped = counts[counts < min_cells_per_group]
if len(dropped) > 0:
    dropped_names = dropped.index.tolist()
    n_dropped_cells = dropped.sum()
    adata = adata[adata.obs[args.groupby].isin(valid_groups)].copy()
    if hasattr(adata.obs[args.groupby], "cat"):
        adata.obs[args.groupby] = adata.obs[args.groupby].cat.remove_unused_categories()
    print(f"Dropped {len(dropped)} group(s) with <{min_cells_per_group} cells: {dropped_names} (n_cells={int(n_dropped_cells)}). Remaining: {adata.shape[0]} cells, {len(valid_groups)} groups.")
else:
    print(f"All groups have >={min_cells_per_group} cells; no filtering.")

######################################################## Main process ########################################################
# Plot order: when groupby is "identity", sort by cell count (high to low); else use categorical/group order
if args.groupby == "identity":
    counts_per_group = adata.obs[args.groupby].value_counts()
    plot_order = counts_per_group.index.tolist()  # already descending by count
    # adata.obs[args.groupby] = pd.Categorical(
    #     adata.obs[args.groupby].astype(str),
    #     categories=plot_order,
    #     ordered=True,
    # )
    adata.obs[args.groupby] = adata.obs[args.groupby].cat.reorder_categories(
    plot_order, ordered=True
    )
    print(f"Plot order for {args.groupby} (by cell count, high to low): {plot_order}")
else:
    if hasattr(adata.obs[args.groupby], "cat"):
        plot_order = list(adata.obs[args.groupby].cat.categories)
    else:
        plot_order = list(group_names)
    print(f"Plot order for {args.groupby}: {plot_order}")

# Fix group (e.g. identity) colors to my_palette for heatmap/dotplot/violin
n_cats = len(plot_order)
adata.uns[args.groupby + "_colors"] = [my_palette[i % len(my_palette)] for i in range(n_cats)]

# run rank_genes_groups and export result table with one sheet per group
sc.tl.rank_genes_groups(adata, args.groupby, method=args.method, use_raw = False, pts=True)
group_names = adata.uns["rank_genes_groups"]["names"].dtype.names

excel_path = os.path.join(output_dir, "rank_genes_groups.xlsx")
with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
    for name in group_names:
        df = sc.get.rank_genes_groups_df(adata, group=name, key='rank_genes_groups')
        sheet_name = str(name)[:31].replace("/", "-").replace("\\", "-").replace("*", "-").replace("?", "-").replace(":", "-").replace("[", "").replace("]", "")
        df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"Saved rank_genes_groups table (one sheet per group) to: {excel_path}")

sc.pl.rank_genes_groups(adata, n_genes=args.n_genes, sharey=False)
fig_path = os.path.join(output_dir, "rank_genes_groups.pdf")
plt.savefig(fig_path, dpi=300, bbox_inches="tight", format="pdf")
plt.close()
print(f"Saved figure to: {fig_path}")

# dotplot (use plot_order so groups follow categorical/defined order)
sc.pl.rank_genes_groups_dotplot(adata, n_genes=args.n_genes, categories_order=plot_order,dendrogram=False)
fig_path = os.path.join(output_dir, "rank_genes_groups_dotplot.pdf")
plt.savefig(fig_path, dpi=300, bbox_inches="tight", format="pdf")
plt.close()
print(f"Saved dotplot to: {fig_path}")

# stacked violin plot (use plot_order so groups follow categorical/defined order)
sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=args.n_genes, categories_order=plot_order,dendrogram=False)
fig_path = os.path.join(output_dir, "rank_genes_groups_stacked_violin.pdf")
plt.savefig(fig_path, dpi=300, bbox_inches="tight", format="pdf")
plt.close()
print(f"Saved stacked violin plot to: {fig_path}")

# heatmap (no row clustering; use plot_order; scale figsize so many groups get readable y-labels)
n_groups = len(plot_order)
heatmap_figsize = (
    max(12, n_groups * 0.55),
    max(10, n_groups * 0.65),
)
sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=args.n_genes,
    dendrogram=False,
    figsize=heatmap_figsize,
)
fig_path = os.path.join(output_dir, "rank_genes_groups_heatmap.pdf")
plt.savefig(fig_path, dpi=300, bbox_inches="tight", format="pdf")
plt.close()
print(f"Saved heatmap to: {fig_path}")

# violin plots: one PDF per groupby category; 2 columns, rows = ceil(n_genes/2)
rank_names = adata.uns["rank_genes_groups"]["names"]

def _sanitize_sheet_name(s):
    return str(s)[:31].replace("/", "-").replace("\\", "-").replace("*", "-").replace("?", "-").replace(":", "-").replace("[", "").replace("]", "")

violin_ncols = 2
violin_dir = os.path.join(output_dir, f"violin_by_{args.groupby}")
os.makedirs(violin_dir, exist_ok=True)
for group in group_names:
    genes_for_group = [g for g in rank_names[group][: args.n_genes].tolist() if g in adata.var_names]
    if not genes_for_group:
        continue
    n_plots = len(genes_for_group)
    nrows = (n_plots + violin_ncols - 1) // violin_ncols
    # Extra row height and bottom margin so rotated x-axis labels are not cut off
    fig, axes = plt.subplots(nrows, violin_ncols, figsize=(6 * violin_ncols, 5 * nrows))
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(bottom=0.18, left=0.08, top=0.92, right=0.96, hspace=0.35, wspace=0.25)
    for idx, gene in enumerate(genes_for_group):
        ax = axes.flat[idx]
        sc.pl.violin(
            adata,
            keys=gene,
            groupby=args.groupby,
            order=plot_order,
            stripplot=True,
            jitter=0.4,
            rotation=90,
            ylabel=f"{gene} Expression level",
            ax=ax,
            show=False,
            inner="box",
        )
    for idx in range(n_plots, axes.size):
        axes.flat[idx].set_visible(False)
    fig_path = os.path.join(violin_dir, f"{_sanitize_sheet_name(group)}.pdf")
    plt.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.2, format="pdf")
    plt.close()
    print(f"Saved violin figure for group '{group}' ({len(genes_for_group)} genes, {nrows}x{violin_ncols}) to: {fig_path}")

print(f"Saved {len(group_names)} violin figures (one per {args.groupby} category) to: {violin_dir}")

# UMAP expression plots: one PDF per groupby category (same layout as violin: 2 cols, one gene per panel)
umap_key = "X_umap" if "X_umap" in adata.obsm else ("umap" if "umap" in adata.obsm else None)
if umap_key is not None:
    umap_dir = os.path.join(output_dir, f"umap_by_{args.groupby}")
    os.makedirs(umap_dir, exist_ok=True)
    umap_ncols = 2
    for group in group_names:
        genes_for_group = [g for g in rank_names[group][: args.n_genes].tolist() if g in adata.var_names]
        if not genes_for_group:
            continue
        n_plots = len(genes_for_group)
        nrows = (n_plots + umap_ncols - 1) // umap_ncols
        fig, axes = plt.subplots(nrows, umap_ncols, figsize=(5 * umap_ncols, 4 * nrows))
        axes = np.atleast_2d(axes)
        fig.subplots_adjust(bottom=0.08, left=0.08, top=0.92, right=0.96, hspace=0.3, wspace=0.25)
        for idx, gene in enumerate(genes_for_group):
            ax = axes.flat[idx]
            sc.pl.umap(adata, color=gene, ax=ax, show=False, title=f"{gene}")
        for idx in range(n_plots, axes.size):
            axes.flat[idx].set_visible(False)
        fig_path = os.path.join(umap_dir, f"{_sanitize_sheet_name(group)}.pdf")
        plt.savefig(fig_path, dpi=300, bbox_inches="tight", pad_inches=0.2, format="pdf")
        plt.close()
        print(f"Saved UMAP figure for group '{group}' ({len(genes_for_group)} genes) to: {fig_path}")
    print(f"Saved {len(group_names)} UMAP figures (one per {args.groupby} category) to: {umap_dir}")
else:
    print("Skip UMAP plots: no 'X_umap' or 'umap' in adata.obsm.")