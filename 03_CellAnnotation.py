# conda activate IntestinalFibroblasts_scRNAseq
######################################################## Import libraries ########################################################
import argparse
import importlib.util
import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

######################################################## Functions ########################################################
from Function_def import plot_umap_by_stage
from Function_def import my_palette
from Function_def import convert_barcode_format
from Function_def import plot_cell_proportion
######################################################## Parse command line arguments ########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Cell annotation: rename cell types from Excel, UMAP visualization, marker dotplot, correlation and proportion plots."
    )
    parser.add_argument(
        "--rename_path",
        type=str,
        default="Data/Cluster_infectiondataset_renamed.xlsx",
        help="Path to Excel file with Barcode, Old_name, New_name (default: Data/Cluster_infectiondataset_renamed.xlsx)",
    )
    parser.add_argument(
        "--parasite_h5ad",
        type=str,
        default="Data/Parasite_Clustered.h5ad",
        help="Path to Parasite h5ad file (default: Data/Parasite_Clustered.h5ad)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="Results/03_CellAnnotation",
        help="Output directory for plots and CSV files (default: Results/03_CellAnnotation)",
    )
    return parser.parse_args()


######################################################## Parse command line and load data ########################################################
# Only assign when argument is non-empty
args = parse_args()

if args.rename_path and str(args.rename_path).strip():
    rename_path = args.rename_path.strip()
    Cluster_infectiondataset_renamed = pd.read_excel(rename_path, sheet_name="active_cluster")
    Cluster_infectiondataset_renamed.columns = ["Barcode", "Old_name", "New_name"]
    Cluster_infectiondataset_renamed = Cluster_infectiondataset_renamed.iloc[1:]
    barcode_to_new_name = dict(zip(Cluster_infectiondataset_renamed["Barcode"], Cluster_infectiondataset_renamed["New_name"]))
else:
    rename_path = None
    Cluster_infectiondataset_renamed = None
    barcode_to_new_name = {}

if args.parasite_h5ad and str(args.parasite_h5ad).strip():
    Parasite = sc.read_h5ad(args.parasite_h5ad.strip())
    if barcode_to_new_name:
        Parasite.obs["New_name"] = Parasite.obs["Unnamed: 0"].map(barcode_to_new_name)
    else:
        Parasite.obs["New_name"] = None
else:
    Parasite = None

if args.output_dir and str(args.output_dir).strip():
    output_dir = args.output_dir.strip()
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory created: {output_dir}")
else:
    output_dir = "Results/03_CellAnnotation"
    os.makedirs(output_dir, exist_ok=True)
    print(f"\nOutput directory created: {output_dir}")

if Parasite is not None:
    available_identities = Parasite.obs["New_name"].dropna().unique()
    print(f"\nNew_name categories: {len(available_identities)}")
    print(f"New_name value counts:")
    print(Parasite.obs["New_name"].value_counts())

######################################################## Main process ########################################################
######################################################## UMAP Visualization for New_name ########################################################

stage_order = ['uninfected_adult', 'early', 'peak', 'repair']
new_name_counts = Parasite.obs['New_name'].value_counts()
new_name_order = [x for x in new_name_counts.index.tolist() if pd.notna(x)]
print(f"\nNew_name order (by cell count):")
for name in new_name_order:
    print(f"  {name}: {new_name_counts[name]} cells")

fig_facet, figs_individual = plot_umap_by_stage(
    Parasite,
    output_dir,
    stage_group_col="stage_grouped",
    identity_col="New_name",
    stage_order=stage_order,
    identity_order=new_name_order,
    color_palette=my_palette,
    figure_prefix="umap_faceted_by_stage_grouped_New_name",
)
plt.close(fig_facet)
for f in figs_individual:
    plt.close(f)
print(f"\nPlots saved to: {output_dir}")
print(f"  - umap_faceted_by_stage_grouped_New_name.pdf / .png")
print(f"  - umap_<stage>_New_name.pdf per stage")

# Export CSV files for Loupe Browser
cell_names = Parasite.obs['CellName'].values
converted_barcodes = [convert_barcode_format(barcode) for barcode in cell_names]

if 'X_umap' in Parasite.obsm:
    umap_df = pd.DataFrame({
        'Barcode': converted_barcodes,
        'UMAP_1': Parasite.obsm['X_umap'][:, 0],
        'UMAP_2': Parasite.obsm['X_umap'][:, 1],
    })
    umap_output_path = os.path.join(output_dir, 'loupe_browser_umap.csv')
    umap_df.to_csv(umap_output_path, index=False)
    print(f"\n1. UMAP CSV file saved to: {umap_output_path}")
    print(f"   Total cells: {len(umap_df)}")

new_name_df = pd.DataFrame({
    'Barcode': converted_barcodes,
    'New_name': Parasite.obs['New_name'].values,
})
new_name_output_path = os.path.join(output_dir, 'loupe_browser_New_name.csv')
new_name_df.to_csv(new_name_output_path, index=False)
print(f"\n2. New_name CSV file saved to: {new_name_output_path}")
print(f"   Total cells: {len(new_name_df)}")

stage_df = pd.DataFrame({
    'Barcode': converted_barcodes,
    'stage_grouped': Parasite.obs['stage_grouped'].values,
})
stage_output_path = os.path.join(output_dir, 'loupe_browser_stage_grouped.csv')
stage_df.to_csv(stage_output_path, index=False)
print(f"\n3. Stage_grouped CSV file saved to: {stage_output_path}")
print(f"   Total cells: {len(stage_df)}")


######################################################## Marker gene dotplot ########################################################

genelist =['Pdpn', 'Pdgfra', 'Igfbp3', 'Fgfr2', 'Pdgfrb', 'Des', 'Acta2', 'Myh11', 'Actg2', 'Cd81', 'Ackr4', 'Cd55', 'Cd34', 'Jam2', 'Ccl11', 'Kit', 'Hand1']

ax = sc.pl.dotplot(Parasite, genelist, groupby='New_name', dendrogram=True, dot_max=1, dot_min=0, standard_scale='var', return_fig=True, use_raw=False)
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Marker_Parasite_dotplot.pdf'))

# Create stacked_violin plot
sc.pl.stacked_violin(Parasite, genelist, groupby='New_name', dendrogram=True, standard_scale='var', show=False, use_raw=False)
# Get the current figure and save it
fig = plt.gcf()
fig.tight_layout()
fig.savefig(os.path.join(output_dir, 'Marker_Parasite_stacked_violin.pdf'), bbox_inches='tight', pad_inches=0.1)
plt.close(fig)

######################################################## Correlation plot########################################################
# /Users/xliu2942/Documents/Projects/IntestinalFibroblasts_scRNAseq/IntestinalFibroblasts_scRNAseq-main/NEW_ParasiteVisualizations.ipynb
sc.tl.dendrogram(Parasite, 'New_name', n_pcs=30)

# Temporarily set rcParams to remove gridlines
with plt.rc_context({'axes.grid': False}):
    sc.pl.correlation_matrix(Parasite, 'New_name', show=False)
    # Get the current figure and adjust layout
    fig = plt.gcf()
    fig.tight_layout()
    # Save the plot as a PDF with tight bounding box to ensure edges are visible
    fig.savefig(os.path.join(output_dir, 'correlation_matrix.pdf'), 
                bbox_inches='tight', 
                pad_inches=0.1,
                dpi=300,
                format='pdf')
    plt.close(fig)

######################################################## cell proportion plot ########################################################
stage_order_proportion = ["uninfected_adult", "early", "peak", "repair"]
stage_colors_proportion = {
    "uninfected_adult": "#A9A9A9",
    "early": "#FFA07A",
    "peak": "#FF6347",
    "repair": "#77DD77",
}
outpath, csv_path = plot_cell_proportion(
    Parasite,
    output_dir,
    stage_order=stage_order_proportion,
    stage_colors=stage_colors_proportion,
    output_filename="representations_parasite.pdf",
    identity_col="New_name",
    stage_col="stage_grouped",
)
print(f"\nCell proportion plot saved to: {outpath}")
if csv_path:
    print(f"Cell proportion table saved to: {csv_path}")