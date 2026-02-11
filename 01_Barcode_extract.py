# conda activate IntestinalFibroblasts_scRNAseq
######################################################## Import libraries ########################################################
import argparse
import scanpy as sc
import anndata as ad
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from matplotlib.lines import Line2D

######################################################## Functions ########################################################
from Function_def import plot_umap_by_stage, convert_barcode_format, my_palette
######################################################## Parse command line arguments ########################################################
# Parse command line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Barcode extract: UMAP visualization and Loupe Browser export for Parasite/Ontogeny h5ad Data.")
    parser.add_argument(
        "--h5ad_file",
        type=str,
        default="Data/Parasite_Clustered.h5ad",
        help="Path to input h5ad file (default: Parasite_Clustered.h5ad)",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default="Results/01_Barcode_extract/Parasite_Clustered",
        help="Output directory for plots and CSV files (default: Parasite_Clustered)",
    )
    parser.add_argument(
        "--ontogeny_h5ad",
        type=str,
        default="Data/OntogenyData_OnlyFBs_Xuan.h5ad",
        help="Path to Ontogeny h5ad file (default: OntogenyData_OnlyFBs_Xuan.h5ad)",
    )
    parser.add_argument(
        "--ontogeny_output_dir",
        type=str,
        default="Results/01_Barcode_extract/Ontogeny_Clustered",
        help="Output directory for Ontogeny plots and CSV files (default: Ontogeny_Clustered)",
    )
    return parser.parse_args()

args = parse_args()
h5ad_file = args.h5ad_file
output_dir = args.output_dir

######################################################## Parasite Dataset ########################################################
# Read h5ad file using scanpy
adata = sc.read_h5ad(h5ad_file)

print(f"Loaded AnnData object with shape: {adata.shape}")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
print(f"\nStage_grouped categories: {adata.obs['stage_grouped'].unique()}")
print(f"\nIdentity value counts:")
print(adata.obs["identity"].value_counts())

# Create output directory
os.makedirs(output_dir, exist_ok=True)
print(f"\nOutput directory created: {output_dir}")

# Parasite: faceted UMAP + separate figures by time point (tableau10 palette)
parasite_stage_order = ['uninfected_adult', 'early', 'peak', 'repair']
parasite_identity_order = [
    'PdgfraloCd55+ FB', 'Cd81+ FB', 'Pdgfrahi FB', 'Pdgfralo FB',
    'Ccl11+ MF', 'Ccl11- MF', 'Pericytes', 'PdgfrahiCd55+ FB',
    'PdgfrahiJam2+ FB', 'Hand1+'
]
fig_facet, figs_individual = plot_umap_by_stage(
    adata,
    output_dir,
    stage_group_col="stage_grouped",
    identity_col="identity",
    stage_order=parasite_stage_order,
    identity_order=parasite_identity_order,
    color_palette=None,
    figure_prefix="umap_faceted_by_stage_grouped_identity",
)
plt.close(fig_facet)
for f in figs_individual:
    plt.close(f)


# Convert barcodes
cell_names = adata.obs['CellName'].values
converted_barcodes = [convert_barcode_format(barcode) for barcode in cell_names]

# Create three separate CSV files for Loupe Browser

# 1. UMAP coordinates file
if 'X_umap' in adata.obsm:
    umap_df = pd.DataFrame({
        'Barcode': converted_barcodes,
        'UMAP_1': adata.obsm['X_umap'][:, 0],
        'UMAP_2': adata.obsm['X_umap'][:, 1],
    })
    umap_output_path = os.path.join(output_dir, 'loupe_browser_umap.csv')
    umap_df.to_csv(umap_output_path, index=False)
    print(f"\n1. UMAP CSV file saved to: {umap_output_path}")
    print(f"   Total cells: {len(umap_df)}")
else:
    print("\n1. UMAP coordinates not found in Data")

# 2. Identity file
identity_df = pd.DataFrame({
    'Barcode': converted_barcodes,
    'identity': adata.obs['identity'].values,
})
identity_output_path = os.path.join(output_dir, 'loupe_browser_identity.csv')
identity_df.to_csv(identity_output_path, index=False)
print(f"\n2. Identity CSV file saved to: {identity_output_path}")
print(f"   Total cells: {len(identity_df)}")

# 3. Stage_grouped file
stage_df = pd.DataFrame({
    'Barcode': converted_barcodes,
    'stage_grouped': adata.obs['stage_grouped'].values,
})
stage_output_path = os.path.join(output_dir, 'loupe_browser_stage_grouped.csv')
stage_df.to_csv(stage_output_path, index=False)
print(f"\n3. Stage_grouped CSV file saved to: {stage_output_path}")
print(f"   Total cells: {len(stage_df)}")

print(f"\nPlots saved to: {output_dir}")
print(f"  - umap_faceted_by_stage_grouped_identity.pdf")
print(f"  - umap_faceted_by_stage_grouped_identity.png")


# have a look test.obs['orig.ident'] and test.obs['Sample_ID'] correspondance
print(adata.obs['orig.ident'].value_counts())
print(adata.obs['Sample_ID'].value_counts())

######################################################## Ontogeny Dataset ########################################################
# Process Ontogeny dataset
Ontogeny = sc.read_h5ad(args.ontogeny_h5ad)
print(f"\n{'='*60}")
print(f"Ontogeny Dataset")
print(f"{'='*60}")
print(f"Loaded Ontogeny object with shape: {Ontogeny.shape}")
print(f"Number of cells: {Ontogeny.n_obs}")
print(f"Number of genes: {Ontogeny.n_vars}")

# Check available columns
print(f"\nAvailable columns in Ontogeny.obs:")
print(Ontogeny.obs.columns.tolist())

if 'stage_grouped' in Ontogeny.obs.columns:
    print(f"\nStage_grouped categories: {Ontogeny.obs['stage_grouped'].unique()}")

if 'identity' in Ontogeny.obs.columns:
    print(f"\nIdentity value counts:")
    print(Ontogeny.obs["identity"].value_counts())
    
    # Filter out identities with 0 cells
    identity_counts_before = Ontogeny.obs["identity"].value_counts()
    identities_with_cells = identity_counts_before[identity_counts_before > 0].index.tolist()
    Ontogeny = Ontogeny[Ontogeny.obs['identity'].isin(identities_with_cells)].copy()
    print(f"\nFiltered out {len(identity_counts_before) - len(identities_with_cells)} identities with 0 cells")
    print(f"Remaining identities: {len(identities_with_cells)}")

# Create output directory for Ontogeny
ontogeny_output_dir = args.ontogeny_output_dir
os.makedirs(ontogeny_output_dir, exist_ok=True)
print(f"\nOutput directory created: {ontogeny_output_dir}")

# Get UMAP coordinates (used for CSV export and by plot function)
ontogeny_umap_coords = Ontogeny.obsm['X_umap']

# Ontogeny: faceted UMAP + separate figures by time point (my_palette, identity order by count)
if 'stage_grouped' in Ontogeny.obs.columns and 'identity' in Ontogeny.obs.columns:
    ontogeny_stage_order = ['embryonic', 'shortly_after_birth', 'after_weaning', 'uninfected_adult']
    ontogeny_identity_counts = Ontogeny.obs['identity'].value_counts()
    ontogeny_identity_values = ontogeny_identity_counts.index.tolist()
    print(f"\nOntogeny identity order (by cell count):")
    for identity in ontogeny_identity_values:
        print(f"  {identity}: {ontogeny_identity_counts[identity]} cells")
    fig_ontogeny_facet, figs_ontogeny_individual = plot_umap_by_stage(
        Ontogeny,
        ontogeny_output_dir,
        stage_group_col="stage_grouped",
        identity_col="identity",
        stage_order=ontogeny_stage_order,
        identity_order=None,
        color_palette=my_palette,
        figure_prefix="umap_faceted_by_stage_grouped_identity",
    )
    plt.close(fig_ontogeny_facet)
    for f in figs_ontogeny_individual:
        plt.close(f)
    print(f"\nOntogeny plots saved to: {ontogeny_output_dir}")
    print(f"  - umap_faceted_by_stage_grouped_identity.pdf")
    print(f"  - umap_faceted_by_stage_grouped_identity.png")
else:
    print("\nWarning: Cannot create Ontogeny UMAP plots - stage_grouped or identity column missing.")

# Export CSV files for Ontogeny dataset
# Get cell names (use CellName if available, otherwise use obs_names)
if 'CellName' in Ontogeny.obs.columns:
    ontogeny_cell_names = Ontogeny.obs['CellName'].values
else:
    ontogeny_cell_names = Ontogeny.obs_names.values

# Convert barcodes
ontogeny_converted_barcodes = [convert_barcode_format(barcode) for barcode in ontogeny_cell_names]

# Create three separate CSV files for Loupe Browser

# 1. UMAP coordinates file
ontogeny_umap_df = pd.DataFrame({
    'Barcode': ontogeny_converted_barcodes,
    'UMAP_1': ontogeny_umap_coords[:, 0],
    'UMAP_2': ontogeny_umap_coords[:, 1],
})
ontogeny_umap_output_path = os.path.join(ontogeny_output_dir, 'loupe_browser_umap.csv')
ontogeny_umap_df.to_csv(ontogeny_umap_output_path, index=False)
print(f"\n1. Ontogeny UMAP CSV file saved to: {ontogeny_umap_output_path}")
print(f"   Total cells: {len(ontogeny_umap_df)}")

# 2. Identity file
ontogeny_identity_df = pd.DataFrame({
    'Barcode': ontogeny_converted_barcodes,
    'identity': Ontogeny.obs['identity'].values,
})
ontogeny_identity_output_path = os.path.join(ontogeny_output_dir, 'loupe_browser_identity.csv')
ontogeny_identity_df.to_csv(ontogeny_identity_output_path, index=False)
print(f"\n2. Ontogeny Identity CSV file saved to: {ontogeny_identity_output_path}")
print(f"   Total cells: {len(ontogeny_identity_df)}")

# 3. Stage_grouped file
ontogeny_stage_df = pd.DataFrame({
    'Barcode': ontogeny_converted_barcodes,
    'stage_grouped': Ontogeny.obs['stage_grouped'].values,
})
ontogeny_stage_output_path = os.path.join(ontogeny_output_dir, 'loupe_browser_stage_grouped.csv')
ontogeny_stage_df.to_csv(ontogeny_stage_output_path, index=False)
print(f"\n3. Ontogeny Stage_grouped CSV file saved to: {ontogeny_stage_output_path}")
print(f"   Total cells: {len(ontogeny_stage_df)}")
# Individual stage plots are created inside plot_umap_by_stage()

