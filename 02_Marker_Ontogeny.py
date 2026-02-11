######################################################## Import libraries ########################################################
# conda activate IntestinalFibroblasts_scRNAseq
import argparse
import scanpy as sc
import anndata as ad
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
######################################################## Functions ########################################################
from Function_def import plot_score

######################################################## Parse command line arguments ########################################################
def parse_args():
    parser = argparse.ArgumentParser(
        description="Marker gene expression and cell cycle analysis for Ontogeny and Parasite datasets."
    )
    parser.add_argument(
        "--ontogeny_h5ad",
        type=str,
        default="Data/OntogenyData_OnlyFBs_Xuan.h5ad",
        help="Path to Ontogeny h5ad file (default: Data/OntogenyData_OnlyFBs_Xuan.h5ad)",
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
        default="Results/02_Marker_Ontogeny",
        help="Output directory for all plots and CSV files (default: Results/02_Marker_Ontogeny)",
    )
    return parser.parse_args()

######################################################## Parse command line and load data ########################################################
# Parse command line and load data (only assign when argument is non-empty)
args = parse_args()
# args.ontogeny_h5ad = "Data/OntogenyData_OnlyFBs_Xuan.h5ad"
# args.parasite_h5ad = "Data/Parasite_Clustered.h5ad"
# args.output_dir = "Results/02_Marker_Ontogeny_test"

if args.ontogeny_h5ad and str(args.ontogeny_h5ad).strip():
    Ontogeny = sc.read_h5ad(args.ontogeny_h5ad.strip())
else:
    Ontogeny = None

if args.output_dir and str(args.output_dir).strip():
    output_dir_base = args.output_dir.strip()
    os.makedirs(output_dir_base, exist_ok=True)
    print(f"\nOutput directory created: {output_dir_base}")
else:
    output_dir_base = "Results/02_Marker_Ontogeny"
    os.makedirs(output_dir_base, exist_ok=True)
    print(f"\nOutput directory created: {output_dir_base}")

if args.parasite_h5ad and str(args.parasite_h5ad).strip():
    Parasite = sc.read_h5ad(args.parasite_h5ad.strip())
else:
    Parasite = None

######################################################## Main process ########################################################
# Marker gene expression for Ontogeny (Figure 6C, starting with the genes mentioned in figure 4b)
# Remove empty category "PdgfraloCD55+ FB" from identity
if 'identity' in Ontogeny.obs.columns:
    Ontogeny = Ontogeny[Ontogeny.obs['identity'] != 'PdgfraloCD55+ FB'].copy()
    print(f"\nFiltered out 'PdgfraloCD55+ FB' category. Remaining cells: {Ontogeny.n_obs}")

# Get unique identity values ordered by cell count (descending)
if 'identity' in Ontogeny.obs.columns:
    ontogeny_identity_counts = Ontogeny.obs['identity'].value_counts()
    # Order by cell count from highest to lowest
    ontogeny_identity_values = ontogeny_identity_counts.index.tolist()
else:
    print("\nWarning: 'identity' column not found in Ontogeny dataset.")
    ontogeny_identity_values = []
    ontogeny_identity_counts = pd.Series()
    

genelist =['Pdpn', 'Pdgfra', 'Igfbp3', 'Fgfr2', 'Pdgfrb', 'Des', 'Acta2', 'Myh11', 'Actg2', 'Cd81', 'Ackr4', 'Cd55', 'Cd34', 'Jam2', 'Ccl11', 'Kit', 'Hand1']

ax = sc.pl.dotplot(Ontogeny, genelist, groupby='identity', dendrogram=True, dot_max=1, dot_min=0, standard_scale='var', return_fig=True, use_raw=False)
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir_base, 'Marker_Ontogeny_dotplot.pdf'))

# Create stacked_violin plot
sc.pl.stacked_violin(Ontogeny, genelist, groupby='identity', dendrogram=True, standard_scale='var', show=False, use_raw=False)
# Get the current figure and save it
fig = plt.gcf()
fig.tight_layout()
fig.savefig(os.path.join(output_dir_base, 'Marker_Ontogeny_stacked_violin.pdf'), bbox_inches='tight', pad_inches=0.1)
plt.close(fig)

########################################################
# Side question: Does parasite infection induce Tnfsf11 expression in intestinal FSCs?
genelist =['Pdpn', 'Pdgfra', 'Igfbp3', 'Fgfr2', 'Pdgfrb', 'Des', 'Acta2', 'Myh11', 'Actg2', 'Cd81', 'Ackr4', 'Cd55', 'Cd34', 'Jam2', 'Ccl11', 'Kit', 'Hand1','Tnfsf11']

ax = sc.pl.dotplot(Parasite, genelist, groupby='identity', dendrogram=True, dot_max=1, dot_min=0, standard_scale='var', return_fig=True, use_raw=False)
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir_base, 'Marker_Parasite_dotplot.pdf'))

########################################################
# Cell cycle score for Ontogeny (Figure 6D)
# genelist from /Users/xliu2942/Documents/Projects/IntestinalFibroblasts_scRNAseq/IntestinalFibroblasts_scRNAseq-main/07_Parasite_Clustering.ipynb
# output_dir is "Cell_Cycle_Score"
output_dir = os.path.join(output_dir_base, 'Cell_Cycle_Score')
os.makedirs(output_dir, exist_ok=True)

G2M_genes_list = ["Hmgb2","Cdk1","Nusap1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
                  "Nuf2","Cks1b","Mki67","Cenpf","Tacc3","Pimreg","Smc4","Ccnb2","Ckap2l",
                  "Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b","Hjurp",
                  "Cdca3","Jpt1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2","Dlgap5",
                  "Cdca2","Cdca8","Ect2", "Kif23","Hmmr","Aurka","Psrc1","Anln","Lbr","Ckap5",
                  "Cenpe","Ctcf", "Nek2","G2e3","Gas2l3", "Cbx5","Cenpa"]
S_genes_list = ["Mcm5","Pcna","Tyms","Fen1","Mcm7","Mcm4","Rrm1","Ung","Gins2","Mcm6","Cdca7",
                "Dtl","Prim1","Uhrf1","Cenpu","Hells","Rfc2","Polr1b","Nasp","Rad51ap1","Gmnn",
                "Wdr76","Slbp","Ccne2","Ubr7","Msh2","Rad51","Rrm2","Cdc45","Cdc6","Exo1","Tipin",
                "Dscc1","Blm","Casp8ap2", "Usp1","Clspn","Pola1","Chaf1b","Mrpl36","E2f8"]
# save the gene list to a excel file (or CSV if openpyxl is not available)
# Create a DataFrame with all genes and their types
all_genes = G2M_genes_list + S_genes_list
gene_types = ['G2M']*len(G2M_genes_list) + ['S']*len(S_genes_list)
gene_list_df = pd.DataFrame({'gene_name': all_genes, 'gene_type': gene_types})
gene_list_df.to_csv(os.path.join(output_dir, 'Cell_Cycle_Score_MarkerGenes_list.csv'), index=False)

sc.tl.score_genes_cell_cycle(Ontogeny, s_genes=S_genes_list, g2m_genes=G2M_genes_list, use_raw=False) 
ax = sc.pl.umap(Ontogeny, color = "phase", use_raw=False, return_fig=True)
ax.tight_layout()
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Ontogeny_cell_cycle_score_phase.pdf'), bbox_inches='tight', pad_inches=0.1)
ax = sc.pl.umap(Ontogeny, color = "G2M_score", use_raw=False, return_fig=True)
ax.tight_layout()
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Ontogeny_cell_cycle_score_G2M_score.pdf'), bbox_inches='tight', pad_inches=0.1)
ax = sc.pl.umap(Ontogeny, color = "S_score", use_raw=False, return_fig=True)
ax.tight_layout()
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Ontogeny_cell_cycle_score_S_score.pdf'), bbox_inches='tight', pad_inches=0.1)
# umap shows "stage_grouped"
ax = sc.pl.umap(Ontogeny, color = "stage_grouped", use_raw=False, return_fig=True)
ax.tight_layout()
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Ontogeny_stage_grouped.pdf'), bbox_inches='tight', pad_inches=0.1)

# G2M_score and S_score faceted by stage_grouped
stage_grouped_order = ['embryonic', 'shortly_after_birth', 'after_weaning', 'uninfected_adult']

# Plot G2M_score faceted by stage_grouped + separate figures per stage
fig_facet, figs_ind = plot_score(
    adata=Ontogeny,
    color_values=Ontogeny.obs['G2M_score'].values,
    output_dir=output_dir,
    color_label='G2M_score',
    stage_order=stage_grouped_order,
    figure_prefix='umap_faceted_by_stage_grouped_G2M_score',
)
if fig_facet is not None:
    plt.close(fig_facet)
for f in figs_ind:
    plt.close(f)

# Plot S_score faceted by stage_grouped + separate figures per stage
fig_facet, figs_ind = plot_score(
    adata=Ontogeny,
    color_values=Ontogeny.obs['S_score'].values,
    output_dir=output_dir,
    color_label='S_score',
    stage_order=stage_grouped_order,
    figure_prefix='umap_faceted_by_stage_grouped_S_score',
)
if fig_facet is not None:
    plt.close(fig_facet)
for f in figs_ind:
    plt.close(f)


########################################################
output_dir = os.path.join(output_dir_base, 'Cdk1_by_stage_and_phase')
os.makedirs(output_dir, exist_ok=True)
print(f"\nOutput directory created: {output_dir}")

# Cdk1 expression by developmental stage and cell cycle phase
# Define the specific order for categories for "stage_grouped"
stage_grouped_order = ['embryonic', 'shortly_after_birth', 'after_weaning', 'uninfected_adult']

# Define the specific order for categories for "phase"
phases_order = ['G1', 'S', 'G2M']

# Determine the number of categories
n_stages = len(stage_grouped_order)
n_phases = len(phases_order)

# Create a figure with subplots - each subplot corresponds to a combination of "stage_grouped" and "phase"
fig, axs = plt.subplots(n_stages, n_phases, figsize=(4*n_phases, 4*n_stages))  # Adjust figsize as needed

for i, stage in enumerate(stage_grouped_order):
    for j, phase in enumerate(phases_order):
        # Subset data for the current combination of stage and phase
        adata_subset = Ontogeny[(Ontogeny.obs['stage_grouped'] == stage) & (Ontogeny.obs['phase'] == phase)].copy()
        
        # Check if there are any samples in the subset
        if not adata_subset.obs.empty:
            # Plot UMAP for the subset on the current subplot
            sc.pl.umap(adata_subset, title=f'{stage} - {phase}', ax=axs[i, j], show=False, use_raw=False, color="Cdk1", size=10)
            axs[i, j].set_xlabel('')
            axs[i, j].set_ylabel('')
        else:
            axs[i, j].axis('off')  # Hide axis if no samples in the subset

# Adjust the layout and add a general title
fig.tight_layout()
fig.suptitle('Cdk1 Expression by Developmental Stage and Cell Cycle Phase', fontsize=16, y=1.02)  # Adjust fontsize as needed

# Save the plot as a PDF
fig.savefig(os.path.join(output_dir, 'Ontogeny_cell_cycle_score_Cdk1_by_stage_and_phase.pdf'), bbox_inches='tight', pad_inches=0.1)

# Single Cdk1 gene expression UMAP
ax = sc.pl.umap(Ontogeny, color = "Cdk1", use_raw=False, return_fig=True)
ax.tight_layout()
# Save the plot as a PDF
ax.savefig(os.path.join(output_dir, 'Ontogeny_Cdk1_expression.pdf'), bbox_inches='tight', pad_inches=0.1)

########################################################
# Cdk1 expression faceted by stage_grouped
# Get Cdk1 expression values (use_raw=False to match the rest of the script)
if 'Cdk1' in Ontogeny.var_names:
    # Get expression from the X matrix (not raw)
    cdk1_idx = list(Ontogeny.var_names).index('Cdk1')
    if hasattr(Ontogeny.X, 'toarray'):
        cdk1_expression = Ontogeny.X[:, cdk1_idx].toarray().flatten()
    else:
        cdk1_expression = Ontogeny.X[:, cdk1_idx]
else:
    print("Warning: Cdk1 gene not found in dataset")
    cdk1_expression = np.zeros(Ontogeny.n_obs)

# Plot Cdk1 expression faceted by stage_grouped + separate figures per stage
fig_facet, figs_ind = plot_score(
    adata=Ontogeny,
    color_values=cdk1_expression,
    output_dir=output_dir,
    color_label='Cdk1 expression',
    stage_order=stage_grouped_order,
    figure_prefix='umap_faceted_by_stage_grouped_Cdk1',
)
if fig_facet is not None:
    plt.close(fig_facet)
for f in figs_ind:
    plt.close(f)

