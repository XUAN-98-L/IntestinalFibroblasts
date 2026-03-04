#conda activate IntestinalFibroblasts_scRNAseq
# Barcode extract for loupe browser export
python Script/01_Barcode_extract.py --h5ad_file Data/Parasite_Clustered.h5ad --output_dir Results/01_Barcode_extract/Parasite_Clustered --ontogeny_h5ad Data/OntogenyData_OnlyFBs_Xuan.h5ad --ontogeny_output_dir Results/01_Barcode_extract/Ontogeny_Clustered

# Marker gene expression and cell cycle analysis for Ontogeny and Parasite datasets
python Script/02_Marker_Ontogeny.py --ontogeny_h5ad Data/OntogenyData_OnlyFBs_Xuan.h5ad --parasite_h5ad Data/Parasite_Clustered.h5ad --output_dir Results/02_Marker_Ontogeny

# Cell annotation for Parasite dataset
python Script/03_CellAnnotation.py --rename_path Data/Cluster_infectiondataset_renamed.xlsx --parasite_h5ad Data/Parasite_Clustered.h5ad --output_dir Results/03_CellAnnotation

# Loupe browser export for Parasite dataset
Rscript Script/04_Loupe.R --input Data/Parasite_Clustered.h5ad --output_dir Results/04_Loupe_Parasite --cluster_col New_name

# Loupe browser export for Ontogeny dataset
Rscript Script/04_Loupe.R --input Data/OntogenyData_OnlyFBs_Xuan.h5ad --output_dir Results/04_Loupe_Ontogeny --cluster_col identity

# Pseudobulk analysis for Parasite dataset
# rm -rf Results/05_Pseudobulk_*
Rscript Script/05_Pseudobulk.R --input Results/03_CellAnnotation/Parasite_annotated.h5ad --output_dir Results/05_Pseudobulk_Parasite/identity_stage_grouped --cluster_col New_name --design "~1+New_name+stage_grouped+New_name:stage_grouped" --iteration_bycelltype FALSE --feature_anno Data/Genomic_Features.tsv

# Pseudobulk analysis for Parasite dataset by celltype
Rscript Script/05_Pseudobulk.R --input Results/03_CellAnnotation/Parasite_annotated.h5ad --output_dir Results/05_Pseudobulk_Parasite/stage_grouped_bycelltype --cluster_col New_name --design "~1+stage_grouped" --iteration_bycelltype TRUE --feature_anno Data/Genomic_Features.tsv

# Pseudobulk analysis for Parasite dataset
Rscript Script/05_Pseudobulk.R --input Results/03_CellAnnotation/Parasite_annotated.h5ad --output_dir Results/05_Pseudobulk_Parasite/stage_grouped --cluster_col New_name --design "~1+stage_grouped" --iteration_bycelltype FALSE --feature_anno Data/Genomic_Features.tsv

# For ontogeny dataset, try "identity" "Ccl11+ MF,Ccl11-Osr1+ MF,Ccl11-Actg2hi MF,Ccl11- MF"
Rscript Script/05_Pseudobulk.R --input Data/OntogenyData_OnlyFBs_Xuan.h5ad --output_dir Results/05_Pseudobulk_Ontogeny/MF_only --cluster_col identity --design "~1+stage_grouped" --iteration_bycelltype FALSE --feature_anno Data/Genomic_Features.tsv --subset "identity" --idents "Ccl11+ MF,Ccl11-Osr1+ MF,Ccl11-Actg2hi MF,Ccl11- MF"

# For ontogeny dataset, try "identity" "Pericytes"
Rscript Script/05_Pseudobulk.R --input Data/OntogenyData_OnlyFBs_Xuan.h5ad --output_dir Results/05_Pseudobulk_Ontogeny/Pericytes_only --cluster_col identity --design "~1+stage_grouped" --iteration_bycelltype FALSE --feature_anno Data/Genomic_Features.tsv --subset "identity" --idents "Pericytes"

# For ontogeny dataset, try "identity" "Colec10+_prolif. FB"
Rscript Script/05_Pseudobulk.R --input Data/OntogenyData_OnlyFBs_Xuan.h5ad --output_dir Results/05_Pseudobulk_Ontogeny/Colec10_prolif_FB_only --cluster_col identity --design "~1+stage_grouped" --iteration_bycelltype FALSE --feature_anno Data/Genomic_Features.tsv --subset "identity" --idents "Colec10+_prolif. FB"

# For celltype subset analysis, compare between shortly_after_birth vs. embryonic, after_weaning vs. shortly_after_birth, 'uninfected_adult vs. after_weaning
Rscript Script/05_Pseudobulk_celltype_subset.R --input Results/05_Pseudobulk_Ontogeny/MF_only --output_dir Results/05_Pseudobulk_Ontogeny/MF_only/pairwise_stage --feature_anno Data/Genomic_Features.tsv --stage_order "embryonic,shortly_after_birth,after_weaning,uninfected_adult"

Rscript Script/05_Pseudobulk_celltype_subset.R --input Results/05_Pseudobulk_Ontogeny/Colec10_prolif_FB_only --output_dir Results/05_Pseudobulk_Ontogeny/Colec10_prolif_FB_only/pairwise_stage --feature_anno Data/Genomic_Features.tsv --stage_order "embryonic,shortly_after_birth,after_weaning,uninfected_adult"

Rscript Script/05_Pseudobulk_celltype_subset.R --input Results/05_Pseudobulk_Ontogeny/Pericytes_only --output_dir Results/05_Pseudobulk_Ontogeny/Pericytes_only/pairwise_stage --feature_anno Data/Genomic_Features.tsv --stage_order "embryonic,shortly_after_birth,after_weaning,uninfected_adult"

# pairwise DEG analysis for Parasite dataset
python Script/06_DEG_rank_genes_groups.py --input Data/OntogenyData_OnlyFBs_Xuan.h5ad --output_dir Results/06_DEG_rank_genes_groups_Ontogeny --groupby "identity" --method "wilcoxon" --n_genes 10

# import sys
# sys.argv = [
#  "script",
#  "--input","Data/OntogenyData_OnlyFBs_Xuan.h5ad",
#  "--output_dir","Results/06_DEG_rank_genes_groups_Ontogeny",
#  "--groupby","identity",
#  "--method","wilcoxon",
#  "--n_genes","10"
# ]
# args = parse_args()