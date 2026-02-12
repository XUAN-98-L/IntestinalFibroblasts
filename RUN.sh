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