# conda activate IntestinalFibroblasts_scRNAseq
#=========================Script Description=================================
# Force reticulate to use Python from conda env (must set before library(reticulate))
python_path <- path.expand("~/.conda/envs/IntestinalFibroblasts_scRNAseq/bin/python")
if (file.exists(python_path)) Sys.setenv(RETICULATE_PYTHON = python_path)
#=========================Loading Packages===================================
#BiocManager::install("beachmat")
#BiocManager::install("DelayedMatrixStats")
#BiocManager::install("SparseArray")
#BiocManager::install("DelayedArray")
#BiocManager::install("SummarizedExperiment")
#BiocManager::install("SingleCellExperiment")
# mamba install -c conda-forge r-base=4.5.2
suppressMessages(library("optparse"))
# mamba install -c conda-forge zlib
# BiocManager::install("HDF5Array")
# BiocManager::install("glmGamPoi")
suppressMessages(library("glmGamPoi"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
#remove.packages("openxlsx")
suppressMessages(library("openxlsx"))
suppressMessages(library("ggrepel"))
suppressMessages(library("reticulate"))
# Error: To use the "deconvolution" method for size factor calculation, you need to install the 'scran' package from Bioconductor. Alternatively, you can use "normed_sum" or "poscounts"to calculate the size factors.
suppressMessages(library("scran"))
#=========================Function Definition================================
# ========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Results/03_CellAnnotation/Parasite_annotated.h5ad",
              help = "Parasite_annotated.h5ad"
  ),make_option(c("--output_dir", "-o"),
                type = "character", default = "Results/05_Pseudobulk",
                help = "output directory path."
  ),make_option(c("--cluster_col", "-c"),
              type = "character", default = "New_name",
              help = "cluster column name."
  ),make_option(c("--design", "-d"),
              type = "character", default = "~1+stage_grouped",
              help = "design formula. default is ~1+stage_grouped, can also be ~1+New_name+stage_grouped+New_name:stage_grouped to identify celltype specific effects"
  ),make_option(c("--iteration_bycelltype"),
              type = "logical", default = FALSE,
              help = "whether to iterate by cell type. default is FALSE"
  ),make_option(c("--feature_anno", "-f"),
              type = "character", default = "Data/Genomic_Features.tsv",
              help = "path to genomic features annotation file. default is /data/gpfs/projects/punim2435/Project/IntestinalFibroblasts/Data/Genomic_Features.tsv"
    ),make_option(c("--subset"),type = "character", default = NULL,
              help = "subset the data by a column name. default is NULL"
    ),make_option(c("--idents"),type = "character", default = NULL,
              help = "Subset by identity. default is NULL"
)
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
#============================================================================
if (is.null(opt$output)) {
  print("NO OUTPUT PATH SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

if (is.null(opt$input)) {
  print("NO INPUT PATH SUPPLIED,current directory will be used!")
  input_dir <- getwd()
} else {
  input_dir <- opt$input
  if (!file.exists(input_dir)) {
    dir.create(input_dir, recursive = T)
  }
}

cluster_col <- opt$cluster_col

if (is.null(opt$design)) {
  print("NO DESIGN FORMULA SUPPLIED,default design formula will be used!")
  design <- "~1+stage_grouped"
} else {
  design <- opt$design
  design <- formula(design)
}

if (is.null(opt$feature_anno)) {
  print("NO FEATURE ANNOTATION FILE SUPPLIED,default feature annotation file will be used!")
  feature_anno <- "Data/Genomic_Features.tsv"
} else {
  feature_anno <- opt$feature_anno
  genomic_features <- read.csv(feature_anno, row.names = NULL, sep="\t")
}

#=========================Loading Data========================================
# need to Get the count matrix
# read h5ad file (uses Python/scanpy via reticulate)
scanpy <- reticulate::import("scanpy")
adata <- scanpy$read_h5ad(opt$input)

if (!is.null(opt$subset) && !is.null(opt$idents)) {
  # opt$subset = "New_name"
  # opt$identity = "Myofibroblast 1 (MF1),Myofibroblast 2 (MF2)"
  # if opt$identity is a comma separated list of values, split it into a vector
  if (grepl(",", opt$idents)) {
    identity_vector <- trimws(strsplit(opt$idents, ",")[[1]])
  } else {
    identity_vector <- trimws(opt$idents)
  }
  # Subset adata in Python (R cannot use .isin() on reticulate objects)
  # Ensure identity_vec is a list in Python (R length-1 vector becomes str otherwise)
  py$subset_col <- opt$subset
  py$identity_vec <- as.list(identity_vector)
  py$adata <- adata
  reticulate::py_run_string("mask = adata.obs[subset_col].astype(str).isin(identity_vec).values; adata = adata[mask].copy()")
  adata <- py$adata
  # save the adata to h5ad file
  adata$write_h5ad(paste0(output_dir, "/adata_subset.h5ad"))
}

genomic_features_unique <- genomic_features %>%
  dplyr::select(c(Symbol, Name)) %>%
  dplyr::distinct()
protein_coding_genes <- unique(subset(genomic_features, Gene.Type == "protein-coding")$Symbol)

#This is to add full names of the genes to the results
gene_metadata <- dplyr::select(genomic_features, c(Name, Symbol)) %>% 
  dplyr::rename("name" = "Symbol",
                "open_name" = "Name") %>%
  distinct()

# extract metadata 
coldata <- adata$obs %>% dplyr::select(CellName, dplyr::all_of(cluster_col), identity, stage, stage_grouped)
rownames(coldata) <- coldata$CellName

# always have cell annotation + stage_grouped in the condition
coldata <- coldata %>%
  dplyr::mutate(condition = paste0(.data[[cluster_col]], "_", stage_grouped))

# Extract the count matrix
cnt <- reticulate::py_to_r(adata$layers[["raw_counts"]])
rownames(cnt) <- coldata$CellName
var_metadata <- adata$var
colnames(cnt) <- rownames(var_metadata)
cnt <- as.data.frame(t(cnt))

#Keep only the protein coding genes
cnt <- subset(cnt, rownames(cnt) %in% protein_coding_genes)

#Keep only the protein coding genes
#cnt <- subset(cnt, rownames(cnt) %in% protein_coding_genes)

# Filter out genes with low counts
cnt <- cnt[rowSums(cnt) > 10, ]

#Generate Pseudobulk Sample
condition_sample_size <- coldata %>%
  dplyr::group_by(condition) %>%
  dplyr::summarise(sample_size = n()) %>%
  dplyr::mutate(subsample_size = ceiling(sample_size/10))

coldata_withsamplesize <- merge(coldata, condition_sample_size, by="condition", all=TRUE)

# Initialize a list to store the data frames from each iteration
pseudobulk_samples <- list()

# Outer loop for each unique condition
for (i in unique(coldata_withsamplesize$condition)) {
  # Inner loop to iterate 10 times for each condition
  for (j in 1:10) {
    # 10-fold subsample sizes were previously calculated in the subsample_size column
    i_subsample_size <- unique(subset(coldata_withsamplesize, condition == i)$subsample_size)
    # Retrieve all the sample IDs that correspond to the condition
    i_sample_pool <- subset(coldata_withsamplesize, condition == i)$CellName
    # Do the subsampling on the retrieved sample IDs (i.e., return a subset of IDs)
    i_subsample_ids <- sample(i_sample_pool, size = i_subsample_size, replace = FALSE)
    # Use these IDs to actually subset the count matrix, ensuring it remains a matrix
    i_cnt <- cnt[,i_subsample_ids, drop = FALSE]
    # Calculate the average count for each gene (i.e., rowwise means). This is 1 pseudobulk replicate.
    i_pseudobulk_sample <- as.integer(rowMeans(i_cnt))

    # Store each replicate in the list with a unique name
    pseudobulk_samples[[paste0(i, "_Rep", j)]] <- i_pseudobulk_sample
  }
}

# Convert the list of vectors to a data frame
final_pseudobulk_df <- as.data.frame(pseudobulk_samples)

# Manually set the column names from the list to retain original characters
colnames(final_pseudobulk_df) <- names(pseudobulk_samples)

# Optional: if you want to set the row names as gene names
# Assuming the row names of 'cnt' are gene names
rownames(final_pseudobulk_df) <- rownames(cnt)


pseudobulk_coldata <- as.data.frame(colnames(final_pseudobulk_df))
colnames(pseudobulk_coldata) <- "SampleName"
rownames(pseudobulk_coldata) <- pseudobulk_coldata$SampleName
# Using stringr to extract the condition from SampleName
pseudobulk_coldata$condition <- stringr::str_extract(pseudobulk_coldata$SampleName, ".*?(?=_Rep)")
pseudobulk_coldata_merged <- merge(pseudobulk_coldata,
                                   dplyr::select(coldata, -c(CellName, stage)),
                                   by = "condition") %>% distinct()
rownames(pseudobulk_coldata_merged) <- pseudobulk_coldata_merged$SampleName

# Reorder the columns of final_pseudobulk_df to match the order in pseudobulk_coldata$SampleName
final_pseudobulk_df <- final_pseudobulk_df[, rownames(pseudobulk_coldata_merged)]

unique(colnames(final_pseudobulk_df) == rownames(pseudobulk_coldata_merged))

# save the final_pseudobulk_df to an RDS file
saveRDS(cnt, file = paste0(output_dir, "/cnt.rds"))
saveRDS(final_pseudobulk_df, file = paste0(output_dir, "/final_pseudobulk_df.rds"))
saveRDS(pseudobulk_coldata_merged, file = paste0(output_dir, "/pseudobulk_coldata_merged.rds"))

# DEX analysis
#design <- formula(paste0("~1+", cluster_col))

# find which column in design formula is used
vars <- all.vars(as.formula(design))

if (cluster_col %in% vars && length(vars) > 1 && opt$iteration_bycelltype == FALSE) {
  # in this case, design = ~1+New_name+stage_grouped+New_name:stage_grouped
  # Create a new workbook
  wb <- createWorkbook()
  fit <- glm_gp(as.matrix(final_pseudobulk_df),
              col_data = pseudobulk_coldata_merged,
              design = design,
              size_factors = "deconvolution",
              subsample = FALSE,
              verbose=TRUE,
              on_disk = FALSE)
  # For each gene, whether the stage change trend is different for different cell types?
  # Is this gene more stage-dependent in some cell types than in others?
  res <- test_de(fit, reduced_design = as.formula(paste0("~1+", cluster_col, "+stage_grouped")),
                            max_lfc = Inf,
                            verbose = TRUE) %>%
        dplyr::select(-c(df1, df2)) %>%
        dplyr::rename("Symbol" = "name") %>%
        merge(genomic_features_unique, by = "Symbol") %>%
        dplyr::mutate(
        pval = case_when(
            pval == 0 ~ .Machine$double.xmin,
            .default = pval),
        adj_pval = case_when(
            adj_pval == 0 ~ .Machine$double.xmin,
            .default = adj_pval)) %>%
        dplyr::arrange(desc(f_statistic))

  # Save the fit and res lists to an RData file
  saveRDS(fit, file = paste0(output_dir, "/DEX_by_", cluster_col, "_stage_grouped_fit.rds"))
  saveRDS(res, file = paste0(output_dir, "/DEX_by_", cluster_col, "_stage_grouped_res.rds"))

  addWorksheet(wb, sheetName = paste0("DEG_by_", cluster_col, "_stage_grouped"))
  writeData(wb, sheet = paste0("DEG_by_", cluster_col, "_stage_grouped"), x = res)
  # Save the workbook to an Excel file
  saveWorkbook(wb, paste0(output_dir, "/DEX_by_", cluster_col, "_stage_grouped.xlsx"), overwrite = TRUE)

} else if (opt$iteration_bycelltype == TRUE) {
    # in this case, deign = ~1+stage_grouped, but we need to iterate by cell type
    # Create a new workbook
    wb <- createWorkbook()
    fit_list <- list()
    res_list <- list()
    # Get unique identities
    unique_identities <- unique(pseudobulk_coldata_merged[[cluster_col]])

    # Loop through each unique identity
    for (i in unique_identities) {
    coldata_subset <- pseudobulk_coldata_merged[pseudobulk_coldata_merged[[cluster_col]] == i, ]
    final_pseudobulk_df_subset <- final_pseudobulk_df[, coldata_subset$SampleName]
    
    fit <- glm_gp(as.matrix(final_pseudobulk_df_subset),
                            col_data = coldata_subset,
                            design = design,
                            size_factors = "deconvolution",
                            subsample = FALSE,
                            verbose = TRUE,
                            on_disk = FALSE)
    fit_list[[i]] <- fit

    res <- test_de(fit, reduced_design = ~1,
                            max_lfc = Inf,
                            verbose = TRUE) %>%
        dplyr::select(-c(df1, df2)) %>%
        dplyr::rename("Symbol" = "name") %>%
        merge(genomic_features_unique, by = "Symbol") %>%
        dplyr::mutate(
        pval = case_when(
            pval == 0 ~ .Machine$double.xmin,
            .default = pval),
        adj_pval = case_when(
            adj_pval == 0 ~ .Machine$double.xmin,
            .default = adj_pval)) %>%
        dplyr::arrange(desc(f_statistic))
    res_list[[i]] <- res

    # Add a new worksheet for the current identity
    addWorksheet(wb, sheetName = substr(as.character(i), 1, 31))
    
    # Write the result data frame to the corresponding worksheet
    writeData(wb, sheet = substr(as.character(i), 1, 31), x = res)
    }
    # Save the workbook to an Excel file
    saveWorkbook(wb, paste0(output_dir, "/DEX_by_stage_grouped_iteration_by_",cluster_col,"_celltype.xlsx"), overwrite = TRUE)
    # Save the fit and res lists to an RData file
    saveRDS(fit_list, file = paste0(output_dir, "/DEX_by_stage_grouped_iteration_by_",cluster_col,"_celltype_fit_list.rds"))
    saveRDS(res_list, file = paste0(output_dir, "/DEX_by_stage_grouped_iteration_by_",cluster_col,"_celltype_res_list.rds"))

} else if (vars == "stage_grouped" && length(vars) == 1 && opt$iteration_bycelltype == FALSE) {
    # in this case, design = ~1+stage_grouped, but without iteration by cell type (general analysis)
    # Create a new workbook
    wb <- createWorkbook()
    
    fit <- glm_gp(as.matrix(final_pseudobulk_df),
                          col_data = pseudobulk_coldata_merged,
                          design = design,
                          size_factors = "deconvolution",
                          subsample = FALSE,
                          verbose = TRUE,
                          on_disk = FALSE)

    res <- test_de(fit, reduced_design = ~1,
                            max_lfc = Inf,
                            verbose = TRUE) %>%
        dplyr::select(-c(df1, df2)) %>%
        dplyr::rename("Symbol" = "name") %>%
        merge(genomic_features_unique, by = "Symbol") %>%
        dplyr::mutate(
        pval = case_when(
            pval == 0 ~ .Machine$double.xmin,
            .default = pval),
        adj_pval = case_when(
            adj_pval == 0 ~ .Machine$double.xmin,
            .default = adj_pval)) %>%
        dplyr::arrange(desc(f_statistic))

    # Save the fit and res lists to an RData file
    saveRDS(fit, file = paste0(output_dir, "/DEX_by_stage_grouped_fit.rds"))
    saveRDS(res, file = paste0(output_dir, "/DEX_by_stage_grouped_res.rds"))
    addWorksheet(wb, sheetName = "DEG_by_stage_grouped")
    writeData(wb, sheet = "DEG_by_stage_grouped", x = res)
    # Save the workbook to an Excel file
    saveWorkbook(wb, paste0(output_dir, "/DEX_by_stage_grouped.xlsx"), overwrite = TRUE)
}



