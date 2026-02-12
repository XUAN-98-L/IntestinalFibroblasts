# conda activate IntestinalFibroblasts_scRNAseq
#=========================Script Description=================================
# Force reticulate to use Python from conda env (must set before library(reticulate))
python_path <- path.expand("~/.conda/envs/IntestinalFibroblasts_scRNAseq/bin/python")
if (file.exists(python_path)) Sys.setenv(RETICULATE_PYTHON = python_path)
#=========================Loading Packages===================================
suppressMessages(library("optparse"))
# https://www.10xgenomics.com/support/software/loupe-browser/latest/tutorials/introduction/lb-louper
suppressMessages(library("loupeR"))
# remotes::install_github("10xGenomics/loupeR")
loupeR::setup()
suppressMessages(library("Seurat"))
suppressMessages(library("tidyverse"))
suppressMessages(library("data.table"))
suppressMessages(library("Matrix"))
suppressMessages(library("openxlsx"))
suppressMessages(library("reticulate"))
#=========================Function Definition================================
# ========================Command Parameters Setting=========================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Results/03_CellAnnotation/Parasite_annotated.h5ad",
              help = "Parasite_annotated.h5ad"
  ),make_option(c("--output_dir", "-o"),
                type = "character", default = "Results/04_Loupe_Parasite",
                help = "output directory path."
  ),make_option(c("--cluster_col", "-c"),
              type = "character", default = "New_name",
              help = "cluster column name."
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

if (is.null(opt$cluster_col)) {
  print("NO CLUSTER COLUMN NAME SUPPLIED, defaulting to New_name!")
  cluster_col <- "New_name"
} else {
  cluster_col <- opt$cluster_col
}

#============================================================================

# before running loupeR, need to Get the count matrix
# read h5ad file (uses Python/scanpy via reticulate)
scanpy <- reticulate::import("scanpy")
adata <- scanpy$read_h5ad(opt$input)

# Get the count matrix (loupeR requires dgCMatrix with cell and gene names)
count_mat <- reticulate::py_to_r(adata$layers[["raw_counts"]])
#if (!inherits(count_mat, "dgCMatrix")) count_mat <- Matrix::Matrix(as.matrix(count_mat), sparse = TRUE)
# AnnData: rows = cells (obs), cols = genes (var) — Index.tolist() gives R vector of names
colnames(count_mat) <- reticulate::py_to_r(adata$var_names$tolist())

# Use Barcode from obs if present, else Unnamed: 0; convert to Loupe-readable format via Python
obs_r <- reticulate::py_to_r(adata$obs)
if ("Barcode" %in% names(obs_r)) {
  raw_barcodes <- as.character(obs_r[["Barcode"]])
} else if ("Unnamed: 0" %in% names(obs_r)) {
  raw_barcodes <- as.character(obs_r[["Unnamed: 0"]])
} else {
  stop("No Barcode or Unnamed: 0 column found in adata$obs!")
}
# Convert barcodes with Function_def.convert_barcode_format (e.g. AAACCCAAGTATCTGC-1_1 -> AAACCCAAGTATCTGC-1)
script_dir <- if (file.exists("Script/Function_def.py")) "Script" else "."
reticulate::py_run_string(paste0("import sys; sys.path.insert(0, '", script_dir, "')"))
reticulate::py_run_string("from Function_def import convert_barcode_format")
py$barcode_list <- raw_barcodes
reticulate::py_run_string("result = [convert_barcode_format(b) for b in barcode_list]")
cell_barcodes <- reticulate::py_to_r(py$result)
rownames(count_mat) <- cell_barcodes

# loupeR expects genes x cells (barcodes as colnames); AnnData has cells x genes → transpose
count_mat <- Matrix::t(count_mat)
if (!inherits(count_mat, "dgCMatrix")) count_mat <- as(count_mat, "dgCMatrix")

# Load the clusters file (one per cell, names = barcodes)
clusters <- factor(reticulate::py_to_r(adata$obs[[cluster_col]]))
#names(clusters) <- cell_barcodes

# extract all other metadata columns
metadata <- adata$obs
metadata <- metadata[, -which(names(metadata) %in% c(cluster_col, "Barcode","Unnamed: 0")), drop = FALSE]
metadata <- as.data.frame(metadata)
rownames(metadata) <- cell_barcodes

# ensure correct order
metadata <- metadata[cell_barcodes, , drop = FALSE]

# keep only columns that are already factors
metadata_factor <- metadata[, sapply(metadata, is.factor), drop = FALSE]
# convert each metadata column to factor (as vectors)
meta_factors <- lapply(metadata_factor, function(x) factor(as.vector(x)))
# combine your main clusters + all metadata factors into one list
cluster_list <- c(list(clusters = clusters), meta_factors)
# set names (cell barcodes) on each element for loupeR; skip elements with wrong length
cluster_list <- lapply(cluster_list, function(x) {
  if (length(x) == length(cell_barcodes)) names(x) <- cell_barcodes
  x
})
# drop any entry that has length 0 so create_loupe does not get invalid clusters
cluster_list <- cluster_list[lengths(cluster_list) == length(cell_barcodes)]

# obsm is dict-like: use bracket access (not $X_umap)
obsm_umap <- adata$obsm[["X_umap"]]
umap_mat <- as.matrix(reticulate::py_to_r(obsm_umap))
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
rownames(umap_mat) <- cell_barcodes
projections <- list(umap = umap_mat)

# Create Loupe object
create_loupe(count_mat = count_mat, clusters = cluster_list, projections = projections, output_name = paste0(output_dir, "/loupe_obj"))