# conda activate /data/gpfs/projects/punim2435/Project/Env/IntestinalFibroblasts
# home dir is full
# /data/gpfs/projects/punim2435/Project/PRO_ACT/Env/IntestinalFibroblasts

# quota -s
# conda create \
# --clone IntestinalFibroblasts_scRNAseq \
# -p /data/gpfs/projects/punim2435/Project/Env/IntestinalFibroblasts

# mamba install r-magick --freeze-installed
# mamba install bioconductor-spatialexperiment

#########Load Packages#########
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("openxlsx"))
suppressMessages(library("GSVA"))
# scran
suppressMessages(library("scran"))
# AnnotationDbi
suppressMessages(library("AnnotationDbi"))
# org.Mm.eg.db
suppressMessages(library("org.Mm.eg.db"))
# SingleCellExperiment
suppressMessages(library("SingleCellExperiment"))
#ReactomePA
suppressMessages(library("ReactomePA"))
suppressMessages(library("clusterProfiler"))
#########Command Parameters#########
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Results/05_Pseudobulk_Ontogeny/stage_grouped/",
              help = "output directory of 05_Pseudobulk.R"
  ),
  make_option(c("--output_dir", "-o"),
              type = "character", default = "Results/07_GSVA_pseudobulk/stage_grouped",
              help = "Output directory for GSVA results."
  ),
  make_option(c("--feature_anno", "-f"),
              type = "character", default = "Data/Genomic_Features.tsv",
              help = "Path to genomic features (for gene names)."
  ),
  make_option(c("--deg_fdr"),
              type = "double", default = 0.05,
              help = "adj_pval threshold for StageGlobal_DEG (default: 0.05)"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
##########################################################
set.seed(123)
# read input RDS
fit <- readRDS(list.files(opt$input, pattern = "fit\\.rds$", full.names = TRUE)[[1]])
deg_res <- readRDS(list.files(opt$input, pattern = "res\\.rds$", full.names = TRUE)[[1]])

if (is.null(opt$output_dir)) {
  print("NO OUTPUT DIRECTORY SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output_dir
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

##########################################################
cnt <- fit[["data"]]@assays@data@listData[["counts"]]
sce <- SingleCellExperiment(assays = list(counts = cnt))
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)
normalized_counts <- logcounts(sce)

#Read the genomic feature table
genomic_features <- read.csv(opt$feature_anno, row.names = NULL, sep="\t")
genomic_features_unique <- genomic_features %>%
  dplyr::select(c(Symbol, Name)) %>%
  dplyr::distinct()
protein_coding_genes <- unique(subset(genomic_features, Gene.Type == "protein-coding")$Symbol)

normalized_counts <- subset(normalized_counts, rownames(normalized_counts) %in% protein_coding_genes)

#This is to add full names of the genes to the results
gene_metadata <- dplyr::select(genomic_features, c(Name, Symbol)) %>% 
  dplyr::rename("name" = "Symbol",
                "open_name" = "Name") %>%
  distinct()

entrezid_mapper <- AnnotationDbi::select(org.Mm.eg.db, keys=keys(org.Mm.eg.db), columns = c("SYMBOL", "ENTREZID")) %>% dplyr::rename("Symbol" = "SYMBOL")

normalized_counts_temp <- as.data.frame(normalized_counts) %>%
  dplyr::mutate(Symbol = rownames(normalized_counts))
normalized_counts_withENTREZID <- merge(normalized_counts_temp, entrezid_mapper, by="Symbol")
rownames(normalized_counts_withENTREZID) <- normalized_counts_withENTREZID$ENTREZID
normalized_counts_withENTREZID <- dplyr::select(normalized_counts_withENTREZID, -c(ENTREZID, Symbol))
rm(normalized_counts_temp)


reactome_gson <- gson_Reactome("mouse")
gobp_gson <- gson_GO(OrgDb = "org.Mm.eg.db", ont = "BP")
gomf_gson <- gson_GO(OrgDb = "org.Mm.eg.db", ont = "MF")

gsid2gene <- reactome_gson@gsid2gene
gsid2name <- reactome_gson@gsid2name
reactome_temp <- merge(gsid2gene, gsid2name, by="gsid")

gsid2gene <- gobp_gson@gsid2gene
gsid2name <- gobp_gson@gsid2name
gobp <- merge(gsid2gene, gsid2name, by="gsid")

gsid2gene <- gomf_gson@gsid2gene
gsid2name <- gomf_gson@gsid2name
gomf <- merge(gsid2gene, gsid2name, by="gsid")

##########################################################
# Build a custom gene set: StageGlobal_DEG (no UP/DOWN)
stage_deg_symbols <- deg_res %>%
  dplyr::filter(!is.na(adj_pval), adj_pval < opt$deg_fdr) %>%
  dplyr::pull(Symbol) %>%
  unique()

stage_deg_entrez <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = stage_deg_symbols,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENTREZID")
) %>%
  dplyr::filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
  dplyr::pull(ENTREZID) %>%
  as.character()

stage_deg_entrez <- intersect(stage_deg_entrez, rownames(normalized_counts_withENTREZID))

custom_name <- paste0("StageGlobal_DEG_FDR", opt$deg_fdr)
StageGlobal_DEG <- tibble::tibble(
  gsid = "CUSTOM",
  gene = stage_deg_entrez,
  name = custom_name
)

##########################################################
# Combine Reactome with custom gene set and run GSVA
reactome <- dplyr::bind_rows(
  reactome_temp %>% dplyr::mutate(gene = as.character(gene)),
  StageGlobal_DEG
)

reactome_list <- split(reactome$gene, f = reactome$name)
reactome_list <- lapply(reactome_list, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))

gobp_list <- split(as.character(gobp$gene), f = gobp$name)
gobp_list <- lapply(gobp_list, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))

gomf_list <- split(as.character(gomf$gene), f = gomf$name)
gomf_list <- lapply(gomf_list, function(x) unique(as.character(x[!is.na(x) & nzchar(x)])))


gsva_object_reactome <- gsvaParam(exprData = as.matrix(normalized_counts_withENTREZID), geneSets = reactome_list, minSize = 10)
gsva_res_reactome <- gsva(gsva_object_reactome, verbose=TRUE)

gsva_object_gobp <- gsvaParam(exprData = as.matrix(normalized_counts_withENTREZID), geneSets = gobp_list, minSize = 10)
gsva_res_gobp <- gsva(gsva_object_gobp, verbose=TRUE)

gsva_object_gomf <- gsvaParam(exprData = as.matrix(normalized_counts_withENTREZID), geneSets = gomf_list, minSize = 10)
gsva_res_gomf <- gsva(gsva_object_gomf, verbose=TRUE)

saveRDS(gsva_res_reactome, file.path(output_dir, "GSVA_reactome.rds"))
saveRDS(gsva_res_gobp, file.path(output_dir, "GSVA_gobp.rds"))
saveRDS(gsva_res_gomf, file.path(output_dir, "GSVA_gomf.rds"))

# save csv
# rownames are the pathway names
gsva_res_reactome = as.data.frame(gsva_res_reactome) %>% 
  dplyr::mutate(pathway = rownames(gsva_res_reactome)) %>% 
  dplyr::select(pathway, everything())
gsva_res_gobp = as.data.frame(gsva_res_gobp) %>%
  dplyr::mutate(pathway = rownames(gsva_res_gobp)) %>%
  dplyr::select(pathway, everything())
gsva_res_gomf = as.data.frame(gsva_res_gomf) %>%
  dplyr::mutate(pathway = rownames(gsva_res_gomf)) %>%
  dplyr::select(pathway, everything())
write.csv(gsva_res_reactome, file.path(output_dir, "GSVA_reactome.csv"), row.names = FALSE)
write.csv(gsva_res_gobp, file.path(output_dir, "GSVA_gobp.csv"), row.names = FALSE)
write.csv(gsva_res_gomf, file.path(output_dir, "GSVA_gomf.csv"), row.names = FALSE)

message("Done. GSVA results in ", output_dir)