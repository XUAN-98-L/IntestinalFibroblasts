# conda activate /data/gpfs/projects/punim2435/Project/Env/IntestinalFibroblasts
# Downstream of 05_Pseudobulk.R (cell-type subset, e.g. MF_only):
# Pairwise stage comparisons: shortly_after_birth vs embryonic,
#   after_weaning vs shortly_after_birth, uninfected_adult vs after_weaning
#=========================Loading Packages===================================
suppressMessages(library("optparse"))
suppressMessages(library("glmGamPoi"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("openxlsx"))
#=========================Command Parameters===================================
option_list <- list(
  make_option(c("--input_dir", "-i"),
              type = "character", default = "Results/05_Pseudobulk_Ontogeny/MF_only",
              help = "Directory with final_pseudobulk_df.rds and pseudobulk_coldata_merged.rds from 05_Pseudobulk.R"
  ),
  make_option(c("--output_dir", "-o"),
              type = "character", default = "Results/05_Pseudobulk_Ontogeny/MF_only_pairwise_stage",
              help = "Output directory for pairwise DEX results."
  ),
  make_option(c("--feature_anno", "-f"),
              type = "character", default = "Data/Genomic_Features.tsv",
              help = "Path to genomic features (for gene names)."
  ),
  make_option(c("--stage_order"),
              type = "character",
              default = "embryonic,shortly_after_birth,after_weaning,uninfected_adult",
              help = "Comma-separated stage order (first = reference). Default: embryonic,shortly_after_birth,after_weaning,uninfected_adult"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir   <- opt$input_dir
output_dir  <- opt$output_dir
feature_anno <- opt$feature_anno
stage_order <- trimws(strsplit(opt$stage_order, ",")[[1]])

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Pairwise comparisons: numerator vs denominator (lfc = numerator - denominator)
comparisons <- list(
  shortly_after_birth_vs_embryonic     = c("shortly_after_birth", "embryonic"),
  after_weaning_vs_shortly_after_birth = c("after_weaning", "shortly_after_birth"),
  uninfected_adult_vs_after_weaning     = c("uninfected_adult", "after_weaning")
)

#=========================Load data===========================================
final_pseudobulk_df    <- readRDS(file.path(input_dir, "final_pseudobulk_df.rds"))
pseudobulk_coldata    <- readRDS(file.path(input_dir, "pseudobulk_coldata_merged.rds"))

genomic_features      <- read.csv(feature_anno, row.names = NULL, sep = "\t")
genomic_features_unique <- genomic_features %>%
  dplyr::select(Symbol, Name) %>%
  dplyr::distinct()

# Ensure stage_grouped is a factor with desired order (first = reference)
pseudobulk_coldata$stage_grouped <- factor(
  as.character(pseudobulk_coldata$stage_grouped),
  levels = stage_order
)
# Drop unused levels if any stage is missing
pseudobulk_coldata$stage_grouped <- droplevels(pseudobulk_coldata$stage_grouped)

# Align count matrix columns to coldata
stopifnot(all(colnames(final_pseudobulk_df) %in% rownames(pseudobulk_coldata)))
pseudobulk_coldata <- pseudobulk_coldata[colnames(final_pseudobulk_df), , drop = FALSE]
stopifnot(identical(colnames(final_pseudobulk_df), rownames(pseudobulk_coldata)))

#=========================Fit full model once=================================
# design ~1+stage_grouped gives coefs: (Intercept), stage_grouped<level2>, ...
fit <- glm_gp(as.matrix(final_pseudobulk_df),
              col_data = pseudobulk_coldata,
              design = ~ 1 + stage_grouped,
              size_factors = "deconvolution",
              subsample = FALSE,
              verbose = TRUE,
              on_disk = FALSE)

# Coefficient names: first level of factor is REFERENCE and has no column.
# So we have Intercept + stage_groupedshortly_after_birth, stage_groupedafter_weaning, ...
reference_level <- stage_order[1]  # e.g. "embryonic" — no coefficient, it's the baseline
coef_names <- colnames(fit$Beta)
stage_coefs <- coef_names[!coef_names %in% c("(Intercept)", "Intercept")]
prefix <- "stage_grouped"
coef_to_level <- setNames(sub(paste0("^", prefix), "", stage_coefs), stage_coefs)

#=========================Helper: build contrast for A vs B===================
# Build a numeric contrast vector: lfc = numerator - denominator.
# fit$Beta has columns in order coef_names; contrast vector same length, 1 at num, -1 at den, 0 else.
# Returns list(contrast = numeric vector for test_de, label = "num vs den").
get_contrast <- function(numerator, denominator) {
  num_coef <- names(coef_to_level)[match(numerator, coef_to_level)]
  if (denominator == reference_level) {
    # A vs reference: contrast = 1 at A's column, 0 elsewhere (omit reference column)
    if (is.na(num_coef)) return(NULL)
    vec <- setNames(rep(0, length(coef_names)), coef_names)
    vec[num_coef] <- 1
    return(list(contrast = unname(vec), label = paste0(numerator, " vs ", denominator)))
  }
  den_coef <- names(coef_to_level)[match(denominator, coef_to_level)]
  if (is.na(num_coef) || is.na(den_coef)) return(NULL)
  # A - B: +1 at num, -1 at den
  vec <- setNames(rep(0, length(coef_names)), coef_names)
  vec[num_coef] <- 1
  vec[den_coef] <- -1
  return(list(contrast = unname(vec), label = paste0(numerator, " vs ", denominator)))
}

#=========================Run pairwise test_de and write Excel================
wb <- createWorkbook()

for (comp_name in names(comparisons)) {
  pair <- comparisons[[comp_name]]
  num <- pair[1]
  den <- pair[2]
  cnt <- get_contrast(num, den)
  if (is.null(cnt)) {
    message("Skipping ", comp_name, ": missing level (levels: ", paste(levels(pseudobulk_coldata$stage_grouped), collapse = ", "), "; reference = ", reference_level, ")")
    next
  }
  message("Testing contrast: ", comp_name, " -> ", cnt$label)

  res <- test_de(fit,
                 contrast = cnt$contrast,
                 max_lfc = Inf,
                 verbose = TRUE) %>%
    dplyr::select(-c(df1, df2)) %>%
    dplyr::rename("Symbol" = "name") %>%
    merge(genomic_features_unique, by = "Symbol", all.x = TRUE) %>%
    dplyr::mutate(
      pval    = ifelse(pval == 0, .Machine$double.xmin, pval),
      adj_pval = ifelse(adj_pval == 0, .Machine$double.xmin, adj_pval)
    ) %>%
    dplyr::arrange(desc(f_statistic))

  sheet_name <- substr(comp_name, 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  writeData(wb, sheet = sheet_name, x = res)

  saveRDS(res, file.path(output_dir, paste0("DEX_", comp_name, "_res.rds")))
}

saveWorkbook(wb, file.path(output_dir, "DEX_pairwise_stage.xlsx"), overwrite = TRUE)
saveRDS(fit, file.path(output_dir, "DEX_pairwise_stage_fit.rds"))
message("Done. Results in ", output_dir)
