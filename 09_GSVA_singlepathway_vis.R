# conda activate /data/gpfs/projects/punim2435/Project/Env/IntestinalFibroblasts
######### Load Packages #########
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("reshape2"))
suppressMessages(library("SummarizedExperiment"))

######### Helper functions (consistent with 08_GSVA_ANOVA.R) #########
short_group_name <- function(x) {
  case_when(
    str_detect(x, "uninfected") ~ "uninfected",
    str_detect(x, "early") ~ "early",
    str_detect(x, "peak") ~ "peak",
    str_detect(x, "repair") ~ "repair",
    str_detect(x, "after_weaning") ~ "after weaning",
    str_detect(x, "shortly_after_birth") ~ "shortly after birth",
    str_detect(x, "embryonic") ~ "embryonic",
    TRUE ~ "other"
  )
}

ordering_factor <- function(x) {
  has_infection <- any(str_detect(x, "early|peak|repair"))
  if (has_infection) {
    case_when(
      str_detect(x, "uninfected") ~ 1L,
      str_detect(x, "early") ~ 2L,
      str_detect(x, "peak") ~ 3L,
      str_detect(x, "repair") ~ 4L,
      TRUE ~ 5L
    )
  } else {
    case_when(
      str_detect(x, "embryonic") ~ 1L,
      str_detect(x, "shortly_after_birth") ~ 2L,
      str_detect(x, "after_weaning") ~ 3L,
      str_detect(x, "uninfected") ~ 4L,
      TRUE ~ 5L
    )
  }
}

prepare_gsva_melted <- function(gsva_matrix, gene_set_prefix = NULL) {
  out <- reshape2::melt(gsva_matrix,
                        varnames = c("GeneSet", "SampleName"),
                        value.name = "GSVAscore") %>%
    dplyr::mutate(
      SampleGroup = str_remove_all(SampleName, "_Rep.+"),
      GeneSet = str_remove_all(GeneSet, ".gene"),
      PlotGroup = short_group_name(SampleGroup),
      identity_simplified = str_extract(SampleGroup, "^[^_]+")
    )
  if (!is.null(gene_set_prefix) && nzchar(gene_set_prefix)) {
    out <- dplyr::mutate(out, GeneSet = str_remove_all(GeneSet, fixed(gene_set_prefix)))
  }
  out %>%
    dplyr::mutate(
      OrderFactor = ordering_factor(SampleGroup),
      SampleGroup = factor(SampleGroup, levels = unique(SampleGroup[order(OrderFactor)]))
    ) %>%
    dplyr::select(-OrderFactor)
}

sanitize_filename <- function(filename) {
  gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", filename))
}

######### Command Parameters #########
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character",
              default = "Results/07_GSVA_pseudobulk/Parasite_stage_grouped/",
              help = "Directory containing GSVA_<method>.rds (e.g. output of 07_GSVA_pseudobulk.R)."),
  make_option(c("--output_dir", "-o"),
              type = "character",
              default = "Results/09_GSVA_singlepathway_vis/Parasite_stage_grouped",
              help = "Output directory for single-pathway GSVA visualization."),
  make_option(c("--pathway", "-p"),
              type = "character",
              default = "Extracellular matrix organization",
              help = "Pathway/GeneSet name to plot (must match GeneSet after cleanup)."),
  make_option(c("--method", "-m"),
              type = "character",
              default = "reactome",
              help = "Which GSVA matrix to plot: reactome | gobp | gomf."),
  make_option(c("--facet_ncol"),
              type = "integer",
              default = 2,
              help = "Number of columns in facet_wrap (default: 2)."),
  make_option(c("--facet_nrow"),
              type = "integer",
              default = 5,
              help = "Number of rows in facet_wrap (default: 5).")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

input_dir <- opt$input
output_dir <- opt$output_dir
if (!file.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

gsva_path <- file.path(input_dir, paste0("GSVA_", opt$method, ".rds"))
if (!file.exists(gsva_path)) {
  stop("Cannot find GSVA RDS: ", gsva_path)
}

gsva_res <- readRDS(gsva_path)
gsva_res <- SummarizedExperiment::assay(gsva_res)

# Melt + annotate (Reactome needs prefix stripping like in 08_GSVA_ANOVA.R)
gene_set_prefix <- if (tolower(opt$method) == "reactome") "REACTOME_" else NULL
gsva_melted <- prepare_gsva_melted(gsva_res, gene_set_prefix = gene_set_prefix)

pathway <- opt$pathway
plotdata <- gsva_melted %>% dplyr::filter(GeneSet == pathway)
if (nrow(plotdata) == 0) {
  available <- sort(unique(gsva_melted$GeneSet))
  hint <- paste(utils::head(available, 30), collapse = "; ")
  stop("Pathway not found after cleanup: '", pathway, "'.\n",
       "Example available GeneSet values (first 30): ", hint)
}

# Order + colors (infection vs ontogeny)
has_infection <- any(plotdata$PlotGroup %in% c("early", "peak", "repair"))
if (has_infection) {
  plot_group_levels <- c("uninfected", "early", "peak", "repair")
  group_colors <- c(
    "uninfected" = "#6baed6",
    "early" = "#fd8d3c",
    "peak" = "#fb6a4a",
    "repair" = "#74c476"
  )
} else {
  plot_group_levels <- c("embryonic", "shortly after birth", "after weaning", "uninfected")
  group_colors <- c(
    "embryonic" = "#9ecae1",
    "shortly after birth" = "#fdae6b",
    "after weaning" = "#fb9a99",
    "uninfected" = "#74c476"
  )
}
plotdata$PlotGroup <- factor(plotdata$PlotGroup, levels = plot_group_levels)

mean_data <- plotdata %>%
  group_by(identity_simplified, PlotGroup) %>%
  summarize(mean_GSVAscore = mean(GSVAscore, na.rm = TRUE), .groups = "drop")

plot_file <- file.path(
  output_dir,
  paste0("GSVA_", opt$method, "_", sanitize_filename(pathway), ".pdf")
)

#pdf(plot_file, width = 2, height = 6)
pdf(plot_file, width = opt$facet_ncol * 2, height = opt$facet_nrow * 2)
print(
  ggplot(plotdata, aes(y = GSVAscore, x = PlotGroup, color = PlotGroup)) +
    geom_point() +
    geom_line(data = mean_data, aes(y = mean_GSVAscore, group = 1), color = "black") +
    scale_color_manual(values = group_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(hjust = 1, size = 5),
      strip.text = element_text(size = 4),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10),
      plot.title = element_text(size = 5)
    ) +
    facet_wrap(~identity_simplified, ncol = opt$facet_ncol, nrow = opt$facet_nrow) +
    ggtitle(pathway) +
    guides(color = "none")
)
dev.off()

# Save the plotted data too
saveRDS(plotdata, file.path(output_dir, paste0("plotdata_", opt$method, "_", sanitize_filename(pathway), ".rds")))
write.csv(plotdata, file.path(output_dir, paste0("plotdata_", opt$method, "_", sanitize_filename(pathway), ".csv")), row.names = FALSE)

message("Done. Plot saved to: ", plot_file)
