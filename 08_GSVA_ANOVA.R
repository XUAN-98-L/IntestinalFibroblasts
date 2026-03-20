# conda activate /data/gpfs/projects/punim2435/Project/Env/IntestinalFibroblasts
#########Load Packages#########
suppressMessages(library("optparse"))
suppressMessages(library("dplyr"))
suppressMessages(library("tidyverse"))
suppressMessages(library("broom"))
suppressMessages(library("openxlsx"))
suppressMessages(library("reshape2"))
########Functions#########
# Function to assign shorter group names based on keywords
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


# Function to create an ordering key for stage levels.
# If early/peak/repair appear in the data: order = uninfected_adult, early, peak, repair.
# Otherwise (ontogeny): order = embryonic, shortly_after_birth, after_weaning, uninfected_adult.
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

# Melt GSVA matrix, add SampleGroup/PlotGroup/identity_simplified, order SampleGroup by stage.
# gene_set_prefix: optional prefix to strip from GeneSet (e.g. "REACTOME_" for Reactome).
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

run_twoway_anova <- function(melted) {
  melted %>%
    dplyr::group_by(GeneSet) %>%
    dplyr::summarize(
      aov_result = list(aov(GSVAscore ~ identity_simplified * PlotGroup, data = cur_data())),
      aov_tidy = list(broom::tidy(aov(GSVAscore ~ identity_simplified * PlotGroup, data = cur_data()))),
      p_value_identity = aov_tidy[[1]] %>% dplyr::filter(term == "identity_simplified") %>% dplyr::pull(p.value),
      p_value_stage = aov_tidy[[1]] %>% dplyr::filter(term == "PlotGroup") %>% dplyr::pull(p.value),
      p_value_interaction = aov_tidy[[1]] %>% dplyr::filter(term == "identity_simplified:PlotGroup") %>% dplyr::pull(p.value),
      .groups = "drop"
    ) %>%
    dplyr::select(-aov_result, -aov_tidy) %>%
    dplyr::mutate(GeneSet = stringr::str_replace_all(GeneSet, "\\.gsid", ""))
}

run_oneway_anova <- function(melted) {
  melted %>%
    dplyr::group_by(GeneSet, identity_simplified) %>%
    dplyr::summarize(
      aov_result = list(aov(GSVAscore ~ SampleGroup, data = cur_data())),
      aov_tidy = list(broom::tidy(aov(GSVAscore ~ SampleGroup, data = cur_data()))),
      p_value_CellTypeVar = aov_tidy[[1]] %>% dplyr::filter(term == "SampleGroup") %>% dplyr::pull(p.value),
      .groups = "drop"
    ) %>%
    dplyr::select(-aov_result, -aov_tidy) %>%
    dplyr::mutate(GeneSet = stringr::str_replace_all(GeneSet, "\\.gsid", ""))
}

run_oneway_anova_single_identity <- function(melted) {
  melted %>%
    dplyr::group_by(GeneSet) %>%
    dplyr::summarize(
      aov_tidy = list(broom::tidy(aov(GSVAscore ~ SampleGroup, data = cur_data()))),
      p_value_SampleGroup = aov_tidy[[1]] %>% dplyr::filter(term == "SampleGroup") %>% dplyr::pull(p.value),
      .groups = "drop"
    ) %>%
    dplyr::select(-aov_tidy) %>%
    dplyr::mutate(GeneSet = stringr::str_replace_all(GeneSet, "\\.gsid", ""))
}
#=========================Command Parameters===================================
option_list <- list(
  make_option(c("--input", "-i"),
              type = "character", default = "Results/07_GSVA_pseudobulk/stage_grouped/",
              help = "output directory of 07_GSVA_pseudobulk.R"
  ),
  make_option(c("--output_dir", "-o"),
              type = "character", default = "Results/08_GSVA_ANOVA/stage_grouped",
              help = "Output directory for GSVA ANOVA results."
  ),
  make_option(c("--subset", "-s"),
              type = "logical", default = FALSE,
              help = "whether the GSVA is a single cell type subset. default is FALSE."
  ),
  make_option(c("--top_n", "-t"),
              type = "integer", default = 10,
              help = "top n pathways to plot. default is 10."
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
##########################################################
if (is.null(opt$output_dir)) {
  print("NO OUTPUT DIRECTORY SUPPLIED,current directory will be used!")
  output_dir <- getwd()
} else {
  output_dir <- opt$output_dir
  if (!file.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }
}

# read input RDS
gsva_res_reactome <- readRDS(paste0(opt$input, "/GSVA_reactome.rds"))
gsva_res_gobp <- readRDS(paste0(opt$input, "/GSVA_gobp.rds"))
gsva_res_gomf <- readRDS(paste0(opt$input, "/GSVA_gomf.rds"))
##########################################################
# Prepare melted data frames for ANOVA
gsva_res_reactome_melted <- prepare_gsva_melted(gsva_res_reactome, gene_set_prefix = "REACTOME_")
saveRDS(gsva_res_reactome_melted, paste0(output_dir, "/gsva_res_reactome_melted.rds"))

gsva_res_gobp_melted <- prepare_gsva_melted(gsva_res_gobp)
saveRDS(gsva_res_gobp_melted, paste0(output_dir, "/gsva_res_gobp_melted.rds"))

gsva_res_gomf_melted <- prepare_gsva_melted(gsva_res_gomf)
saveRDS(gsva_res_gomf_melted, paste0(output_dir, "/gsva_res_gomf_melted.rds"))


### Perform ANOVA analysis ###
# If subset = TRUE: single cell type; one-way ANOVA per pathway (GSVA ~ SampleGroup) to test whether pathway activity varies across stage within that cell type.
# If subset = FALSE: multiple cell types; two-way ANOVA per pathway (GSVA ~ identity * stage) to test main effects and interaction (pathways that change with stage and/or differ by cell type).

if (opt$subset) {
  message("Mode: one-way ANOVA per pathway — pathway activity vs SampleGroup (stage) within the subset cell type.")
  anova_oneway_reactome <- run_oneway_anova_single_identity(gsva_res_reactome_melted)
  anova_oneway_gobp    <- run_oneway_anova_single_identity(gsva_res_gobp_melted)
  anova_oneway_gomf    <- run_oneway_anova_single_identity(gsva_res_gomf_melted)

  wb <- createWorkbook()
  addWorksheet(wb, sheetName = "reactome")
  writeData(wb, sheet = "reactome", x = anova_oneway_reactome)
  addWorksheet(wb, sheetName = "gobp")
  writeData(wb, sheet = "gobp", x = anova_oneway_gobp)
  addWorksheet(wb, sheetName = "gomf")
  writeData(wb, sheet = "gomf", x = anova_oneway_gomf)
  saveWorkbook(wb, file.path(output_dir, "ANOVA_oneway_subset.xlsx"), overwrite = TRUE)
  message("Saved ", file.path(output_dir, "ANOVA_oneway_subset.xlsx"))
} else if (opt$subset == FALSE) {
  message("Mode: two-way ANOVA per pathway — main effects of identity and stage plus interaction.")
  anova_twoway_reactome <- run_twoway_anova(gsva_res_reactome_melted) %>% dplyr::arrange(p_value_interaction)
  anova_twoway_gobp    <- run_twoway_anova(gsva_res_gobp_melted) %>% dplyr::arrange(p_value_interaction)
  anova_twoway_gomf    <- run_twoway_anova(gsva_res_gomf_melted) %>% dplyr::arrange(p_value_interaction)

  wb_twoway <- createWorkbook()
  addWorksheet(wb_twoway, sheetName = "reactome")
  writeData(wb_twoway, sheet = "reactome", x = anova_twoway_reactome)
  addWorksheet(wb_twoway, sheetName = "gobp")
  writeData(wb_twoway, sheet = "gobp", x = anova_twoway_gobp)
  addWorksheet(wb_twoway, sheetName = "gomf")
  writeData(wb_twoway, sheet = "gomf", x = anova_twoway_gomf)
  saveWorkbook(wb_twoway, file.path(output_dir, "ANOVA_TwoWay.xlsx"), overwrite = TRUE)
  message("Saved ", file.path(output_dir, "ANOVA_TwoWay.xlsx"))
}

### GSVA PLOTS ###
plot_dir <- file.path(output_dir, "GSVA_plots_reactome")
if (!file.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

sanitize_filename <- function(filename) {
  gsub("[^A-Za-z0-9_]", "", gsub(" ", "_", filename))
}

# Choose plotting order + colors depending on whether infection stages exist
has_infection <- any(gsva_res_reactome_melted$PlotGroup %in% c("early", "peak", "repair"))
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

# sort anova_twoway_reactome based on p_value_stage
# top n pathways that are significantly different between stages
anova_twoway_reactome_topn <- anova_twoway_reactome %>% arrange(p_value_stage) %>% head(opt$top_n)
anova_twoway_reactome_topn <- c(anova_twoway_reactome_topn$GeneSet, "Extracellular matrix organization")

gene_sets <- unique(anova_twoway_reactome_topn)
for (gene_set in gene_sets) {
  tryCatch({
    plotdata <- subset(gsva_res_reactome_melted, GeneSet == gene_set)
    plotdata$PlotGroup <- factor(plotdata$PlotGroup, levels = plot_group_levels)

    mean_data <- plotdata %>%
      group_by(identity_simplified, PlotGroup) %>%
      summarize(mean_GSVAscore = mean(GSVAscore, na.rm = TRUE), .groups = "drop")

    output_file <- file.path(
      plot_dir,
      paste0("GSVA_reactome_", sanitize_filename(gene_set), ".pdf")
    )

    pdf(output_file, width = 4, height = 6)
    print(
      ggplot(plotdata, aes(y = GSVAscore, x = PlotGroup, color = PlotGroup)) +
        geom_point() +
        geom_line(
          data = mean_data,
          aes(y = mean_GSVAscore, group = 1),
          color = "black"
        ) +
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
        facet_wrap(~identity_simplified, ncol = 4, nrow = 5) +
        ggtitle(gene_set) +
        guides(color = "none")
    )
    dev.off()
  }, error = function(e) {
    message("Error generating GSVA plot for gene set: ", gene_set)
    message("Error message: ", e$message)
    try(dev.off(), silent = TRUE)
  })
}

message("Done. GSVA plots in ", plot_dir)

message("Done. ANOVA results in ", output_dir)