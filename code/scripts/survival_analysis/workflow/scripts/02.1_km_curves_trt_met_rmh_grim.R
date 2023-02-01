# @created: 20 Oct 21
# @modified: 07 Dec 22
# @authors: Yoann Pradat
#
# Draw kaplan-meier curves with stratification based on metastatic sites and treatments received.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

source("workflow/functions/utils.R")
source("workflow/functions/km_curves.R")

# functions ============================================================================================================

main <- function(args){
  # load tables
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  df_bio <- load_bio(args$cohort)
  df_cln <- df_cln %>% filter(!is.na(Survival_Time))

  # transform survival columns
  df_cln <- df_cln %>% mutate(Survival_Time=Survival_Time/30, Survival_Status=ifelse(Survival_Status=="Deceased", 1, 0))

  # select tumor types
  if (tolower(args$tumor_types)!="all"){
    df_cln <- df_cln %>% filter(Tumor_Type %in% toupper(args$tumor_types))
  }

  # draw plots
  cols <- c("Drugs_Before_Biopsy", "Metastatic_Sites", "RMH_Score", "GRIM_Score")
  subtitles <- c("# of drugs received", "# of metastatic sites", "RMH score", "GRIM score")

  if (args$tumor_types=="All"){
    cols <- c(cols, "Tumor_Type")
    subtitles <- c(subtitles, "Tumor types")
  }

  for (i in 1:length(cols)){
    col_strata <- cols[[i]]

    # filter out patients with NA
    if (col_strata=="Drugs_Before_Biopsy"){
      df_sub <- df_cln %>% filter(!is.na(Systematic_Treatment_Before_Biopsy))
    } else {
      df_sub <- df_cln %>% filter(!is.na(.data[[col_strata]]))
    }
    
    if (col_strata=="Tumor_Type"){
      tt_keep <- c("LUAD", "PRAD", "BRCA", "PAAD", "BLCA")
      df_sub_a <- df_sub %>% filter(.data[[col_strata]]  %in% tt_keep)
      df_sub_b <- df_sub %>% mutate(Tumor_Type="All")
      df_sub <-  bind_rows(df_sub_a, df_sub_b)
    }
    
    if (!grepl("Score|Tumor_Type", col_strata)){
      df_sub <- add_col_count(df_sub, col_strata)
      col_strata <- paste("Count", col_strata, sep="_")
    }

    # aggregate strata that are too small
    if (!grepl("Tumor_Type", col_strata)){
      df_sub <- aggregate_small_strata(df_sub, col_strata, max_strata=args$max_strata)
      tables.height <- 0.3
      palette <- "lancet"
    } else {
      # ordered alphabetically All BLCA BRCA LUAD PAAD PRAD 
      colors_gl <- load_colors(sheet="Global")
      colors_tt <- load_colors(sheet="Project_TCGA_More_Clear")
      palette <- c(setNames(colors_gl[["prism"]], "All"), unlist(colors_tt[sort(tt_keep)]))
      tables.height <- 0.35
    }

    # # risk table height
    # if (grepl("Score", col_strata)){
    #   tables.height <- 0.25
    # } else {
    #   n_strata <- length(unique(df_sub[[col_strata]]))
    #   tables.height <- 0.1 + (n_strata/8)*0.2
    # }

    title <- subtitles[[i]]
    if (title=="Tumor types") title <- NULL
    fit <- eval(parse(text=(paste("survfit(Surv(Survival_Time, Survival_Status) ~", col_strata, ", data=df_sub)"))))
    draw_km_curves(df_sub, fit, output=args$outputs[[i]], width=args$widths[[i]], height=args$heights[[i]],
                   title=title, xlim=88, break.time.by=12, tables.height=tables.height,  break.y.by=0.5,
                   palette=palette, legend.labs=sort(unique(df_sub[[col_strata]])), pval.coord=c(12,0.95), fontsize=3.5)
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Kaplan-Meier curves for counts of metastases and drugs.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument("--tumor_types", type="character", nargs="+", default="All", 
                      help="List of tumor types to be considered.")
  parser$add_argument("--max_strata", type="integer", default=4,
                      help="Maximum number of strata. Small-size strata are grouped with nearby strata.")
  parser$add_argument("--widths", type="double", nargs="+", default=c(5,5,5,5,5),
                      help="Width of output plots in inches.")
  parser$add_argument("--heights", type="double", nargs="+", default=c(3,3,3,3,3.75),
                      help="Height of output plots in inches.")
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plots.",
    default=c("../../../results/survival_analysis/km_curves/All/trt_met_rmh_grim/km_curves_trt_All.pdf",
              "../../../results/survival_analysis/km_curves/All/trt_met_rmh_grim/km_curves_met_All.pdf",
              "../../../results/survival_analysis/km_curves/All/trt_met_rmh_grim/km_curves_rmh_All.pdf",
              "../../../results/survival_analysis/km_curves/All/trt_met_rmh_grim/km_curves_grim_All.pdf",
              "../../../results/figures_paper/F1e.pdf"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  for (name in names(args)){
    if (!is.null(args[[name]]) && length(args[[name]])==1){
      if (args[[name]]=="None"){
        args[[name]] <- NULL
      }
    }
  }

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
