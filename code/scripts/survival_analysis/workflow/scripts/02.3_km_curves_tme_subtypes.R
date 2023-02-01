# @created: 23 Nov 21
# @modified: 27 Jun 22
# @authors: Yoann Pradat
#
# Draw Kaplan-Meier curves where tumors are stratified according to the presence or absence of some genomic alterations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

source("workflow/functions/km_curves.R")

# functions ============================================================================================================

main <- function(args){
  # load tables
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  df_cln <- df_cln %>% filter(!is.na(Survival_Time)) %>% filter(!is.na(Sample_Id_RNA_T))
  df_tme <- load_table(args$tme_subtypes) %>% rename(Sample_Id=`...1`)

  # select tumor types
  if (tolower(args$tumor_types)!="all"){
    df_cln <- df_cln %>% filter(Tumor_Type %in% toupper(args$tumor_types))
  }

  # transform survival columns
  df_cln <- df_cln %>%
    mutate(Survival_Time=Survival_Time/30, Survival_Status=ifelse(Survival_Status=="Deceased", 1, 0)) %>%
    rename(Sample_Id=Sample_Id_DNA_T)

  # merge
  df_cln <- left_join(df_cln, df_tme, by=c("Sample_Id_RNA_T"="Sample_Id"))
  col_strata <- "Label"
  df_sub <- df_cln %>% filter(!is.na(.data[[col_strata]]))

  title <- paste("Immune subtypes", "-", paste0(args$tumor_types, collapse=" & "))
  fit <- eval(parse(text=(paste("survfit(Surv(Survival_Time, Survival_Status) ~", col_strata, ", data=df_sub)"))))
  draw_km_curves(df_sub, fit, output=args$output, width=args$width, height=args$height,
                 title=title, xlim=88, break.time.by=12, tables.height=0.3, conf.int=F, palette="lancet",
                 legend.labs=sort(unique(df_sub[[col_strata]])), pval.coord=c(12,0.95), ylab=NULL)

  # special curves for PRAD and BLCA
  if (args$tumor_types=="PRAD"){
    tme <- "D"
  } else if (args$tumor_types=="BLCA"){
    tme <- "F"
  } else {
    tme <- NULL
  }

  if (!is.null(tme)){
    df_sub <- df_sub %>% mutate(!!col_strata:=ifelse(.data[[col_strata]]==tme, tme, paste("non", tme)))
    title <- paste("Immune subtypes", "-", paste0(args$tumor_types, collapse=" & "))
    fit <- eval(parse(text=(paste("survfit(Surv(Survival_Time, Survival_Status) ~", col_strata, ", data=df_sub)"))))
    draw_km_curves(df_sub, fit, output=gsub(".pdf", paste0("_", tme, ".pdf"), args$output), width=args$width,
                   height=args$height, title=title, xlim=88, break.time.by=12, tables.height=0.25, conf.int=F,
                   palette="lancet", legend.labs=sort(unique(df_sub[[col_strata]])), pval.coord=c(12,0.95), ylab=NULL)
  }

}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--tumor_types", type="character", nargs="+", default="PRAD", 
                      help="List of tumor types to be considered.")
  parser$add_argument("--tme_subtypes", type="character", help="Path to table of TME subtypes.",
                    default="../../../results/immuno_analysis/prism/tables/mfp_subtypes_predicted_LogisticRegression.tsv")
  parser$add_argument("--width", type="double", default=8, help="Width of output plot in inches.")
  parser$add_argument("--height", type="double", default=6, help="Height of output plot in inches.")
  parser$add_argument("--output", type="character", help="Path to output folder.",
                      default="../../../results/survival_analysis/km_curves/PRAD/by_subtype/km_curves_by_tme_subtypes.pdf")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
