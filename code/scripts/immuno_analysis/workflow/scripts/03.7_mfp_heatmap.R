# @created: 06 Aug 21
# @modified: 15 Feb 22
# @authors: Yoann Pradat
#
# Using computed signature scores and predicted subtypes, draw a heatmap showing the values of each signature score for 
# each sample. Samples are grouped by subtype.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))

source("workflow/functions/drawer.R")

# functions ============================================================================================================

load_annotation <- function(input_cln, input_bio=NULL){
  if (!is.null(input_bio)){
    df_cln <- read_tsv(input_cln, progress=F, show_col_types=F)
    df_bio <- read_tsv(input_bio, progress=F, show_col_types=F)

    cols_cln <- c("Subject_Id", "Sample_Id_RNA_T", "Project_TCGA_More", "Date_Death", "Age_At_Biopsy",
                  "Histological_Type", "Survival_Status", "OS")
    cols_cln <- intersect(names(df_cln), cols_cln)
    cols_bio <- c("Sample_Id", "Biopsy_Type", "Biopsy_Site", "Biopsy_Subsite")
    cols_bio <- intersect(names(df_bio), cols_bio)
    df_ann_cln <- df_cln %>% select(all_of(cols_cln)) %>%
      rename(Sample_Id=Sample_Id_RNA_T) %>% filter(!is.na(Sample_Id)) %>% distinct()
    df_ann_bio <- df_bio %>% select(all_of(cols_bio)) %>% distinct()
    df_ann <- left_join(df_ann_cln, df_ann_bio, by="Sample_Id")
    if ("Survival_Status" %in% names(df_cln)){
      df_ann <- df_ann %>% mutate(Alive=ifelse(Survival_Status=="Deceased", "No", "Yes"))
    }
    if ("OS" %in% names(df_cln)){
      df_ann <- df_ann %>% mutate(Alive=ifelse(OS==1, "No", "Yes"))
    }
  } else {
    df_ann = load_table(input_cln)
  }

  df_ann
}

add_subtypes <- function(df_ann, input_subtypes){
  if (!"MFP" %in% names(df_ann)){
    df_mfp <- read_tsv(input_subtypes, progress=F, show_col_types=F) %>% rename(Sample_Id=`...1`)
    df_mfp <- df_mfp %>% rename(MFP=Label)
    df_ann <- left_join(df_mfp, df_ann, on="Sample_Id")
  }

  df_ann
}

standardize_colnames <- function(df_ann){
  if ("Project_TCGA_More" %in% names(df_ann)){
    df_ann <- df_ann %>% rename(Tumor_Type=Project_TCGA_More)
  } else if ("TCGA_project" %in% names(df_ann)){
    df_ann <- df_ann %>% rename(Tumor_Type=TCGA_project)
  }

  df_ann
}

get_top_annotation <- function(df_annot){
  annots <- list()
  col <- list()

  if ("Alive" %in% names(df_annot)){
    annots[["Alive"]] <- df_annot$Alive
    col[["Alive"]] <- c("No"="#5e15cd", "Yes"="#64e6a6")
  }

  if ("Biopsy_Site" %in% names(df_annot)){
    annots[["Biopsy_Site"]] <- df_annot$Biopsy_Site
    col[["Biopsy_Site"]] <- unlist(load_colors("Biopsy_Site"))
  }


  if ("Proba_D" %in% names(df_annot)){
    probas <- df_annot[,c("Proba_D", "Proba_F", "Proba_IE", "Proba_IE/F" )]
    annots[["Probas"]] <- anno_barplot(probas, axis=F, gp=gpar(fill = 2:5, lwd=0), bar_width=1, border=F)
  }

  params <- list(col=col, simple_anno_size=unit(3, "mm"))
  do.call("HeatmapAnnotation", c(annots, params))
}


draw_heatmap <- function(df_score, df_annot, filepath){
  pdf(
      file    = filepath,
      width   = max(8, min(20, 4 + 0.10*ncol(df_score))),
      height  = 6,
      onefile = FALSE
  )
  plot.new()

  if (ncol(df_score) < 20){
    row_title_gp <- gpar(fontsize=8)
    column_title_gp <- gpar(fontsize=6, fontface="bold")
  } else if (ncol(df_score) < 40) {
    row_title_gp <- gpar(fontsize=8)
    column_title_gp <- gpar(fontsize=8, fontface="bold")
  } else {
    row_title_gp <- gpar(fontsize=10)
    column_title_gp <- gpar(fontsize=10, fontface="bold")
  }

  top_annotation <- get_top_annotation(df_annot)
  draw_heatmap_mfp(df_score=df_score, df_annot=df_annot, use_bagaev_params=T, use_raster=F,
                   row_title_gp=row_title_gp, column_title_gp=column_title_gp,
                   top_annotation=top_annotation)
  dev.off()
}


main <- function(args){
  df_score_all <- read_tsv(args$input_signatures, progress=F, show_col_types=F) %>% rename(Signature=`...1`)
  df_annot_all <- load_annotation(args$input_cln, args$input_bio)
  df_annot_all <- add_subtypes(df_annot_all, args$input_subtypes)
  df_annot_all <- standardize_colnames(df_annot_all)
  dir.create(args$output, showWarnings=F, recursive=T)

  for (tumor_type in unique(df_annot_all$Tumor_Type)){
    cat(paste("-drawing heatmap for tumor type", tumor_type, "...\n"))
    df_annot <- df_annot_all %>% filter(Tumor_Type==tumor_type) %>% arrange(MFP)
    df_score <- df_score_all %>% select(all_of(c("Signature", df_annot$Sample_Id))) %>%
      column_to_rownames(var="Signature")

    filepath <- file.path(args$output, paste0("heatmap_mfp_subtypes_", tumor_type, ".pdf"))
    draw_heatmap(df_score, df_annot, filepath)
    cat("ok!\n")
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw heatmap showing matrix of signature scores organized by immune subtypes.')
  parser$add_argument('--input_signatures', type="character",
                      default="../../../results/immuno_analysis/tcga/tables/signatures_mfp_model_preprocessed.tsv",
                      help=paste("Path to input signature scores table."))
  parser$add_argument('--input_cln', type="character", help="Table of samples annotations.",
                      default="../../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv")
  parser$add_argument('--input_bio', type="character", help="Table of samples annotations.",
                      default="../../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv")
  parser$add_argument('--input_subtypes', type="character", help='Path to table of predicted subtypes.',
                    default="../../../results/immuno_analysis/tcga/tables/mfp_subtypes_predicted_LogisticRegression.tsv")
  parser$add_argument('--output', default="../../../results/immuno_analysis/tcga/plots",
                      type="character", help='Path to output folder.')
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)

  main(args)
}
