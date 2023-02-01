
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))

#' Draw a heatmap of immune subtypes as Figure 2.A of PMID: 34019806
#'
#' @param df_score A data.frame containing sample scores in columns. Column names should be sample names while row names
#'   should be signature names.
#' @param df_annot A data.frame containing annotations of samples. It must contain at least 1 column 'Sample_Id' with
#'   values matching the column names of df_score.
#' @param use_bagaev_params A boolean set to TRUE by default. If TRUE, the row order, row split groups and column split
#'   groups are obtained from the function `params_bagaev_signatures`.
draw_heatmap_mfp <- function(df_score, df_annot=NULL, use_bagaev_params=T, ...){
  stopifnot(setequal(colnames(df_score), df_annot$Sample_Id))
  df_annot <- df_annot %>% arrange(match(Subject_Id, colnames(df_score)))
  heatmap_legend_at <- c(floor(min(df_score)), 0, ceiling(max(df_score)))

  dots <- list(...)
  params <- list(matrix=as.matrix(df_score),
                 show_row_dend=F, 
                 show_column_dend=F,
                 show_column_names=F,
                 show_row_names=T,
                 row_names_gp=gpar(fontsize=8),
                 row_names_side="left", 
                 use_raster=T,
                 row_title_rot=0,
                 row_title_side="right", 
                 row_title_gp=gpar(fontsize=10),
                 cluster_rows=F,
                 cluster_columns=T,
                 cluster_row_slices=F,
                 cluster_column_slices=F,
                 column_title_gp=gpar(fontsize=10, fontface="bold"),
                 heatmap_legend_param=list(legend_direction="horizontal", title="z-score", at=heatmap_legend_at),
                 raster_by_magick=T)

  if (use_bagaev_params){
    params_bagaev <- params_bagaev_signatures()
    params$row_order <- params_bagaev$row_order
    df_score <- df_score[params$row_order,]
    params$matrix <- as.matrix(df_score)

    params$row_split <- params_bagaev$row_split
    params$row_title_gp <- gpar(fontsize=10, col=params_bagaev$row_split_title_col, fontface="bold")
    params$column_split <- factor(df_annot$MFP, levels=c("D", "F", "IE/F", "IE"))

    column_split_mins <- list("D"=5,"F"=5,"IE/F"=10,"IE"=10)
    for (old_name in names(params_bagaev$column_split_names)){
      new_name <- params_bagaev$column_split_names[[old_name]]
      if (sum(params$column_split==old_name)>column_split_mins[[old_name]]){
        levels(params$column_split)[levels(params$column_split)==old_name] <- new_name
      }
    }
  }

  params <- modifyList(params, dots)
  ht <- do.call('Heatmap', c(params, dots[!names(dots) %in% names(params)]))
  draw(ht, heatmap_legend_side="bottom", padding=unit(c(2, 2, 10, 2), "mm"))
}

params_bagaev_signatures <- function(){
  df_row_annot <- read_tsv("resources/bagaev_2021/signatures/gene_signatures_categories.tsv",
                           progress=F, show_col_types=F)
  df_row_annot$Category <- sapply(df_row_annot$Category, function(s) str_replace(s, ", ", ",\n"))
  column_split_names <- list("IE"="Immune-Enriched,\nNon-Fibrotic",
                             "IE/F"="Immune-Enriched,\nFibrotic",
                             "F"="Fibrotic",
                             "D"="Depleted")
  row_split_title_col <- df_row_annot %>% distinct(Category, Color) %>% pull(Color)
  row_split=factor(df_row_annot$Category, levels=df_row_annot %>% distinct(Category) %>% pull(Category))

  list(row_order=df_row_annot$Signature, row_split=row_split, row_split_title_col=row_split_title_col,
       column_split_names=column_split_names)
}
