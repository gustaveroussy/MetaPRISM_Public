suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

get_cln_sample <- function(df_cln, sample_type){
  col_sid <- paste0("Sample_Id_", sample_type)
  cols <- c("Tumor_Type", col_sid, "Subject_Id")
  df_cln %>%
    filter(grepl(sample_type, Sample_Type)) %>%
    select(all_of(cols)) %>%
    rename(Sample_Id=.data[[col_sid]])
}


get_bio_sample <- function(df_bio, sample_type, algos){
  analyte <- unlist(strsplit(sample_type, "_"))[[1]]
  stopifnot(analyte %in% c("RNA", "DNA"))
  cols_algos <- paste(analyte, algos, sep="_")
  cols_algos <- intersect(cols_algos, colnames(df_bio))
  df_bio %>%
    filter_at(vars(all_of(cols_algos)), all_vars(.==1))
}


get_alg_sample <- function(df_cln, df_bio, sample_type, algos){
  df_cln_sam <- get_cln_sample(df_cln, sample_type)
  df_bio_sam <- get_bio_sample(df_bio, sample_type, algos)
  df_cln_sam %>% filter(Sample_Id %in% df_bio_sam$Sample_Id)
}


