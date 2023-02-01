suppressPackageStartupMessages(library(UltiSig))
suppressPackageStartupMessages(library(readr))


load_signature_profiles <- function(method=NULL, basis_name=NULL, denovo_name=NULL, denovo_dim=NULL){
  if (is.null(method)){
    # load published sets of signatures
    df_W <- eval(parse(text=paste0("signatures_", basis_name)))
  }  else {
    # load de novo extracted sets of signatures
    if (method=="sigprofilerjulia"){
      folder <- file.path(folder_results, method)
      subfolder <- denovo_name
      filename <- paste0("processes_", denovo_dim)
      filepath <- file.path(folder, subfolder, filename)
      df_W <- read.table(args$count_mut, sep="\t", header=T, row.names=1, as.is=T, check.names=F)
      colnames(df_W) <- sapply(1:denovo_dim, function(x) paste("signature", x, sep="_"))
      if (grepl("sbs_96", denovo_name, ignore.case=T)){
        rownames(df_W) <- mutations_order(mode="sbs_96")
      } else if (grepl("sbs_192", denovo_name, ignore.case=T)){
        rownames(df_W) <- mutations_order(mode="sbs_192")
      } else if (grepl("dbs_78", denovo_name, ignore.case=T)){
        rownames(df_W) <- mutations_order(mode="dbs_78")
      }
    }
  }

  df_W 
}

load_signature_counts <- function(method, basis_name=NULL, cohort_name=NULL, sparse_name=NULL, denovo_name=NULL,
                                  denovo_dim=NULL){
  folder <- file.path(folder_results, method)
  
  if (method=="sigprofilerjulia"){
    subfolder <- denovo_name
    filename <- paste0("exposures_", denovo_dim)
    filepath <- file.path(folder, subfolder, filename)
    df_H <- read_tsv(filepath, progress=F, show_col_types=F)
    df_H[,"signature"] <- sapply(1:denovo_dim, function(x) paste("signature", x, sep="_"))
    df_H <- df_H %>% column_to_rownames(var="signature")
  } else {
    if (is.null(sparse)){
      filename <- paste0("signature_counts_", basis_name, "_", cohort_name, ".csv")
    } else {
      filename <- paste0("signature_counts_", basis_name, "_", cohort_name, "_sparse_", sparse_name, ".csv")
    }
    filepath <- file.path(folder, filename)
    df_H <- read_tsv(filepath, progress=F, show_col_types=F)
    df_H <- df_H %>% column_to_rownames(var=names(df_H)[1])
  }

  df_H
}
