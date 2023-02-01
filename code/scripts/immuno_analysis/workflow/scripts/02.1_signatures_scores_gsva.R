# @created: 26 Jul 21
# @modified: 15 Feb 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(GSVA))
suppressPackageStartupMessages(library(GSEABase))

# functions ============================================================================================================

remove_duplicates_summary <- function(df_sum, bio_data){
  has_duplicates <- length(unique(df_sum[["Sample_Id"]])) < nrow(df_sum)
  is_tcga <- sum(grepl("^TCGA", df_sum[["Col_Name"]]))==nrow(df_sum) 
  if (has_duplicates){
    if (is_tcga){
      df_bio <- read_tsv(bio_data, progress=F, show_col_types=F)
      df_bio <- df_bio %>% filter(Sample_Type %in% c("RNA_T", "RNA_N"))
      cols_rna <- colnames(df_bio)[grepl("^RNA", colnames(df_bio))]
      df_bio[["RNA_All"]] <- rowSums(df_bio[cols_rna])
      df_sum <- left_join(df_sum, df_bio[,c("Aliquot_Id", "RNA_All")], by="Aliquot_Id")
      df_sum <- df_sum %>% arrange(Sample_Id, desc(RNA_All))
      df_sum <- df_sum[!duplicated(df_sum$Sample_Id),]
      df_sum$RNA_All <- NULL
    } else {
      print("-ERROR: unexpected duplicated ids in summary table")
    }
  }
  df_sum
}


load_expression <- function(expression_data, expression_summary, bio_data){
  df_rna <- read_tsv(expression_data, progress=F, show_col_types=F)
  df_sum <- read_tsv(expression_summary, progress=F, show_col_types=F)
  df_sum <- remove_duplicates_summary(df_sum, bio_data)

  cols_keep <- c(c("ensembl_gene_id", "geneName", "geneID"), df_sum[["Col_Name"]])
  cols_keep <- intersect(cols_keep, colnames(df_rna))
  df_rna <- df_rna[,cols_keep]
  oldnames <- df_sum$Col_Name
  newnames <- df_sum$Sample_Id
  newnames <- newnames[which(oldnames %in% colnames(df_rna))]
  oldnames <- oldnames[which(oldnames %in% colnames(df_rna))]
  df_rna <- df_rna %>% 
    rename_with(~ newnames[which(oldnames == .x)], .cols = all_of(oldnames))

  df_rna
}

connect_gene_ids_to_symbols <- function(df_rna, df_genes){
  if ("ensembl_gene_id" %in% colnames(df_rna)){
    df_rna <- df_rna %>% dplyr::rename(Ensembl_Gene_Id_Version=ensembl_gene_id)
  } else {
    df_rna <- df_rna %>% dplyr::rename(Ensembl_Gene_Id_Version=geneID)
  }
  df_genes <- df_genes[,c("Ensembl_Gene_Id_Version", "Hugo_Symbol")]
  df_rna <- left_join(df_rna, df_genes, by="Ensembl_Gene_Id_Version")

  df_rna
}

select_genes <- function(df_rna, df_genes, genesets, on="Ensembl_Gene_Id_Version"){
  set_symbols <- unique(unlist(geneIds(genesets)))

  if (on=="Ensembl_Gene_Id"){
    df_rna["Ensembl_Gene_Id"] <- unlist(lapply(df_rna$Ensembl_Gene_Id_Version,
                                               function(x) unlist(strsplit(x, "\\."))[1]))
  }

  mask_id <- df_rna[[on]] %in% df_genes[[on]]
  mask_symbols <- df_rna[["Hugo_Symbol"]] %in% set_symbols
  df_rna[mask_id | mask_symbols,]
}

clear_gene_ids <- function(df_rna){
  df_rna$Ensembl_Gene_Id <- NULL
  df_rna$Ensembl_Gene_Id_Version <- NULL
  df_rna$geneID <- NULL
  df_rna$geneName <- NULL

  df_rna
}

check_genesets_symbols <- function(df_rna, genesets){
  set_symbols <- unique(unlist(geneIds(genesets)))
  rna_symbols <- unique(df_rna$Hugo_Symbol)

  if (length(setdiff(set_symbols, rna_symbols))){
    cat("-WARNING: gene symbols of some inputs signatures are not in the input expression table\n")
    cat(paste0("\t", paste(paste0("-", setdiff(set_symbols, rna_symbols)), collapse="\n\t"), "\n"))
  } else{
    cat("-INFO: all gene symbols from gene sets found in expression table\n")
  }
}

summarize_duplicated_symbols <- function(df_rna){
  cat("-summarizing duplicated gene symbols expression by mean ... ")
  mask_duplicated <- duplicated(df_rna$Hugo_Symbol)
  duplicated_symbols <- df_rna$Hugo_Symbol[mask_duplicated]
  symbols <- unique(df_rna$Hugo_Symbol)
  samples <- colnames(df_rna)[!colnames(df_rna) %in% "Hugo_Symbol"]
  rna_data <- list()
  if (length(duplicated_symbols) > 0){
    for (symbol in duplicated_symbols){
      rna_data_symbol <- list(colMeans(df_rna[df_rna$Hugo_Symbol==symbol,samples]))
      rna_data_symbol <- setNames(rna_data_symbol, symbol)
      rna_data <- c(rna_data, rna_data_symbol)
    }

    df_rna_duplicated <- data.frame(rna_data)
    df_rna_duplicated <- t(df_rna_duplicated)
    df_rna_not_duplicated <- df_rna %>% 
      filter(!Hugo_Symbol %in% duplicated_symbols) %>%
      column_to_rownames(var="Hugo_Symbol")

    df_rna <- rbind(df_rna_not_duplicated, df_rna_duplicated)
    df_rna <- df_rna[symbols,samples]
  } else {
    df_rna <- df_rna %>% column_to_rownames(var="Hugo_Symbol")
  }
  cat("ok!\n")

  df_rna
}

compute_signature_scores <- function(df_rna, genesets, args){
  cat("-running gsva function on expression table ...\n")
  df_scores <- gsva(expr=as.matrix(df_rna), gset.idx.list=genesets,
                    method=args$gsva_method, ssgsea.norm=args$gsva_ssgsea_norm, tau=args$gsva_tau)

  as.data.frame(df_scores)
}

# run ==================================================================================================================

main <- function(args){
  df_rna <- load_expression(args$expression_data, args$expression_summary, args$bio_data)
  genesets <- getGmt(args$signatures, geneIdType=SymbolIdentifier())
  df_genes_table <- read_tsv(args$gene_table, progress=F, show_col_types=F)
  df_genes_selection <- read_tsv(args$gene_selection, progress=F, show_col_types=F)

  df_rna <- connect_gene_ids_to_symbols(df_rna, df_genes_table)
  df_rna <- select_genes(df_rna, df_genes_selection, genesets=genesets, on="Ensembl_Gene_Id")
  check_genesets_symbols(df_rna, genesets)
  df_rna <- clear_gene_ids(df_rna)
  df_rna <- summarize_duplicated_symbols(df_rna)

  df_scores <- compute_signature_scores(df_rna, genesets, args)

  # save
  cat("-saving matrix of signatures scores ... ")
  dir.create(dirname(args$output), showWarnings=F, recursive=T)
  write.table(df_scores, args$output, quote=F, row.names=T, sep="\t")
  cat("ok!\n")
  cat(paste("--file saved at", args$output, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Compute signature scores using GSVA R implementation.')
  parser$add_argument('--expression_data', type="character", help="Path to input expression data table.",
                      default="../../../data/tcga/rna/kallisto-tximport/quantification_genes_10656_tpm_full.txt.gz")
  parser$add_argument('--expression_summary', type="character", help="Path to input expression summary table.",
                      default="../../../data/tcga/rna/summary/kallisto-tximport_quantification_genes_tpm_full.tsv")
  parser$add_argument('--bio_data', type="character", help="Path to biospecimen data table.",
                        default="../../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv")
  parser$add_argument('--gene_table', type="character", help='Path to gene table linking gene ids to gene symbols.',
                      default="../../../results/immuno_analysis/gene_tables/gencode_v27_updated.tsv")
  parser$add_argument('--gene_selection', type="character", help="Path to gene selection table.",
                default="../../../results/immuno_analysis/gene_tables/gencode_v27_coding_no-mt_no-histones_updated.tsv")
  parser$add_argument('--signatures', type="character",
                      default="resources/bagaev_2021/signatures/gene_signatures.gmt",
                      help='Path to gene sets file in GMT format.')
  parser$add_argument('--gsva_method', type="character", default="ssgsea",
                      help='Value for method argument of gsva function.')
  parser$add_argument('--gsva_ssgsea_norm', type="character", default="False",
                      help='Value for ssgsea.norm argument of gsva function.')
  parser$add_argument('--gsva_tau', type="double", default=0.25,
                      help='Value for tau argument of gsva function.')
  parser$add_argument('--output',type="character", help='Path to output table of signature scores.')
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)

  main(args)
}
