# @created: 21 Feb 22
# @modified: 21 Feb 22
# @authors: Yoann Pradat
#
# Annotate a table of fusions using fusions lists from different sources and using gene annotations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

load_fusions_chimerkb <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "chimerkb_4.0.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% rename(Gene_1=H_gene, Chr_1=H_chr, Strand_1=H_strand, Breakpoint_1=H_position,
                      Gene_2=T_Gene, Chr_2=T_chr, Strand_2=T_strand, Breakpoint_2=T_position) %>%
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      mutate(Chr_1=gsub("^chr", "", Chr_1), Chr_2=gsub("^chr", "", Chr_2)) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_chimerseq_normal <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "chimerseq_normal_4.0_curated.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% rename(Gene_1=H_gene, Chr_1=H_chr, Strand_1=H_strand, Breakpoint_1=H_position,
                      Gene_2=T_gene, Chr_2=T_chr, Strand_2=T_strand, Breakpoint_2=T_position) %>%
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      mutate(Chr_1=gsub("^chr", "", Chr_1), Chr_2=gsub("^chr", "", Chr_2)) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_chitars_cancer <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "chitars_5.0/chitar_cancer_fusion.tsv")
  cols_keep <- c("Fusion_Id")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_chitars_all_human <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "chitars_5.0/all_human_ChiTaRS_coord.csv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% rename(Gene_1=gene1, Chr_1=chromosome1, Strand_1=strand1, Breakpoint_1=`genomic-stop1(eloc)`,
                      Gene_2=gene2, Chr_2=chromosome2, Strand_2=strand2, Breakpoint_2=`genomic-start2(sloc)`) %>%
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_cosmic <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "cosmic_V95_fusions_curated_enhanced.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_chitars_tic <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "Chitars-TIC.tsv")
  cols_keep <- c("Fusion_Id")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_tic <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "tic_3.3.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Gene_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>%  
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_babiceanu <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "babiceanu_recurrent_normal.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% rename(Gene_1=up_gene, Chr_1=up_chr, Strand_1=up_strand, Breakpoint_1=up_Genome_pos,
                      Gene_2=dw_gene, Chr_2=dw_chr, Strand_2=dw_strand, Breakpoint_2=dw_Genome_pos) %>%
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      mutate(Chr_1=gsub("^chr", "", Chr_1), Chr_2=gsub("^chr", "", Chr_2)) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_gtex <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "gtex_v6.tsv")
  cols_keep <- c("Fusion_Id", "Gene_1", "Chr_1", "Strand_1", "Breakpoint_1",
                 "Gene_2", "Chr_2", "Strand_2", "Breakpoint_2")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% rename(Gene_1=GeneName1, Chr_1=chr1, Strand_1=strand1, Breakpoint_1=Breakpoint1,
                      Gene_2=GeneName2, Chr_2=chr2, Strand_2=strand2, Breakpoint_2=Breakpoint2) %>%
      unite(Fusion_Id, Gene_1, Gene_2, sep="--", remove=F) %>%
      select(all_of(cols_keep)) %>% distinct()

  df
}


load_fusions_whitelist <- function(fusions_lists){
  filepath <- file.path(fusions_lists, "whitelist_clean.tsv")
  cols_keep <- c("Fusion_Id")
  df <- read_tsv(filepath, progress=F, show_col_types=F)
  df <- df %>% select(all_of(cols_keep)) %>% distinct()

  df
}


annotate_using_fusion_list <- function(df_fus, fusions_lists, name, label=NULL){
  df_lst <- switch(name,
                   "babiceanu"=load_fusions_babiceanu(fusions_lists),
                   "chimerkb"=load_fusions_chimerkb(fusions_lists),
                   "chimerseq_normal"=load_fusions_chimerseq_normal(fusions_lists),
                   "chitars_cancer"=load_fusions_chitars_cancer(fusions_lists),
                   "chitars_all_human"=load_fusions_chitars_all_human(fusions_lists),
                   "chitars_tic"=load_fusions_chitars_tic(fusions_lists),
                   "cosmic"=load_fusions_cosmic(fusions_lists),
                   "gtex"=load_fusions_gtex(fusions_lists),
                   "tic"=load_fusions_tic(fusions_lists),
                   "whitelist"=load_fusions_whitelist(fusions_lists))

  # annotate only at gene-pair level (disregard breakpoint)
  if (is.null(label)) label <- name
  col_call <- paste0(label,"_AnnotC")
  df_lst <- df_lst %>% select(Fusion_Id) %>% distinct() %>% mutate(!!col_call:=1)

  left_join(df_fus, df_lst, by="Fusion_Id")
}


match_col_classes <- function(df1, df2) {
  sharedColNames <- names(df1)[names(df1) %in% names(df2)]
  sharedColTypes <- sapply(df1[,sharedColNames], class)

  for (n in sharedColNames) {
    class(df2[[n]]) <- sharedColTypes[n]
  }

  df2
}


add_gene_annotation <- function(df_fus, df_gen, col_gen, col_sym_fus, col_sym_gen, col_gid_fus=NULL, col_gid_gen=NULL){
  cat(paste("-adding gene annotation", col_gen, "for", col_sym_fus, "..."))
  df_fus <- df_fus %>% mutate(!!col_gen:=NA) %>% mutate(Row_Order=1:nrow(df_fus))

  if (!is.null(col_gid_fus) & !is.null(col_gid_gen)){
    df_fus_a <- df_fus %>% filter(!is.na(.data[[col_gen]]))
    df_fus_b <- df_fus %>% filter(is.na(.data[[col_gen]]))
    df_fus_b[[col_gen]] <- NULL

    join_on <- setNames(col_gid_gen, col_gid_fus)
    df_gen_gid <- df_gen %>% select(all_of(c(col_gid_gen, col_gen))) %>% distinct() %>% 
      group_by(.data[[col_gid_gen]]) %>% summarize(!!col_gen:=paste0(.data[[col_gen]], collapse="|"), .groups="drop")
    df_fus_b <- left_join(df_fus_b, df_gen_gid, by=join_on)
    df_fus_a <- match_col_classes(df_fus_b, df_fus_a)
    df_fus <- bind_rows(df_fus_a, df_fus_b)
  }

  df_fus_a <- df_fus %>% filter(!is.na(.data[[col_gen]]))
  df_fus_b <- df_fus %>% filter(is.na(.data[[col_gen]]))
  df_fus_b[[col_gen]] <- NULL

  join_on <- setNames(col_sym_gen, col_sym_fus)
  df_gen_sym <- df_gen %>% select(all_of(c(col_sym_gen, col_gen))) %>% distinct() %>% 
    group_by(.data[[col_sym_gen]]) %>% summarize(!!col_gen:=paste0(.data[[col_gen]], collapse="|"), .groups="drop")
  df_fus_b <- left_join(df_fus_b, df_gen_sym, by=join_on)
  df_fus_a <- match_col_classes(df_fus_b, df_fus_a)
  df_fus <- bind_rows(df_fus_a, df_fus_b)

  cat("done!\n")
  df_fus %>% arrange(Row_Order) %>% select(-Row_Order)
}


aggregate_fusions_annotations <- function(df_fus, col_id){
  cols_call <- colnames(df_fus)[grepl("_AnnotC$", colnames(df_fus))]
  dfs_fus_call <- list()
  for (col_call in cols_call){
    dfs_fus_call[[col_call]] <- df_fus %>% filter(.data[[col_call]]==1) %>% select(all_of(col_id)) %>%
      distinct() %>% mutate(Caller=gsub("_AnnotC$", "", col_call))
  }
  df_fus_call <- bind_rows(dfs_fus_call)
  df_fus_call <- df_fus_call %>% group_by_at(all_of(col_id)) %>%
    summarize(Annotations_Custom=paste0(Caller, collapse="|"), .groups="drop")

  left_join(df_fus %>% select(-all_of(cols_call)), df_fus_call, by=col_id)
}


main <- function(args){
  # load fusions
  df_fus <- read_tsv(args$input, progress=F, show_col_types=F)

  # add [Label]_AnnotC columns from different fusion lists
  fusions_lists <- list(babiceanu="Babiceanu_Normal",
                        chimerkb="ChimerKB_v4.0",
                        chimerseq_normal="ChimerSeq_Normal_v4.0",
                        chitars_cancer="Chitars_Cancer_v5.0",
                        chitars_all_human="Chitars_All_Human_v5.0",
                        # chitars_tic="Chitars_v5.0_TIC_v3.3",
                        cosmic="COSMIC_curated_v95",
                        gtex="GTEX_V6_NAR_2020",
                        # whitelist="Whitelist_Custom", 
                        tic="TIC_v3.3")

  for (name in names(fusions_lists)){
    label <- fusions_lists[[name]]
    cat(paste("-adding annotations from", label, "..."))
    df_fus <- annotate_using_fusion_list(df_fus, args$fusions_lists, name=name, label=label)
    cat("done!\n")
  }

  # add ONE_PARTNER_IS_DRIVER_AnnotC column
  df_drv <- read_tsv(args$drivers, progress=F, show_col_types=F) %>% filter(Inclusion_Level==1)
  label <- "ONE_PARTNER_IS_DRIVER"
  col_gen <- "Inclusion_Level"
  col_gen_1 <- paste0(col_gen, "_Gene_1")
  col_gen_2 <- paste0(col_gen, "_Gene_2")
  df_fus <- add_gene_annotation(df_fus, df_gen=df_drv, col_gen=col_gen, col_sym_fus="Gene_1", col_sym_gen="Hugo_Symbol")
  df_fus <- df_fus %>% rename(!!col_gen_1:=.data[[col_gen]])
  df_fus <- add_gene_annotation(df_fus, df_gen=df_drv, col_gen=col_gen, col_sym_fus="Gene_2", col_sym_gen="Hugo_Symbol")
  df_fus <- df_fus %>% rename(!!col_gen_2:=.data[[col_gen]])
  df_fus <- df_fus %>% mutate(!!paste0(label, "_AnnotC"):=ifelse(.data[[col_gen_1]]==1 | .data[[col_gen_2]]==1,1,0)) %>%
    select(-all_of(c(col_gen_1, col_gen_2)))

  # add gencode gene type
  df_gcd <- read_tsv(args$gencode, progress=F, show_col_types=F)
  col_gen <- "gene_type"
  col_gen_1 <- "Gene_Type_1"
  col_gen_2 <- "Gene_Type_2"
  df_fus <- add_gene_annotation(df_fus, df_gen=df_gcd, col_gen=col_gen, col_sym_fus="Gene_1", col_sym_gen="gene_name",
                                col_gid_fus="Gene_Id_1", col_gid_gen="Ensembl_Gene_ID_Gencode")
  df_fus <- df_fus %>% rename(!!col_gen_1:=.data[[col_gen]])
  df_fus <- add_gene_annotation(df_fus, df_gen=df_gcd, col_gen=col_gen, col_sym_fus="Gene_2", col_sym_gen="gene_name",
                                col_gid_fus="Gene_Id_2", col_gid_gen="Ensembl_Gene_ID_Gencode")
  df_fus <- df_fus %>% rename(!!col_gen_2:=.data[[col_gen]])

  # aggregate all new annotations into a new column "Annotations_Custom"
  df_fus <- aggregate_fusions_annotations(df_fus, col_id="Fusion_Id")

  # save
  write.table(df_fus, gsub(".gz$", "", args$output), sep="\t", row.names=F, quote=F)
  if (grepl(".gz", args$output)){
    system(paste("gzip", gsub(".gz$", "", args$output)))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Custom annotator of fusions.')
  parser$add_argument('--input', type="character", help='Path to input table of fusions.',
                      default="../../../data/tcga_validation/rna/fusions/tcga_validation_annotated_FusionAnnotator.tsv.gz")
  parser$add_argument('--fusions_lists', type="character", help='Path to fusions_lists folder containing fusions lists.',
                      default="resources/fusions_lists")
  parser$add_argument('--gencode', type="character", help='Path to gencode table.',
                      default="../../../data/resources/gencode/gencode.v27.annotation_genes.tsv")
  parser$add_argument('--drivers', type="character", help='Path to table of driver genes.',
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--output", type="character", help="Path where output fusion table will be saved.")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
