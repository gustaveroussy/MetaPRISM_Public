# @created: 02 Nov 21
# @modified: 27 Dec 22
# @authors: Yoann Pradat
#
# Assemble multiple tables containing potential features into one large table. Models will later be trained
# on the subparts of this large table.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pryr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

source("../common/functions/model_utils.R")

# functions ============================================================================================================

print_size_object <- function(obj){
  obj_name <- deparse(substitute(obj))
  obj_size <- object_size(obj)
  if (obj_size > 1e9){
    print_size <- round(obj_size/1e9, 2)
    unit_size <- "GB"
  } else {
    print_size <- round(obj_size/1e6, 2)
    unit_size <- "MB"
  }
  cat(paste("-INFO MEMORY: size of", obj_name, "is", print_size, unit_size, "\n"))
}


print_total_memory <- function(){
  mem_size <- mem_used()
  if (mem_size > 1e9){
    print_size <- round(mem_size/1e9, 2)
    unit_size <- "GB"
  } else {
    print_size <- round(mem_size/1e6, 2)
    unit_size <- "MB"
  }
  cat(paste("-INFO MEMORY: total memory usage is", print_size, unit_size, "\n\n"))
}


check_size <- function(df_dat, n_row_start, name){
  if (nrow(df_dat)>n_row_start){
    ids_duplicated <- unique(df_dat[duplicated(df_dat$Subject_Id),"Subject_Id",drop=T])
    list_ids <- paste0("\t", paste0(ids_duplicated, collapse="\n\t"))
    stop(paste("the following ids are duplicated after running", name, ":\n", list_ids))
  }
}


add_cln_brca <- function(filepath, dfs){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  cols_ids <- c("Subject_Id")
  cols_bin <- c("HER2_Status", "HR_Status", "ER_Status", "PR_Status")
  cols_bin <- intersect(cols_bin, colnames(df_dat))
  cols_cat <- c("Subtype_Ihc")
  cols_dat <- c(cols_ids, cols_bin, cols_cat)
  df_dat <- df_dat[,cols_dat]
  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding CLN breast cancer specific variables from table", basename(filepath), dat_size, "..."))

  df_cov_bin <- data.frame(Covariate=cols_bin, Plot_Name=gsub("_", " ", cols_bin), Nature="Binary",
                           Reference_Level="Negative", Class_Lvl_2="Tumor")
  df_cov_cat <- data.frame(Covariate=cols_cat, Plot_Name=gsub("_Ihc", " IHC", cols_cat), Nature="Discrete_Exclusive",
                           Redundancy=paste0(cols_bin, collapse="&"), Reference_Level="HER2_plus",
                           Class_Lvl_2="Tumor_Class")
  df_cov <- bind_rows(df_cov_bin, df_cov_cat) %>%
    mutate(Class_Lvl_1="Clinical", Class_Lvl_3="At_Diagnosis")

  dfs$dat <- left_join(dfs$dat, df_dat, by="Subject_Id")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_cln_brca")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_dna_cna_chr_arm <- function(filepath, dfs, cohort){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_id <- "Sample_Id_DNA_P"
  df_dat <- df_dat %>% unite(!!col_id, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # select samples
  df_sam <- load_table(paste0("../../../data/", cohort, "/wes/somatic_cna/sample_list.tsv"))
  df_sam <- df_sam %>% unite(!!col_id, Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_", remove=F)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
  }
  df_dat <- df_dat %>% filter(.data[[col_id]] %in% df_sam[[col_id]])

  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding CNA events per chromosome arm from table", basename(filepath), dat_size, "..."))

  # get list of all analyzed pairs
  df_ids <- df_dat %>% distinct(Sample_Id_DNA_P)

  # spread losses into one table and gains into another table
  df_dat_loss <- df_dat %>% filter(copy_number %in% c(-1))
  df_dat_loss <- df_dat_loss %>% select(Sample_Id_DNA_P, arm) %>% mutate(value=1) %>% spread(arm, value)
  df_dat_loss <- bind_rows(df_dat_loss, df_ids %>% filter(!Sample_Id_DNA_P %in% df_dat_loss$Sample_Id_DNA_P))
  df_dat_loss <- df_dat_loss %>% replace(is.na(.), 0) %>% column_to_rownames(var="Sample_Id_DNA_P")

  df_dat_gain <- df_dat %>% filter(copy_number %in% c(1, 2))
  df_dat_gain <- df_dat_gain %>% select(Sample_Id_DNA_P, arm) %>% mutate(value=1) %>% spread(arm, value)
  df_dat_gain <- bind_rows(df_dat_gain, df_ids %>% filter(!Sample_Id_DNA_P %in% df_dat_gain$Sample_Id_DNA_P))
  df_dat_gain <- df_dat_gain %>% replace(is.na(.), 0) %>% column_to_rownames(var="Sample_Id_DNA_P")

  # order columns
  cols_chr_arm_p <- paste0(1:23, "p")
  cols_chr_arm_q <- paste0(1:23, "q")
  cols_chr_arm <- c(rbind(cols_chr_arm_p, cols_chr_arm_q))
  cols_chr_arm_gain <- intersect(cols_chr_arm, colnames(df_dat_gain))
  cols_chr_arm_loss <- intersect(cols_chr_arm, colnames(df_dat_loss))
  df_dat_loss <- df_dat_loss[,cols_chr_arm_loss]
  df_dat_gain <- df_dat_gain[,cols_chr_arm_gain]

  colnames(df_dat_loss) <- paste0(colnames(df_dat_loss), "_Loss")
  colnames(df_dat_gain) <- paste0(colnames(df_dat_gain), "_Gain")

  # data table
  df_dat <- bind_cols(df_dat_loss, df_dat_gain) %>% rownames_to_column(var="Sample_Id_DNA_P")

  # covs table
  cols_bin <- setdiff(colnames(df_dat), c("Sample_Id_DNA_P"))
  df_cov <- data.frame(Covariate=cols_bin, Plot_Name=gsub("_", " ", cols_bin), Nature="Binary") %>%
    mutate(Class_Lvl_1="DNA", Class_Lvl_2="CNA_Chr_Arm")

  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_DNA_P")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_cna_chr_arm")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_dna_cna_per_gene <- function(filepath, dfs, df_gen, cohort){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_id <- "Sample_Id_DNA_P"
  df_dat <- df_dat %>% unite(!!col_id, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # select samples
  df_sam <- load_table(paste0("../../../data/", cohort, "/wes/somatic_cna/sample_list.tsv"))
  df_sam <- df_sam %>% unite(!!col_id, Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_", remove=F)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
  }
  df_dat <- df_dat %>% filter(.data[[col_id]] %in% df_sam[[col_id]])

  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding CNA events per gene from table", basename(filepath), dat_size, "..."))

  df_dat <- df_dat %>% filter(Hugo_Symbol %in% df_gen$Hugo_Symbol)

  # loss per gene
  df_dat_loss <- df_dat %>% filter(Copy_Number==-2) %>% select(Sample_Id_DNA_P, Hugo_Symbol, Copy_Number) %>%
    spread(Hugo_Symbol, Copy_Number) %>% column_to_rownames(var="Sample_Id_DNA_P")
  df_dat_loss[df_dat_loss==-2] <- 1
  df_dat_loss <- df_dat_loss %>% replace(is.na(.), 0)
  colnames(df_dat_loss) <- paste0(colnames(df_dat_loss), "_Loss")
  df_dat_loss <- df_dat_loss %>% rownames_to_column(var="Sample_Id_DNA_P")

  # gain per gene
  df_dat_gain <- df_dat %>% filter(Copy_Number==2) %>% select(Sample_Id_DNA_P, Hugo_Symbol, Copy_Number) %>%
    spread(Hugo_Symbol, Copy_Number) %>% column_to_rownames(var="Sample_Id_DNA_P")
  df_dat_gain[df_dat_gain==2] <- 1
  df_dat_gain <- df_dat_gain %>% replace(is.na(.), 0)
  colnames(df_dat_gain) <- paste0(colnames(df_dat_gain), "_Gain")
  df_dat_gain <- df_dat_gain %>% rownames_to_column(var="Sample_Id_DNA_P")

  # data table
  df_dat <- full_join(df_dat_loss, df_dat_gain, by="Sample_Id_DNA_P")
  df_dat <- df_dat %>% replace(is.na(.), 0)

  # covs table
  cols_bin <- setdiff(colnames(df_dat), c("Sample_Id_DNA_P"))
  df_cov <- data.frame(Covariate=cols_bin, Plot_Name=gsub("_", " ", cols_bin), Nature="Binary") %>%
    mutate(Class_Lvl_1="DNA", Class_Lvl_2="CNA_Per_Gene", Class_Lvl_3="All")

  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_DNA_P")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_cna_per_gene")
  cat(" done!\n")

  print_size_object(dfs)
  dfs
}


add_dna_mut_signatures <- function(filepath, dfs){
  df_sig_count <- load_table(filepath) %>% column_to_rownames(var="Signature")
  if (is.null(df_sig_count)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  dat_size <- paste0("(", nrow(df_sig_count), ",", ncol(df_sig_count), ")")
  cat(paste("-INFO: adding MUT signatures from table", basename(filepath), dat_size, "..."))

  df_sig_perct <- df_sig_count %>% mutate_all(~(./sum(.) %>% as.vector))

  df_sig_count <- df_sig_count %>% t()
  df_sig_perct <- df_sig_perct %>% t()

  plot_names_count <- colnames(df_sig_count)
  plot_names_perct <- colnames(df_sig_perct)
  colnames(df_sig_count) <- paste0(colnames(df_sig_count), "_Count")
  colnames(df_sig_perct) <- paste0(colnames(df_sig_perct), "_Perct")
  df_dat <- cbind(df_sig_count, df_sig_perct)

  df_cov_count <- data.frame(Covariate=colnames(df_sig_count), Plot_Name=plot_names_count) %>% 
    mutate(Nature="Continuous", Class_Lvl_1="DNA", Class_Lvl_2="Mutational_Signatures", Class_Lvl_3="Count")
  df_cov_perct <- data.frame(Covariate=colnames(df_sig_perct), Plot_Name=plot_names_perct) %>% 
    mutate(Nature="Continuous", Class_Lvl_1="DNA", Class_Lvl_2="Mutational_Signatures", Class_Lvl_3="Percent")
  df_dat <- df_dat %>% as.data.frame() %>% rownames_to_column(var="Sample_Id_DNA_T")
  df_cov <- rbind(df_cov_count, df_cov_perct)

  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_DNA_T")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_mut_signatures")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_counts <- function(filepaths, dfs, agg="gene", data_type="DNA", evt_type="Alteration", cohort="prism"){
  dfs_data <- lapply(filepaths, load_table)
  if (length(dfs_data)==0) return(dfs)
  n_row_start <- nrow(dfs$dat)

  if (agg %in% c("pathway")){
    pattern_1 <- paste0("(?<=count_by_", agg, "_)[a-z0-9A-Z\\_\\-]+(?=_", data_type, ")")
    level_names_1 <- sapply(filepaths, function(s) str_extract(s, pattern_1), USE.NAMES=F)
  } else {
    level_names_1 <- agg
  }

  pattern_2 <- paste0("(?<=_", data_type, "_)[a-z0-9A-Z\\_\\-]+(?=_", cohort, ")")
  level_names_2 <- sapply(filepaths, function(s) str_extract(s, pattern_2), USE.NAMES=F)

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_psb <- "Sample_Id_DNA_P"

  for (i in 1:length(filepaths)){
    df_dat_cnt <- dfs_data[[i]]
    level_1 <- level_names_1[[i]]
    level_2 <- level_names_2[[i]]
    
    dat_size <- paste0("(", nrow(df_dat_cnt), ",", ncol(df_dat_cnt), ")")
    cat(paste("-INFO: adding", data_type, evt_type, "counts per", agg, "from table", basename(filepaths[[i]]), 
              dat_size, "..."))

    if (all(c(col_tsb, col_nsb) %in% colnames(df_dat_cnt))){
      df_dat_cnt <- df_dat_cnt %>% unite(!!col_psb, all_of(c(col_tsb, col_nsb)), sep="_vs_")
    }

    col_row <- intersect(c(col_psb, "Tumor_Sample_Barcode", "Sample_Id", "Subject_Id"), colnames(df_dat_cnt))[1]
    df_dat_cnt <- df_dat_cnt %>% column_to_rownames(var=col_row)
    df_dat_sts <- df_dat_cnt %>% mutate_if(is.numeric, function(x) as.integer(as.logical(x))) 

    plot_names <- colnames(df_dat_cnt)
    colnames(df_dat_cnt) <- paste0(colnames(df_dat_cnt), "_", data_type, "_count_", level_2)
    colnames(df_dat_sts) <- paste0(colnames(df_dat_sts), "_", data_type, "_status_", level_2)

    df_cov_cnt <- data.frame(Covariate=colnames(df_dat_cnt), Plot_Name=plot_names) %>% 
      mutate(Nature="Continuous", Class_Lvl_1=data_type,
             Class_Lvl_2=paste0(evt_type, "_Counts_", str_to_title(agg)), 
             Class_Lvl_3=paste0(level_1, "_", level_2))

    df_cov_sts <- data.frame(Covariate=colnames(df_dat_sts), Plot_Name=plot_names) %>% 
      mutate(Nature="Binary", Class_Lvl_1=data_type,
             Class_Lvl_2=paste0(evt_type, "_Status_", str_to_title(agg)), 
             Class_Lvl_3=paste0(level_1, "_", level_2))

    if (col_row!="Subject_Id" & col_row!=col_psb){
      if (data_type=="RNA"){
        col_row <- "Sample_Id_RNA_T"
      } else {
        col_row <- "Sample_Id_DNA_T"
      }
    }

    df_dat_cnt <- df_dat_cnt %>% as.data.frame() %>% rownames_to_column(var=col_row)
    df_dat_sts <- df_dat_sts %>% as.data.frame() %>% rownames_to_column(var=col_row)

    dfs$dat <- left_join(dfs$dat, df_dat_cnt, by=col_row)
    dfs$cov <- bind_rows(dfs$cov, df_cov_cnt)
    
    if (!agg %in% c("total")){
      dfs$dat <- left_join(dfs$dat, df_dat_sts, by=col_row)
      dfs$cov <- bind_rows(dfs$cov, df_cov_sts)
    }

    cat(" done!\n")
  }
  check_size(dfs$dat, n_row_start, paste("add ", data_type, evt_type, "counts per", agg))

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_dna_mut_counts_total <- function(filepath, dfs, level_2="All"){
  df_dat_cnt <- load_table(filepath)
  if (is.null(df_dat_cnt)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  dat_size <- paste0("(", nrow(df_dat_cnt), ",", ncol(df_dat_cnt), ")")
  cat(paste("-INFO: adding MUT counts total from table", basename(filepath), dat_size, "..."))

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_par <- "Sample_Id_DNA_P"

  stopifnot(all(c(col_tsb, col_nsb)%in% colnames(df_dat_cnt)))

  df_dat_cnt <- df_dat_cnt %>% unite(!!col_par, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
  df_dat_cnt <- df_dat_cnt[,c(col_par,"Total")] %>% column_to_rownames(var=col_par)
  df_dat_ori <- dfs$dat

  # add high tmb status (>10 mt/mb) and very high tmb status (>100 mt/mb)
  filepath_bed <- "../../../data/resources/target_files/all_targets_intersect_padded_10n.bed"
  df_bed <- read_tsv(filepath_bed, col_names=F, show_col_types=F, progress=F)
  target_size <- sum(df_bed$X3-df_bed$X2)/1e6

  df_dat_cnt <- df_dat_cnt %>% rename(TMB="Total") %>% mutate(TMB=TMB/target_size) %>% 
    mutate(High_TMB=ifelse(TMB >= 100, 2, ifelse(TMB >= 10, 1, 0)))
  plot_names <- paste(c("TMB", "High TMB"), level_2)
  colnames(df_dat_cnt) <- paste0(colnames(df_dat_cnt), "_DNA_count_", level_2)

  df_cov <- data.frame(Covariate=colnames(df_dat_cnt), Plot_Name=plot_names) %>% 
    mutate(Nature=c("Continuous", "Discrete_Ordered"),
           Class_Lvl_1="DNA",
           Class_Lvl_2=c("Mutation_Counts_Total", "Mutation_Status_Total"), 
           Class_Lvl_3=level_2)

  df_dat_cnt <- df_dat_cnt %>% rownames_to_column(var=col_par)
  df_dat_ori <- left_join(df_dat_ori, df_dat_cnt, by=col_par)

  dfs$dat <- df_dat_ori
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_mut_counts_total")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_dna_cna_summary_stats <- function(filepath, dfs, cohort){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_id <- "Sample_Id_DNA_P"
  df_dat <- df_dat %>% unite(!!col_id, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # select samples
  df_sam <- load_table(paste0("../../../data/", cohort, "/wes/somatic_cna/sample_list.tsv"))
  df_sam <- df_sam %>% unite(!!col_id, Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_", remove=F)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
  }
  df_dat <- df_dat %>% filter(.data[[col_id]] %in% df_sam[[col_id]])

  colnames(df_dat) <- gsub(":", "_", colnames(df_dat))
  df_dat <- df_dat %>% rename(Whole_Genome_Duplication=WGD, Aneuploidy_Score=aneuploidy_score,
                              Average_ploidy=Ploidy, Nb_Large_Scale_Transitions=lst)
  # select columns
  cols_keep <- c(col_id, col_tsb, col_nsb, "Average_ploidy", "Whole_Genome_Duplication",
                 "LOSS", "LOSS_Deletion", "LOSS_LOH_cnLOH", "GAIN", "GAIN_HL_amplification", "GAIN_LL_ML_amplification")
  df_dat <- df_dat %>% select(all_of(cols_keep))

  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding CNA summary statistics from table", basename(filepath), dat_size, "..."))

  cols_ids <- c(col_id, col_tsb, col_nsb)
  cols_bin <- c("Whole_Genome_Duplication")
  cols_con <- setdiff(colnames(df_dat), c(cols_ids, cols_bin))

  df_cov_bin <- data.frame(Covariate=cols_bin, Plot_Name=gsub("_", " ", cols_bin), Nature="Binary")
  df_cov_con <- data.frame(Covariate=cols_con, Plot_Name=gsub("_", " ", cols_con), Nature="Continuous")
  df_cov <- bind_rows(df_cov_bin, df_cov_con) %>%
    mutate(Class_Lvl_1="DNA", Class_Lvl_2="CNA_Summary_Statistics")
  df_cov[df_cov$Covariate=="LOSS", "Redundancy"] <- paste0(cols_con[grepl("^LOSS_", cols_con)], collapse="&")
  df_cov[df_cov$Covariate=="GAIN", "Redundancy"] <- paste0(cols_con[grepl("^GAIN_", cols_con)], collapse="&")

  dfs$dat <- left_join(dfs$dat, df_dat, by=col_id)
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_cna_summary_stats")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_dna_msi_status_score <- function(filepath, dfs, cohort){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  col_tid <- "Tumor_Sample_Id"
  col_nid <- "Normal_Sample_Id"
  col_id <- "Sample_Id_DNA_P"
  df_dat <- df_dat %>% unite(!!col_id, all_of(c(col_tid, col_nid)), sep="_vs_", remove=F)

  # select samples
  df_sam <- load_table(paste0("../../../data/", cohort, "/wes/somatic_msi/sample_list.tsv"))
  df_sam <- df_sam %>% unite(!!col_id, Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_", remove=F)
  if ("QC" %in% colnames(df_sam)){
    df_sam <- df_sam %>% filter(grepl("PASS|TBD", QC))
  }
  df_dat <- df_dat %>% filter(.data[[col_id]] %in% df_sam[[col_id]])

  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding DNA MSI status and score from table", basename(filepath), dat_size, "..."))

  # rename
  df_dat <- df_dat %>% rename(MS_Stepwise_Difference=Stepwise_Difference, MS_Status=Status)

  cols_ids <- c("Sample_Id_DNA_P", "Sample_Id_DNA_T", "Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode",
                "Tumor_Sample_Id", "Normal_Sample_Id")
  cols_dis <- c("MS_Status")
  cols_con <- setdiff(colnames(df_dat), c(cols_ids, cols_dis))

  df_cov_dis <- data.frame(Covariate=cols_dis, Plot_Name=gsub("_", " ", cols_dis), Nature="Discrete_Exclusive")
  df_cov_con <- data.frame(Covariate=cols_con, Plot_Name=gsub("_", " ", cols_con), Nature="Continuous")
  df_cov <- bind_rows(df_cov_dis, df_cov_con) %>%
    mutate(Class_Lvl_1="DNA", Class_Lvl_2="MSI")
  df_cov[df_cov$Covariate=="MS_Status", "Redundancy"] <- "MS_Stepwise_Difference"

  dfs$dat <- left_join(dfs$dat, df_dat, by=col_id)
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_dna_msi_status_score")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_rna_gex_tf_signatures <- function(filepath, dfs){
  df_dat <- load_table(filepath) 
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  df_dat <- df_dat %>% rename(Signature=`...1`) %>%
    column_to_rownames(var="Signature") %>% t()

  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding RNA transcription factor signatures from table", basename(filepath), dat_size, "..."))

  plot_names <- paste(colnames(df_dat))
  colnames(df_dat) <- paste0(colnames(df_dat), "_", "TF_Score")

  df_cov <- data.frame(Covariate=colnames(df_dat), Plot_Name=plot_names) %>% 
    mutate(Nature="Continuous", Class_Lvl_1="RNA", Class_Lvl_2="Expression_Signature", Class_Lvl_3="TF")
  df_dat <- df_dat %>% as.data.frame() %>% rownames_to_column(var="Sample_Id_RNA_T")

  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_RNA_T")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_rna_gex_tf_signatures")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_rna_gex_ssgsea_pathway <- function(filepaths, dfs, cohort="prism"){
  dfs_data <- lapply(filepaths, load_table)
  if (length(dfs_data)==0) return(dfs)
  n_row_start <- nrow(dfs$dat)

  pattern <- paste0("(?<=ssgsea_by_pathway_)[a-z0-9A-Z\\_\\-]+(?=_", cohort, ")")

  levels <- sapply(filepaths, function(s) str_extract(s, pattern), USE.NAMES=F)
  dfs_data <- setNames(dfs_data, levels)
  filepaths <- setNames(filepaths, levels)

  for (level in levels){
    filepath <- filepaths[[level]]
    df_dat <- dfs_data[[level]] %>% column_to_rownames(var="Geneset") %>% t()
    dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
    cat(paste("-INFO: adding RNA ssgsea signatures from table", filepath, dat_size, "..."))

    plot_names <- colnames(df_dat)
    colnames(df_dat) <- paste0(colnames(df_dat),"_", level)
    df_cov <- data.frame(Covariate=colnames(df_dat), Plot_Name=plot_names) %>% 
      mutate(Nature="Continuous", Class_Lvl_1="RNA", Class_Lvl_2="Expression_Signature", Class_Lvl_3=level)
    df_dat <- df_dat %>% as.data.frame() %>% rownames_to_column(var="Sample_Id_RNA_T")

    dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_RNA_T")
    dfs$cov <- bind_rows(dfs$cov, df_cov)

    cat(" done!\n")
  }
  check_size(dfs$dat, n_row_start, "add_rna_gex_ssgsea_pathway")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_rna_gex_tme_bagaev <- function(filepath, dfs){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  if (!"Sample_Id" %in% colnames(df_dat)){
    df_dat <- df_dat %>% rename(Sample_Id="...1")
  }
  df_dat <- df_dat %>% column_to_rownames(var="Sample_Id")
  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding RNA ssgsea Bagaev signatures from table", basename(filepath), dat_size, "..."))

  level <- "TME_Bagaev"
  plot_names <- colnames(df_dat)
  plot_names[plot_names=="Label"] <- "TME"
  plot_names[grepl("^Proba_", plot_names)] <- paste("prob. TME", 
                                                    gsub("Proba_", "", plot_names[grepl("^Proba_", plot_names)]))
  colnames(df_dat) <- paste0(level, "_", colnames(df_dat))
  df_cov <- data.frame(Covariate=colnames(df_dat), Plot_Name=plot_names) %>% 
    mutate(Nature=ifelse(grepl("Label", Covariate), "Discrete_Exclusive", "Continuous"),
           Class_Lvl_1="RNA", Class_Lvl_2="Immune_TME", 
           Class_Lvl_3=ifelse(grepl("Label", Covariate), paste0(level, "_Label"), paste0(level, "_Proba")))
  df_dat <- df_dat %>% as.data.frame() %>% rownames_to_column(var="Sample_Id_RNA_T")


  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_RNA_T")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_rna_gex_tme_bagaev")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


add_rna_fus_counts_total <- function(filepath, dfs){
  df_dat <- load_table(filepath)
  if (is.null(df_dat)) return(dfs)
  n_row_start <- nrow(dfs$dat)

  df_dat <- df_dat %>% rename(Sample_Id_RNA_T=Sample_Id, TFB_RNA_count_All=Total)
  dat_size <- paste0("(", nrow(df_dat), ",", ncol(df_dat), ")")
  cat(paste("-INFO: adding RNA total fusion counts from table", basename(filepath), dat_size, "..."))

  df_cov <- data.frame(Covariate="TFB_RNA_count_All", Plot_Name="TFB All") %>% 
    mutate(Nature="Continuous", Class_Lvl_1="RNA", Class_Lvl_2="Fusion_Count_Total", Class_Lvl_3="All")

  dfs$dat <- left_join(dfs$dat, df_dat, by="Sample_Id_RNA_T")
  dfs$cov <- bind_rows(dfs$cov, df_cov)
  check_size(dfs$dat, n_row_start, "add_rna_fus_counts_total")
  cat(" done!\n")

  print_size_object(dfs)
  print_total_memory()
  dfs
}


main <- function(args){
  # load data and covariates table for clinical data
  df_cov <- load_table(args$cln_covariates)
  df_cln <- load_table(args$cln)
  df_gen <- load_table(args$driver_genes) %>% filter(Inclusion_Level==1)

  # add pair DNA id
  col_tsb <- "Sample_Id_DNA_T"
  col_nsb <- "Sample_Id_DNA_N"
  col_par <- "Sample_Id_DNA_P"
  df_cln <- df_cln %>% unite(!!col_par, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # link data and covs into a list for convenience
  dat_size <- paste0("(", nrow(df_cln), ",", ncol(df_cln), ")")
  cat(paste("-INFO: initialzing tables with clinical table", basename(args$cln), dat_size, "\n"))
  dfs <- list(dat=df_cln, cov=df_cov)
  print_size_object(dfs)

  # add clinical BRCA specific covariates
  dfs <- add_cln_brca(args$cln_brca, dfs)

  # add CNA status per chromosome arm (loss/gain)
  dfs <- add_dna_cna_chr_arm(args$dna_cna_chr_arm, dfs, cohort=args$cohort)

  # add CNA status per gene for "driver" genes (loss/gain)
  dfs <- add_dna_cna_per_gene(args$dna_cna_per_gene, dfs, df_gen, cohort=args$cohort)

  # add CNA summary statistics
  dfs <- add_dna_cna_summary_stats(args$dna_cna_summary_stats, dfs, cohort=args$cohort)

  # add DNA mutational signatures 
  dfs <- add_dna_mut_signatures(args$dna_mut_signatures, dfs)

  # add mutation counts from DNA only
  dfs <- add_counts(args$dna_mut_counts_gene, dfs, agg="gene", evt_type="Mutation", data_type="DNA", cohort=args$cohort)
  dfs <- add_counts(args$dna_mut_counts_pathway, dfs, agg="pathway", evt_type="Mutation", data_type="DNA",
                    cohort=args$cohort)
  dfs <- add_dna_mut_counts_total(args$dna_mut_counts_total, dfs, level_2="All")

  # add alterations counts from DNA only
  dfs <- add_counts(args$dna_alt_counts_gene, dfs, agg="gene", evt_type="Alteration", data_type="DNA",
                    cohort=args$cohort)
  dfs <- add_counts(args$dna_alt_counts_pathway, dfs, agg="pathway", evt_type="Alteration", data_type="DNA",
                    cohort=args$cohort)
  dfs <- add_counts(args$dna_alt_counts_total, dfs, agg="total", evt_type="Alteration", data_type="DNA",
                    cohort=args$cohort)

  # add alterations counts from DNA and RNA
  dfs <- add_counts(args$dna_rna_alt_counts_gene, dfs, agg="gene", evt_type="Alteration", data_type="DNA_RNA")
  dfs <- add_counts(args$dna_rna_alt_counts_pathway, dfs, agg="pathway", evt_type="Alteration", data_type="DNA_RNA")
  dfs <- add_counts(args$dna_rna_alt_counts_total, dfs, agg="total", evt_type="Alteration", data_type="DNA_RNA")

  # add DNA msi score and status
  dfs <- add_dna_msi_status_score(args$dna_msi_status_score, dfs, cohort=args$cohort)

  # add misc RNA signatures
  dfs <- add_rna_gex_tf_signatures(args$rna_gex_tf_signatures, dfs)
  dfs <- add_rna_gex_tme_bagaev(args$rna_gex_tme_bagaev, dfs)
  dfs <- add_rna_gex_ssgsea_pathway(args$rna_gex_ssgsea_pathway, dfs, cohort=args$cohort)

  # add fusion burden
  dfs <- add_rna_fus_counts_total(args$rna_fus_counts_total, dfs)

  # save
  save_table(dfs$dat, args$output_dat, sep="\t")
  save_table(dfs$cov, args$output_cov, sep="\t")
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe distribution of survival times.')
  parser$add_argument("--cohort", nargs="?", type="character", help="Name of the cohort.",
                      default="prism")
  parser$add_argument("--cln", nargs="?", type="character", help="Path to clinical table.",
                      default="../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv")
  parser$add_argument("--cln_brca", nargs="?", type="character", help="Path to clinical table for breast cancers.",
                      default="../../../data/prism/clinical/curated_other/cln_prism_brca_curated.tsv")
  parser$add_argument("--cln_covariates", nargs="?", type="character", default="resources/table_clinical_covariates.xlsx",
                      help="Table defining classes for each of the variable in the clinical table.")
  parser$add_argument("--driver_genes", nargs="?", type="character", help="Path to table of driver genes.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--dna_cna_chr_arm", nargs="?", type="character",
                      help="Path to table of copy-number alterations per chromosome arm.",
                      default="../../../data/prism/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz")
  parser$add_argument("--dna_cna_per_gene", nargs="?", type="character",
                      help="Path to table of copy-number alterations per gene",
                      default="../../../data/prism/wes/somatic_cna/somatic_calls.tsv.gz")
  parser$add_argument("--dna_cna_summary_stats", nargs="?", type="character",
                      help="Path to table of summary statistics of somatic cna events.",
                      default="../../../data/prism/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz")
  parser$add_argument("--dna_mut_signatures", nargs="?", type="character", help="Path to table of mutational signatures",
                      default="../../../results/mutational_signatures/projection_known_signatures/MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_prism.tsv")
  parser$add_argument("--dna_mut_counts_gene", nargs="?", type="character",
                      help="Path(s) to table(s) of mutation counts per gene.",
              default="../../../results/somatic_mutations/count/count_by_gene_oncokb_cosmic_DNA_non_synonymous_prism.tsv")
  parser$add_argument("--dna_mut_counts_pathway", nargs="*", type="character",
                      help="Path(s) to table(s) of mutation counts per gene.",
      default=c("../../../results/somatic_mutations/count/count_by_pathway_cpad_kegg_DNA_non_synonymous_prism.tsv",
                "../../../results/somatic_mutations/count/count_by_pathway_sanchez_vega_DNA_non_synonymous_prism.tsv",
                "../../../results/somatic_mutations/count/count_by_pathway_msigdb_hallmarks_DNA_non_synonymous_prism.tsv"))
  parser$add_argument("--dna_mut_counts_total", nargs="?", type="character",
                      help="Path(s) to table(s) of total mutation counts.",
                      default="../../../results/somatic_mutations/count/count_total_DNA_all_prism.tsv")
  parser$add_argument("--dna_alt_counts_gene", nargs="?", type="character",
                      help="Path(s) to table(s) of DNA alteration counts per gene.",
                      default="../../../results/combined_alterations/count/count_by_gene_DNA_annotated_prism.tsv")
  parser$add_argument("--dna_alt_counts_pathway", nargs="*", type="character",
                      help="Path(s) to table(s) of DNA alteration counts per pathway.",
  default=c("../../../results/combined_alterations/count/count_by_pathway_cpad_kegg_DNA_annotated_prism.tsv",
            "../../../results/combined_alterations/count/count_by_pathway_sanchez_vega_DNA_annotated_prism.tsv",
            "../../../results/combined_alterations/count/count_by_pathway_msigdb_hallmarks_DNA_annotated_prism.tsv"))
  parser$add_argument("--dna_alt_counts_total", nargs="?", type="character",
                      help="Path(s) to table(s) of total DNA alteration counts.",
                      default="../../../results/combined_alterations/count/count_total_DNA_annotated_prism.tsv")
  parser$add_argument("--dna_rna_alt_counts_gene", nargs="?", type="character",
                      help="Path(s) to table(s) of DNA_RNA alteration counts per gene.",
                      default="../../../results/combined_alterations/count/count_by_gene_DNA_RNA_annotated_prism.tsv")
  parser$add_argument("--dna_rna_alt_counts_pathway", nargs="*", type="character",
                      help="Path(s) to table(s) of DNA_RNA alteration counts per pathway.",
  default=c("../../../results/combined_alterations/count/count_by_pathway_cpad_kegg_DNA_RNA_annotated_prism.tsv",
            "../../../results/combined_alterations/count/count_by_pathway_sanchez_vega_DNA_RNA_annotated_prism.tsv",
            "../../../results/combined_alterations/count/count_by_pathway_msigdb_hallmarks_DNA_RNA_annotated_prism.tsv"))
  parser$add_argument("--dna_rna_alt_counts_total", nargs="?", type="character",
                      help="Path(s) to table(s) of total DNA_RNA alteration counts.",
                      default="../../../results/combined_alterations/count/count_total_DNA_RNA_annotated_prism.tsv")
  parser$add_argument("--dna_msi_status_score", nargs="?", type="character",
                      help="Path to table of MSI status and score.",
                      default="../../../data/prism/wes/somatic_msi/somatic_msi.tsv")
  parser$add_argument("--rna_gex_tf_signatures", nargs="?", type="character",
                      help="Path to table of TF signature acitivity.",
                      default="../../../results/rna_general/dorothea/tf_oncogenic_prism.tsv")
  parser$add_argument("--rna_gex_tme_bagaev", nargs="?", type="character", 
                      help="Path(s) to table(s) of membership probability to each of Bagaev's immune TME classes.",
                      default="../../../results/immuno_analysis/prism/tables/mfp_subtypes_predicted_LogisticRegression.tsv")
  parser$add_argument("--rna_gex_ssgsea_pathway", nargs="*", type="character", 
                      help="Path(s) to table(s) of ssGSEA per pathway.",
                      default=c("../../../results/rna_general/ssgsea/ssgsea_by_pathway_immune_bagaev_prism.tsv.gz",
                                "../../../results/rna_general/ssgsea/ssgsea_by_pathway_msigdb_H_prism.tsv.gz",
                                "../../../results/rna_general/ssgsea/ssgsea_by_pathway_msigdb_C2_Reactome_prism.tsv.gz",
                                "../../../results/rna_general/ssgsea/ssgsea_by_pathway_msigdb_C6_prism.tsv.gz",
                                "../../../results/rna_general/ssgsea/ssgsea_by_pathway_msigdb_C8_prism.tsv.gz"))
  parser$add_argument("--rna_fus_counts_total", nargs="?", type="character",
                      help="Path(s) to table(s) of total fusion counts.",
                      default="../../../results/fusions_analysis/count/count_total_RNA_all_prism.tsv")
  parser$add_argument("--output_dat", nargs="?", type="character", help="Path to output data table.",
                      default="../../../results/survival_analysis/data_prism/all_features/data.tsv.gz")
  parser$add_argument("--output_cov", nargs="?", type="character", help="Path to output covariates table.",
                      default="../../../results/survival_analysis/data_prism/all_features/covs.tsv.gz")
  parser$add_argument("--log", nargs="?", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  for (name_arg in names(args)){
    arg <- args[[name_arg]]
    if (is.character(arg) & length(arg)==1){
      if (arg=="null"){
        args[[name_arg]] <- NULL
      }
    } 
  }

  print(args)
  main(args)
}
