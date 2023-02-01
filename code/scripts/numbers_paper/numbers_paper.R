# @created: 28 Jun 2022
# @modified: 23 Aug 2022
# @author: Yoann Pradat
# 
#     CentraleSupelec
#     MICS laboratory
#     9 rue Juliot Curie, Gif-Sur-Yvette, 91190 France
# 
#     Institut Gustave Roussy
#     Prism Center
#     114 rue Edouard Vaillant, Villejuif, 94800 France
# 
# Compute numbers and percentages inside the text.

suppressMessages(library(argparse))
suppressMessages(library(dplyr))
suppressMessages(library(forcats))
suppressMessages(library(foreach))
suppressMessages(library(rprism))
suppressMessages(library(readr))
suppressMessages(library(readxl))
suppressMessages(library(tibble))
suppressMessages(library(tidyr))
suppressMessages(library(stringr))
suppressMessages(library(Exact))

# functions ============================================================================================================

select_by_sample_type <- function(df, cohort, sample_type="DNA_T|RNA_T", col_id="Subject_Id", cols_cln=c()){
  df_cln <- load_cln(cohort)
  df_cln <- df_cln %>% rename(Tumor_Type=Project_TCGA_More)
  mask <- grepl(sample_type, df_cln$Sample_Type)
  df_cln_st <- df_cln[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," subjects with ", sample_type, "\n"))

  if (col_id=="Sample_Id"){
    dfs_cln_sub <- list()
    for (st in unlist(str_split(sample_type, "\\|"))){
      col_id_st <- paste0("Sample_Id_", st)
      df_cln_st_sub <- df_cln_st %>% filter(!is.na(.data[[col_id_st]]))
      df_cln_st_sub <- df_cln_st_sub %>% rename(Sample_Id=.data[[col_id_st]])
      dfs_cln_sub[[st]] <- df_cln_st_sub
    }
    df_cln_sub <- bind_rows(dfs_cln_sub)
  } else if (col_id=="DNA_P"){
      cols_sid <- c("Sample_Id_DNA_T", "Sample_Id_DNA_N")
      df_cln_sub <- df_cln_st
      df_cln_sub <- df_cln_st %>% unite(!!col_id, all_of(c(cols_sid)), sep="_vs_", remove=F)
  } else {
    df_cln_sub <- df_cln_st
  }

  mask <- df[[col_id]] %in% df_cln_sub[[col_id]]
  df_st <- df[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," lines from subjects with ", sample_type, "\n"))

  # add columns
  if (length(cols_cln)>0){
    df_st <- left_join(df_st, df_cln_sub[,c(col_id, cols_cln)], by=col_id)
  }

  df_st
}


select_by_selections <- function(df, filepath_sam, name, col_id, col_id_sam="Sample_Id"){
  df_sam <- load_table(filepath_sam)
  df_sam <- df_sam[df_sam[[paste0("Use_", name)]]==1,]
  mask <- df[[col_id]] %in% df_sam[[col_id_sam]]
  df <- df[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," lines from selected samples\n"))
  df
}


get_margins <- function(df_cln){
  margins <- df_cln %>% group_by(Project_TCGA_More) %>% summarize(n=n())
  margins <- margins %>% column_to_rownames(var="Project_TCGA_More") 

  margins
}


get_counts <- function(df_alt, margins, col="Res_Level_Simple", genes=c()){
  # filter
  ## keep only alterations with Tiers
  df_alt_nna <- df_alt %>% filter(!is.na(.data[[col]]))

  ## keep only some genes (if any)
  if (length(genes)>0){
    mask <- Reduce(function(m1, m2) m1 | m2, lapply(genes, function(gene) grepl(gene, df_alt_nna$Hugo_Symbol)))
    df_alt_nna <- df_alt_nna %>% filter(mask)
  }

  df_count <- df_alt_nna %>% select(Subject_Id, Tumor_Type, .data[[col]]) %>% distinct() %>%
    group_by(Tumor_Type, .data[[col]]) %>% summarize(Count=n(), .groups="keep")

  tts <- rownames(margins)
  tts_no_alt <- setdiff(tts, df_count$Tumor_Type)
  for (tt in tts_no_alt){
    df_count <- bind_rows(df_count, tibble(Tumor_Type=tt, !!col:="Tier1", Count=0))
  }

  df_count <- df_count %>% spread(Tumor_Type, Count) %>% replace(is.na(.), 0) %>% column_to_rownames(var=col)
  df_count <- df_count[,tts]

  df_count
}


get_pvals_cmh <- function(df_count_a, df_count_b, margins_a, margins_b, n_cores=1){
  stopifnot(identical(rownames(df_count_a), rownames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), colnames(df_count_b)))
  stopifnot(identical(colnames(df_count_a), rownames(margins_a)))
  stopifnot(identical(colnames(df_count_b), rownames(margins_b)))

  # can only compare across common strats
  common_strats <- (margins_a > 0) & (margins_b > 0)
  dfc_count_a <- df_count_a[,common_strats,drop=F]
  dfc_count_b <- df_count_b[,common_strats,drop=F]
  marginsc_a <- margins_a[common_strats,,drop=F]
  marginsc_b <- margins_b[common_strats,,drop=F]

  n_row <- nrow(dfc_count_a)
  n_strat <- ncol(dfc_count_a)

  pvals <- foreach(i = 1:n_row, .packages=c("Exact"), .combine=c) %do% {
    c_a <- dfc_count_a[i,]
    c_b <- dfc_count_b[i,]
    
    if (sum(c_a)==0 & sum(c_b)==0){
      pval <- -999
    } else {
      x  <- sapply(1:n_strat, function(j) 
                   {matrix(c(c_a[j], marginsc_a[j,] - c_a[j], c_b[j], marginsc_b[j,] - c_b[j]), byrow=T, nrow=2)},
                   simplify="array")

      table2x2 <- function(j){
        c(c_a[j], c_b[j], marginsc_a[j,] - c_a[j],  marginsc_b[j,] - c_b[j])
      }

      x <- array(Reduce(c, sapply(1:n_strat, table2x2)), dim=c(2,2,n_strat))

      # NOTE: more work is needed here.
      # The exact CMH test is a conditional test that assumes fixed row and column margins.
      # Is the asymptotic CMH test (chi-squared distributed statistic) a conditional test?
      # Let us choose exact=F for now.
      #
      # correct=T applies Yates' correction.
      #
      # References https://doi.org/10.3390/stats1010008 and https://onlinelibrary.wiley.com/doi/abs/10.1111/anzs.12215
      # may help in choosing a test that is conditional only on fixed row margins.
      out <- mantelhaen.test(x=x, alternative="two.sided", correct=T, exact=F)

      pval <- out$p.value
      common_or <- out$estimate
      cat(paste("-for", rownames(dfc_count_a)[i], "commonds OR estimate is:", common_or, "and 95% CI is:",
                paste0("[", out$conf.int[1], ";", out$conf.int[2], "]\n")))
    }

    pval
  }

  df_pvals <- data.frame("p.value"=pvals, row.names=rownames(df_count_a))
  df_pvals
}


main <- function(args){
  cohorts <- c("prism", "met500", "tcga")

  # clinical tables
  dfs_cln <- setNames(lapply(cohorts, load_cln), cohorts)
  dfs_cln <- lapply(cohorts, function(cohort)
    select_by_sample_type(df=dfs_cln[[cohort]], cohort=cohort, sample_type="DNA_T\\|RNA_T", col_id="Subject_Id"))
  dfs_cln <- setNames(dfs_cln, cohorts)
  dfs_cln <- lapply(cohorts, function(cohort)
    select_by_selections(df=dfs_cln[[cohort]], filepath_sam=eval(parse(text=paste0("args$CA_selection_", cohort))),
                         name="heatmap_all", col_id="Subject_Id", col_id_sam="Subject_Id"))
  dfs_cln <- setNames(dfs_cln, cohorts)

  # alteration tables
  dfs_alt <- lapply(cohorts, function(cohort) eval(parse(text=paste0("load_table(args$CA_alt_table_", cohort, ")"))))
  dfs_alt <- setNames(dfs_alt, cohorts)
  dfs_alt <- lapply(cohorts, function(cohort)
    select_by_selections(df=dfs_alt[[cohort]], filepath_sam=eval(parse(text=paste0("args$CA_selection_", cohort))),
                         name="heatmap_all", col_id="Sample_Id", col_id_sam="Sample_Id"))
  dfs_alt <- setNames(dfs_alt, cohorts)

  # PRISM vs. TCGA
  cohort_a <- "prism"
  cohort_b <- "tcga"
  margins_a <- get_margins(dfs_cln[[cohort_a]])
  margins_b <- get_margins(dfs_cln[[cohort_b]])

  df_count_a <- get_counts(dfs_alt[[cohort_a]], margins_a, col="Res_Level_Simple")
  df_count_b <- get_counts(dfs_alt[[cohort_b]], margins_b, col="Res_Level_Simple")

  get_pvals_cmh(df_count_a, df_count_b, margins_a, margins_b, n_cores=1)

  # MET500 vs. TCGA
  cohort_a <- "met500"
  cohort_b <- "tcga"
  margins_a <- get_margins(dfs_cln[[cohort_a]])
  margins_b <- get_margins(dfs_cln[[cohort_b]])

  df_count_a <- get_counts(dfs_alt[[cohort_a]], margins_a, col="Res_Level_Simple")
  df_count_b <- get_counts(dfs_alt[[cohort_b]], margins_b, col="Res_Level_Simple")

  get_pvals_cmh(df_count_a, df_count_b, margins_a, margins_b, n_cores=1)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Compute numbers and percentages in the text.')
  parser$add_argument('--CA_alt_table_prism', type="character", help='Table of aggregated alterations from CA.',
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv")
  parser$add_argument('--CA_alt_table_met500', type="character", help='Table of aggregated alterations from CA.',
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_met500.tsv")
  parser$add_argument('--CA_alt_table_tcga', type="character", help='Table of aggregated alterations from CA.',
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_tcga.tsv")
  parser$add_argument('--CA_selection_prism', type="character", help='Table of sample selections from CA.',
                      default="../../../results/combined_alterations/selection/selection_samples_prism.tsv")
  parser$add_argument('--CA_selection_met500', type="character", help='Table of sample selections from CA.',
                      default="../../../results/combined_alterations/selection/selection_samples_met500.tsv")
  parser$add_argument('--CA_selection_tcga', type="character", help='Table of sample selections from CA.',
                      default="../../../results/combined_alterations/selection/selection_samples_tcga.tsv")
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")


  main(args)
}
