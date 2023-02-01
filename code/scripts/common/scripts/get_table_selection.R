# @created: 02 Aug 21
# @modified: 14 Dec 22
# @authors: Yoann Pradat
#
# Create a table with tumor types in rows, data sets in columns and additional binary columns indicating
# which tumor types should be used in which analysis.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(yaml))

# functions ============================================================================================================

get_cln_sample <- function(df_cln, sample_type){
  if (grepl("__", sample_type)){
    cols_sid <- paste0("Sample_Id_", unlist(strsplit(sample_type, "__")))
  } else {
    cols_sid <- paste0("Sample_Id_", sample_type)
  }
  masks <- lapply(unlist(strsplit(sample_type,"__")), function(st) grepl(st, df_cln$Sample_Type))
  mask <- Reduce(function(m1,m2) m1 & m2, masks)
  df_cln <- df_cln[mask,]

  if (length(cols_sid)==1){
    df_cln <- df_cln %>% rename(Sample_Id=.data[[cols_sid[1]]])
  } else {
    df_cln <- df_cln %>% unite("Sample_Id", all_of(cols_sid), sep="_vs_", remove=T)
  }
    
  cols <- c("Tumor_Type", "Sample_Id", "Subject_Id")
  df_cln %>% select(all_of(cols))
}


get_algo_sample <- function(df_bio, sample_type, algos, cohort){
  analyte <- unlist(strsplit(sample_type, "_"))[[1]]
  if (analyte=="DNA"){
    analyte <- "WES"
  }

  dfs_sam <- list()
  for (algo in algos){
    fp_sam <- file.path(cohort, tolower(analyte), algo, "sample_list.tsv")
    df_sam <- load_from_data(fp_sam)
    if (("Tumor_Sample_Id" %in% colnames(df_sam)) & ("Normal_Sample_Id" %in% colnames(df_sam))){
      cols_sid <- c("Tumor_Sample_Id", "Normal_Sample_Id")
      df_sam <- df_sam %>% unite("Sample_Id", all_of(cols_sid), sep="_vs_", remove=F)
    }

    if ("QC" %in% colnames(df_sam)){
      df_sam <- df_sam %>% filter(QC %in% c("PASS", "TBD"))
    }

    dfs_sam[[algo]] <- df_sam %>% select(Sample_Id)
  }

  df_sam <- Reduce(inner_join, dfs_sam)
  df_sam
}


get_alg_sample <- function(df_cln, df_bio, sample_type, algos, cohort){
  df_cln_sam <- get_cln_sample(df_cln, sample_type)
  df_bio_sam <- get_algo_sample(df_bio, sample_type, algos, cohort)
  df_cln_sam %>% filter(Sample_Id %in% df_bio_sam$Sample_Id)
}


get_samples_one <- function(df_cln, df_bio, sample_types, algos, cohort){
  if (is.null(sample_types)){
    df_sam <- df_cln %>% select(c(Tumor_Type, Subject_Id))
  } else {
    dfs_cln_sam <- lapply(sample_types, function(x) get_alg_sample(df_cln, df_bio, x, algos[[x]], cohort))
    dfs_cln_sam <- setNames(dfs_cln_sam, sample_types)

    if (length(sample_types) > 1){
      df_cln_sam <- Reduce(function(a, b) inner_join(dfs_cln_sam[[a]], dfs_cln_sam[[b]],
                                                     by=c("Subject_Id", "Tumor_Type"), 
                                                     suffix=paste0("_", c(a,b))), names(dfs_cln_sam))

      df_sam <- tibble()
      for (sample_type in sample_types){
        df_sam_st <- df_cln_sam %>% rename(Sample_Id=.data[[paste0("Sample_Id", "_", sample_type)]]) %>%
          select(Tumor_Type, Subject_Id, Sample_Id) %>% mutate(Sample_Type=sample_type)
        df_sam <- bind_rows(df_sam, df_sam_st)
      }
    } else {
      df_sam <- dfs_cln_sam[[sample_types]]
    }
  }

  df_sam
}


get_samples <- function(df_cln, df_bio, params, cohort){
  sample_types <- params$sample_types
  algos <- params$algos
  is_list_sample_types <- class(sample_types)=="list"
  is_list_algos <- class(algos[[names(algos)[1]]])=="list"

  if (is_list_sample_types & is_list_algos){
    stopifnot(setequal(names(sample_types), names(algos)))
    dfs_sam <- list()
    for (name in names(sample_types)){
      dfs_sam[[name]] <- get_samples_one(df_cln, df_bio, sample_types[[name]], algos[[name]], cohort)
    }
  } else if (is_list_sample_types){
    dfs_sam <- list()
    for (name in names(sample_types)){
      dfs_sam[[name]] <- get_samples_one(df_cln, df_bio, sample_types[[name]], algos, cohort)
    }
  } else if (is_list_algos){
    dfs_sam <- list()
    for (name in names(sample_types)){
      dfs_sam[[name]] <- get_samples_one(df_cln, df_bio, sample_types[[name]], algos, cohort)
    }
  } else {
    names_tt_keep_min_sizes <- names(params$tt_keep_min_sizes)
    names_tt_keeps <- names(params$tt_keeps)
    names_tt_drops <- names(params$tt_drops)
    names <- c(names_tt_keep_min_sizes, names_tt_keeps, names_tt_drops)
    names <- unique(names)

    if (is.null(names)) names <- "Use"

    dfs_sam <- list()
    for (name in names){
      dfs_sam[[name]] <- get_samples_one(df_cln, df_bio, sample_types, algos, cohort)
    }
  }

  dfs_sam
}


counts_per_tumor_type_one <- function(df_sam){
  df_sam %>% select(Subject_Id, Tumor_Type) %>% distinct() %>% group_by(Tumor_Type) %>% summarize(n=n())
}


counts_per_tumor_type <- function(dfs_sam){
  is_df_dfs_sam <- is.data.frame(dfs_sam)
  
  if (is_df_dfs_sam){
    dfs_cnt <- counts_per_tumor_type_one(dfs_sam)
  } else {
    dfs_cnt <- lapply(dfs_sam, counts_per_tumor_type_one)
  }

  dfs_cnt
}


init_samples_table <- function(df_cln, params){
  is_list_sample_types <- class(params$sample_types)=="list"

  if (is_list_sample_types){
    sample_types <- unique(unlist(params$sample_types))
  } else {
    sample_types <- params$sample_types
  }
  
  dfs_sam <- lapply(sample_types, function(sample_type) get_cln_sample(df_cln, sample_type))
  dfs_sam <- setNames(dfs_sam, sample_types)

  df_sam_ini <- Reduce(bind_rows, dfs_sam) %>% arrange(Subject_Id)
  df_sam_ini
}


aggregate_counts <- function(dfs_cnt){
  for (cohort in names(dfs_cnt)){
    dfs_cnt_cohort <- dfs_cnt[[cohort]]
    is_df_dfs_cnt_cohort <- is.data.frame(dfs_cnt_cohort)

    if (is_df_dfs_cnt_cohort){
      dfs_cnt[[cohort]] <- dfs_cnt[[cohort]] %>% rename(!!toupper(cohort):=n)
    } else {
      for (name in names(dfs_cnt_cohort)){
        dfs_cnt[[cohort]][[name]] <- dfs_cnt[[cohort]][[name]] %>% rename(!!paste0(toupper(cohort), "_", name):=n)
      }
      dfs_cnt[[cohort]] <-  dfs_cnt[[cohort]] %>% Reduce(function(df1,df2) full_join(df1,df2,by="Tumor_Type"), .)
    }
  }

  dfs_cnt %>%
    Reduce(function(df1,df2) full_join(df1,df2,by="Tumor_Type"), .) %>%
    arrange(Tumor_Type)
}


aggregate_samples <- function(dfs_sam){
  is_df_dfs_sam <- is.data.frame(dfs_sam)

  if (is_df_dfs_sam) {
    cols_sel <- intersect(colnames(dfs_sam), c("Subject_Id", "Sample_Id"))
    df_sam_agg <- dfs_sam %>% select(all_of(cols_sel)) %>% mutate(Use=1)
  } else {
    uses <- names(dfs_sam)
    uses_sub_only <- c()
    for (use in uses){
      cols_sel <- intersect(colnames(dfs_sam[[use]]), c("Subject_Id", "Sample_Id"))
      if (setequal(cols_sel, "Subject_Id")){
        uses_sub_only <- c(uses_sub_only, use)
      }
    }
    uses <- c(setdiff(uses, uses_sub_only), uses_sub_only)

    dfs_sam_agg <- list()
    cols_use <- c()
    for (use in uses){
      col_use <- paste0("Use_", use)
      cols_use <- c(cols_use, col_use)
      cols_sel <- intersect(colnames(dfs_sam[[use]]), c("Subject_Id", "Sample_Id"))
      dfs_sam_agg[[use]] <- dfs_sam[[use]] %>% select(all_of(cols_sel)) %>% mutate(!!col_use:=1)
    }

    df_sam_agg <- Reduce(function(df1,df2) full_join(df1, df2, by=intersect(colnames(df1), colnames(df2))), dfs_sam_agg)
    if ("Sample_Id" %in% colnames(df_sam_agg)){
      mask_sample <- !is.na(df_sam_agg$Sample_Id)
      for (col_use in cols_use){
        df_sam_agg[mask_sample & is.na(df_sam_agg[[col_use]]), col_use] <- 0
      }
    }
  }

  df_sam_agg
}


annotate_usage <- function(df_cnt, tt_keep_min_sizes, tt_keeps, tt_drops){
  usages <- names(tt_keep_min_sizes)
  stopifnot(usages==names(tt_keeps))
  stopifnot(usages==names(tt_drops))

  for (usage in usages){
    tt_keep_min_size <- tt_keep_min_sizes[[usage]]
    tt_keep <- tt_keeps[[usage]]
    tt_drop <- tt_drops[[usage]]
    col_usage <- paste0("Use_", usage)
    if (!is.null(tt_keep)) tt_keep <- paste(tt_keep, collapse="|")
    if (!is.null(tt_drop)) tt_drop <- paste(tt_drop, collapse="|")
    if (!is.null(tt_keep_min_size)){
      if ("PRISM" %in% colnames(df_cnt)){
        col_cnt <- "PRISM"
      } else {
        col_cnt <- paste0("PRISM", "_", usage)
      }
      df_cnt <- df_cnt %>% 
        mutate(!!col_usage:=ifelse(.data[[col_cnt]] >= tt_keep_min_size, 1, 0))
    } else {
      df_cnt <- df_cnt %>% 
        mutate(!!col_usage:=1)
    }
    if (!is.null(tt_drop)) df_cnt <- df_cnt %>%
      mutate(!!col_usage:=ifelse(.data[[col_usage]]==0, 0, ifelse(grepl(tt_drop, Tumor_Type), 0, 1)))
    if (!is.null(tt_keep)) df_cnt <- df_cnt %>%
      mutate(!!col_usage:=ifelse(.data[[col_usage]]==0, 0, ifelse(grepl(tt_keep, Tumor_Type), 1, 0)))
  }

  df_cnt
}


select_samples_tumor_type <- function(dfs_sam, cohort, df_cnt){
  is_df_dfs_sam <- is.data.frame(dfs_sam)

  if (is_df_dfs_sam){
    col_use <- "Use"
    tumor_types_use <- df_cnt %>% filter(.data[[col_use]]==1) %>% pull(var="Tumor_Type")
    dfs_sam <- dfs_sam %>% filter(Tumor_Type %in% tumor_types_use)
  } else {
    for (use in names(dfs_sam)){
      col_use <- paste0("Use_", use)
      tumor_types_use <- df_cnt %>% filter(.data[[col_use]]==1) %>% pull(var="Tumor_Type")
      dfs_sam[[use]] <- dfs_sam[[use]] %>% filter(Tumor_Type %in% tumor_types_use)
    }
  }

  dfs_sam
}


join_ini_and_sel <- function(df_ini, df_sel){
  if ("Subject_Id" %in% colnames(df_sel)){
    df_sel_a <- df_sel %>% filter(!is.na(Sample_Id))
    df_sel_b <- df_sel %>% filter(is.na(Sample_Id)) %>% select(-Sample_Id)
    df_ini_a <- df_ini %>% filter(!Subject_Id %in% df_sel_b$Subject_Id)
    df_ini_b <- df_ini %>% filter(Subject_Id %in% df_sel_b$Subject_Id)
    df_ini_a <- left_join(df_ini_a, df_sel_a)
    df_ini_b <- left_join(df_ini_b, df_sel_b)
    df_join <- bind_rows(df_ini_a, df_ini_b) %>% arrange(Subject_Id)
  } else {
    df_join <- left_join(df_ini, df_sel) %>% arrange(Subject_Id)
  }

  df_join <- df_join %>% replace(is.na(.), 0)
  df_join
}


main <- function(args){
  # load data
  if (length(args$cohorts)==1){
    dfs_cln <- setNames(list(load_cln(study=args$cohorts)), args$cohorts)
    dfs_bio <- setNames(list(load_bio(study=args$cohorts)), args$cohorts)
  } else {
    dfs_cln <- sapply(args$cohorts, function(cohort) load_cln(study=cohort))
    dfs_bio <- sapply(args$cohorts, function(cohort) load_bio(study=cohort))
  }

  # load params
  params <- read_yaml(args$config)
  sections <- unlist(strsplit(args$section, "/"))
  for (section in sections){
    params <- params[[section]]
  }

  # rename Project_TCGA_More to Tumor_Type
  dfs_cln <- lapply(dfs_cln, function(df) df %>% rename(Tumor_Type=.data[[params$col_tt]]))

  # init samples table
  dfs_sam_ini <- lapply(dfs_cln, function(df_cln) init_samples_table(df_cln, params))

  # get tables of sample/subject selections
  dfs_sam_sel <- lapply(args$cohorts, function(x) get_samples(df_cln=dfs_cln[[x]], df_bio=dfs_bio[[x]], params, x))
  dfs_sam_sel <- setNames(dfs_sam_sel, args$cohorts)

  # get table of tumor type counts
  dfs_cnt <- lapply(dfs_sam_sel, function(dfs_sam) counts_per_tumor_type(dfs_sam=dfs_sam))
  df_cnt <- aggregate_counts(dfs_cnt)
  df_cnt <- annotate_usage(df_cnt, tt_keep_min_sizes=params$tt_keep_min_sizes, tt_keeps=params$tt_keeps,
                           tt_drops=params$tt_drops)
  df_cnt <- df_cnt %>% replace(is.na(.), 0)

  # subset tables of samples selection on selected tumor types
  dfs_sam_sel <- lapply(args$cohorts, function(cohort) select_samples_tumor_type(dfs_sam=dfs_sam_sel[[cohort]],
                                                                                 cohort=cohort, df_cnt=df_cnt))
  dfs_sam_sel <- setNames(dfs_sam_sel, args$cohorts)

  # merge table of samples initialized and selected
  dfs_sam_sel <- lapply(dfs_sam_sel, aggregate_samples)
  dfs_sam <- lapply(args$cohorts, function(x) join_ini_and_sel(df_ini=dfs_sam_ini[[x]], df_sel=dfs_sam_sel[[x]]))
  dfs_sam <- setNames(dfs_sam, args$cohorts)

  # save
  write.table(df_cnt, args$output_cnt, quote=F, row.names=F, sep="\t")
  cat(paste("-file saved at", args$output_cnt, "\n"))

  for (i in 1:length(args$cohorts)){
    cohort <- args$cohort[i]
    output <- args$output_sam[i]
    df_sam <- dfs_sam[[cohort]]
    mask_sample <- !is.na(df_sam$Sample_Id)
    cols_use <- colnames(df_sam)[grepl("^Use_", colnames(df_sam))]
    for (col_use in cols_use){
      df_sam[mask_sample & is.na(df_sam[[col_use]]), col_use] <- 0
    }

    write.table(df_sam, output, quote=F, row.names=F, sep="\t")
    cat(paste("-file saved at", output, "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description=paste("Build table of sample counts per tumor type along with a binary",
                                             "indicator for each analysis."))
  parser$add_argument("--cohorts", type="character", nargs="+", default=c("tcga", "prism", "met500"),
                      help="Names of the cohorts.")
  parser$add_argument("--config", type="character", default="config/config.yaml",
                      help=paste("Path to the config file."))
  parser$add_argument("--section", type="character", default="selection",
                      help=paste("Name of the section in the config file config/config.yaml. This section must contain",
                                 "the keys: 'sample_types', 'tt_keep_min_sizes', 'tt_keeps', 'tt_drops'. You can use",
                                 "/ separators if the section is nested."))
  parser$add_argument("--output_cnt", type="character", help="Path to output table of tumor types counts and selection.")
  parser$add_argument("--output_sam", nargs="+", type="character", 
                      help="Path to output tables of samples selection for each cohort.")
  parser$add_argument("--log", type="character", help="Path to log file.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
