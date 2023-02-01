# @created: 21 Sep 21
# @modified: 30 Dec 22
# @authors: Yoann Pradat
#
# Prepare a table with of mutations counts per sample. One column is one sample while one line is one of the mutation
# types considered (for sbs 96, there are 96 mutation types).

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(UltiSig))

# functions ============================================================================================================

get_count_table <- function(df_maf, mutation_mode){
  cat("--building matrix of somatic mutation counts \n")
  df_count <- UltiSig::count_table(df_maf=as.data.frame(df_maf),
                                   col_subject_id="Tumor_Sample_Barcode",
                                   remove_duplicates=T,
                                   remove_non_std_chr=T,
                                   keep_all_var=T, 
                                   mutation_mode=mutation_mode,
                                   genome_build="hg19",
                                   package="sigminer")
  cat("--done!\n")
  df_count
}

select_min_count <- function(df_count, min_count){
  df_count_sub <- df_count[,colSums(df_count) >= min_count]
  cat(paste0("--selected ", ncol(df_count_sub), "/",  ncol(df_count), " tumors with more than ", min_count,
               " mutations\n"))
  df_count_sub
}

main <- function(args){
  # message
  cat("-computing mutation counts matrix ...\n")

  # load maf file
  cat("--loading maf file ...")
  df_maf <- rprism::load_table(args$mut, header_prefix="##", progress=F, show_col_types=F)
  cat("done!\n")

  # load maf sum file
  cat("--loading maf summmary file ...")
  df_maf_sum <- rprism::load_table(args$mut_sum, progress=F, show_col_types=F)
  cat("done!\n")

  # select samples
  cat("--selecting pairs of tumor normal ...")
  df_maf <- rprism::preprocess_wes_mut(df_mut=df_maf, cohort=args$cohort, cols_bio=c("Subject_Id"), select_pairs=T)
  cat("done!\n")

  if (args$split_by_vaf){
    if (!"t_tumor_f" %in% colnames(df_maf)){
      df_maf <- df_maf %>% mutate(t_tumor_f=t_alt_count/t_depth)
    }

    df_vaf_high <- df_maf %>% group_by(Subject_Id) %>% summarize(High_Vaf=quantile(t_tumor_f, 2/3))
    df_maf <- left_join(df_maf, df_vaf_high, by="Subject_Id")
    df_vaf_low <- df_maf %>% group_by(Subject_Id) %>% summarize(Low_Vaf=quantile(t_tumor_f, 1/3))
    df_maf <- left_join(df_maf, df_vaf_low, by="Subject_Id")

    df_maf_high <- df_maf %>% filter(t_tumor_f>=High_Vaf)
    df_maf_low <- df_maf %>% filter(t_tumor_f<=Low_Vaf)

    df_count_1 <- get_count_table(df_maf_high, toupper(args$mode))
    df_count_1 <- select_min_count(df_count_1, floor(1/3*args$min_mut))
    df_count_2 <- get_count_table(df_maf_low, toupper(args$mode))
    df_count_2 <- select_min_count(df_count_2, floor(1/3*args$min_mut))

    common_samples <- intersect(colnames(df_count_1), colnames(df_count_2))
    cat(paste("-df_count_1 and df_count_2 have", length(common_samples), "samples in common\n"))
    df_count_1 <- df_count_1[,common_samples]
    df_count_2 <- df_count_2[,common_samples]
  } else {
    df_count_1 <- get_count_table(df_maf, toupper(args$mode))
    df_count_2 <- select_min_count(df_count_1, args$min_mut)
  }

  df_count_1 <- cbind(data.frame(Context=rownames(df_count_1)), df_count_1)
  df_count_2 <- cbind(data.frame(Context=rownames(df_count_2)), df_count_2)

  # save
  write.table(df_count_1, args$output_1,  sep="\t", quote=F, row.names=F)
  cat(paste("-file saved at", args$output_1, "\n"))

  write.table(df_count_2, args$output_2,  sep="\t", quote=F, row.names=F)
  cat(paste("-file saved at", args$output_2, "\n"))
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build mutation count matrices')
  parser$add_argument("--cohort", type="character", help="Cohort name", default="prism")
  parser$add_argument("--mut", type="character", help="Path to mutations table relative.",
                      default="../../../data/prism/wes/somatic_maf/somatic_calls.maf.gz")
  parser$add_argument("--mut_sum", type="character", help="Path to summary mutations table.",
                      default="../../../data/prism/wes/summary/somatic_maf.tsv")
  parser$add_argument("--mode", type="character",
                      help="Mode defining mutation categories. 'SBS96' is the most standard mode.")
  parser$add_argument("--min_mut", type="integer", default=50,
                      help="Minimum number of mutations required for a sample to be included in output_min_mut file.")
  parser$add_argument("--split_by_vaf", action="store_true", default=F,
                      help="If used, --output_1 and --output_2 will be high vaf and low vaf respectively.")
  parser$add_argument("--output_1", type="character",
                      help="Path to output count matrix of all mutations or high vaf with --min_mut.")
  parser$add_argument("--output_2", type="character",
                      help="Path to output count matrix of samples with --min_mut or low vaf with --min_mut.")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
