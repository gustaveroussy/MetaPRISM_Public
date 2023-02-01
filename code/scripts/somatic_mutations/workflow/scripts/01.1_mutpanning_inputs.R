# @modified: 04 Jul 21
# @modified: 28 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

load_mutations_and_select_samples <- function(cohort, samples, mutations){
  # load mutations
  cat("-loading mutations...\n")
  df_mut <- load_table(mutations, header_prefix="##", guess_max=1e4)
  df_mut <- preprocess_wes_mut(df_mut, cohort,
                               cols_cln=c("Subject_Id", "Project_TCGA_More"), select_pairs=T)
  df_mut$Tumor_Type <- df_mut$Project_TCGA_More
  df_mut <- df_mut %>% unite("Sample_Id", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)

  # select samples
  cat("-selecting samples ...\n")
  df_sam <- load_table(samples)
  df_sam <- df_sam %>% filter(Use_mutpanning==1)
  df_mut <- df_mut %>% filter(Sample_Id %in% df_sam$Sample_Id)
  n_sub_sam <- df_sam %>% distinct(Subject_Id) %>% nrow()
  n_sub_mut <- df_mut %>% distinct(Subject_Id) %>% nrow()
  cat(paste0("-INFO: ", n_sub_sam, " subjects selected", ";", n_sub_mut, " with at least 1 mutation.\n"))

  df_mut
}

main <- function(args){
  df_mut <- load_mutations_and_select_samples(cohort=args$cohort, samples=args$samples, mutations=args$mutations)
  if (args$tumor_type!="PanCancer") df_mut <- df_mut %>% filter(Tumor_Type==args$tumor_type)

  df_ann <-  df_mut %>% distinct(Tumor_Sample_Barcode, Tumor_Type) %>%
    rename(Sample=Tumor_Sample_Barcode, Cohort=Tumor_Type)

  if (!file.exists(args$output_ann) | args$overwrite){
    write.table(df_ann, args$output_ann, quote=F, row.names=F, sep="\t")
    cat(paste("-input samples annotations saved at", args$output_ann))
  }
  if (!file.exists(args$output_maf) | args$overwrite){
    write.table(df_mut, args$output_maf, quote=F, row.names=F, sep="\t")
    cat(paste("-input mutations saved at", args$output_maf, "\n"))
  }
}

# parameters ===========================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build table of filtered fusions.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort. Used for loading bio,cln and sum data.",
                      default="prism")
  parser$add_argument("--samples", type="character", help="Path to input table of selected samples.",
                    default="../../../results/somatic_mutations/selection/selection_samples_prism.tsv")
  parser$add_argument("--mutations", type="character", help="Path to input mutations table.",
                    default="../../../data/prism/wes/somatic_maf/somatic_calls.maf.gz")
  parser$add_argument("--tumor_type", type="character", help="The tumor type analyzed.", default="LUAD")
  parser$add_argument("--output_ann", type="character", help="Path to output samples annotations file.",
                      default="../../../results/somatic_mutations/mutpanning/prism/LUAD_inputs/samples_LUAD.txt")
  parser$add_argument("--output_maf", type="character", help="Path to output mutations file.",
                      default="../../../results/somatic_mutations/mutpanning/prism/LUAD_inputs/mutations_LUAD.maf")
  parser$add_argument("-w", "--overwrite", action="store_true", default=F,
                      help=paste("Specify this option to overwrite existing inputs. Be aware that overwriting input",
                                 "files triggers rerun of downstream analyses by snakemake."))
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)

  main(args)
}
