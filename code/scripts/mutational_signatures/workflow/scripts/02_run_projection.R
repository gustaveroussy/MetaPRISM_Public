# @created: 22 Sep 21
# @modified: 22 Sep 21
# @authors: Yoann Pradat
#
# Project matrix of mutation counts onto known basis of mutational signatures.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(UltiSig))
source("workflow/functions/load.R")

save_projection <- function(projection, M, filepath_H, filepath_R){
  W <- projection$W # basis of signatures
  H <- projection$H # signatures weights in tumors
  R <- projection$R # R = WH is the reconstructed mutation profiles

  H_save <- cbind(Signature=rownames(H), H)
  write.table(H_save, filepath_H, sep="\t", row.names=F, quote=F)

  R_save <- cbind(Context=rownames(R), R)
  write.table(R_save, filepath_R, sep="\t", row.names=F, quote=F)
}

main <- function(args){
  # select
  df_M <- read.table(args$count_mut, sep="\t", header=T, as.is=T, check.names=F)
  df_M <- df_M %>% column_to_rownames(var="Context")
  df_W <- load_signature_profiles(basis_name=args$basis_name)

  # project
  cat("--running the projection...\n")

  if (args$package=="deconstructSigs"){
    projection <- UltiSig::project_profiles(M=df_M, W=df_W, package=args$package,
                                            signature_cutoff=0, n_cores=args$n_cores)
  } else {
    projection <- UltiSig::project_profiles(M=df_M, W=df_W, package=args$package)
  }

  cat("--done!\n")

  # save
  save_projection(projection, M=df_M, filepath_H=args$filepath_H, filepath_R=args$filepath_R)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Run projection methods')
  parser$add_argument("-c", "--count_mut", type="character", help="Path to the matrix of mutation counts.",
    default="../../../results/mutational_signatures/counts_mutations/counts_mutations_sbs_96_min_mut_prism.tsv")
  parser$add_argument("-b", "--basis_name", type="character", help="Name of the set of known signatures.",
                      default="cosmic_sbs_96_v3.2")
  parser$add_argument("-k", "--package", type="character", help="Name of the method/package employed.",
                      default="MutationalPatterns")
  parser$add_argument("--filepath_H", type="character",
                      help="Path where output signature count matrix will be saved.")
  parser$add_argument("--filepath_R", type="character",
                      help="Path where output reconstructed mutation count matrix will be saved.")
  parser$add_argument("-n", "--n_cores", type="integer", default=1,
                      help="Number of cores for the projection method. Used only if package='deconstructSigs'.")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
