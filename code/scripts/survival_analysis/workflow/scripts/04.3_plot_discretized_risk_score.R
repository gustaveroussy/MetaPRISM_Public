# @created: 24 Jun 22
# @modified: 05 Aug 22
# @authors: Yoann Pradat
#
# Uses the coefficients estimates pooled across imputations and repeats in order to compute the "risk score" (linear
# predictors from Cox model) of each individual patient.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(survminer))

source("workflow/functions/km_curves.R")
source("workflow/functions/metrics.R")
source("../common/functions/model_utils.R")

# functions ============================================================================================================

compute_risk_score <- function(df_dat, df_model, i=NULL){
  covs <- df_model$Covariate
  betas <- df_model$Coefficient
  X <- as.matrix(df_dat[,covs])
  risks <- X %*% betas
  if (is.null(i)){
    col_risk <- "Risk_Score"
  } else {
    col_risk <- paste0("Risk_Score_", i)
  }
  df_out <- df_dat %>% select(Subject_Id, Survival_Status, Survival_Time)
  df_out <- bind_cols(df_out, tibble(!!col_risk:=as.vector(risks)))
  df_out
}


decode_covariable_names_list <- function(dfs_dat, df_cov){
  for (i in 1:length(dfs_dat)){
    out <- decode_covariable_names(dfs_dat[[i]], df_cov)
    dfs_dat[[i]] <- out$df_dat
  }

  dfs_dat
}

decode_covariable_names_model <- function(df_model, df_cov){
  df_model <- df_model %>% rename(Code=Covariate)
  df_model <- left_join(df_model, df_cov[,c("Code", "Covariate")], by="Code") %>% select(-Code)

  df_model
}

apply_min_max_from_model <- function(dfs_dat, df_cov_dat, df_cov_model){
  for (i in 1:length(dfs_dat)){
    df_dat <- dfs_dat[[i]]
    for (cov in df_cov_model$Covariate){
      s_dat <- df_dat[[cov]]
      min_cov_dat <- df_cov_dat[df_cov_dat$Covariate==cov, "Min", drop=T]
      max_cov_dat <- df_cov_dat[df_cov_dat$Covariate==cov, "Max", drop=T]
      min_cov_model <- df_cov_model[df_cov_model$Covariate==cov, "Min", drop=T]
      max_cov_model <- df_cov_model[df_cov_model$Covariate==cov, "Max", drop=T]
      s_dat_undo <- s_dat * (max_cov_dat - min_cov_dat) + min_cov_dat
      s_dat_model <- (s_dat_undo - min_cov_model)/(max_cov_model - min_cov_model)
      df_dat[[cov]] <- s_dat
    }
    dfs_dat[[i]] <- df_dat
  }

  dfs_dat
}


main <- function(args){
  # get risk scores for each imputed dataset
  dfs_dat <- lapply(args$input_dats, load_table)
  df_model <- load_table(args$input_model)
  df_cov_model <- load_table(args$input_cov_model)
  df_cov_dat <- load_table(args$input_cov_dats)

  # translate from code to Covariate_Name
  dfs_dat <- decode_covariable_names_list(dfs_dat, df_cov_dat)
  df_model <- decode_covariable_names_model(df_model, df_cov_model)

  # apply min-max transformation from the model
  dfs_dat <- apply_min_max_from_model(dfs_dat, df_cov_dat, df_cov_model)
  
  # load model manually
  df_model <- list("Age_At_Biopsy"=0.77,
                   "Subtype_Ihc_HER2_minus_HR_plus"=-0.99,
                   "Subtype_Ihc_TNBC"=3.19,
                   "Burden_Annotated_No_Level_DNA_count_annotated"=3.79,
                   "Burden_Annotated_Tier1_DNA_count_annotated"=1.55,
                   "Burden_Annotated_Tier2_DNA_count_annotated"=-0.76,
                   "Burden_Annotated_Tier3_DNA_count_annotated"=0.47,
                   "Average_ploidy"=-0.73,
                   "LOSS_Deletion"=-2.82,
                   "LOSS_LOH_cnLOH"=0.37,
                   "GAIN_HL_amplification"=3.56,
                   "GAIN_LL_ML_amplification"=-2.69,
                   "MS_Stepwise_Difference"=2.35,
                   "High_TMB_DNA_count_All"=-0.55,
                   "TFB_RNA_count_All"=0.78,
                   "TME_Bagaev_Label_F"=-0.55,
                   "TME_Bagaev_Label_IE"=-1.92,
                   "TME_Bagaev_Label_IE/F"=-2.98)
  df_model <- tibble(Covariate=names(df_model), Coefficient=as.vector(unlist(df_model)))

  # compute risk score
  dfs_risk <- lapply(1:length(args$input_dats), function(i) compute_risk_score(dfs_dat[[i]], df_model, i))
  df_risk <- Reduce(left_join, dfs_risk)
  df_risk$Survival_Time <- df_risk$Survival_Time/30

  # compute average risk score
  if (args$name_cohort=="prism"){
    lim <- 6
  } else {
    lim <- 60
  }
  col_lim <- paste(">",lim,"month survival")

  cols_risk <- colnames(df_risk)[grepl("^Risk_Score", colnames(df_risk))]
  df_risk[, "Risk score"] <- rowMeans(df_risk[,cols_risk])
  df_risk[,col_lim] <- ifelse(df_risk$Survival_Time < lim & df_risk$Survival_Status==1, "No", "Yes")
  counts <- table(df_risk[,col_lim])
  df_risk[,col_lim] <- ifelse(df_risk$Survival_Time < lim & df_risk$Survival_Status==1, 
                                           paste0("< ", lim, "m (n=", counts["No"], ")"),
                                           paste0("> ", lim, "m (n=", counts["Yes"], ")"))
  
  # for plot titles
  title_part <- paste(unlist(str_split(args$name_samples, "_"))[[1]], toupper(args$name_cohort))

  # draw density plot < 6 months and > 6 months
  theme_set(theme_minimal())

  p <- ggplot(
    df_risk, 
    aes(x=`Risk score`, y=.data[[col_lim]], fill = stat(x))
    ) +
    geom_density_ridges_gradient(scale=1.4, size=0.6, rel_min_height=0.05, jittered_points=T, point_shape="|",
                                 point_size=2, point_alpha=1, alpha=1,
                                 position=position_points_jitter(width=0.05, height=0)) +
    scale_fill_gradient(low="#7FD8BE", high="#FCAB64", name="Risk score") +
    labs(title = paste(args$name_model, '-', title_part)) +
    ylab("Survival") +
    geom_vline(xintercept=c(0), linetype="dashed", color="#ec0000") +
    annotate("text", x=1, y=0.5, label="High risk", color="#ec0000", hjust=0) +
    annotate("text", x=-1, y=0.5, label="Low risk", color="#00468b", hjust=1) +
    theme(legend.position="none",
          title = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_blank()) +
    scale_y_discrete(expand = c(1, 0.5))

  # ggsave(args$output_densities, p, device="pdf", width=12, height=6, units="cm")
  ggsave(args$output_densities, p, device="pdf", width=8, height=6, units="cm")

  # draw survival curves according to risk
  if (args$name_cohort=="prism"){
    pval.coord <- c(12,0.95)
  } else {
    pval.coord <- c(12,0.5)
  }

  df_risk[,"High_Risk"] <- ifelse(df_risk$`Risk score` > 0, "High risk", "Low risk")
  col_strata <- "High_Risk"
  fit <- eval(parse(text=(paste("survfit(Surv(Survival_Time, Survival_Status) ~", col_strata, ", data=df_risk)"))))
  draw_km_curves(df_risk, fit, output=args$output_km_curves, width=3.25, height=2.25,
                 title=paste(args$name_model, "-", title_part), xlim=88, break.time.by=24, tables.height=0.25, 
                 palette=c("#ec0000", "#00468b"), legend.labs=rev(unique(df_risk[[col_strata]])), pval.coord=pval.coord,
                 ylab=NULL)

  # save table
  write.table(df_risk, gsub(".gz", "", args$output_risks), sep="\t", row.names=F, quote=F)
  system(paste("gzip -f",gsub(".gz", "", args$output_risks)))
  cat(paste("-table of risks saved at", args$output_risks, "\n"))

  # compute c-index TCGA
  taus <- c(6, 12, 24, 48, 72)
  weight.ij <- matrix(1, nrow(df_risk), nrow(df_risk))
  weight.it <- matrix(1, nrow(df_risk), length(taus))
  ipcw <- list("weight.ij"=weight.ij, "weight.it"=weight.it)
  surv <- Surv(df_risk$Survival_Time, df_risk$Survival_Status)
  risk <- df_risk$`Risk score`

  for (tau in taus){
    c_index <- c_index_ipcw(surv, risk, ipcw, tau)
    cat("c-index at", tau, ":", c_index[["s.index"]], "\n")
  }
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder_data_test <- "../../../results/survival_analysis/data_tcga/sub_features/BRCA_dna_and_rna/cln_dna_1_rna_1"
  default_folder_data_model <- "../../../results/survival_analysis/data_prism/sub_features/BRCA_dna_and_rna/cln_dna_1_rna_1"
  default_folder_models <- "../../../results/survival_analysis/models_prism/BRCA_dna_and_rna/cln_dna_1_rna_1"
  default_folder_plots <- "../../../results/survival_analysis/plots_prism/BRCA_dna_and_rna/cln_dna_1_rna_1"

  parser <- ArgumentParser(description='Prepare tables with sub features selected.')
  parser$add_argument("--name_cohort", type="character", help="Cohort name.", default="tcga")
  parser$add_argument("--name_samples", type="character", help="Samples selection.", default="BRCA_dna_rna")
  parser$add_argument("--input_dats", type="character", nargs="+", help="Path(s) to input data table.",
                     #default=file.path(default_folder_data_test, "processed", paste0("data.complete_cases.tsv.gz")))
                     default=file.path(default_folder_data_test, "processed", paste0("data.imputed_", 1:5, ".tsv.gz")))
  parser$add_argument("--input_model", type="character", nargs="+", help="Path(s) to input data table.",
                    default=file.path(default_folder_models, "none_coxph_standard", "covs.pooled_ax_rep.tsv.gz"))
  parser$add_argument("--input_cov_model", type="character", help="Path to input covariates table.",
                    default=file.path(default_folder_data_model, "processed", "covs.final.tsv.gz"))
  parser$add_argument("--input_cov_dats", type="character", help="Path to input covariates table.",
                    default=file.path(default_folder_data_test, "processed", "covs.final.tsv.gz"))
  parser$add_argument("--name_model", type="character", help="Name of the model as shown on plots.",
                    default="M7bis")
  parser$add_argument("--output_densities", type="character", help="Path to output plot of risk score distribution.",
                    default=file.path(default_folder_plots, "densities_6_month_none_coxph_standard_tcga.pdf"))
  parser$add_argument("--output_km_curves", type="character", help="Path to output plot of km curves.",
                    default=file.path(default_folder_plots, "km_curves_risk_none_coxph_standard_tcga.pdf"))
  parser$add_argument("--output_risks", type="character", help="Path to output table of predicted risk scores.",
                    default=file.path(default_folder_models, "none_coxph_standard", "predicted_risk_scores_tcga.tsv.gz"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
