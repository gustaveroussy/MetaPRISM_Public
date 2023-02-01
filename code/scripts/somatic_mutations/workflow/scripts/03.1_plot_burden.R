# @created: 29 Sep 21
# @modified: 13 Dec 22
# @authors: Yoann Pradat
#
# Draw a plot showing the mutational burden distribution per tumor type and per cohort.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# functions ============================================================================================================
load_mutations_and_select_samples <- function(args, select_samples=T){
  n_cohorts <- length(args$cohorts)
  dfs_mut <- list()
  for (i in 1:n_cohorts){
    cohort <- args$cohorts[[i]]
    input_table <- args$mutations[[i]]
    input_sample <- args$samples[[i]]
    cat(paste("-processing",  cohort, "...\n"))

    # load mutations
    cat("-loading mutations...\n")
    df_mut <- load_table(input_table, guess_max=1e4)
    df_mut <- df_mut %>% rename(Tumor_Sample_Barcode=Tumor_Sample_Id, Matched_Norm_Sample_Barcode=Normal_Sample_Id)

    # for tcga, rename aliquot_id to barcodes
    df_mut <- preprocess_wes_mut(df_mut, cohort, cols_cln=c("Subject_Id", "Project_TCGA_More"), select_pairs=T)
    df_mut <- df_mut %>% rename(Tumor_Type=Project_TCGA_More) %>%
      unite("Sample_Id", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)

    # select samples
    cat("-selecting samples ...\n")
    df_sam <- load_table(input_sample)

    if (select_samples){
      df_sam <- df_sam %>% filter(Use_plotting_sel==1)
    } else {
      df_sam <- df_sam %>% filter(Use_plotting_all==1)
    }

    df_mut <- df_mut %>% filter(Sample_Id %in% df_sam$Sample_Id)
    cat("-number of unique subjects", df_mut %>% distinct(Subject_Id) %>% nrow(), "\n")

    dfs_mut[[cohort]] <- df_mut
  }

  dfs_mut
}


prepare_burden_table <- function(df_evt, df_sam, norm_factor=1, mode="raw"){
  if (mode=="raw"){
    # count events per tumor
    df_burden <- df_evt %>%
      group_by(Subject_Id, Tumor_Type) %>%
      summarize(Burden=n()/norm_factor, .groups="keep")
  } else if (mode=="precomputed") {
    df_burden <- df_evt %>% mutate(Burden=N_Mutations/norm_factor) %>%
      select(Subject_Id, Tumor_Type, Burden)
  } else {
    stop(paste("error! unrecognized value", mode, "for mode"))
  }

  # add samples with no detected event
  df_burden_mis <- df_sam %>%
    filter(!Sample_Id %in% df_evt$Sample_Id) %>%
    select(Subject_Id, Tumor_Type) %>%
    mutate(Burden=0)

  rbind(df_burden, df_burden_mis)
}


add_pvals <- function(dfs, cohort_a, cohort_b, col_pop, tts){
  pvals <- list()
  for (tt in tts){
    dfs_a <- dfs[[cohort_a]]
    dfs_b <- dfs[[cohort_b]]
    if (tt %in% names(dfs_a) & tt %in% names(dfs_b)){
      pop_a <- dfs_a[[tt]][[col_pop]]
      pop_b <- dfs_b[[tt]][[col_pop]]
      test <- wilcox.test(x=pop_a, y=pop_b, alternative="two.sided")
      pvals[[tt]] <- test$p.value
    }
  }
  pvals
}


convert_pval_to_symbol <- function(pval){
  if (is.null(pval)){
    pattern <- "" 
  } else {
    if (pval < 0.001){
      pattern <- "***"
    } else if (pval < 0.01) {
      pattern <- "**"
    } else if (pval < 0.05) {
      pattern <- "*"
    } else {
      pattern <- ""
    }
  }
  pattern
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

add_fold_change <- function(fig, x, y, cohort, tt, dfs_plot, cohorts2colors){
  fold_change <- median(dfs_plot[[cohort]][[tt]]$Burden)/median(dfs_plot[["tcga"]][[tt]]$Burden)
  fold_change_text <- paste0(specify_decimal(fold_change, 1), "x")
  fig <- fig %>% add_annotations(x=x, y=y, text=fold_change_text, textangle=0, xanchor="centered", yanchor="bottom",
                                 xref="x", yref="paper", showarrow=F,
                                 font=list(size=20, family="Helvetica", color=cohorts2colors[[cohort]]))

  fig
}


main <- function(args){
  # load tumor counts
  df_counts_tt <- load_table(args$counts)

  # load mutations
  dfs_mut <- load_mutations_and_select_samples(args, select_samples=!args$use_all_samples)

  # load samples
  dfs_sam <- setNames(lapply(args$samples, load_table), args$cohorts)

  # select samples
  if (args$use_all_samples){
    tt_keep <- df_counts_tt %>% filter(Use_plotting_all==1) %>% pull(Tumor_Type)
    dfs_sam <- lapply(dfs_sam, function(df) df%>% filter(Use_plotting_all==1))
  } else {
    tt_keep <- df_counts_tt %>% filter(Use_plotting_sel==1) %>% pull(Tumor_Type)
    dfs_sam <- lapply(dfs_sam, function(df) df%>% filter(Use_plotting_sel==1))
  }

  # get target size
  df_bed <- read.table(args$target_bed)
  df_bed <- df_bed %>% mutate(Size=V3-V2+1)
  target_size <- sum(df_bed$Size)/1e6 

  # get tables for plot
  dfs_burden <- lapply(args$cohorts, function(x) prepare_burden_table(df_evt=dfs_mut[[x]], df_sam=dfs_sam[[x]], 
                                                                      norm_factor=target_size, mode="precomputed"))
  dfs_burden <- setNames(dfs_burden, args$cohorts)

  # for each cohort, show only tumor types with sufficient data
  for (cohort in args$cohorts){
    sid_keep <- dfs_sam[[cohort]] %>% filter(Tumor_Type %in% tt_keep) %>%
      group_by(Tumor_Type) %>% filter(n()>=args$min_counts) %>% pull(Subject_Id)
    dfs_burden[[cohort]] <- dfs_burden[[cohort]] %>% filter(Subject_Id %in% sid_keep)
  }

  # prepare plot
  offset <- 0.1
  tt_order <- compute_groups_order_burden_plot(df=bind_rows(dfs_burden),
                                               col_burden="Burden",
                                               col_groups="Tumor_Type")

  dfs_plot <- lapply(args$cohorts, function(x) compute_coordinates_burden_plot(dfs_burden[[x]], groups_order=tt_order,
                                                                               col_burden="Burden", 
                                                                               col_groups="Tumor_Type",
                                                                               offset=offset))
  dfs_plot <- setNames(dfs_plot, args$cohorts)


  # add pvalues from tests comparing prism to tcga and met500 to tcga
  cohorts2colors <- load_colors(sheet="Global")[args$cohorts]

  pvals_p_vs_t <- add_pvals(dfs=dfs_plot, cohort_a="prism", cohort_b="tcga", col_pop="Burden", tts=tt_keep)
  pvals_m_vs_t <- add_pvals(dfs=dfs_plot, cohort_a="met500", cohort_b="tcga", col_pop="Burden", tts=tt_keep)
  pvals <- c(pvals_p_vs_t, pvals_m_vs_t)
  qvals <- p.adjust(pvals, method="fdr")

  qvals_p_vs_t <- qvals[1:length(pvals_p_vs_t)]
  qvals_m_vs_t <- qvals[(length(pvals_p_vs_t)+1):(length(pvals_p_vs_t)+length(pvals_m_vs_t))]

  tt_old2new <- list()
  for (tt in tt_keep){
    if (tt %in% names(qvals_p_vs_t)){
      qval_p_vs_t <- qvals_p_vs_t[[tt]]
    } else {
      qval_p_vs_t <- NULL
    }
    if (tt %in% names(qvals_m_vs_t)){
      qval_m_vs_t <- qvals_m_vs_t[[tt]]
    } else {
      qval_m_vs_t <- NULL
    }

    pattern_p_vs_t <- convert_pval_to_symbol(qval_p_vs_t)
    pattern_m_vs_t <- convert_pval_to_symbol(qval_m_vs_t)

    sup_p_vs_t <- paste0("<span style='color:", cohorts2colors[["prism"]], "'>", pattern_p_vs_t, "</span>")
    sup_m_vs_t <- paste0("<span style='color:", cohorts2colors[["met500"]], "'>", pattern_m_vs_t, "</span>")

    if (pattern_p_vs_t=="" & pattern_m_vs_t==""){
      sup_tt <- ""
    } else if (pattern_p_vs_t!="" & pattern_m_vs_t!=""){
      sup_tt <- paste0(sup_p_vs_t,",",sup_m_vs_t)
    } else if (pattern_p_vs_t!=""){
      sup_tt <- sup_p_vs_t
    } else {
      sup_tt <- sup_m_vs_t
    }

    tt_old2new[[tt]] <- paste0('<sup>', sup_tt, '</sup>', tt)
  }

  # replace tumor types names
  tt_order_new <- sapply(tt_order, function(x) tt_old2new[[x]], USE.NAMES=F)

  dfs_plot_new <- list()
  for (cohort in args$cohorts){
    dfs_plot_new_cohort <- list()
    for (tt in tt_keep){
      if (tt %in% names(dfs_plot[[cohort]])){
        df_plot <- dfs_plot[[cohort]][[tt]]
        df_plot$Tumor_Type <- tt_old2new[[tt]]
        dfs_plot_new_cohort[[tt_old2new[[tt]]]] <- df_plot
      }
    }
    dfs_plot_new[[cohort]] <- dfs_plot_new_cohort
  }

  # draw plot
  if (args$use_all_samples){
    lwd_axes <- 2
  } else {
    lwd_axes <- 1
  }

  if (args$show_pvals){
    fig <- draw_burden_plot(dfs=dfs_plot_new, groups=tt_order_new, stacks=args$cohorts, yrange=c(10^-2, 10^3),
                            ytitle="TMB (mut/Mb)", ytitlefont=list(size=20, family="Helvetica", color="black"),
                            yname="Mutational burden", margin=list(t=args$margin_top), marker_size=4, marker_opacity=1,
                            stacks2colors=cohorts2colors, offset=offset, min_median_size=10, lwd_axes=lwd_axes,
                            xlabelsfont=list(size=20, family="Helvetica", color="black"),
                            ylabelsfont=list(size=20, family="Helvetica", color="black"),
                            legendfont=list(size=20, family="Helvetica", color="black"))
  } else {
    fig <- draw_burden_plot(dfs=dfs_plot, groups=tt_order, stacks=args$cohorts, yrange=c(10^-2, 10^3),
                            ytitle="TMB (mut/Mb)", ytitlefont=list(size=20, family="Helvetica", color="black"),
                            yname="Mutational burden", margin=list(t=args$margin_top), marker_size=4, marker_opacity=1,
                            stacks2colors=cohorts2colors, offset=offset, min_median_size=10, lwd_axes=lwd_axes,
                            xlabelsfont=list(size=20, family="Helvetica", color="black"),
                            ylabelsfont=list(size=20, family="Helvetica", color="black"),
                            legendfont=list(size=20, family="Helvetica", color="black"))
  }

  # add fold change vs TCGA when the difference is significant
  j <- 1
  for (tt in tt_order){
    y_gap <- 0.15
    y_cur <- 0.8

    for (cohort in c("prism", "met500")){
      if (cohort=="prism"){
        qvals <- qvals_p_vs_t
      } else if (cohort=="met500"){
        qvals <- qvals_m_vs_t
      }

      if (tt %in% names(qvals)){
        qval <- qvals[[tt]]
      } else {
        qval <- 1
      }

      if (qval < 0.05){
        fig <- add_fold_change(fig=fig, x=j+0.5, y=y_cur, cohort=cohort, tt=tt, dfs_plot=dfs_plot, 
                               cohorts2colors=cohorts2colors)
        y_cur <- y_cur - y_gap
      }
    }
    j <- j+1
  }

  # save
  save_image(fig, args$output, width=args$output_width, height=args$output_height)
  cat(paste("-plot saved at", args$output, "\n"))

  # save
  if (!is.null(args$output_paper)){
    save_image(fig, args$output_paper, width=args$output_width, height=args$output_height)
    cat(paste("-plot saved at", args$output_paper, "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Prepare tables for the fusions burden plots.')
  parser$add_argument("--cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                  default=c("../../../results/somatic_mutations/selection/selection_samples_prism.tsv",
                            "../../../results/somatic_mutations/selection/selection_samples_met500.tsv",
                            "../../../results/somatic_mutations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--mutations", nargs="+", help="Paths to input mutation tables.",
                      default=c("../../../data/prism/wes/summary/somatic_maf.tsv",
                                "../../../data/met500/wes/summary/somatic_maf.tsv",
                                "../../../data/tcga/wes/summary/somatic_maf.tsv"))
  parser$add_argument("--target_bed", type="character", help="Path to target bed file.",
                      default="../../../data/resources/target_files/all_targets_intersect_padded_10n.bed")
  parser$add_argument("--show_pvals", action="store_true", help="Specify to show pvalues next to tumor types.",
                      default=F)
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/somatic_mutations/selection/selection_tumor_types.tsv")
  parser$add_argument("--min_counts", type="integer", default=10,
                      help="Only tumor types with at least this number of tumors will be drawn")
  parser$add_argument("--use_all_samples", action="store_true", default=F, help="If specified, all samples are used.")
  parser$add_argument("--output_width", type="integer", default=650, help="Width of output plot in pixels.")
  parser$add_argument("--output_height", type="integer", default=350, help="Width of output plot in pixels.")
  parser$add_argument("--margin_top", type="integer", default=120, help="Top margin for labels.")
  parser$add_argument("--output", type="character", help="Path where output plot will be saved.",
                      default="../../../results/somatic_mutations/burden/cumulative_scatter_plot_prism_tumor_types.pdf")
  parser$add_argument("--output_paper", type="character", help="Path where output plot will be saved.",
                      default=NULL)
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
