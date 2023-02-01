# @created: 28 Sep 21
# @modified: 02 Jun 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(xlsx))

# functions ============================================================================================================

load_cnas_and_select_samples <- function(args, select_samples=T){
  n_cohorts <- length(args$cohorts)
  dfs_cna <- list()

  col_pair <- "DNA_P"
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_tid <- "Sample_Id_DNA_T"
  col_nid <- "Sample_Id_DNA_N"
  col_sub <- "Subject_Id"
  col_tt <- "Project_TCGA_More"

  for (i in 1:n_cohorts){
    cohort <- args$cohorts[[i]]
    input_table <- args$cna_chr[[i]]
    input_sample <- args$samples[[i]]
    cat(paste("-processing",  cohort, "...\n"))

    # load cnas
    cat("-loading cnas\n")
    df_cna <- load_table(input_table)
    df_cna <- df_cna %>% unite(!!col_pair, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

    # select in-design samples
    df_cln <- load_cln(study=cohort)
    df_cln <- df_cln %>% unite(!!col_pair, all_of(c(col_tid, col_nid)), sep="_vs_", remove=F)

    mask_in <- df_cna[[col_pair]] %in% df_cln[[col_pair]]
    df_cna <- df_cna[mask_in,]
    cat("-selected", paste0(sum(mask_in), "/", length(mask_in)) , "alterations from in-design DNA samples ...\n")
    df_cna <- left_join(df_cna, df_cln[c(col_sub, col_tt, col_pair)], by=col_pair)

    if (select_samples){
      # select samples
      df_sam <- load_table(input_sample)
      df_sam <- df_sam %>% filter(Use_selections==1)
      mask_in <- df_cna[[col_pair]] %in% df_sam$Sample_Id
      df_cna <- df_cna[mask_in,]
      cat("-selected", paste0(sum(mask_in), "/", length(mask_in)) , "alterations from selected DNA samples ...\n")
      cat("-number of unique subjects", df_cna %>% distinct(Subject_Id) %>% nrow(), "\n")
    }

    # add Alteration
    # removing copy_number=-2, i.e homozyguous deletions of chr arms
    df_cna <- df_cna %>% filter(copy_number %in% c(-1, 1, 2)) %>%
      mutate(Alteration=ifelse(copy_number < 0, "loss", "gain")) %>%
      rename(Chromosome_Arm=arm)

    # format cna
    cols_idx <- c(col_tsb, col_nsb, col_sub, col_tt, "Chromosome_Arm", "Alteration")
    df_cna <- df_cna %>% select(all_of(cols_idx)) %>% rename(Tumor_Type=.data[[col_tt]])

    dfs_cna[[cohort]] <- df_cna
  }

  dfs_cna
}


select_non_wgd_samples <- function(dfs_cna, dfs_wgd){
  col_pai <- "DNA_P"
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  cohorts <- names(dfs_cna)
  
  for (cohort in cohorts){
    df_cna <- dfs_cna[[cohort]] %>% unite(!!col_pai, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
    df_wgd <- dfs_wgd[[cohort]] %>% unite(!!col_pai, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F) %>%
      filter(WGD > 0)
    cna_pai <- unique(df_cna[[col_pai]])
    wgd_pai <- unique(df_wgd[[col_pai]])
    cna_non_wgd_pai <- setdiff(cna_pai, wgd_pai)
    mask_in <- df_cna[[col_pai]] %in% cna_non_wgd_pai
    cat("-INFO: selected", sum(mask_in), "lines from", paste0(length(cna_non_wgd_pai), "/", length(cna_pai)),
        "without WGD from cohort", cohort, "\n")
    dfs_cna[[cohort]] <- df_cna[mask_in,]
  }

  dfs_cna
}


check_min_counts <- function(min_counts, cohorts){
  if (length(min_counts)==1){
    min_counts <- rep(min_counts, length.out=length(cohorts))
  } else {
    stopifnot(length(min_counts)==length(cohorts))
  }

  setNames(rep(min_counts, length.out=length(cohorts)), cohorts)
}


select_alterations <- function(dfs_cna, col_alt, alterations){
  if (alterations=="losses"){
    dfs_cna <- lapply(dfs_cna, function(df) df %>% filter(.data[[col_alt]]=="loss"))
  } else if (alterations=="gains"){
    dfs_cna <- lapply(dfs_cna, function(df) df %>% filter(.data[[col_alt]]=="gain"))
  } else {
    stop("-unsupported value", alterations, "for --alterations option.")
  }

  dfs_cna
}


get_tables_for_plot <- function(dfs_alt, df_tt_count, cohorts_keep, evt_keep, tt_keep, col_evt, col_tt, col_stack,
                                cohorts_a=c("prism", "met500"), cohorts_b=c("tcga", "tcga"), min_counts_evt=NULL,
                                max_evt=NULL, n_cores=4){
  dfs_plot <- list()
  
  # For each cohort, compute tables of event counts/percentages aggregated by Tumor_Type
  for (cohort in cohorts_keep){
    df_alt_count_cohort <- dfs_alt[[cohort]] %>% filter(.data[[col_evt]] %in% evt_keep) %>%
      group_by(.data[[col_tt]], .data[[col_evt]]) %>% summarize(Count=n(), .groups="keep")
    df_tt_count_cohort <- df_tt_count %>% rename(Count=.data[[paste0(toupper(cohort), "_", "selections")]])

    dfs_tabs <- get_tables_for_facetted_heatmap_barplots(df_alt_count_cohort, df_tt_count_cohort, tt_keep, col_evt, col_tt)
    dfs_plot[[cohort]] <- dfs_tabs
  }

  # Select list of genes for the plot. This step is important as it conditions how many tests will be performed
  in_plot_evt <- select_in_plot_evt(dfs_plot, min_counts_evt, col_evt, max_evt=max_evt, max_evt_cohort="prism")
  dfs_plot <- lapply(dfs_plot, function(dfs) select_vals_in_dfs(dfs, col_evt, in_plot_evt))

  # differential testing per tumor type
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b, test="fisher", col_evt=col_evt,
                               in_plot_evt=in_plot_evt, tt_keep=tt_keep, n_cores=n_cores, suffix="")

  # differential testing across tumor types 
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b, test="cmh", col_evt=col_evt,
                               in_plot_evt=in_plot_evt, tt_keep=tt_keep, n_cores=n_cores, suffix="_strats")

  dfs_plot
}
 

save_tables_xlsx <- function(dfs_plot, col_tt, filepath){
  df_desc <- data.frame(Supplementary=c("Tables for facetted heatmap"))
  write.xlsx(df_desc, file=filepath, sheetName="description", row.names=FALSE)

  for (name_a in names(dfs_plot)){
    for (name_b in names(dfs_plot[[name_a]])){
      sheet_name <- paste(name_a, name_b, sep="_")
      df_name <- dfs_plot[[name_a]][[name_b]]

      if (!is.null(df_name)){
        if (name_b=="count_col") df_name <- df_name %>% rownames_to_column(var=col_tt)
        write.xlsx(as.data.frame(df_name), file=filepath, sheetName=sheet_name, row.names=FALSE, append=T)
      }
    }
  }

  cat(paste("-INFO: excel workbook saved at", filepath), "\n")
}


order_dfs <- function(dfs, col_var, order_var=NULL){
  if (is.null(order_var)){
    order_var <- dfs$count_row %>% arrange(Count) %>% pull(col_var) 
  } else {
    order_var <- intersect(order_var, dfs$count_row %>% arrange(Count) %>% pull(col_var))
  }

  for (name in names(dfs)){
    if (col_var %in% colnames(dfs[[name]])){
      dfs[[name]] <- dfs[[name]] %>% arrange(match(.data[[col_var]], order_var))
    }
  }

  dfs
}


main <- function(args){
  # load alterations and tt counts
  dfs_cna <- load_cnas_and_select_samples(args)

  # load wgd status and select only samples without WGD
  dfs_wgd <- lapply(args$cna_wgd, load_table)
  dfs_wgd <- setNames(dfs_wgd, args$cohorts)
  dfs_cna <- select_non_wgd_samples(dfs_cna, dfs_wgd)

  # exceptionnaly, modify counts of TCGA (2 reasons, 480 and WGD) and of PRISM (WGD) and MET500 (WGD)
  df_tt_count_met500 <- dfs_cna$met500 %>% select(Subject_Id, Tumor_Type) %>% distinct() %>% group_by(Tumor_Type) %>%
    summarize(MET500_selections=n())
  df_tt_count_prism <- dfs_cna$prism %>% select(Subject_Id, Tumor_Type) %>% distinct() %>% group_by(Tumor_Type) %>%
    summarize(PRISM_selections=n())
  df_tt_count_tcga <- dfs_cna$tcga %>% select(Subject_Id, Tumor_Type) %>% distinct() %>% group_by(Tumor_Type) %>%
    summarize(TCGA_selections=n())
  dfs_tt_count <- list(df_tt_count_met500, df_tt_count_prism, df_tt_count_tcga)
  df_tt_count <- Reduce(function(df1, df2) left_join(df1, df2, by="Tumor_Type"), dfs_tt_count)
  df_tt_count$Use_selections <-1 

  # read input min counts
  min_counts_evt <- check_min_counts(args$min_counts_evt, args$cohorts)

  # define util variable names
  col_tt <- "Tumor_Type"
  col_evt <- "Chromosome_Arm"
  col_alt <- "Alteration"

  # select alterations
  dfs_cna <- select_alterations(dfs_cna=dfs_cna, col_alt=col_alt, alterations=args$alterations)

  # select data for columns
  cohorts_keep <- c("prism", "met500", "tcga")
  df_tt_count <- df_tt_count %>% filter(Use_selections==1)
  tt_keep <- df_tt_count %>% pull(.data[[col_tt]])
  evt_keep <- Reduce(union, lapply(dfs_cna, function(df) unique(df[[col_evt]])))

  dfs_plot <- get_tables_for_plot(dfs_alt=dfs_cna, df_tt_count=df_tt_count, cohorts_keep=cohorts_keep,
                                  evt_keep=evt_keep, tt_keep=tt_keep, col_evt=col_evt, col_tt=col_tt,
                                  col_stack=NULL, cohorts_a=c("prism", "met500"), cohorts_b=c("tcga", "tcga"), 
                                  min_counts_evt=min_counts_evt, n_cores=args$n_cores)

  # save all tables in an excel workbook
  save_tables_xlsx(dfs_plot=dfs_plot, col_tt=col_tt, filepath=args$output_tables)

  # aggregate tables across cohorts to prepare for the plot
  cohorts_keep_order <- c("prism", "met500", "tcga")
  df_tt_count_ori <- load_table(args$counts)
  tt_keep_order <- df_tt_count_ori %>% filter(Use_selections==1) %>% arrange(desc(PRISM_selections)) %>%
    pull(.data[[col_tt]])

  dfs_plot_agg <- dfs_plot$prism

  for (name_comp in names(dfs_plot)[grepl("_vs_", names(dfs_plot))]){
    for (subname in names(dfs_plot[[name_comp]])){
      new_name <- paste(name_comp, subname, sep="_")
      dfs_plot_agg[[new_name]] <- dfs_plot[[name_comp]][[subname]]
    }
  }

  for (name in names(dfs_plot_agg)){
    df_for_plot <- dfs_plot_agg[[name]]
    colnames_df <- colnames(df_for_plot)
    if (length(intersect(colnames_df, tt_keep_order))>0){
      colnames_df_other <- setdiff(colnames_df, tt_keep_order)
      df_for_plot <- df_for_plot[c(colnames_df_other, tt_keep_order)]
    }
    dfs_plot_agg[[name]] <- df_for_plot
  }

  # order rows for the heatmap
  order_var <- rev(c(rbind(paste0(seq(1,23), "p"), paste0(seq(1,23), "q"))))
  dfs_plot_agg <- order_dfs(dfs_plot_agg, col_var=col_evt, order_var=order_var)

  # define arguments for the function draw_facetted_heatmap_barplots_2
  cohorts2colors <- load_colors(sheet="Global")[c("prism", "met500")]
  tumors2colors <- load_colors(sheet="Project_TCGA_More")

  if (args$alterations=="gains"){
    colors_palette_heatmap <- "Reds"
  } else {
    colors_palette_heatmap <- "Blues"
  }


  names2colors <- list(met500_vs_tcga_qvals_strats=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals_strats=cohorts2colors[["prism"]],
                       met500_vs_tcga_qvals=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals=cohorts2colors[["prism"]],
                       colnames=tumors2colors)

  names2plots <- list(heatmap="percent",
                      heatmap_hover="count",
                      count_col="count_col",
                      heatmap_pvals=c("prism_vs_tcga_qvals", "met500_vs_tcga_qvals"),
                      left_pvals=c("prism_vs_tcga_qvals_strats", "met500_vs_tcga_qvals_strats"),
                      middle_right_bar="count_row")

  fonts <- list(x_tick_heatmap=list(size=14, color="black", family="Helvetica"),
                y_tick_heatmap=list(size=10, color="black", family="Helvetica"),
                x_tick_row_bar=list(size=12, color="black", family="Helvetica"),
                z_heatmap=list(size=12, color="black", family="Helvetica"),
                cohort_size=list(size=14, color="black", family="Helvetica"),
                legend=list(size=12, color="black", family="Helvetica"))

  # draw figure
  figs <- draw_facetted_heatmap_barplots_1(dfs=dfs_plot_agg, 
                                           col_var=col_evt,
                                           names2plots=names2plots,
                                           names2colors=names2colors,
                                           colors_palette_heatmap=colors_palette_heatmap,
                                           alpha_left=0.1, 
                                           alpha_heatmap=0.1,
                                           col_stack=NULL,
                                           width_one=args$output_plot_width_one,
                                           height_one=args$output_plot_height_one,
                                           fonts=fonts,
                                           stacks2colors=NULL,
                                           add_colors_names=F,
                                           x_tick_heatmap_side="top",
                                           x_title_row_bar="Fraction of samples",
                                           add_cohorts_size=T)

  # assemble figures
  fig_left <- figs[["heatmap"]]
  fig_right <- figs[["middle_right_bar"]]
  fig <- subplot(fig_left, fig_right, margin=0.03, widths=c(0.7, 0.3))
  fig <- fig %>% layout(plot_bgcolor='rgba(0,0,0,0)')

  # save
  width <- args$output_plot_width_one*ncol(dfs_plot_agg$percent)
  height <- args$output_plot_height_one*nrow(dfs_plot_agg$percent)
  orca(fig, args$output_plot, width=width, height=height)
  cat(paste("-plot saved at", args$output_plot, "\n"))

  if (!is.null(args$output_plot_paper)){
    orca(fig, args$output_plot_paper, width=width, height=height)
    cat(paste("-plot saved at", args$output_plot_paper, "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw heatmap of chromosome arm events.')
  parser$add_argument("--cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/somatic_cnas/selection/selection_samples_prism.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_met500.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--cna_chr", nargs="+", help="Paths to tables of CNAs per chromsome arm.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz"))
  parser$add_argument("--cna_wgd", nargs="+", help="Paths to tables of CNA summary statistics.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz"))
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/somatic_cnas/selection/selection_tumor_types.tsv")
  parser$add_argument("--min_counts_evt", type="integer", nargs="+", default=c(1,1,1),
                      help="Criteria for building the list of genes in the plot.")
  parser$add_argument("--alterations", type="character", default="gains",
                      help="Type of alterations considered in the plot. Choose gains or losses.")
  parser$add_argument("--n_cores", type="integer", default=4, help="Number of cores for parallel computations.")
  parser$add_argument("--output_tables", type="character",
        help="Path to where output tables are saved.",
        default="../../../results/somatic_cnas/heatmap_chr_arm/tables_losses.xlsx")
  parser$add_argument("--output_plot_width_one", type="integer", default=35, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot_height_one", type="integer", default=15, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot", type="character",
                      help="Path to where output facetted heatmap plot is saved.",
                      default="../../../results/somatic_cnas/heatmap_chr_arm/heatmap_gains.pdf")
  parser$add_argument("--output_plot_paper", type="character",
                      help="Path to where output facetted heatmap plot is saved.",
                      default=NULL)
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
