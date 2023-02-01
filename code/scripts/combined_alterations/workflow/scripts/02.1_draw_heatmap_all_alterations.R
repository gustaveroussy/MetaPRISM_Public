# @created: 10 Mar 22
# @modified: 22 Aug 22
# @authors: Yoann Pradat
#
# Draw heatmap like Figure 7 of Sanchez-Vega et al. 2018

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(xlsx))

# functions ============================================================================================================

select_samples <- function(df_agg, df_sam){
  cat("-selecting samples ...\n")
  df_sam <- df_sam %>% filter(Use_heatmap_all==1)
  n_row_bef <- nrow(df_agg)
  df_agg <- df_agg %>% filter(Sample_Id %in% df_sam$Sample_Id)
  n_row_aft <- nrow(df_agg)
  cat("-INFO: selected", paste0(n_row_aft, "/", n_row_bef), "lines by selecting samples\n")
  cat("-INFO: number of unique subjects after selection:", df_agg %>% distinct(Subject_Id) %>% nrow(), "\n")
  df_agg
}


get_alt_keep <- function(df_alt, cohorts_keep, col_lvl, col_gen, col_alt, col_id, min_counts, tiers_keep=c(), 
                         force_complete_gene=F){
  alt_keep <- c()

  for (cohort in cohorts_keep){
    df_count <- df_alt %>% filter(Cohort==cohort) %>%
      group_by_at(all_of(c(col_id, col_alt))) %>%
      summarize(Count=1, .groups="keep") %>%
      filter(!grepl(" Oncogenic, no Level$", .data[[col_alt]])) %>%
      group_by_at(all_of(c(col_alt))) %>%
      summarize(Count=sum(Count), .groups="drop")

    alt_keep_cohort <- df_count %>% filter(Count >= min_counts[[cohort]]) %>%
      pull(var=col_alt)
    alt_keep <- union(alt_keep, alt_keep_cohort)
  }

  for (tier in tiers_keep){
    alt_keep_tier <- df_alt %>% select(all_of(c(col_alt, col_lvl))) %>%
      distinct() %>% filter(.data[[col_lvl]]==tier) %>%
      pull(var=col_alt)

    alt_keep <- union(alt_keep, alt_keep_tier)
  }

  if (force_complete_gene){
    genes_alt_keep <- df_alt %>% filter(.data[[col_alt]] %in% alt_keep) %>%
      select(.data[[col_gen]]) %>% distinct() %>% pull(var=col_gen)

    alt_keep_genes <- df_alt %>% select(all_of(c(col_alt, col_gen))) %>%
      distinct() %>% filter(.data[[col_gen]] %in% genes_alt_keep) %>%
      filter(!grepl(" Oncogenic, no Level$", .data[[col_alt]])) %>%
      pull(var=col_alt)

    alt_keep <- union(alt_keep, alt_keep_genes)
  }

  sort(alt_keep)
}


get_tables_for_plot_cohort <- function(df_alt_cohort, df_tt_count_cohort, alt_keep, tt_keep, col_lvl, col_alt, col_tt, 
                                       col_id){
  # count for each tumor type
  df_count_tt <- df_tt_count_cohort %>%
    filter(.data[[col_tt]] %in% tt_keep) %>%
    select(.data[[col_tt]], Count) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var=col_tt)

  # counts alt
  # - consider only tumor types in tt_keep
  # - consider only alterations in alt_keep
  # - make sure to count each alteration only once per patient
  df_count_alt <- df_alt_cohort %>%
    filter(.data[[col_tt]] %in% tt_keep) %>%
    filter(.data[[col_alt]] %in% alt_keep) %>%
    group_by_at(all_of(c(col_id, col_alt, col_tt))) %>%
    summarize(Count=1, .groups="drop") %>%
    group_by_at(all_of(c(col_tt, col_alt))) %>%
    summarize(Count=sum(Count), .groups="keep") %>% 
    spread(.data[[col_tt]], Count) %>%
    replace(is.na(.), 0) %>%
    ungroup(.data[[col_alt]])


  # levels alt
  # - consider only tumor types in tt_keep
  # - consider only alterations in alt_keep
  df_level_alt <- df_alt_cohort %>%
    filter(.data[[col_tt]] %in% tt_keep) %>%
    filter(.data[[col_alt]] %in% alt_keep) %>%
    select(all_of(c(col_alt, col_tt, col_lvl))) %>%
    distinct() %>%
    spread(.data[[col_tt]], .data[[col_lvl]])


  # if some alterations are missing, add rows with 0
  alt_keep_mis <- setdiff(alt_keep, df_count_alt[[col_alt]])
  for (alt_mis in alt_keep_mis){
    df_count_alt <- bind_rows(df_count_alt, tibble(!!col_alt:=alt_mis)) %>% replace(is.na(.), 0)
    df_level_alt <- bind_rows(df_level_alt, tibble(!!col_alt:=alt_mis))
  }

  # if some tumor types are missing, add columns with 0
  tt_keep_mis <- setdiff(tt_keep, colnames(df_count_alt))
  for (tt_mis in tt_keep_mis){
    df_count_alt <- df_count_alt %>% mutate(!!tt_mis:=0)
    df_level_alt <- df_level_alt %>% mutate(!!tt_mis:=NA)
  }

  # set columns and rows to correct order
  df_count_alt <- df_count_alt %>% column_to_rownames(var=col_alt)
  df_count_alt <- df_count_alt[alt_keep, tt_keep] %>% rownames_to_column(var=col_alt)
  df_level_alt <- df_level_alt %>% column_to_rownames(var=col_alt)
  df_level_alt <- df_level_alt[alt_keep, tt_keep] %>% rownames_to_column(var=col_alt)

  # percent alt
  df_percent_alt <- df_count_alt %>% mutate(across(c(all_of(tt_keep)), function(x)  x/df_count_tt[cur_column(),]))

  # total count alt
  df_tot_alt <- df_count_alt %>%
    mutate(Count=rowSums(across(where(is.numeric))), Percent=Count/sum(df_count_tt$Count)) %>%
    select(c(.data[[col_alt]], Count, Percent))

  # count of samples having at least one tier1, one tier2 or one tier3 events.
  df_count_best_lvl <- df_alt_cohort %>% select(all_of(c(col_id, col_tt, col_lvl))) %>% distinct() %>% 
    replace(is.na(.), "Tier4") %>%
    group_by_at(all_of(col_id)) %>% slice_min(.data[[col_lvl]]) %>% ungroup() %>% 
    group_by_at(all_of(c(col_tt, col_lvl))) %>% summarize(Count=n(), .groups="drop") %>%
    replace(.=="Tier4", "No level")

  # percent of samples having at least one tier1, one tier2 or one tier3 events.
  df_percent_best_lvl <- df_count_best_lvl %>% spread(Tumor_Type, Count) %>% replace(is.na(.), 0) %>%
    mutate(across(c(all_of(tt_keep)), function(x)  x/df_count_tt[cur_column(),]))

  percent_no_oncogenic_evt <- 1 - colSums(df_percent_best_lvl %>% column_to_rownames(var=col_lvl))
  df_percent_no_oncogenic_evt <-  bind_cols(tibble(!!col_lvl:="No event"), 
                                              data.frame(as.list(percent_no_oncogenic_evt)))
  df_percent_best_lvl <- bind_rows(df_percent_best_lvl, df_percent_no_oncogenic_evt)

  # count of samples stratified by the number of actionable alterations in each sample
  df_count_number <- df_alt_cohort %>% filter(!is.na(.data[[col_lvl]])) %>%
    select(all_of(c(col_id, col_tt, col_alt))) %>% distinct() %>%
    group_by_at(all_of(c(col_id, col_tt))) %>% summarize(N_Alteration=n(), .groups="drop") %>%
    mutate(N_Alteration=ifelse(N_Alteration>=4, ">=4", as.character(N_Alteration))) %>%
    group_by_at(all_of(c(col_tt, "N_Alteration"))) %>% summarize(Count=n(), .groups="drop")

  # add missing tumor types
  for (tt_mis in tt_keep_mis){
    df_count_number <- bind_rows(df_count_number, tibble(!!col_tt:=tt_mis, N_Alteration="1", Count=0))
  }

  # percent of samples stratified by the number of actionable alterations in each sample
  df_percent_number <- df_count_number %>% spread(Tumor_Type, Count) %>% replace(is.na(.), 0) %>%
    mutate(across(c(all_of(tt_keep)), function(x)  x/df_count_tt[cur_column(),]))

  # add missing numbers of alterations
  nalt_keep <- c("1", "2", "3", ">=4")
  nalt_keep_mis <- setdiff(nalt_keep, df_percent_number[["N_Alteration"]])
  if (length(nalt_keep_mis)>0)
    df_percent_number <- bind_rows(df_percent_number, tibble(N_Alteration=nalt_keep_mis)) %>% replace(is.na(.), 0)
  df_percent_number <- df_percent_number %>% arrange(match(.data[["N_Alteration"]], nalt_keep))

  percent_no_oncogenic_evt <- 1 - colSums(df_percent_number %>% column_to_rownames(var="N_Alteration"))
  df_percent_no_oncogenic_evt <-  bind_cols(tibble(N_Alteration="0"), 
                                              data.frame(as.list(percent_no_oncogenic_evt)))
  df_percent_number <- bind_rows(df_percent_number, df_percent_no_oncogenic_evt)

  list(count_col=df_count_tt,
       count_row=df_tot_alt,
       count=df_count_alt,
       level=df_level_alt,
       percent=df_percent_alt,
       percent_col_top=df_percent_best_lvl,
       percent_col_bot=df_percent_number)
}


get_tables_for_plot <- function(df_alt, df_tt_count, cohorts_keep, alt_keep, tt_keep, col_lvl, col_alt, col_tt, col_id, 
                                cohorts_a=c("prism", "met500"), cohorts_b=c("tcga", "tcga"), n_cores=4){
  dfs_plot <- list()

  for (cohort in cohorts_keep){
    df_alt_cohort <- df_alt %>% filter(Cohort==cohort)
    df_tt_count_cohort <- df_tt_count %>% rename(Count=.data[[paste0(toupper(cohort), "_heatmap_all")]]) %>%
      select(all_of(col_tt), Count)
    dfs_plot[[cohort]] <- get_tables_for_plot_cohort(df_alt_cohort=df_alt_cohort, 
                                                     df_tt_count_cohort=df_tt_count_cohort,
                                                     alt_keep=alt_keep, tt_keep=tt_keep,
                                                     col_lvl=col_lvl, col_alt=col_alt,
                                                     col_tt=col_tt, col_id=col_id)
  }

  # differential testing per tumor type
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b,
                               test="fisher", col_evt=col_alt, in_plot_evt=alt_keep, tt_keep=tt_keep,
                               n_cores=n_cores, suffix="")

  # differential testing across tumor types (stratified test)
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b,
                               test="cmh", col_evt=col_alt, in_plot_evt=alt_keep, tt_keep=tt_keep,
                               n_cores=n_cores, suffix="_strats")

  dfs_plot
}


save_tables_xlsx <- function(dfs_plot, col_tt, filepath){
  df_desc <- data.frame(Supplementary=c("Tables for facetted heatmap"))
  write.xlsx(df_desc, file=filepath, sheetName="description", row.names=FALSE)

  for (name_a in names(dfs_plot)){
    for (name_b in names(dfs_plot[[name_a]])){
      sheet_name <- paste(name_a, name_b, sep="_")
      df_name <- dfs_plot[[name_a]][[name_b]]

      if (name_b=="count_col") df_name <- df_name %>% rownames_to_column(var=col_tt)
      write.xlsx(as.data.frame(df_name), file=filepath, sheetName=sheet_name, row.names=FALSE, append=T)
    }
  }

  cat(paste("-INFO: excel workbook saved at", filepath), "\n")
}


aggregate_tables_across_cohorts <- function(dfs, cohorts_keep_order, tt_keep_order, table_name, how="col"){
  dfs_table <- list()
  
  if (how=="col"){
    col_alt <- unique(unlist(sapply(dfs, function(df) colnames(df[[table_name]])[1])))
    stopifnot(length(col_alt)==1)
  }

  for (cohort in cohorts_keep_order){
    if (how=="col") df_tab <- dfs[[cohort]][[table_name]] %>% column_to_rownames(var=col_alt)
    if (how=="row") df_tab <- dfs[[cohort]][[table_name]] %>% t()
    df_tab <- as.data.frame(df_tab[,tt_keep_order,drop=F])
    dfs_table[[cohort]] <- df_tab
  }
  
  df_table <- data.frame(row.names=rownames(dfs_table[[1]]))
  for (tt in tt_keep_order){
    for (cohort in cohorts_keep_order){
      df_table_cohort_tt <- dfs_table[[cohort]][tt] %>% rename(!!paste(tt, "-", toupper(cohort)):=.data[[tt]])
      df_table <- bind_cols(df_table, df_table_cohort_tt)
    }
  }

  if (how=="col") df_table <- df_table %>% rownames_to_column(var=col_alt)
  if (how=="row") df_table <- t(df_table)

  df_table
}


get_tables_for_plot_agg <- function(dfs_plot, cohorts_keep_order, tt_keep_order){
  cat("-INFO: aggregating table across cohorts\n")
  dfs_plot_agg <- list()
  table_names <- lapply(cohorts_keep_order, function(cohort) names(dfs_plot[[cohort]]))
  table_names <- table_names[[1]]
  stopifnot(all(sapply(cohorts_keep_order, function(cohort) setequal(names(dfs_plot[[cohort]]), table_names))))

  table_names <- setdiff(table_names, c("count_row"))

  for (table_name in table_names){
    cat(paste("--processing", table_name,"\n"))

    if (table_name=="count_col"){
      how <- "row"
    } else {
      how <- "col"
    }

    dfs_plot_agg[[table_name]] <- aggregate_tables_across_cohorts(dfs=dfs_plot, 
                                                                  cohorts_keep_order=cohorts_keep_order,
                                                                  tt_keep_order=tt_keep_order,
                                                                  table_name=table_name,
                                                                  how=how)
  }

  # add pvalues from fisher test to dfs_plot_agg
  dfs_qvals <- list("prism"=list("qvals"=dfs_plot[["prism_vs_tcga"]][["qvals"]]),
                    "tcga"=list("qvals"=dfs_plot[["prism_vs_tcga"]][["qvals"]] %>% 
                                mutate_if(is.numeric, function(x) 1)),
                    "met500"=list("qvals"=dfs_plot[["met500_vs_tcga"]][["qvals"]] %>% 
                                  mutate_if(is.numeric, function(x) 1)))

  dfs_plot_agg[["prism_vs_tcga_qvals"]] <- aggregate_tables_across_cohorts(dfs=dfs_qvals,
                                                                           cohorts_keep_order=cohorts_keep_order,
                                                                           tt_keep_order=tt_keep_order,
                                                                           table_name="qvals",
                                                                           how="col")
  cat("--processed prism_vs_tcga_qvals\n")


  dfs_qvals <- list("prism"=list("qvals"=dfs_plot[["prism_vs_tcga"]][["qvals"]] %>%
                                mutate_if(is.numeric, function(x) 1)),
                    "tcga"=list("qvals"=dfs_plot[["prism_vs_tcga"]][["qvals"]] %>% 
                                mutate_if(is.numeric, function(x) 1)),
                    "met500"=list("qvals"=dfs_plot[["met500_vs_tcga"]][["qvals"]]))

  dfs_plot_agg[["met500_vs_tcga_qvals"]] <- aggregate_tables_across_cohorts(dfs=dfs_qvals,
                                                                           cohorts_keep_order=cohorts_keep_order,
                                                                           tt_keep_order=tt_keep_order,
                                                                           table_name="qvals",
                                                                           how="col")
  cat("--processed met500_vs_tcga_qvals\n")

  # add pvalues from stratified test to dfs_plot_agg
  dfs_plot_agg[["prism_vs_tcga_qvals_strats"]] <- dfs_plot[["prism_vs_tcga"]][["qvals_strats"]]
  cat("--processed prism_vs_tcga_qvals_strats\n")
  dfs_plot_agg[["met500_vs_tcga_qvals_strats"]] <- dfs_plot[["met500_vs_tcga"]][["qvals_strats"]]
  cat("--processed met500_vs_tcga_qvals_strats\n")

  dfs_plot_agg
}


main <- function(args){
  # load alterations and tt counts
  df_alt <- load_table(args$alterations)
  df_tt_count <- load_table(args$counts)
  df_tt_count <- df_tt_count %>% filter(Use_heatmap_all==1)

  # select samples (only in-design samples from patients with both DNA & RNA) and for selected tumor types
  dfs_sam <- setNames(lapply(args$samples, function(x) load_table(x)), args$cohorts)
  df_sam <- Reduce(bind_rows, dfs_sam)
  df_alt <- select_samples(df_alt, df_sam)

  # define parameters specific to direction
  if (args$direction=="sensitivity"){
    col_lvl <- "Sen_Level_Simple"
    min_counts <- list("tcga"=100, "prism"=12, "met500"=8)

    levels2colors <- list(Tier1="#00af56",
                          Tier2="#6fbfdd",
                          Tier3="#d4a8d0",
                          "No level"="#ffe0c3",
                          "No event"="#dadada")
    numbers2colors <- list(">=4"="#00691f",
                           "3"="#009c36",
                           "2"="#36be61",
                           "1"="#c8deb4",
                           "0"="#f2f3f8")

    legend_titles <- c("ESCAT<br>levels", "Actionable<br>alterations")
  } else if (args$direction=="resistance"){
    col_lvl <- "Res_Level_Simple"
    min_counts <- list("tcga"=10, "prism"=2, "met500"=2)

    levels2colors <- list(Tier1="#DD2D4A",
                          Tier2="#ffb700",
                          Tier3="#e0aaff",
                          "No level"="#ffe0c3",
                          "No event"="#dadada")
    numbers2colors <- list(">=4"="#9d0208",
                           "3"="#d00000",
                           "2"="#e85d04",
                           "1"="#faa307",
                           "0"="#f2f3f8")

    legend_titles <- c("ESCAT<br>levels", "Resistance<br>alterations")
  } else {
    stop(paste("-ERROR: unsupported value", args$direction, "for --direction parameter."))
  }

  # define util variable names
  col_gen <- paste0("Hugo_Symbol_", col_lvl)
  col_alt <- paste0("Alteration_", col_lvl)
  col_tt <- "Tumor_Type"
  col_id <- "Subject_Id"

  # for alterations that are not alreay equalqto gene symbol, add gene symbol as prefix
  mask <- df_alt[[col_alt]] != df_alt[[col_gen]]
  df_alt[mask,col_alt] <- df_alt %>% filter(mask) %>%
    unite("Unite", all_of(c(col_gen, col_alt)), sep=" ") %>% pull(var="Unite")

  # select data for columns
  cohorts_keep <- c("prism", "met500", "tcga")
  tt_keep <- df_tt_count %>% pull(.data[[col_tt]])
  alt_keep <- get_alt_keep(df_alt=df_alt, cohorts_keep=cohorts_keep, col_lvl=col_lvl, col_gen=col_gen,
                           col_alt=col_alt, col_id=col_id, min_counts=min_counts, tiers_keep=c("Tier1"),
                           force_complete_gene=T)

  # get tables for plot split per cohort
  dfs_plot <- get_tables_for_plot(df_alt=df_alt, df_tt_count=df_tt_count, cohorts_keep=cohorts_keep,
                                  alt_keep=alt_keep, tt_keep=tt_keep, col_lvl=col_lvl, col_alt=col_alt,
                                  col_tt=col_tt, col_id=col_id, cohorts_a=c("prism", "met500"),
                                  cohorts_b=c("tcga", "tcga"), n_cores=args$n_cores)

  # save all tables in an excel workbook
  save_tables_xlsx(dfs_plot=dfs_plot, col_tt=col_tt, filepath=args$output_tables)

  # aggregate tables across cohorts to prepare for the plot
  cohorts_keep_order <- c("prism", "met500", "tcga")
  tt_keep_order <- df_tt_count %>% arrange(desc(PRISM_heatmap_all)) %>% pull(.data[[col_tt]])
  dfs_plot_agg <- get_tables_for_plot_agg(dfs_plot=dfs_plot, cohorts_keep_order=cohorts_keep_order, 
                                          tt_keep_order=tt_keep_order)

  # define arguments for the function draw_facetted_heatmap_barplots_2
  cohorts2colors <- load_colors(sheet="Global")[c("prism", "met500")]
  tumors2colors <- load_colors(sheet="Project_TCGA_More")

  names2colors <- list(met500_vs_tcga_qvals_strats=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals_strats=cohorts2colors[["prism"]],
                       met500_vs_tcga_qvals=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals=cohorts2colors[["prism"]],
                       heatmap=levels2colors[1:3],
                       barplot_top=levels2colors,
                       barplot_bot=numbers2colors,
                       colnames=tumors2colors)

  names2plots <- list(heatmap="percent",
                      heatmap_color="level",
                      heatmap_hover="count",
                      heatmap_pvals=c("prism_vs_tcga_qvals", "met500_vs_tcga_qvals"),
                      left_pvals=c("prism_vs_tcga_qvals_strats", "met500_vs_tcga_qvals_strats"),
                      barplot_top="percent_col_top",
                      barplot_bot="percent_col_bot",
                      count_col="count_col")

  fonts <- list(x_tick_heatmap=list(size=10),
                y_tick_heatmap=list(size=10),
                x_tick_row_bar=list(size=10),
                y_title_row_bar=list(size=9),
                legend=list(size=10),
                legend_label=list(size=10),
                legend_title=list(size=12))

  # size ajdustments
  if (args$direction=="resistance"){
    margin <- 0.015
    legend_x_to_y_ratio <- 0.15
  } else {
    margin <- 0.0075
    legend_x_to_y_ratio <- 0.15
  }

  # call the main function draw_facetted_heatmap_barplots_2
  figs <- draw_facetted_heatmap_barplots_2(dfs=dfs_plot_agg,
                                           col_var=col_alt,
                                           names2plots=names2plots,
                                           names2colors=names2colors,
                                           linewidth=2,
                                           width_one=15,
                                           height_one=15,
                                           alpha_left=0.1, 
                                           alpha_heatmap=0.1,
                                           width_edges_heatmap=0.3,
                                           legend_titles=legend_titles,
                                           legend_x_to_y_ratio=legend_x_to_y_ratio,
                                           fonts=fonts)


  height_barplot <- args$output_plot_barplot_height/args$output_plot_height
  fig <- subplot(figs[["barplot_top"]], figs[["barplot_bot"]], figs[["heatmap"]], nrows=3, margin=margin,
                 heights=c(height_barplot/2, height_barplot/2, 1-height_barplot), titleY=T)
  fig <- fig %>% layout(margin=list(l=250))

  # save
  save_image(fig, args$output_plot, width=args$output_plot_width, height=args$output_plot_height)
  cat(paste("-plot saved at", args$output_plot, "\n"))

  if (!is.null(args$output_plot_paper)){
    save_image(fig, args$output_plot_paper, width=args$output_plot_width, height=args$output_plot_height)
    cat(paste("-plot saved at", args$output_plot_paper, "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw heatmap of somatic mutations-hit genes.')
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/combined_alterations/selection/selection_samples_met500.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_prism.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--alterations", help="Paths to input table of alterations.",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations.tsv")
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/combined_alterations/selection/selection_tumor_types.tsv")
  parser$add_argument("--direction", help="Choose either resistance or sensitivity.",
                      default="resistance")
  parser$add_argument("--n_cores", type="integer", default=4, help="Number of cores for parallel computations.")
  parser$add_argument("--output_tables", type="character", help="Path to where output tables are saved.",
                      default="../../../results/combined_alterations/heatmap_all/tables_heatmap_resistance.xlsx")
  parser$add_argument("--output_plot_barplot_height", type="integer", default=200,
                      help="Height of combined barplots above the heatmap in pixels.")
  parser$add_argument("--output_plot_width", type="integer", default=900, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot_height", type="integer", default=650, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot", type="character",
                      help="Path to where output facetted heatmap plot is saved.",
                      default="../../../results/combined_alterations/heatmap_all/heatmap_resistance.pdf")
  parser$add_argument("--output_plot_paper", type="character",
                      help="Path to where output facetted heatmap plot is saved.", default=NULL)
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
