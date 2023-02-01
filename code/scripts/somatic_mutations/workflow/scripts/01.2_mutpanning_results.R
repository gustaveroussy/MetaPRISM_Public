# @created: 28 Dec 22
# @modified: 29 Dec 22
# @authors: Yoann Pradat
#
# Extra heatmap for describing the activity of mutational signatures across tumor types and across cohorts.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tableExtra))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))

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

  # select columns
  cols_mut <- c("Hugo_Symbol", "Variant_Classification", "Sample_Id", "Subject_Id", "Tumor_Type", "Start_Position",
                "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
  df_mut[,cols_mut]
}


load_mutpanning_results <- function(mpg, tumor_types, value="FDR"){
  dfs_mpg <- list()
  for (tumor_type in tumor_types){
    fp_mpg <- file.path(mpg, tumor_type, paste0("SignificanceFiltered/Significance", tumor_type, ".txt"))
    if (file.exists(fp_mpg)){
      df_mpg <- load_table(fp_mpg)
      df_mpg <- df_mpg %>% group_by(Name) %>% summarize(!!value:=min(.data[[value]])) %>% arrange(.data[[value]])
      df_mpg <- df_mpg %>% select(Name, all_of(value)) %>%
        rename(!!tumor_type:={{value}})
      dfs_mpg[[tumor_type]] <- df_mpg
    }
  }
  Reduce(function(df1,df2) full_join(df1, df2, by="Name"), dfs_mpg)
}


get_signif_genes <- function(df_mpg, thresh=1e-1){
  tts <- setdiff(names(df_mpg), c("Name"))
  df_mpg[,"Min"] <- sapply(1:nrow(df_mpg), function(i) min(df_mpg[i, tts], na.rm=T))
  df_mpg_thresh <- df_mpg %>% filter(Min<=thresh)
  if (nrow(df_mpg_thresh)==0){
    genes <- c()
  } else {
    genes <- df_mpg_thresh$Name
  }
  genes
}


compute_scale_color <- function(df_mut, df_mpg, df_tt, cohort, gen_keep, col_gen="Hugo_Symbol", tt_keep,
                                col_tt="Tumor_Type", col_sid="Subject_Id"){ 
  # percentage of mutations in each gene
  df_cnt <- df_mut %>% filter(.data[[col_gen]] %in% gen_keep) %>%
    distinct_at(all_of(c(col_sid, col_gen, col_tt))) %>% group_by_at(all_of(c(col_gen, col_tt))) %>%
    summarize(Count=n(), .groups="drop")
  df_cnt <- left_join(df_cnt, df_tt, how="left", by=col_tt)
  df_pct <- df_cnt %>% mutate(Percent=Count/Count_Total) %>% select(-Count)
  df_pct <- df_pct %>% spread(.data[[col_gen]], Percent) %>% rename(n=Count_Total)

  # add missing tumor types, missing genes, and fill 0s
  df_mis_tt <- tibble(!!col_tt:=setdiff(tt_keep, df_pct[[col_tt]]))
  if (nrow(df_mis_tt)>0) df_pct <- bind_rows(df_pct, df_mis_tt)

  df_mis_gen <- tibble(.rows=nrow(df_pct))
  for (gen in setdiff(gen_keep, colnames(df_pct))){
    df_mis_gen[,gen] <- NA
  }
  if (ncol(df_mis_gen)>0) df_pct <- bind_cols(df_pct, df_mis_gen)
  df_pct <- df_pct %>% replace(is.na(.), 0)

  # dscale for plot
  dscale <- df_pct %>%
    mutate(!!col_tt:=paste0(.data[[col_tt]], paste(" -", toupper(cohort))))

  # pivot mutpanning results table
  df_col <- df_mpg %>% filter(Name %in% gen_keep) %>% column_to_rownames(var="Name") %>% t()
  df_col <- as.data.frame(df_col)
  df_col[df_col <= 1e-10] <- 1e-10+1e-11
  df_col <- -log10(df_col)
  df_col <- bind_cols(tibble(!!col_tt:=rownames(df_col)), as_tibble(df_col))

  # add missing tts and fill 0s
  df_mis_tt <- tibble(!!col_tt:=setdiff(tt_keep, df_col[[col_tt]]))
  if (nrow(df_mis_tt)>0) df_col <- bind_rows(df_col, df_mis_tt)

  df_mis_gen <- tibble(.rows=nrow(df_col))
  for (gen in setdiff(gen_keep, colnames(df_col))){
    df_mis_gen[,gen] <- NA
  }
  if (ncol(df_mis_gen)>0) df_col <- bind_cols(df_col, df_mis_gen)

  df_col <- bind_rows(df_col, tibble(!!col_tt:=setdiff(tt_keep, df_col[[col_tt]])))
  df_col <- df_col %>% replace(is.na(.), 0)

  # dcolor for plot
  dcolor <- df_col %>%
    mutate(!!col_tt:=paste0(.data[[col_tt]], paste(" -", toupper(cohort))))

  list(dscale=dscale, dcolor=dcolor)
}


get_plot_theme <- function(dscale, base_size=12, core_size=5, scale_breaks=seq(from=0, to=1, by=0.1),
                           color_palette=NULL, color_bg=NULL,  color_breaks=NULL){
  # graphical parameters

  if (is.null(color_palette)){
    color_palette <- c("#ffc651", "#ffa759", "#ff8962", "#ff6b6b", "#cc6999", "#9968c8", "#6767f8", "#4459ce",
                       "#224ba5", "#013d7c")
  }

  if (is.null(color_bg)){
    color_bg <- c("#f8f9fa", "#e9ecef")
  }

  if (is.null(color_breaks)){
    color_breaks <- c(0, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 1e6)
  }
  
  # background colors
  fill_bg <- rep(NA, length.out=ncol(dscale))

  if (class(color_bg)=="list"){
    for (name in names(color_bg)){
      fill_bg[grep(name, colnames(dscale), ignore.case=T)] <- color_bg[[name]]
    }
  } else {
    fill_bg <- color_bg
  }
  

  theme <- ttheme_awesome(base_size=base_size,
                          rep_mode="col",
                          core_size=core_size, 
                          scale_breaks=scale_breaks,
                          color_palette=color_palette, 
                          color_breaks=color_breaks, 
                          legend_position="top_right",
                          core=list(bg_params=list(fill=fill_bg)))

  theme
}


get_plot_data <- function(cohorts, dfs_mut, dfs_tt, dfs_mpg, tt_keep, gen_keep){
  # one table of color scale for each cohort
  col_tt <- "Tumor_Type"
  col_gen <- "Hugo_Symbol"
  col_sid <- "Subject_Id"
  scales_colors <- list()
  for (cohort in cohorts){
    scale_color <- compute_scale_color(df_mut=dfs_mut[[cohort]], df_mpg=dfs_mpg[[cohort]], cohort=cohort,
                                       df_tt=dfs_tt[[cohort]], gen_keep=gen_keep,
                                       col_gen=col_gen, tt_keep=tt_keep, col_tt=col_tt, col_sid=col_sid)
    scales_colors[[cohort]] <- scale_color
  }

  # merge dscale
  dscale <- do.call(rbind, lapply(cohorts, function(cohort) scales_colors[[cohort]]$dscale))
  dscale <- dscale %>% rowwise() %>% 
    mutate(TT = unlist(str_split(.data[[col_tt]], " - "))[1]) %>% 
    mutate(Cohort = unlist(str_split(.data[[col_tt]], " - "))[2])
  dscale$Cohort <- factor(dscale$Cohort, levels=toupper(cohorts))
  dscale <- dscale %>% arrange(TT, Cohort) %>% select(-Cohort, -TT)
  cols_more <- list("n="=dscale$n)
  dscale$n <- NULL
  dscale <- column_to_rownames(.data=dscale, var=col_tt)
  dscale <- t(as.matrix(dscale))

  # merge dcolor
  dcolor <- do.call(rbind, lapply(cohorts, function(cohort) scales_colors[[cohort]]$dcolor))
  dcolor <- dcolor %>% rowwise() %>% 
    mutate(TT = unlist(str_split(.data[[col_tt]], " - "))[1]) %>% 
    mutate(Cohort = unlist(str_split(.data[[col_tt]], " - "))[2])
  dcolor$Cohort <- factor(dcolor$Cohort, levels=toupper(cohorts))
  dcolor <- dcolor %>% arrange(TT, Cohort) %>% select(-Cohort, -TT)
  dcolor <- column_to_rownames(.data=dcolor, var=col_tt)
  dcolor <- t(as.matrix(dcolor))

  # order dscale and dcolor in same order as tt_keep
  cols_order <- c()
  for (tt in tt_keep){
    cols_order <- c(cols_order, colnames(dscale)[grepl(paste0("^", tt, " -"), colnames(dscale))])
  }

  index <- match(cols_order, colnames(dscale))
  for (name in names(cols_more)){
    cols_more[[name]] <- cols_more[[name]][index]
  }

  dscale <- dscale[gen_keep,cols_order]
  dcolor <- dcolor[gen_keep,cols_order]

  # get theme
  scale_breaks <- seq(from=0, to=1, by=0.1)
  color_palette <- c("#dee2e6", "#ffa759", "#ff8962", "#ff6b6b", "#cc6999", "#9968c8", "#6767f8", "#4459ce", "#224ba5",
                     "#013d7c")
  color_breaks <- c(0:10)
  color_bg <- c("#f8f9fa", "#e9ecef")

  base_size <- 12
  core_size <- 5

  theme <- get_plot_theme(dscale, base_size=base_size, scale_breaks=scale_breaks, core_size=core_size,
                          color_palette=color_palette, color_bg=color_bg, color_breaks=color_breaks)

  list(dscale=dscale, dcolor=dcolor, cols_more=cols_more, rows_more=NULL, theme=theme)
}


main <- function(args){
  # load mutations data
  n_coh <- length(args$cohorts)
  dfs_mut <- lapply(1:n_coh, function(i) load_mutations_and_select_samples(cohort=args$cohorts[i],
                                                                           samples=args$sams[i],
                                                                           mutations=args$muts[i]))
  dfs_mut <- setNames(dfs_mut, args$cohorts)

  # load table of counts
  df_counts_tt <- load_table(args$counts)

  # load tumor counts
  df_counts_tt <- load_table(args$counts)
  tt_keep <- df_counts_tt %>% filter(Use_mutpanning==1) %>% 
    arrange(desc(.data[[paste0(toupper(args$cohorts[1]), "_mutpanning")]])) %>%
    pull(Tumor_Type)
  df_counts_tt <- df_counts_tt %>% filter(Tumor_Type %in% tt_keep)
  dfs_tt <- lapply(args$cohorts, function(cohort) {
                     Cohort <- paste0(toupper(cohort), "_mutpanning")
                     df_cohort <- df_counts_tt %>% select(Tumor_Type, Cohort) %>% 
                       rename(Count_Total=Cohort)
                     df_cohort})
  dfs_tt <- setNames(dfs_tt, args$cohorts)

  # load samples 
  dfs_sam <- lapply(args$sams, function(sam) load_table(sam))
  dfs_sam <- setNames(dfs_sam, args$cohorts)

  # load mutpanning results
  dfs_mpg <- lapply(args$mpgs, function(mpg) load_mutpanning_results(mpg, tt_keep))
  dfs_mpg <- setNames(dfs_mpg, args$cohorts)

  # get list of genes to include in plot
  gen_keep <- Reduce(union, lapply(dfs_mpg, function(df_mpg) get_signif_genes(df_mpg, thresh=1e-1)))

  # sort them by frequency
  df_mut_all <- bind_rows(dfs_mut) %>% filter(Hugo_Symbol %in% gen_keep)
  df_mut_all <- df_mut_all %>% group_by(Hugo_Symbol) %>% summarize(Count=n()) %>% arrange(desc(Count))
  gen_keep <- intersect(df_mut_all$Hugo_Symbol, gen_keep)

  # plot data
  plot_data <- get_plot_data(cohorts=args$cohorts, dfs_mut=dfs_mut, dfs_tt=dfs_tt, dfs_mpg=dfs_mpg,
                             tt_keep=tt_keep, gen_keep=gen_keep)

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$output_plot,
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of mutated patients",
                   dcolor_title_legend="-Log10 of mutpanning FDR")

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$output_plot_paper,
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of mutated patients",
                   dcolor_title_legend="-Log10 of mutpanning FDR")


  # table results
  dfs_res <- lapply(args$cohorts, function(cohort) {
                      tt_mis <- setdiff(tt_keep, colnames(dfs_mpg[[cohort]]))
                      for (tt in tt_mis) dfs_mpg[[cohort]][[tt]] <- NA
                      old_names <- c("Name", tt_keep)
                      new_names <- c("Hugo_Symbol", paste(tt_keep, "-",  toupper(cohort)))
                      dfs_mpg[[cohort]] %>% rename_at(vars(old_names), ~ new_names)})
  df_res <- Reduce(function(df1, df2) full_join(df1, df2, by="Hugo_Symbol"), dfs_res)
  df_res <- df_res %>% filter(Hugo_Symbol %in% gen_keep)

  # order cols and rows
  cols_all <- lapply(args$cohorts, function(cohort) paste(tt_keep, "-",  toupper(cohort)))
  cols_ord <- c()
  for (i in 1:length(cols_all[[1]])){
    cols_ord <- c(cols_ord, sapply(cols_all, function(cols) cols[[i]]))
  }

  df_res <- df_res[match(gen_keep, df_res$Hugo_Symbol),c("Hugo_Symbol", cols_ord)]
  write.table(df_res, args$output_data, quote=F, row.names=F, sep="\t")
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Extra heatmap showing fdr values mutpanning per gene and  tumor type.')
  parser$add_argument("--cohorts", type="character", nargs="+", default=c("prism", "met500", "tcga"),
                      help="Cohort names.")
  parser$add_argument("--counts", type="character", help="Tumor types.",
                      default="../../../results/somatic_mutations/selection/selection_tumor_types.tsv")
  parser$add_argument("--sams", nargs="+", help="Paths to sample tables.",
                  default=c("../../../results/somatic_mutations/selection/selection_samples_prism.tsv",
                            "../../../results/somatic_mutations/selection/selection_samples_met500.tsv",
                            "../../../results/somatic_mutations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--muts", type="character", nargs="+", help="Paths to mutations tables.",
                      default=c("../../../data/prism/wes/somatic_maf/somatic_calls.maf.gz",
                                "../../../data/met500/wes/somatic_maf/somatic_calls.maf.gz",
                                "../../../data/tcga/wes/somatic_maf/somatic_calls.maf.gz"))
  parser$add_argument("--mpgs", type="character", nargs="+", help="Paths to folders with mutpanning results.",
                      default=c("../../../results/somatic_mutations/mutpanning/prism",
                                "../../../results/somatic_mutations/mutpanning/met500",
                                "../../../results/somatic_mutations/mutpanning/tcga"))
  parser$add_argument("--output_plot", type="character", help="Path to output plot.",
                      default="../../../results/somatic_mutations/mutpanning/table_extra_mutpanning.pdf")
  parser$add_argument("--output_plot_paper", type="character", help="Path to output plot.",
                      default="../../../results/figures_paper/FS12.pdf")
  parser$add_argument("--output_data", type="character", help="Path to output plot.",
                      default="../../../results/somatic_mutations/mutpanning/results_mutpanning.tsv")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
