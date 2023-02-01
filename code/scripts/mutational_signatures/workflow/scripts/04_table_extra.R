# @created: 03 Feb 22
# @modified: 27 May 22
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

compute_scale_color <- function(df_sig, df_cln, cohort, target_size, tt_keep, col_tt="Project_TCGA_More",
                                col_sid="Subject_Id"){ 
  df_sig <- df_sig %>% column_to_rownames(var="Signature")
  df_sig <- as.data.frame(t(df_sig))
  df_sig[,col_sid] <- rownames(df_sig)
  df_sig <- as_tibble(df_sig)

  # add tumor type (Project_TCGA_More) and select
  df_sig <- left_join(df_sig, df_cln[,c(col_sid, col_tt)])
  df_sig <- df_sig %>%
    filter(.data[[col_tt]] %in% tt_keep)

  # add tumor type (Project_TCGA_More) and select
  # df_sig <- df_sig %>%
  #   filter(.data[[col_tt]] %in% tt_keep) %>%
  #   group_by(.data[[col_tt]]) %>%
  #   mutate(n=n())
  # df_sig <- df_sig[df_sig$n >= min_tt_size,]
  # df_sig$n <- NULL

  # dscale for plot
  dscale <- df_sig %>%
    group_by(.data[[col_tt]]) %>%
    mutate(n=n()) %>%
    summarize_at(vars(-{{col_sid}}), ~sum(.x>0)) %>%
    mutate_at(vars(-{{col_tt}},-n), ~./n) %>%
    mutate(!!col_tt:=paste0(.data[[col_tt]], paste(" -", toupper(cohort))))

  # dcolor for plot
  dcolor <- df_sig %>%
    group_by(.data[[col_tt]]) %>%
    summarize_at(vars(-{{col_sid}}), ~median(.[.!=0]*1e6/target_size)) %>%
    replace(is.na(.),0) %>%
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

compute_pvals_rowwise <- function(df_table_a, df_table_b, test="mann-whitney-u", cores=1){
  stopifnot(identical(rownames(df_table_a), rownames(df_table_b)))
  n_row <- nrow(df_table_a)

  cl <- makeCluster(cores)
  registerDoParallel(cl)

  pvals <- foreach(i = 1:n_row, .packages=c("Exact"), .combine=c) %dopar% {
    sample_a <- as.numeric(df_table_a[i,,drop=T])
    sample_b <- as.numeric(df_table_b[i,,drop=T])
    if (length(unique(sample_a)) > 1 & length(unique(sample_b)) > 1){
      if (test=="mann-whitney-u"){
        out <- wilcox.test(x=sample_a, y=sample_b,
                           alternative="two.sided",
                           exact=F,
                           correct=T,
                           paired=F)
        pval <- out$p.value
      } else {
        stop("Only 'mann-whitney-u' test is available'")
      }
    } else {
      pval <- -999
    }
    pval
  }
  stopCluster(cl)

  data.frame(pval=pvals, row.names=rownames(df_table_a))
}

subset_table_tumor_type <- function(df_table, df_cln, col_tt, tt_val, col_names="Sample_Id_DNA_T"){
  cols_names <- colnames(df_table)
  tumor_types <- df_cln[match(cols_names, df_cln[[col_names]]), col_tt, drop=T]
  df_table_tt <- df_table[,cols_names[grepl(tt_val, tumor_types)]]
  df_table_tt
}

compute_pvals_nonparametric <- function(dfs_table, dfs_cln, cohort_a, cohort_b, col_tt, tt_vals, cores=1){
  df_pvals <- NULL
  for (tt_val in tt_vals){
    df_table_a <- subset_table_tumor_type(dfs_table[[cohort_a]], dfs_cln[[cohort_a]],  col_tt=col_tt, tt_val=tt_val)
    df_table_b <- subset_table_tumor_type(dfs_table[[cohort_b]], dfs_cln[[cohort_b]],  col_tt=col_tt, tt_val=tt_val)
    pvals_tt <- compute_pvals_rowwise(df_table_a, df_table_b, cores=cores)
    if (is.null(df_pvals)){
      df_pvals <- pvals_tt
    } else {
      df_pvals <- cbind(df_pvals, pvals_tt)
    }
  }
  colnames(df_pvals) <- paste(tt_vals, "-", toupper(cohort_b))

  df_pvals
}


get_plot_data <- function(cohorts, dfs_sig, dfs_cln, df_aetiologies, target_size, tt_keep, show_all_signatures,
                          first_half=F, second_half=F){
  # one table of color scale for each cohort
  col_tt <- "Project_TCGA_More"
  col_sid <- "Sample_Id_DNA_T"
  scales_colors <- list()
  for (cohort in cohorts){
    scale_color <- compute_scale_color(df_sig=dfs_sig[[cohort]], df_cln=dfs_cln[[cohort]], cohort=cohort,
                                       target_size=target_size, tt_keep=tt_keep, col_tt=col_tt, col_sid=col_sid)
    scales_colors[[cohort]] <- scale_color
  }

  # merge dscale
  dscale <- do.call(rbind, lapply(cohorts, function(cohort) scales_colors[[cohort]]$dscale))
  dscale <- dscale %>% rowwise() %>% 
    mutate(Tumor_Type = unlist(str_split(.data[[col_tt]], " - "))[1]) %>% 
    mutate(Cohort = unlist(str_split(.data[[col_tt]], " - "))[2])
  dscale$Cohort <- factor(dscale$Cohort, levels=toupper(cohorts))
  dscale <- dscale %>% arrange(Tumor_Type, Cohort) %>% select(-Cohort, -Tumor_Type)
  cols_more <- list("n="=dscale$n)
  dscale$n <- NULL
  dscale <- column_to_rownames(.data=dscale, var=col_tt)
  dscale <- t(as.matrix(dscale))
  df_aetiologies <- df_aetiologies %>% column_to_rownames(var="Signature")

  # merge dcolor
  dcolor <- do.call(rbind, lapply(cohorts, function(cohort) scales_colors[[cohort]]$dcolor))
  dcolor <- dcolor %>% rowwise() %>% 
    mutate(Tumor_Type = unlist(str_split(.data[[col_tt]], " - "))[1]) %>% 
    mutate(Cohort = unlist(str_split(.data[[col_tt]], " - "))[2])
  dcolor$Cohort <- factor(dcolor$Cohort, levels=toupper(cohorts))
  dcolor <- dcolor %>% arrange(Tumor_Type, Cohort) %>% select(-Cohort, -Tumor_Type)
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

  dscale <- dscale[,cols_order]
  dcolor <- dcolor[,cols_order]

  # pvalues
  dfs_pvals <- list()

  # prism vs tcga
  cohort_a <- "tcga"
  cohort_b <- "prism"
  dfs_pvals[[cohort_b]] <- compute_pvals_nonparametric(dfs_table=dfs_sig, dfs_cln=dfs_cln,
                                                       cohort_a=cohort_a, cohort_b=cohort_b,
                                                       col_tt=col_tt, tt_vals=tt_keep, cores=4)

  # met500 vs tcga
  cohort_a <- "tcga"
  cohort_b <- "met500"
  dfs_pvals[[cohort_b]] <- compute_pvals_nonparametric(dfs_table=dfs_sig, dfs_cln=dfs_cln,
                                                       cohort_a=cohort_a, cohort_b=cohort_b,
                                                       col_tt=col_tt, tt_vals=tt_keep, cores=4)

  # correct for multiple testing
  dfs_qvals <- correct_multiple_testing(dfs_pvals, cohorts=c("prism", "met500"), method="fdr")

  # build dframes from dfs_qvals
  dframes <- lapply(dfs_qvals, function(df) {m <- as.matrix(df); m[m>=0.1] <- 9; m[abs(m)<=0.1] <- 1; m[m==-1] <- 0;
                    rownames(m) <- rownames(dscale); m})
  colors_frames <- load_colors(sheet="Global")[c("prism", "met500")]

  # select relevant rows dscale, dcolor
  if (!show_all_signatures){
    if (first_half){
      rows_keep <- rownames(dscale)[1:round(4*nrow(dscale)/5)]
    } else if (second_half) {
      rows_keep <- rownames(dscale)[(round(4*nrow(dscale)/5)+1):nrow(dscale)]
    } else {
      rows_keep <- union(rownames(dscale)[rowSums(dframes$prism)>0],rownames(dscale)[rowSums(dframes$met500)>0])
      rows_keep <- intersect(rownames(dscale),rows_keep)
    }
    dscale <- dscale[rows_keep,]
    dcolor <- dcolor[rows_keep,]
    dframes <- lapply(dframes, function(d) d[rows_keep,])
  }

  # add rows more
  rows_more <- list("Aetiology"=df_aetiologies[rownames(dscale), "Aetiology"])

  # get theme
  scale_breaks <- seq(from=0, to=1, by=0.1)
  color_palette <- c("#ffc651", "#ffa759", "#ff8962", "#ff6b6b", "#cc6999", "#9968c8", "#6767f8", "#4459ce", "#224ba5",
                     "#013d7c")
  color_breaks <- c(0, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10, 25, 1e6)
  color_bg <- c("#f8f9fa", "#e9ecef")
  base_size <- 12
  core_size <- 5

  theme <- get_plot_theme(dscale, base_size=base_size, scale_breaks=scale_breaks, core_size=core_size,
                          color_palette=color_palette, color_bg=color_bg, color_breaks=color_breaks)

  list(dscale=dscale, dcolor=dcolor, cols_more=cols_more, rows_more=rows_more, theme=theme, dframes=dframes,
       colors_frames=colors_frames)
}
  

main <- function(args){
  # load clinical data
  dfs_cln <- lapply(args$clns, function(cln) load_table(cln))
  dfs_cln <- setNames(dfs_cln, args$cohorts)

  # load signatures data
  dfs_sig <- lapply(args$sigs, function(sig) load_table(sig))
  dfs_sig <- setNames(dfs_sig, args$cohorts)

  # load samples 
  dfs_sam <- lapply(args$sams, function(sam) load_table(sam))
  dfs_sam <- setNames(dfs_sam, args$cohorts)

  # get target size
  df_bed <- read.table(args$target_bed)
  df_bed <- df_bed %>% mutate(Size=V3-V2+1)
  target_size <- sum(df_bed$Size)

  # load signatures aetiologies
  df_aetiologies <- load_table(args$aetiologies)

  # load tumor counts
  df_counts_tt <- load_table(args$counts)

  # if use_all_samples = F, select samples using dfs_sam tables
  if (args$use_all_samples){
    tt_keep <- df_counts_tt %>% arrange(desc(.data[[toupper(args$cohorts[1])]])) %>% pull(Tumor_Type)
    tt_keep <- tt_keep[!grepl("MISC - Not_TCGA|Unknown_Primary|LAML", tt_keep)]
  } else {
    for (cohort in args$cohorts){
      df_sig <- dfs_sig[[cohort]]
      df_sam <- dfs_sam[[cohort]] %>% filter(Use_plotting==1)
      df_sam <- df_sam %>% separate(Sample_Id, c("DNA_T", "DNA_N"), sep="_vs_", remove=F)

      sam_keep <- intersect(df_sam$DNA_T, colnames(df_sig))
      df_sig <- df_sig %>% select(all_of(c("Signature", sam_keep)))

      dfs_sig[[cohort]] <- df_sig
    }
    tt_keep <- df_counts_tt %>% filter(Use_plotting==1) %>% 
      arrange(desc(.data[[paste0(toupper(args$cohorts[1]), "_plotting")]])) %>%
      pull(Tumor_Type)
  }

  # plot figure 1 with all signatures
  plot_data <- get_plot_data(cohorts=args$cohorts, dfs_sig=dfs_sig, dfs_cln=dfs_cln, df_aetiologies=df_aetiologies,
                             target_size=target_size, tt_keep=tt_keep, show_all_signatures=T)

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[1],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors with the signature",
                   dcolor_title_legend="Median mut/Mb due to signature",
                   dframes=plot_data$dframes, colors_frames=plot_data$colors_frames)

  # plot figure 2 with only first half of signatures
  plot_data <- get_plot_data(cohorts=args$cohorts, dfs_sig=dfs_sig, dfs_cln=dfs_cln, df_aetiologies=df_aetiologies,
                             target_size=target_size, tt_keep=tt_keep, show_all_signatures=F, first_half=T)


  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs_paper[1],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors with the signature",
                   dcolor_title_legend="Median mut/Mb due to signature",
                   dframes=plot_data$dframes, colors_frames=plot_data$colors_frames)


  # plot figure 2 with only second half of signatures
  plot_data <- get_plot_data(cohorts=args$cohorts, dfs_sig=dfs_sig, dfs_cln=dfs_cln, df_aetiologies=df_aetiologies,
                             target_size=target_size, tt_keep=tt_keep, show_all_signatures=F, second_half=T)

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs_paper[2],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors with the signature",
                   dcolor_title_legend="Median mut/Mb due to signature",
                   dframes=plot_data$dframes, colors_frames=plot_data$colors_frames)

  # plot figure 3 with only platinum signatures
  plot_data <- get_plot_data(cohorts=args$cohorts, dfs_sig=dfs_sig, dfs_cln=dfs_cln, df_aetiologies=df_aetiologies,
                             target_size=target_size, tt_keep=tt_keep, show_all_signatures=T)

  rows_keep <- c("SBS31", "SBS35")
  cols_keep_old <- paste(tt_keep, "-", "PRISM")
  cols_keep_new <- paste(tt_keep)
  index_col_keep <- match(cols_keep_old, colnames(plot_data$dscale))
  index_row_keep <- match(rows_keep, rownames(plot_data$dscale))

  plot_data$dscale <- plot_data$dscale[rows_keep,cols_keep_old]
  colnames(plot_data$dscale) <- cols_keep_new

  plot_data$dcolor <- plot_data$dcolor[rows_keep,cols_keep_old]
  colnames(plot_data$dcolor) <- cols_keep_new

  plot_data$cols_more[["n="]] <- plot_data$cols_more[["n="]][index_col_keep]
  plot_data$rows_more[["Aetiology"]] <- plot_data$rows_more[["n="]][index_row_keep]
  plot_data$dframes <- list(prism=plot_data$dframes[["prism"]][rows_keep,cols_keep_old])
  colnames(plot_data$dframes[["prism"]]) <- cols_keep_new

  plot_data$theme$legend$show <- F

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs[2],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors with the signature",
                   dcolor_title_legend="Median mut/Mb due to signature",
                   dframes=plot_data$dframes, colors_frames=plot_data$colors_frames)

  draw_table_extra(dscale=plot_data$dscale, theme=plot_data$theme, output=args$outputs_paper[3],
                   dcolor=plot_data$dcolor, cols_more=plot_data$cols_more,
                   rows_more=plot_data$rows_more, margin_x=unit(3, "inches"),
                   dscale_title_legend="Proportion of tumors with the signature",
                   dcolor_title_legend="Median mut/Mb due to signature",
                   dframes=plot_data$dframes, colors_frames=plot_data$colors_frames)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder <- "../../../results/mutational_signatures/projection_known_signatures/MutationalPatterns"
  basename_sig <- "counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia"

  parser <- ArgumentParser(description='Extra heatmap showing activity of mutational signatures per tumor type.')
  parser$add_argument("--cohorts", type="character", nargs="+", default=c("prism", "met500", "tcga"),
                      help="Cohort names.")
  parser$add_argument("--clns", type="character", nargs="+", help="Path to clinical tables.",
                      default=c("../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv",
                                "../../../data/met500/clinical/curated/cln_met500_in_design_curated.tsv",
                                "../../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv"))
  parser$add_argument("--sigs", type="character", nargs="+", help="Path to signature tables.",
    default=c(file.path(default_folder, paste(basename_sig, "prism.tsv", sep="_")),
              file.path(default_folder, paste(basename_sig, "met500.tsv", sep="_")),
              file.path(default_folder, paste(basename_sig, "tcga.tsv", sep="_"))))
  parser$add_argument("--sams", nargs="+", help="Paths to sample tables.",
                  default=c("../../../results/mutational_signatures/selection/selection_samples_prism.tsv",
                            "../../../results/mutational_signatures/selection/selection_samples_met500.tsv",
                            "../../../results/mutational_signatures/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/mutational_signatures/selection/selection_tumor_types.tsv")
  parser$add_argument("--target_bed", type="character", help="Path to target bed file.",
                      default="../../../data/resources/target_files/all_targets_intersect_padded_10n.bed")
  parser$add_argument("--aetiologies", type="character", help="Path to table of signature aetiologies.",
                      default="resources/signatures_cosmic/signatures_cosmic_aetiologies_v3.2.tsv")
  parser$add_argument("--use_all_samples", action="store_true", default=F, help="If specified, all samples are used.")
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plots.",
  default=c(file.path(default_folder,"table_extra_signatures_cosmic_sbs_96_min_mut_v3.2_sparse_sigprofilerjulia_full.pdf"),
            file.path(default_folder,"table_extra_signatures_cosmic_sbs_96_min_mut_v3.2_sparse_sigprofilerjulia_plat.pdf")))
  parser$add_argument("--outputs_paper", type="character", nargs="+", help="Path to output plots for paper.",
    default=c("../../../results/figures_paper/FS4a.pdf",
              "../../../results/figures_paper/FS4b.pdf",
              "../../../results/figures_paper/F2b_top.pdf"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
