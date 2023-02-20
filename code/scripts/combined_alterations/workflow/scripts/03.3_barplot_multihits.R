# @created: 30 Mar 22
# @modified: 31 Dec 22
# @authors: Yoann Pradat
#
# Draw beautiful stacked bar plots.

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


load_oncogenes_or_tumorsuppressor <- function(gen_table, evt_type){
  df_genes <- load_table(gen_table)

  if (evt_type=="oncogenes"){
    genes <- df_genes %>%
      filter(Is_Oncogene=="Yes" & Inclusion_Level==1) %>%
      pull(Hugo_Symbol)
  } else if (evt_type=="tumorsuppressors"){
    genes <- df_genes %>%
      filter(Is_Tumor_Suppressor_Gene=="Yes" & Inclusion_Level==1) %>%
      pull(Hugo_Symbol)
  } else {
    stop("- Unsupported value of --evt_type. Choose 'oncogenes' or 'tumorsuppressors.")
  }


  genes
}


aggregate_tables_across_cohorts <- function(dfs, col_alt, cohorts_order, tt_keep_order, table_name, how="col"){
  dfs_table <- list()
  for (cohort in cohorts_order){
    if (how=="col") df_tab <- dfs[[cohort]][[table_name]] %>% column_to_rownames(var=col_alt)
    if (how=="row") df_tab <- dfs[[cohort]][[table_name]] %>% t()
    for (tt_keep_mis in  setdiff(tt_keep_order, colnames(df_tab))){
      df_tab[,tt_keep_mis] <- 0
    }
    df_tab <- as.data.frame(df_tab[,tt_keep_order,drop=F])
    dfs_table[[cohort]] <- df_tab
  }
  
  df_table <- data.frame(row.names=rownames(dfs_table[[1]]))
  for (tt in tt_keep_order){
    for (cohort in cohorts_order){
      df_table_cohort_tt <- dfs_table[[cohort]][tt] %>% rename(!!paste(tt, "-", toupper(cohort)):=.data[[tt]])
      df_table <- bind_cols(df_table, df_table_cohort_tt)
    }
  }

  if (how=="col") df_table <- df_table %>% rownames_to_column(var=col_alt)
  if (how=="row") df_table <- t(df_table)

  df_table
}


load_cnas_loh <- function(df_alt, cna_tables, cohorts){
  cna_keep <- c("cnLOH", "LOH")
  dfs_cna_loh <- list()

  for (i in 1:length(cohorts)){
    cohort <- cohorts[i]
    cna_table <- cna_tables[i]
    cat(paste("-loading CNAs for", cohort, "..."))
    df_cln <- load_cln(study=cohort)
    df_cna <- load_table(cna_table)
    cat("done!\n")
    
    # select in-design samples
    df_cln <- df_cln %>% filter(!is.na(Sample_Id_DNA_T))
    df_cln <- df_cln %>% unite(DNA_P, Sample_Id_DNA_T, Sample_Id_DNA_N, remove=F)
    df_cna <- load_table(cna_table) %>% unite(DNA_P, Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, remove=F)
    mask <- df_cna$DNA_P %in% df_cln$DNA_P
    cat(paste0("-INFO: selected ", sum(mask), "/", length(mask), " lines from selection of in-design DNA_P\n"))
    df_cna <- df_cna %>% filter(mask)

    df_cna_loh <- df_cna %>% filter(Hugo_Symbol %in% df_alt$Hugo_Symbol) %>%
      filter(Copy_Number_More %in% cna_keep) %>% rename(Sample_Id=DNA_P) %>%
      select(Sample_Id, Hugo_Symbol) %>%  distinct() %>% 
      mutate(Alteration_Category="LOH", Alteration="LOH", Cohort=cohort)

    df_cln <- df_cln %>% rename(Sample_Id=DNA_P, Tumor_Type=Project_TCGA_More) %>%
      select(Sample_Id, Subject_Id, Tumor_Type, MSKCC_Oncotree, Civic_Disease)
    df_cna_loh <- left_join(df_cna_loh, df_cln, by="Sample_Id")

    dfs_cna_loh[[cohort]] <- df_cna_loh 
  }

  bind_rows(dfs_cna_loh)
}


get_tabs_count_percent <- function(df_alt_cohort, col_alt, col_id, col_gen, alt_keep){
  # count
  df_alt_cohort_a <- df_alt_cohort %>% filter(.data[[col_alt]]=="Mutation")
  df_alt_cohort_b <- df_alt_cohort %>% filter(.data[[col_alt]]!="Mutation")

  df_alt_cohort_a <- df_alt_cohort_a %>% select(all_of(c(col_id, col_alt, col_gen))) %>% 
      group_by_at(all_of(c(col_id, col_gen))) %>%
      summarize(!!col_alt:=paste0(sort(.data[[col_alt]]), collapse=" & ")) %>%
      mutate(!!col_alt:=ifelse(.data[[col_alt]]=="Mutation", .data[[col_alt]], "Multiple mutations"))

  df_alt_cohort_b <-  df_alt_cohort_b %>% select(all_of(c(col_id, col_alt, col_gen)))
  df_alt_cohort <- bind_rows(df_alt_cohort_a, df_alt_cohort_b)

  df_count <- df_alt_cohort %>% select(all_of(c(col_id, col_alt, col_gen))) %>% 
    group_by_at(all_of(c(col_id, col_gen))) %>%
    summarize(!!col_alt:=paste0(sort(.data[[col_alt]]), collapse=" & ")) %>% 
    filter(.data[[col_alt]] %in% alt_keep) %>%
    group_by_at(all_of(c(col_gen, col_alt))) %>% summarize(Count=n(), .groups="drop")

  # percent 
  df_percent <- df_count %>% spread(.data[[col_gen]], Count) %>% replace(is.na(.), 0) %>% 
    column_to_rownames(var=col_alt)
  df_percent <- t(t(df_percent) / rowSums(t(df_percent)))
  df_percent <- as.data.frame(df_percent) %>% rownames_to_column(var=col_alt)

  # add missing alt_keep if any
  alt_keep_mis <- setdiff(alt_keep, df_percent[[col_alt]])
  if (length(alt_keep_mis)>0)
    df_percent <- bind_rows(df_percent, tibble(!!col_alt:=alt_keep_mis)) %>% replace(is.na(.), 0)
  df_percent <- df_percent %>% filter(.data[[col_alt]] %in% alt_keep) %>% arrange(match(.data[[col_alt]], alt_keep))
  
  list(count=df_count, percent=df_percent) 
}


main <- function(args){
  # load alterations
  df_alt <- load_table(args$alt_table)

  # select samples (only in-design samples from patients with both DNA & RNA) and for selected tumor types
  dfs_sam <- setNames(lapply(args$sam_tables, function(x) load_table(x)), args$cohorts)
  df_sam <- Reduce(bind_rows, dfs_sam)
  df_alt <- select_samples(df_alt, df_sam)

  # load cnas to get LOH events
  df_cna_loh <- load_cnas_loh(df_alt=df_alt, cna_tables=args$cna_tables, cohorts=args$cohorts)
  
  # concatenate alterations and LOH events
  df_alt <- bind_rows(df_alt, df_cna_loh)

  col_gen <- "Hugo_Symbol"
  col_alt <- "Alteration_Category"
  col_id <- "Subject_Id"

  # select genes to keep
  genes_keep <- load_oncogenes_or_tumorsuppressor(args$gen_table, args$gene_type)
  df_alt <- df_alt %>% filter(Hugo_Symbol %in% genes_keep)

  # simplify alterations
  df_alt <- df_alt %>%
    mutate(Alteration_Category=ifelse(Alteration_Category %in% c("Del","Ins","Mut"), "Mutation",  Alteration_Category))

  # select alterations to keep
  if (args$gene_type=="oncogenes"){
    alt_keep <- c("Amplification & Multiple mutations", "Amplification & Mutation", "Multiple mutations",
                  "Amplification", "Mutation")
    alt_selc <- alt_keep
  } else if (args$gene_type=="tumorsuppressors"){
    alt_keep <- c("LOH & Mutation", "LOH & Multiple mutations", "Multiple mutations", "Deletion", "Mutation")
    alt_selc <- union(alt_keep, c("LOH"))
  }

  df_alt <- df_alt %>% filter(Alteration_Category %in% alt_selc)

  # draw per tumor type
  df_alt <- bind_rows(df_alt, df_alt %>% mutate(Tumor_Type="All"))
  n_tumor_types <- length(args$tumor_types)
  df_counts <- load_table(args$counts)

  for (i in 1:n_tumor_types){
    tumor_type <- args$tumor_types[i]
    output_table <- args$output_tables[i]
    output_plot <- args$output_plots[i]
    output_paper <- args$output_papers[i]

    # For each cohort, compute tables of event counts/percentages aggregated by Tumor_Type
    cohorts <- args$cohorts
    dfs_plot <- lapply(cohorts, function(cohort) 
                       get_tabs_count_percent(df_alt %>% filter(Tumor_Type==tumor_type, Cohort==cohort),
                                              col_alt=col_alt, col_id=col_id, col_gen=col_gen, alt_keep=alt_keep))
    dfs_plot <- setNames(dfs_plot, cohorts)

    # aggregate tables across cohorts to prepare for the plot
    dfs_for_plot <- list()
    cohorts_order <- args$cohorts
    cohorts_order <- intersect(cohorts_order, unique(df_alt$Cohort))

    # apply threshold selection
    if (tumor_type == "All"){
      size_tumor_type_p <- df_counts %>% filter(Use_heatmap_all==1) %>% pull(var="PRISM_heatmap_all") %>% sum()
      size_tumor_type_m <- df_counts %>% filter(Use_heatmap_all==1) %>% pull(var="MET500_heatmap_all") %>% sum()
      size_tumor_type_t <- df_counts %>% filter(Use_heatmap_all==1) %>% pull(var="TCGA_heatmap_all") %>% sum()
      min_count_p <- max(5, round(0.05*size_tumor_type_p))
      min_count_m <- max(5, round(0.05*size_tumor_type_m))
      min_count_t <- max(5, round(0.05*size_tumor_type_t))
    } else {
      size_tumor_type_p <- df_counts %>% filter(Tumor_Type==tumor_type) %>% pull(var="PRISM_heatmap_all")
      size_tumor_type_m <- df_counts %>% filter(Tumor_Type==tumor_type) %>% pull(var="MET500_heatmap_all")
      size_tumor_type_t <- df_counts %>% filter(Tumor_Type==tumor_type) %>% pull(var="TCGA_heatmap_all")
      min_count_p <- max(5, round(0.10*size_tumor_type_p))
      min_count_m <- max(5, round(0.10*size_tumor_type_m))
      min_count_t <- max(5, round(0.10*size_tumor_type_t))
    }

    genes_keep_order_p <- dfs_plot$prism$count %>% group_by_at(all_of(col_gen)) %>% summarize(Count=sum(Count)) %>%
      arrange(desc(Count)) %>% filter(Count >= min_count_p) %>% pull(var=col_gen)
    genes_keep_order_m <- dfs_plot$met500$count %>% group_by_at(all_of(col_gen)) %>% summarize(Count=sum(Count)) %>%
      arrange(desc(Count)) %>% filter(Count >= min_count_m) %>% pull(var=col_gen)
    genes_keep_order_t <- dfs_plot$tcga$count %>% group_by_at(all_of(col_gen)) %>% summarize(Count=sum(Count)) %>%
      arrange(desc(Count)) %>% filter(Count >= min_count_t) %>% pull(var=col_gen)
    genes_keep_order <- union(union(genes_keep_order_p, genes_keep_order_m), genes_keep_order_t)

    table_names <- c("percent")
    
    for (j in 1:length(table_names)){
      dfs_for_plot[[table_names[[j]]]] <- aggregate_tables_across_cohorts(dfs=dfs_plot, col_alt=col_alt, 
                                                                          cohorts_order=cohorts_order,
                                                                          tt_keep_order=genes_keep_order,
                                                                          table_name=table_names[[j]],
                                                                          how="col")
    }

    # save table
    df_percent <- dfs_for_plot$percent
    write.table(df_percent, output_table, sep="\t", quote=F, row.names=F)

    # get colors for the plot
    if (args$gene_type=="oncogenes"){
      alt2colors <- list("Amplification & Multiple mutations"="#00af56",
                         "Amplification & Mutation"="#6fbfdd",
                         "Multiple mutations"="#d4a8d0",
                         "Amplification"="#ffe0c3",
                         "Mutation"="#dadada")
    } else {
      alt2colors <- list("LOH & Multiple mutations"="#00af56",
                         "LOH & Mutation"="#6fbfdd",
                         "Multiple mutations"="#d4a8d0",
                         "Deletion"="#ffe0c3",
                         "Mutation"="#dadada")
    }

    cohorts2colors <- load_colors(sheet="Global")[c("prism", "met500", "tcga")]
    names(cohorts2colors) <- as.vector(sapply(names(cohorts2colors), toupper))

    # draw the barplot
    dfs <- dfs_for_plot
    fonts <- list(x_tick_heatmap=list(size=8),
                  y_tick_heatmap=list(size=8),
                  x_tick_row_bar=list(size=7),
                  legend=list(size=6),
                  y_title_row_bar=list(size=9),
                  legend_label=list(size=7),
                  legend_title=list(size=12))
    linewidth <- 2

    # Shared x and y vectors
    x <- colnames(dfs$percent[,2:ncol(dfs$percent)]) 
    x_ann <- rprismtools:::add_colors_names(x, cohorts2colors, pos=2)

    # width, height and legend
    if (tumor_type=="All"){
      if (args$gene_type=="tumorsuppressors"){
        legend_x <- c(-0.45, -0.45)
        left_margin <- 160
      } else {
        legend_x <- c(-0.8, -0.8)
        left_margin <- 220
      }
      width <- left_margin + length(genes_keep_order) * args$output_width_one
      height <- args$output_height
    } else {
      legend_x=c(-10, -10)
      width <- 70 + length(genes_keep_order) * args$output_width_one
      height <- args$output_height
    }

    fig <- rprismtools:::draw_barplot_above_htmp(df_plot=dfs$percent, stacks2colors=alt2colors,
                                                 fonts=fonts, linewidth=linewidth,
                                                 x=x, x_ann=x_ann, xaxis_showticklabels=T,
                                                 legend_title="Events",
                                                 legend_x=legend_x, legend_x_to_y_ratio=0.3, legend_y_gap=0.06)


    if (tumor_type=="All"){
      fig <- fig %>% layout(margin=list(l=left_margin))
    }
    save_image(fig, output_plot, width=width, height=height)
    cat(paste("-plot saved at", output_plot, "\n"))

    if (!is.null(output_paper)){
      save_image(fig, output_paper, width=width, height=height)
      cat(paste("-plot saved at", output_paper, "\n"))
    }
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw heatmap of somatic mutations-hit genes.')
  parser$add_argument("--cohorts", nargs="+", help="Cohort names.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--alt_table", help="Path to input table of alterations.",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations.tsv")
  parser$add_argument("--cna_tables", nargs="+", help="Paths to input tables of CNA calls per gene.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls_filters.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls_filters.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls_filters.tsv.gz"))
  parser$add_argument("--sam_tables", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/combined_alterations/selection/selection_samples_prism.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_met500.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--gen_table", help="Path to table of genes.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--gene_type", help="Choose either oncogenes or tumorsuppressors",
                      default="tumorsuppressors")
  parser$add_argument("--counts", type="character", help="Table of tumor type counts.",
                      default="../../../results/combined_alterations/selection/selection_tumor_types.tsv")
  parser$add_argument("--tumor_types", type="character", nargs="+", help="Tumor types drawn.",
                      default=c("All", "BRCA", "BLCA", "LUAD", "PAAD", "PRAD"))
  parser$add_argument("--output_width_one", type="integer", default=30, help="Width of plot in pixels.")
  parser$add_argument("--output_height", type="integer", default=200, help="Width of plot in pixels.")
  parser$add_argument("--output_tables", type="character", help="Paths to where tables are saved.",
                  default=c("../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_All.tsv",
                            "../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_BRCA.tsv",
                            "../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_BLCA.tsv",
                            "../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_LUAD.tsv",
                            "../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_PAAD.tsv",
                            "../../../results/combined_alterations/other_plots/table_multihits_tumorsuppressors_PRAD.tsv"))
  parser$add_argument("--output_plots", type="character", nargs="+", help="Paths to where plots are saved.",
                default=c("../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_All.pdf",
                          "../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_BRCA.pdf",
                          "../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_BLCA.pdf",
                          "../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_LUAD.pdf",
                          "../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_PAAD.pdf",
                          "../../../results/combined_alterations/other_plots/barplot_multihits_tumorsuppressors_PRAD.pdf"))
  parser$add_argument("--output_papers", type="character", nargs="+", help="Path to where plot is saved.", default=NULL,
                default=c("../../results/figures_paper/FS9_10_tumorsuppressors_All.pdf",
                          "../../results/figures_paper/FS9_10_tumorsuppressors_BRCA.pdf"
                          "../../results/figures_paper/FS9_10_tumorsuppressors_BLCA.pdf"
                          "../../results/figures_paper/FS9_10_tumorsuppressors_LUAD.pdf"
                          "../../results/figures_paper/FS9_10_tumorsuppressors_PAAD.pdf"
                          "../../results/figures_paper/FS9_10_tumorsuppressors_PRAD.pdf"))
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
