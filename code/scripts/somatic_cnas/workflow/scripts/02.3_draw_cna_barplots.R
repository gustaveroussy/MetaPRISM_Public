# @created: 29 Aug 22
# @modified: 02 Jan 23
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# functions ============================================================================================================

select_by_sample_type <- function(df, cohort, sample_type="DNA_T|RNA_T", col_id="Subject_Id", cols_cln=c()){
  df_cln <- load_cln(cohort)
  df_cln <- df_cln %>% rename(Tumor_Type=Project_TCGA_More)
  mask <- grepl(sample_type, df_cln$Sample_Type)
  df_cln_st <- df_cln[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," subjects with ", sample_type, "\n"))

  if (col_id=="Sample_Id"){
    dfs_cln_sub <- list()
    for (st in unlist(str_split(sample_type, "\\|"))){
      col_id_st <- paste0("Sample_Id_", st)
      df_cln_st_sub <- df_cln_st %>% filter(!is.na(.data[[col_id_st]]))
      df_cln_st_sub <- df_cln_st_sub %>% rename(Sample_Id=.data[[col_id_st]])
      dfs_cln_sub[[st]] <- df_cln_st_sub
    }
    df_cln_sub <- bind_rows(dfs_cln_sub)
  } else if (col_id=="DNA_P"){
      cols_sid <- c("Sample_Id_DNA_T", "Sample_Id_DNA_N")
      df_cln_sub <- df_cln_st
      df_cln_sub <- df_cln_st %>% unite(!!col_id, all_of(c(cols_sid)), sep="_vs_", remove=F)
  } else {
    df_cln_sub <- df_cln_st
  }

  mask <- df[[col_id]] %in% df_cln_sub[[col_id]]
  df_st <- df[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," lines from subjects with ", sample_type, "\n"))

  # add columns
  if (length(cols_cln)>0){
    df_st <- left_join(df_st, df_cln_sub[,c(col_id, cols_cln)], by=col_id)
  }

  df_st
}


select_by_selections <- function(df, filepath_sam, name, col_id, col_id_sam="Sample_Id"){
  df_sam <- load_table(filepath_sam)
  df_sam <- df_sam[df_sam[[paste0("Use_", name)]]==1,]
  mask <- df[[col_id]] %in% df_sam[[col_id_sam]]
  df <- df[mask,]
  cat(paste0("-INFO: selected ", sum(mask), "/", length(mask) ," lines from selected samples\n"))
  df
}


main <- function(args){
  # load cna calls
  dfs_cln <- lapply(args$cln_tables, load_table)
  dfs_cln <- setNames(dfs_cln, args$cohorts)
  dfs_cna <- lapply(args$cna_tables, load_table)
  dfs_cna <- setNames(dfs_cna, args$cohorts)

  # add DNA_P
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  dfs_cna$prism <- dfs_cna$prism %>% unite(DNA_P, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
  dfs_cna$met500 <- dfs_cna$met500 %>% unite(DNA_P, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
  dfs_cna$tcga <- dfs_cna$tcga %>% unite(DNA_P, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)

  # select by sample type and in-design
  samples <- setNames(args$samples, args$cohorts)
  dfs_cln <- lapply(args$cohorts, function(cohort)
    select_by_sample_type(df=dfs_cln[[cohort]], cohort=cohort, sample_type="DNA_T", col_id="Subject_Id"))
  dfs_cln <- setNames(dfs_cln, args$cohorts)
  dfs_cln <- lapply(args$cohorts, function(cohort)
    select_by_selections(df=dfs_cln[[cohort]], filepath_sam=samples[[cohort]], name="selections", col_id="Subject_Id", 
                         col_id_sam="Subject_Id"))
  dfs_cln <- setNames(dfs_cln, args$cohorts)

  # select by sample-type and selections
  cols_id <- list(prism="DNA_P", met500="DNA_P", tcga="DNA_P")
  dfs_cna <- lapply(args$cohorts, function(cohort)
    select_by_sample_type(df=dfs_cna[[cohort]], cohort=cohort, sample_type="DNA_T", col_id=cols_id[[cohort]], 
                          cols_cln=c("Subject_Id", "Tumor_Type")))
  dfs_cna <- setNames(dfs_cna, args$cohorts)

  cols_id <- list(prism="DNA_P", met500="DNA_P", tcga="DNA_P")
  dfs_cna <- lapply(args$cohorts, function(cohort)
    select_by_selections(df=dfs_cna[[cohort]], filepath_sam=samples[[cohort]], name="selections",
                         col_id=cols_id[[cohort]], col_id_sam="Sample_Id"))
  dfs_cna <- setNames(dfs_cna, args$cohorts)

  # prepare tables
  dfs_cna <- lapply(args$cohorts, function(cohort) dfs_cna[[cohort]] %>% mutate(Cohort=cohort))
  dfs_cna <- setNames(dfs_cna, args$cohorts)
  df_plot <- bind_rows(dfs_cna) %>% filter(Copy_Number %in% c(-2,2)) %>%
    select(Subject_Id, Cohort, Tumor_Type, Hugo_Symbol, Chromosome, Copy_Number) %>% distinct()
  n_pats <- lapply(args$cohorts, function(cohort) dfs_cln[[cohort]] %>% nrow())
  n_pats <- setNames(n_pats, args$cohorts)
  df_n_pats <- data.frame(Cohort=names(n_pats), Size=unlist(n_pats))

  # amplified genes ====================================================================================================
  df_plot_hlg <- df_plot %>% filter(Copy_Number==2)

  # select genes shown in the plot using prism
  n_pat_p <- n_pats$prism
  df_cna_p_hlg <- dfs_cna$prism %>% filter(Copy_Number==2)
  prop_hlg_p <- df_cna_p_hlg %>% group_by(Hugo_Symbol, Chromosome) %>% summarize(Prop=n()/n_pat_p, .groups="drop") %>%
    arrange(desc(Prop))

  # select only genes mutated in at least 2.5% of prism
  # prop_hlg_p_plot <- prop_hlg_p %>% filter(Prop >= 0.025)

  # simplify the selection taking into account related genes
  related_genes <- list()
  related_genes[["CCND1"]] <- c("FGF3", "PPFIA1", "FADD", "ANO1", "FGF19", "FGF4", "CCND1", "ORAOV1", "MYEOV", "CTTN",
                                "SHANK2", "TPCN2", "MRGPRD", "MRGPRF", "IGHMBP2", "MRPL21", "CPT1A", "MTL5", "GAL",
                                "PPP6R3", "DHCR7", "NADSYN1", "LRP5")
  related_genes[["AR"]] <- c("AR", "OPHN1")
  related_genes[["19q13 genes"]] <- c("NCCRP1", "LRFN1", "LGALS13", "ZFP36", "MED29", "PAF1", "GMFG", "EID2B", "EID2",
                                      "PLEKHG2", "SAMD4B", "RPS16", "SUPT5H", "IFNL1", "SYCN", "IFNL2", "TIMM50", "IFNL3",
                                      "DLL3", "LEUTX", "FCGBP", "LGALS14", "LGALS16", "FBL", "DYRK1B", "CLC")
  related_genes[["EGFR"]] <- c("EGFR", "VOPP1", "LANCL2", "CCT6A", "GBAS", "MRPS17", "CHCHD2", "VSTM2A", "NUPR1L",
                               "PHKG1", "SEPT14", "SUMF2", "ZNF713", "PSPH", "SEC61G")
  related_genes[["KRAS"]] <- c("LYRM5", "KRAS", "CASC1", "LRMP", "C12orf77", "BCAT1", "MDM2", "CPSF6", "CPM", "IFLTD1")
  related_genes[["MYC"]] <- c("MYC", "POU5F1B")

  prop_hlg_p_plot <- prop_hlg_p %>% filter(Hugo_Symbol %in% Reduce(union, related_genes))

  for (name_gene_group in names(related_genes)){
    genes_group <- related_genes[[name_gene_group]]

    if (name_gene_group %in% prop_hlg_p_plot$Hugo_Symbol){
      genes_group <- setdiff(genes_group, name_gene_group)
      prop_hlg_p_plot <- prop_hlg_p_plot %>% filter(!Hugo_Symbol %in% genes_group)
      df_plot_hlg <- df_plot_hlg %>% filter(!Hugo_Symbol %in% genes_group)
    } else {
      mean_prop_group <- prop_hlg_p_plot %>% filter(Hugo_Symbol %in% genes_group) %>% {mean(.$Prop)}
      prop_hlg_p_plot <- bind_rows(prop_hlg_p_plot %>% filter(!Hugo_Symbol %in% genes_group),
                                   tibble(Hugo_Symbol=name_gene_group, Prop=mean_prop_group))
      df_plot_hlg <- df_plot_hlg %>%
        mutate(Hugo_Symbol=ifelse(Hugo_Symbol %in% genes_group, name_gene_group, Hugo_Symbol)) %>% distinct()
    }
  }

  # apply the selection
  df_plot_hlg <- df_plot_hlg %>% filter(Hugo_Symbol %in% prop_hlg_p_plot$Hugo_Symbol)
  df_plot_hlg <- df_plot_hlg %>% group_by(Cohort, Hugo_Symbol) %>% summarize(Count=n())
  df_plot_hlg <- left_join(df_plot_hlg, df_n_pats, by="Cohort")
  df_plot_hlg <- df_plot_hlg %>% mutate(Prop=Count/Size)
  
  # choose order
  genes_order <- df_plot_hlg %>% filter(Cohort=="prism") %>% arrange(desc(Prop)) %>% pull(var="Hugo_Symbol")
  df_plot_hlg$Cohort <- factor(df_plot_hlg$Cohort, levels=args$cohorts)
  df_plot_hlg$Hugo_Symbol <- factor(df_plot_hlg$Hugo_Symbol, levels=genes_order)

  # deleted genes ======================================================================================================
  df_plot_hd <- df_plot %>% filter(Copy_Number==-2)

  # select genes shown in the plot using prism
  n_pat_p <- n_pats$prism
  df_cna_p_hd <- dfs_cna$prism %>% filter(Copy_Number==-2)
  prop_hd_p <- df_cna_p_hd %>% group_by(Hugo_Symbol, Chromosome) %>% summarize(Prop=n()/n_pat_p, .groups="drop") %>%
    arrange(desc(Prop))

  # select only genes mutated in at least 2.5% of prism
  # prop_hd_p_plot <- prop_hd_p %>% filter(Prop >= 0.025)

  # simplify the selection taking into account related genes
  related_genes <- list()
  related_genes[["CDKN2A"]] <- c("CDKN2B", "CDKN2A", "C9orf53", "MTAP", "DMRTA1", "IFNE", "IFNA1", "ELAVL2", "IFNA8",
                                 "IFNA5", "IFNA6", "KLHL9", "IFNA13", "IFNA2", "IFNA7", "IFNA14", "IFNA4", "IFNA16",
                                 "IFNA21", "IFNA10", "IFNA17", "IFNW1", "IZUMO3", "IFNB1", "PTPLAD2", "FOCAD", "TUSC1",
                                 "MLLT3")
  related_genes[["PTEN"]] <- c("PTEN", "RNLS", "KLLN", "ATAD1", "LIPJ", "PAPSS2")
  related_genes[["FAM160A / LGALS9C"]] <- c("FAM106A", "LGALS9C")
  related_genes[["KIR genes"]] <- c("KIR3DL1", "KIR2DL4", "KIR2DL3", "KIR3DL3", "KIR2DL1")
  related_genes[["RHD / RSRP1"]] <- c("RHD", "RSRP1", "C1orf63")

  prop_hd_p_plot <- prop_hd_p %>% filter(Hugo_Symbol %in% Reduce(union, related_genes))

  for (name_gene_group in names(related_genes)){
    genes_group <- related_genes[[name_gene_group]]

    if (name_gene_group %in% prop_hd_p_plot$Hugo_Symbol){
      genes_group <- setdiff(genes_group, name_gene_group)
      prop_hd_p_plot <- prop_hd_p_plot %>% filter(!Hugo_Symbol %in% genes_group)
      df_plot_hd <- df_plot_hd %>% filter(!Hugo_Symbol %in% genes_group)
    } else {
      mean_prop_group <- prop_hd_p_plot %>% filter(Hugo_Symbol %in% genes_group) %>% {mean(.$Prop)}
      prop_hd_p_plot <- bind_rows(prop_hd_p_plot %>% filter(!Hugo_Symbol %in% genes_group),
                                   tibble(Hugo_Symbol=name_gene_group, Prop=mean_prop_group))
      df_plot_hd <- df_plot_hd %>%
        mutate(Hugo_Symbol=ifelse(Hugo_Symbol %in% genes_group, name_gene_group, Hugo_Symbol)) %>% distinct()
    }
  }

  # apply the selection
  df_plot_hd <- df_plot_hd %>% filter(Hugo_Symbol %in% prop_hd_p_plot$Hugo_Symbol)
  df_plot_hd <- df_plot_hd %>% group_by(Cohort, Hugo_Symbol) %>% summarize(Count=n())
  df_plot_hd <- left_join(df_plot_hd, df_n_pats, by="Cohort")
  df_plot_hd <- df_plot_hd %>% mutate(Prop=Count/Size)

  # choose order
  genes_order <- df_plot_hd %>% filter(Cohort=="prism") %>% arrange(desc(Prop)) %>% pull(var="Hugo_Symbol")
  df_plot_hd$Cohort <- factor(df_plot_hd$Cohort, levels=args$cohorts)
  df_plot_hd$Hugo_Symbol <- factor(df_plot_hd$Hugo_Symbol, levels=genes_order)

  # plots ==============================================================================================================
  
  # plot parameters
  cohorts2colors <- load_colors("Global")[args$cohorts]
  cohorts2colors_pale <- lapply(cohorts2colors, function(col) paste0(col, "40"))

  # left  barplot 
  p1 <- ggbarplot(df_plot_hlg %>% filter(Cohort=="prism"), x="Hugo_Symbol", y="Prop", color="Cohort", fill="Cohort",
                  position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(size=9, family="Helvetica", face="plain", color="black", angle=90, 
                                     hjust=1, vjust=0.5),
          axis.text.y = element_text(size=9, family="Helvetica", face="plain", color="black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=0.1, color="black"),
          axis.line.x = element_line(size=0.1, color="black"),
          axis.line.y = element_line(size=0.1, color="black"),
          title = element_text(size=11, family="Helvetica", face="plain", color="black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=9, family="Helvetica", face="plain", color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size=9, family="Helvetica", face="plain", color="black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "none",
          legend.key.size = unit(0.4, "cm"),
          plot.margin=unit(c(1,0.1,0.1,0.1), "cm")) +
    labs(title="Focal amplifications", y="Fraction of samples", x=NULL) +
    scale_color_manual(breaks = names(cohorts2colors),
                       values = unlist(cohorts2colors)) +
    scale_fill_manual(breaks = names(cohorts2colors_pale),
                       values = unlist(cohorts2colors_pale)) +
    scale_x_discrete(position="bottom") +
    scale_y_continuous(limits=c(0, 0.16), breaks=c(0, 0.04, 0.08, 0.12, 0.16), labels=c("0%", "4%", "8%", "12%", "16%"))

  # right  barplot 
  p2 <- ggbarplot(df_plot_hd %>% filter(Cohort=="prism"), x="Hugo_Symbol", y="Prop", color="Cohort", fill="Cohort",
                  position = position_dodge(0.9)) +
    theme(axis.text.x = element_text(size=9, family="Helvetica", face="plain", color="black", angle=90, 
                                     hjust=1, vjust=0.5),
          axis.text.y = element_text(size=9, family="Helvetica", face="plain", color="black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=0.1, color="black"),
          axis.line.x = element_line(size=0.1, color="black"),
          axis.line.y = element_line(size=0.1, color="black"),
          title = element_text(size=11, family="Helvetica", face="plain", color="black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=9, family="Helvetica", face="plain", color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size=9, family="Helvetica", face="plain", color="black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "none",
          legend.key.size = unit(0.4, "cm"),
          plot.margin=unit(c(1,0.1,0.1,0.1), "cm")) +
    labs(title="Focal deletions", y="Fraction of samples", x=NULL) +
    scale_color_manual(breaks = names(cohorts2colors),
                       values = unlist(cohorts2colors)) +
    scale_fill_manual(breaks = names(cohorts2colors_pale),
                       values = unlist(cohorts2colors_pale)) +
    scale_x_discrete(position="bottom") +
    scale_y_continuous(limits=c(0, 0.16), breaks=c(0, 0.04, 0.08, 0.12, 0.16), labels=c("0%", "4%", "8%", "12%", "16%"))

  p <- ggarrange(p1, p2, ncol=2, nrow=1, widths=c(0.5, 0.5), align="h")
  ggsave(args$output_plot, p, width=14, height=8, units="cm")
  cat(paste("-plot saved at", args$output_plot, "\n"))
  ggsave(args$output_plot_paper, p, width=14, height=8, units="cm")
  cat(paste("-plot saved at", args$output_plot_paper, "\n"))
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw barplots of CNA frequency.')
  parser$add_argument("--cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/somatic_cnas/selection/selection_samples_prism.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_met500.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--cln_tables", nargs="+", help="Paths to tables of clinical attributes.",
                      default=c("../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv",
                                "../../../data/met500/clinical/curated/cln_met500_in_design_curated.tsv",
                                "../../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv"))
  parser$add_argument("--cna_tables", nargs="+", help="Paths to tables of CNA calls.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls.tsv.gz"))
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/somatic_cnas/selection/selection_tumor_types.tsv")
  parser$add_argument("--output_plot", type="character", help="Path to output plots.",
                      default="../../../results/somatic_cnas/other_plots/barplots_main_genes.pdf")
  parser$add_argument("--output_plot_paper", type="character",
                      default="../../../results/figures_paper/FS11d-e.pdf")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
