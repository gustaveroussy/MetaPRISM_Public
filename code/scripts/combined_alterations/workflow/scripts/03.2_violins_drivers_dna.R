# @created: 30 Mar 22
# @modified: 20 Aug 22
# @authors: Yoann Pradat
#
# Draw violin plots showing age and drug count distribution per tumor type. 

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggtext))

# functions ============================================================================================================

main <- function(args){
  # load data
  dfs_cln <- lapply(args$cohorts, load_cln)
  dfs_cln <- setNames(dfs_cln, args$cohorts)
  dfs_cln <- lapply(args$cohorts, function(x) dfs_cln[[x]] %>% mutate(Cohort=toupper(x)))
  dfs_cln <- setNames(dfs_cln, args$cohorts)
  df_cln <- bind_rows(dfs_cln) %>% rename(Tumor_Type=Project_TCGA_More)
  df_cln$Cohort <- factor(df_cln$Cohort, level=toupper(args$cohorts))

  dfs_alt <- lapply(args$alt_tables, load_table)
  dfs_alt <- setNames(dfs_alt, args$cohorts)
  dfs_alt <- lapply(args$cohorts, function(x) dfs_alt[[x]] %>% mutate(Cohort=toupper(x)))
  dfs_alt <- setNames(dfs_alt, args$cohorts)
  df_alt <- bind_rows(dfs_alt)

  dfs_sam <- lapply(args$samples, load_table)
  dfs_sam <- setNames(dfs_sam, args$cohorts)
  dfs_sam <- lapply(args$cohorts, function(x) dfs_sam[[x]] %>% mutate(Cohort=toupper(x)))
  dfs_sam <- setNames(dfs_sam, args$cohorts)
  df_sam <- bind_rows(dfs_sam)

  # select same samples as for heatmap dna
  df_sam <- df_sam %>% filter(Use_heatmap_dna==1)
  df_cln <- df_cln %>% unite("Sample_Id", Sample_Id_DNA_T, Sample_Id_DNA_N, sep="_vs_", remove=F)
  df_cln <- df_cln %>% filter(Sample_Id %in% df_sam$Sample_Id)
  df_alt <- df_alt %>% filter(Sample_Id %in% df_sam$Sample_Id)

  # select only mutations, small indels, and CNAs
  cat_mut_ind_cna  <- c("Mut", "Ins", "Del", "Amplification", "Deletion")
  df_alt <- df_alt %>% filter(Alteration_Category %in% cat_mut_ind_cna)

  # add number of drivers to clinical table
  df_cnt_alt <- df_alt %>% group_by(Subject_Id) %>% summarize(Count_Alt=n())
  df_cln <- left_join(df_cln, df_cnt_alt, by="Subject_Id")
  df_cln[is.na(df_cln["Count_Alt"]), "Count_Alt"] <- 0

  # add additional tumor type correspond to full cohort
  # df_cln <- bind_rows(df_cln, df_cln %>% mutate(Tumor_Type:="All"))

  # set tumor type as factor which chosen order
  df_count <- load_table(args$counts)
  df_count <- df_count %>% arrange(desc(PRISM_heatmap_dna))
  tt_order <- intersect(df_count$Tumor_Type, unique(df_cln$Tumor_Type))
  # tt_order <- c(tt_order, "All")
  df_cln$Tumor_Type <- factor(df_cln$Tumor_Type, levels=tt_order)
  df_cln$Cohort <- factor(df_cln$Cohort, levels=toupper(args$cohorts))

  # prepare boxplots
  cohorts2colors <- load_colors("Global")[args$cohorts]
  cohorts2colors_clear <- load_colors("Global")[paste0(args$cohorts, "_clear")]
  names(cohorts2colors) <- toupper(names(cohorts2colors))
  names(cohorts2colors_clear) <- toupper(names(cohorts2colors))
  cohorts2colors_pale <- lapply(cohorts2colors, function(col) paste0(col, "40"))
  cohorts2black <- lapply(cohorts2colors, function(col) "black")

  # boxplot
  p <- ggboxplot(df_cln, x="Tumor_Type", y="Count_Alt", color="Cohort", fill="Cohort",
                  outlier.shape=NA, lwd=0.3, coef=0) +
    theme(axis.text.x = element_text(size=10, family="Helvetica", face="plain", color="black", angle=90, 
                                     hjust=1, vjust=0.5),
          axis.text.y = element_text(size=10, family="Helvetica", face="plain", color="black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(size=0.25, color="black"),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=10, family="Helvetica", face="plain", color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = c(0.875, 1.1),
          legend.key.size = unit(0.4, "cm"),
          plot.margin=unit(c(1,0.1,0.1,0.1), "cm")) +
    labs(title=NULL, y="# WES-derived driver events", x=NULL) +
    scale_fill_manual(breaks = names(cohorts2colors_clear),
                       values = unlist(cohorts2colors_clear)) +
    scale_color_manual(breaks = names(cohorts2colors),
                       values = unlist(cohorts2colors)) +
    scale_x_discrete(position="bottom") +
    scale_y_continuous(limits=c(0, 15), breaks=c(1,5,9,13), labels=c(1,5,9,13))

  # wilcoxon adjusted tests
  stat.test <- df_cln %>% group_by(Tumor_Type) %>%
    wilcox_test(Count_Alt ~ Cohort, ref.group="TCGA") %>%
    adjust_pvalue(method="BH") %>% add_significance("p.adj")
  stat.test <- stat.test %>% add_xy_position(x="Tumor_Type", dodge=0.8, fun="median_iqr")
  p <- p + stat_pvalue_manual(stat.test, label="p.adj.signif", tip.length=0, remove.bracket=F, hide.ns=T, size=3)

  ggsave(args$output, p, width=3.75, height=2.5, units="in", dpi=600)
  cat(paste("-INFO: file saved at", args$output, "\n"))

  if (!is.null(args$output_paper)){
    ggsave(args$output_paper, p, width=3.75, height=2.5, units="in", dpi=600)
    cat(paste("-INFO: file saved at", args$output_paper, "\n"))
  }
}

# run ==================================================================================================================
 
if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw violin plots showing age and drug count distribution per tumor type.')
  parser$add_argument("--cohorts", type="character", nargs="+", help="Name of the cohort.", 
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--counts_tt", type="character", help="Path to table of counts per tumor type",
                      default="../../../results/combined_alterations/selection/selection_tumor_types.tsv")
  parser$add_argument("--alt_tables", type="character", nargs="+", help="Path to table of alt_table",
                      default=c("../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv",
                                "../../../results/combined_alterations/alterations/aggregated_alterations_met500.tsv",
                                "../../../results/combined_alterations/alterations/aggregated_alterations_tcga.tsv"))
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/combined_alterations/selection/selection_samples_met500.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_prism.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_tcga.tsv"))
  parser$add_argument('--output', type="character", help='Path to output plot.',
                  default="../../../results/combined_alterations/other_plots/violins_drivers_dna.pdf")
  parser$add_argument('--output_paper', type="character", help='Path to output plot for paper.', default=NULL)
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
