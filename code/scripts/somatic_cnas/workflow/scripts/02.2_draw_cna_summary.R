# @created: 17 Aug 22
# @modified: 31 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(rstatix))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# functions ============================================================================================================

select_samples_and_add_cln_data <- function(df_dat, df_cln, df_sam, name, cohort,
                                            cols_cln=c("Subject_Id", "Project_TCGA_More")){
  cat(paste0("-INFO: processing ", cohort, "...\n"))

  df_dat <- df_dat %>% unite("DNA_P", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_", remove=F)

  # select samples in tables selection
  df_sam <- df_sam[df_sam[[paste0("Use_", name)]]==1,]
  mask_sam <- df_dat$DNA_P %in% df_sam$Sample_Id
  cat(paste0("-INFO: selected ", sum(mask_sam), "/", length(mask_sam) ," lines from tumor/normal pairs in sam.\n"))
  df_dat <- df_dat[mask_sam,]

  # add clinical data
  if (!is.null(df_cln) & !is.null(cols_cln)){
    df_cln <- df_cln %>% unite("DNA_P", Sample_Id_DNA_T, Sample_Id_DNA_N, sep="_vs_", remove=F)
    df_dat <- left_join(df_dat, df_cln[,c("DNA_P", cols_cln)], by="DNA_P")
  }

  df_dat
}


get_pvalues_bars <- function(cohorts, df_cna, col_x, levels_x, bracket_width=0.267){
  # wilcoxon adjusted tests
  dfs_margins <- lapply(cohorts, function(x) df_cna %>% filter(Cohort==toupper(x)) %>% group_by_at(col_x) %>% 
                        summarize(Count=n()) %>% column_to_rownames(var=col_x))
  dfs_margins <- setNames(dfs_margins, cohorts)

  dfs_counts <- lapply(cohorts, function(x) df_cna %>% filter(Cohort==toupper(x), WGD>0) %>% 
                       group_by_at(col_x) %>% summarize(Count=n()) %>% mutate(Type="WGD") %>% 
                       spread({{col_x}}, Count))
  dfs_counts <- setNames(dfs_counts, cohorts)
  dfs <- lapply(cohorts, function(x) list(count=dfs_counts[[x]], count_col=dfs_margins[[x]]))
  dfs <- setNames(dfs, cohorts)

  # prism vs tcga
  pvals_p_vs_t <- compute_pvals_table(dfs, cohort_a="prism", cohort_b="tcga", row_col="Type", row_names=c("WGD"),
                                      col_names=levels_x, test="fisher")
  df_pvals_p_vs_t <- pvals_p_vs_t %>% gather({{col_x}}, p) %>% 
    mutate(`.y.`="Fraction", group1="TCGA", group2="PRISM") %>% as_tibble()
  df_pvals_p_vs_t <- left_join(df_pvals_p_vs_t,
                               dfs_margins$tcga %>% rownames_to_column(var=col_x) %>% rename(n1=Count), by=col_x)
  df_pvals_p_vs_t <- left_join(df_pvals_p_vs_t,
                               dfs_margins$prism %>% rownames_to_column(var=col_x) %>% rename(n2=Count), by=col_x)

  # met500 vs tcga
  pvals_m_vs_t <- compute_pvals_table(dfs, cohort_a="met500", cohort_b="tcga", row_col="Type", row_names=c("WGD"),
                                      col_names=levels_x, test="fisher")
  df_pvals_m_vs_t <- pvals_m_vs_t %>% gather({{col_x}}, p) %>% 
    mutate(`.y.`="Fraction", group1="TCGA", group2="MET500") %>% as_tibble()
  df_pvals_m_vs_t <- left_join(df_pvals_m_vs_t,
                               dfs_margins$tcga %>% rownames_to_column(var=col_x) %>% rename(n1=Count), by=col_x)
  df_pvals_m_vs_t <- left_join(df_pvals_m_vs_t,
                               dfs_margins$met500 %>% rownames_to_column(var=col_x) %>% rename(n2=Count), by=col_x)

  # adjust pvalues
  df_pvals <- bind_rows(df_pvals_p_vs_t, df_pvals_m_vs_t) %>% mutate(p=abs(p)) %>% 
    adjust_pvalue(method="BH") %>% add_significance("p.adj")
  
  # add y.position
  dfs_cna_wgd <- lapply(cohorts, function(x) df_cna %>% filter(Cohort==toupper(x)) %>% group_by_at(col_x) %>%
                        summarize(!!paste0("Fraction_", x):=mean(WGD)))
  dfs_cna_wgd <- setNames(dfs_cna_wgd, cohorts)
  df_ypos_p_vs_t <- full_join(dfs_cna_wgd$tcga %>% mutate(group1="TCGA"),
                              dfs_cna_wgd$prism %>% mutate(group2="PRISM"), by=col_x) %>%
                    mutate(y.position=pmax(Fraction_tcga, Fraction_prism)) %>% 
                    select(-Fraction_tcga, -Fraction_prism)
  df_ypos_p_vs_t$groups <- lapply(1:nrow(df_ypos_p_vs_t), function(i) c("TCGA", "PRISM"))
  df_ypos_m_vs_t <- full_join(dfs_cna_wgd$tcga %>% mutate(group1="TCGA"),
                              dfs_cna_wgd$met500 %>% mutate(group2="MET500"), by=col_x) %>%
                    mutate(y.position=pmax(Fraction_tcga, Fraction_met500)) %>%
                    select(-Fraction_tcga, -Fraction_met500)
  df_ypos_m_vs_t$groups <- lapply(1:nrow(df_ypos_m_vs_t), function(i) c("TCGA", "MET500"))
  df_ypos <- bind_rows(df_ypos_p_vs_t, df_ypos_m_vs_t)
  df_pvals <- left_join(df_pvals, df_ypos, by=c(col_x, "group1", "group2"))

  # add x, xmin, xmax
  df_pvals[[col_x]] <- factor(df_pvals[[col_x]], levels=rev(levels_x))
  df_pvals$x <- as.numeric(df_pvals[[col_x]])
  df_pvals$xmin <- df_pvals$x - bracket_width
  df_pvals$xmax <- df_pvals$x + ifelse(df_pvals$group2=="MET500", 0, 1) * bracket_width 

  df_pvals
}


select_samples_from_wgd <- function(dfs_cna, dfs_wgd, wgd_value=0){
  col_pai <- "DNA_P"
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  cohorts <- names(dfs_cna)
  
  for (cohort in cohorts){
    df_cna <- dfs_cna[[cohort]] %>% unite(!!col_pai, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F)
    df_wgd <- dfs_wgd[[cohort]] %>% unite(!!col_pai, all_of(c(col_tsb, col_nsb)), sep="_vs_", remove=F) %>%
      filter(WGD %in% wgd_value)
    cna_pai <- unique(df_cna[[col_pai]])
    wgd_pai <- unique(df_wgd[[col_pai]])
    cna_non_wgd_pai <- setdiff(cna_pai, wgd_pai)
    mask_in <- df_cna[[col_pai]] %in% cna_non_wgd_pai
    cat("-INFO: selected", sum(mask_in), "lines from", paste0(length(cna_non_wgd_pai), "/", length(cna_pai)),
        paste0("with WGD=", wgd_value), "from cohort", cohort, "\n")
    dfs_cna[[cohort]] <- df_cna[mask_in,]
  }

  dfs_cna
}


draw_density_plot <- function(cohorts, dfs_cna){
  # aggregate
  df_cna <- bind_rows(dfs_cna) %>%
    select(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, Ploidy, WGD, Cohort, Tumor_Type)

  # add cohort with ploidy status for plot
  df_cna <- df_cna %>% mutate(Subcohort=ifelse(WGD==0, paste(Cohort, "diploid"), paste(Cohort, "polyploid")))

  # order cohorts
  levels_cohort <- rev(toupper(cohorts))
  levels_subcohort <- Reduce(union, lapply(levels_cohort, function(x) paste(x, c("diploid", "polyploid"))))
  df_cna$Cohort <- factor(df_cna$Cohort, levels=levels_cohort)
  df_cna$Subcohort <- factor(df_cna$Subcohort, levels=levels_subcohort)

  # colors
  colors_global <- load_colors(sheet="Global")

  # draw density plot
  palette_fill <- list(`TCGA diploid`=colors_global$tcga_clear,
                       `TCGA polyploid`=colors_global$tcga,
                       `MET500 diploid`=colors_global$met500_clear,
                       `MET500 polyploid`=colors_global$met500,
                       `PRISM diploid`=colors_global$prism_clear,
                       `PRISM polyploid`=colors_global$prism)

  palette_color <- list(`TCGA diploid`=colors_global$tcga,
                       `TCGA polyploid`=colors_global$tcga,
                       `MET500 diploid`=colors_global$met500,
                       `MET500 polyploid`=colors_global$met500,
                       `PRISM diploid`=colors_global$prism,
                       `PRISM polyploid`=colors_global$prism)

  plot_den <- ggdensity(df_cna, x="Ploidy", color="Subcohort", fill="Subcohort", alpha=1) +
    fill_palette(palette=unlist(palette_fill)) +
    color_palette(palette=unlist(palette_color)) +
    coord_cartesian(xlim=c(0.75,6.25), ylim=c(0,NA)) + 
    scale_x_continuous(breaks=seq(1,6,1)) +
    scale_y_continuous(expand=c(0.01,0)) +
    labs(title=NULL, y=NULL, x="Genome ploidy") +
    theme(legend.title=element_blank(),
          axis.ticks.x = element_line(size=0.25, color="black"),
          axis.ticks.y=element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_text(family="Helvetica", colour="black", size=8),
          axis.text.x=element_text(family="Helvetica", colour="black", size=8),
          legend.key.size = unit(0.5, "lines"),
          legend.text=element_text(family="Helvetica", colour="black", size=7),
          legend.direction="vertical",
          legend.position=c(0, 0.825),
          panel.background=element_rect(fill = "transparent",colour=NA),
          plot.background=element_rect(fill = "transparent",colour=NA),
          plot.margin=unit(c(0.1,1.3,0.1,1.3), "cm")) +
        guides(shape = guide_legend(override.aes = list(size = 0.5)),
               color = guide_legend(override.aes = list(size = 0.5)))

  # data right boxplot
  df_cna_wgd <- df_cna %>% group_by(Cohort, Tumor_Type) %>% summarize(Fraction=mean(WGD), .groups="drop")

  # order tumor types
  levels_tumor_type <- df_cna %>% filter(Cohort=="PRISM") %>% group_by(Tumor_Type) %>% mutate(N=n()) %>%
    distinct(Tumor_Type, N) %>% arrange(desc(N)) %>% pull(Tumor_Type)
  df_cna_wgd$Tumor_Type <- factor(df_cna_wgd$Tumor_Type, levels=rev(levels_tumor_type))

  # draw right boxplot
  palette_fill <- list(TCGA=colors_global$tcga,
                       MET500=colors_global$met500,
                       PRISM=colors_global$prism) 

  palette_color <- list(TCGA=colors_global$tcga,
                       MET500=colors_global$met500,
                       PRISM=colors_global$prism)

  # comute pvalues
  df_pvals <- get_pvalues_bars(cohorts=args$cohorts, df_cna=df_cna, col_x="Tumor_Type", levels_x=levels_tumor_type)

  plot_bar <- ggbarplot(df_cna_wgd, x="Tumor_Type", y="Fraction", color="Cohort", fill="Cohort", width=0.2,
                  position=position_dodge(0.9)) +
    theme(axis.text.x = element_text(size=6, family="Helvetica", face="plain", color="black"),
          axis.text.y = element_text(size=6, family="Helvetica", face="plain", color="black", hjust=1),
          axis.ticks.x = element_line(size=0.25, color="black"),
          axis.ticks.y = element_blank(), 
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_blank(),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=6, family="Helvetica", face="plain", color="black"),
          legend.position = "none") +
    labs(title=NULL, y="WGD frequency", x=NULL) +
    fill_palette(palette=unlist(palette_fill)) +
    color_palette(palette=unlist(palette_color)) +
    scale_y_continuous(expand=c(0.002,0), breaks=c(0, 0.2, 0.4, 0.6, 0.8), labels=c("0", "20", "40", "60", "80"))


  plot_bar <- plot_bar + stat_pvalue_manual(df_pvals, label="p.adj.signif", tip.length=0, remove.bracket=F,
                                            hide.ns=T, size=2, coord.flip=T, bracket.nudge.y=0.15) + coord_flip()

  # combine and save
  plot <- plot_den + annotation_custom(grob=ggplotGrob(plot_bar), xmin=3.75, xmax=7.75, ymin=0.5, ymax=4.5)
}


draw_inverted_boxplots <- function(df_dat, col_y1, col_y2, ylab1, ylab2, col_x, col_z, colors_z, colors_z_clear,
                                   scale_y1_breaks, scale_y1_labels,
                                   scale_y2_breaks, scale_y2_labels,
                                   z_ref=NULL, top_margin=0, legend_position="none"){

  quantiles <- df_dat %>% group_by_at(col_x) %>%
    summarize(Max_1=quantile(.data[[col_y1]], 0.75), Max_2=quantile(.data[[col_y2]], 0.75))
  quantile_max <- max(max(quantiles$Max_1, quantiles$Max_2))
  quantile_max <- round(1.25*quantile_max, max(0, -floor(log10(quantile_max))))

  # left boxplot
  p1 <- ggboxplot(df_dat, x=col_x, y=col_y1, color=col_z, fill=col_z,
                  outlier.shape=NA, lwd=0.3, coef=0) +
    theme(axis.text.x = element_text(size=8, family="Helvetica", face="plain", color="black"),
          axis.text.y = element_blank(),
          axis.ticks.x = element_line(size=0.25, color="black"),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.title = element_blank(),
          legend.text = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = legend_position, 
          legend.key.size = unit(0.4, "cm"),
          plot.margin=unit(c(top_margin,0.1,0,0.4), "cm")) +
    labs(title=NULL, y=ylab1, x=NULL) +
    scale_color_manual(breaks = names(colors_z),
                       values = unlist(colors_z)) +
    scale_fill_manual(breaks = names(colors_z_clear),
                       values = unlist(colors_z_clear)) +
    scale_x_discrete(position="top") +
    scale_y_reverse(limits=c(quantile_max, NA), expand=c(0,0), breaks=scale_y1_breaks, labels=scale_y1_labels)


  # wilcoxon adjusted tests
  stat.test.1 <- df_dat %>% group_by_at(col_x) %>%
    wilcox_test(eval(parse(text=paste0("`", col_y1, "` ~ ", "`", col_z, "`"))), ref=z_ref) %>%
    adjust_pvalue(method="BH") %>% add_significance("p.adj")
  stat.test.1 <- stat.test.1 %>% add_xy_position(x=col_x, dodge=0.8, fun="median_iqr") %>%
    mutate(y.position=-y.position)

  p1 <- p1 + stat_pvalue_manual(stat.test.1, label="p.adj.signif", tip.length=0, remove.bracket=F, hide.ns=T, 
                                size=2, coord.flip=T, bracket.nudge.y=-0.1, vjust=1.5) + coord_flip()


  # right boxplot
  p2 <- ggboxplot(df_dat, x=col_x, y=col_y2, color=col_z, fill=col_z,
                  outlier.shape=NA, lwd=0.3, coef=0) +
    theme(axis.text.x = element_text(size=8, family="Helvetica", face="plain", color="black"),
          axis.text.y = element_text(size=8, family="Helvetica", face="plain", color="black", hjust=0.5),
          axis.ticks.x = element_line(size=0.25, color="black"),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size=0.25, color="black"),
          axis.line.y = element_line(size=0.25, color="black"),
          panel.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.title = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.text = element_text(size=8, family="Helvetica", face="plain", color="black"),
          legend.background = element_blank(),
          legend.key = element_blank(),
          legend.position = "none",
          legend.key.size = unit(0.4, "cm"),
          plot.margin=unit(c(top_margin,0.5,0,0), "cm")) +
    labs(title=NULL, y=ylab2, x=NULL) +
    scale_color_manual(breaks = names(colors_z),
                       values = unlist(colors_z)) +
    scale_fill_manual(breaks = names(colors_z_clear),
                       values = unlist(colors_z_clear)) +
    scale_y_continuous(limits=c(NA, quantile_max), expand=c(0,0), breaks=scale_y2_breaks, labels=scale_y2_labels)

  # wilcoxon adjusted tests
  stat.test.2 <- df_dat %>% group_by_at(col_x) %>%
    wilcox_test(eval(parse(text=paste0("`", col_y2, "` ~ ", "`", col_z, "`"))), ref=z_ref) %>%
    adjust_pvalue(method="BH") %>% add_significance("p.adj")
  stat.test.2 <- stat.test.2 %>% add_xy_position(x=col_x, dodge=0.8, step.increase=0.025, fun="median_iqr")

  p2 <- p2 + stat_pvalue_manual(stat.test.2, label="p.adj.signif", tip.length=0, remove.bracket=F, hide.ns=T, 
                                size=2, coord.flip=T, bracket.nudge.y=0.1) + coord_flip()

  ggarrange(p1, p2, ncol=2, nrow=1, widths=c(0.45, 0.55))
}


draw_inverted_boxplots_non_wgd <- function(cohorts, dfs_cna){
  # aggregate
  df_cna <- bind_rows(dfs_cna)

  # select tumors without WGD
  df_cna <- df_cna %>% filter(WGD==0)

  # order cohorts
  levels_cohort <- rev(toupper(cohorts))
  df_cna$Cohort <- factor(df_cna$Cohort, levels=levels_cohort)

  # order tumor types
  levels_tumor_type <- df_cna %>% filter(Cohort=="PRISM") %>% group_by(Tumor_Type) %>% mutate(N=n()) %>%
    distinct(Tumor_Type, N) %>% arrange(desc(N)) %>% pull(Tumor_Type)
  df_cna$Tumor_Type <- factor(df_cna$Tumor_Type, levels=rev(levels_tumor_type))

  # prepare boxplots
  cohorts2colors <- load_colors(sheet="Global")[cohorts]
  cohorts2colors_clear <- load_colors("Global")[paste0(args$cohorts, "_clear")]
  names(cohorts2colors) <- toupper(names(cohorts2colors))
  names(cohorts2colors_clear) <- toupper(names(cohorts2colors))

  # draw
  draw_inverted_boxplots(df_dat=df_cna, col_y1="GAIN:LL_ML_amplification", col_y2="LOSS:LOH_cnLOH",
                         ylab1="Genome fraction covered by\nLLG or MLG",
                         ylab2="Genome fraction covered by\nLOH or cn-LOH", col_x="Tumor_Type", col_z="Cohort",
                         colors_z=cohorts2colors, colors_z_clear=cohorts2colors_clear,
                         scale_y1_breaks=c(0.2, 0.4, 0.6), scale_y1_labels=c("20", "40", "60"),
                         scale_y2_breaks=c(0.2, 0.4, 0.6), scale_y2_labels=c("20", "40", "60"),
                         z_ref="TCGA", top_margin=0, legend_position=c(0.2, 0.5))
}


draw_inverted_boxplots_wgd <- function(df_cna, col_y1, col_y2, ylab1, ylab2, scale_y1_breaks, scale_y1_labels,
                                       scale_y2_breaks, scale_y2_labels){
  # add cohort with ploidy status for plot
  df_cna <- df_cna %>% mutate(Subcohort=ifelse(WGD==0, paste(Cohort, "diploid"), paste(Cohort, "polyploid")))

  # order subcohorts
  levels_subcohort <- Reduce(union, lapply(unique(df_cna$Cohort), function(x) paste(x, c("diploid", "polyploid"))))
  df_cna$Subcohort <- factor(df_cna$Subcohort, levels=rev(levels_subcohort))

  # order tumor types
  levels_tumor_type <- df_cna %>% group_by(Tumor_Type) %>% mutate(N=n()) %>%
    distinct(Tumor_Type, N) %>% arrange(desc(N)) %>% pull(Tumor_Type)
  df_cna$Tumor_Type <- factor(df_cna$Tumor_Type, levels=rev(levels_tumor_type))

  # prepare boxplots
  subcohorts2colors <- setNames(list("#ffb4a2", "#9d4edd"), levels_subcohort)
  subcohorts2colors_clear <- setNames(list("#ffcdb2", "#e0aaff"), levels_subcohort)

  # draw
  draw_inverted_boxplots(df_dat=df_cna, col_y1=col_y1, col_y2=col_y2,
                         ylab1=ylab1,ylab2=ylab2, col_x="Tumor_Type", col_z="Subcohort",
                         colors_z=subcohorts2colors, colors_z_clear=subcohorts2colors_clear,
                         scale_y1_breaks=scale_y1_breaks, scale_y1_labels=scale_y1_labels,
                         scale_y2_breaks=scale_y2_breaks, scale_y2_labels=scale_y2_labels,
                         z_ref=NULL, top_margin=1, legend_position=c(0.35, 1.25))
}


draw_inverted_boxplots_wgd_b <- function(df_cna){
  # add cohort with ploidy status for plot
  df_cna <- df_cna %>% mutate(Subcohort=ifelse(WGD==0, paste(Cohort, "diploid"), paste(Cohort, "polyploid")))

  # order subcohorts
  levels_subcohort <- Reduce(union, lapply(unique(df_cna$Cohort), function(x) paste(x, c("diploid", "polyploid"))))
  df_cna$Subcohort <- factor(df_cna$Subcohort, levels=rev(levels_subcohort))

  # order tumor types
  levels_tumor_type <- df_cna %>% group_by(Tumor_Type) %>% mutate(N=n()) %>%
    distinct(Tumor_Type, N) %>% arrange(desc(N)) %>% pull(Tumor_Type)
  df_cna$Tumor_Type <- factor(df_cna$Tumor_Type, levels=rev(levels_tumor_type))

  # prepare boxplots
  subcohorts2colors <- setNames(list("#ffb4a2", "#9d4edd"), levels_subcohort)
  subcohorts2colors_clear <- setNames(list("#ffcdb2", "#e0aaff"), levels_subcohort)

  # draw
  draw_inverted_boxplots(df_dat=df_cna, col_y1="Count_og", col_y2="Count_tsg",
                         ylab1="Number of oncogenes\naltered ",
                         ylab2="Number of tumorsuppressors\naltered", col_x="Tumor_Type", col_z="Subcohort",
                         colors_z=subcohorts2colors, colors_z_clear=subcohorts2colors_clear,
                         scale_y1_breaks=c(1,2,4), scale_y1_labels=c("1", "2", "4"),
                         scale_y2_breaks=c(1,2,4), scale_y2_labels=c("1", "2", "4"),
                         z_ref=NULL, top_margin=1, legend_position=c(0.35, 1.25))
}



add_number_genes_mut <- function(df_cna, df_mut, df_gen, type="og"){
  if (type=="og"){
    genes <- df_gen %>% filter(Inclusion_Level==1, Is_Oncogene=="Yes") %>% pull(var="Hugo_Symbol")
  } else if (type=="tsg"){
    genes <- df_gen %>% filter(Inclusion_Level==1, Is_Oncogene!="Yes", Is_Tumor_Suppressor_Gene=="Yes") %>%
      pull(var="Hugo_Symbol")
  }

  col_count <- paste0("Count_", type)
  df_cnt <- df_mut %>% filter(Hugo_Symbol %in% genes) %>%
    group_by(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>%
    summarize(!!col_count:=length(unique(Hugo_Symbol)))
  left_join(df_cna, df_cnt) %>% mutate(!!col_count:=ifelse(is.na(.data[[col_count]]), 0, .data[[col_count]]))
}


main <- function(args){
  # load data
  dfs_sam_cna <- setNames(lapply(args$sam_cna_tables, load_table), args$cohorts)
  dfs_sam_mut <- setNames(lapply(args$sam_mut_tables, load_table), args$cohorts)
  dfs_cna <- setNames(lapply(args$cna_tables, load_table), args$cohorts)
  dfs_cln <- setNames(lapply(args$cln_tables, load_table), args$cohorts)
  dfs_mut <- setNames(lapply(args$mut_tables, load_table), args$cohorts)
  dfs_arm <- setNames(lapply(args$arm_tables, load_table), args$cohorts)

  # select samples
  dfs_cna <- lapply(args$cohorts, function(x) select_samples_and_add_cln_data(df_dat=dfs_cna[[x]], df_cln=dfs_cln[[x]],
                                                                              df_sam=dfs_sam_cna[[x]], name="selections",
                                                                              cohort=x))
  dfs_cna <- setNames(dfs_cna, args$cohorts)

  dfs_arm <- lapply(args$cohorts, function(x) select_samples_and_add_cln_data(df_dat=dfs_arm[[x]], df_cln=dfs_cln[[x]],
                                                                              df_sam=dfs_sam_cna[[x]], name="selections",
                                                                              cohort=x))
  dfs_arm <- setNames(dfs_arm, args$cohorts)

  # polishing
  dfs_cna <- lapply(dfs_cna, function(df_cna) df_cna %>% rename(Tumor_Type=Project_TCGA_More))
  dfs_cna <- lapply(args$cohorts, function(x) dfs_cna[[x]] %>% mutate(Cohort=toupper(x)))
  dfs_cna <- setNames(dfs_cna, args$cohorts)

  # draw plot with ploidy densities + barplot ==========================================================================
  p <- draw_density_plot(cohorts=args$cohorts, dfs_cna=dfs_cna)
  ggsave(plot=p, width=3.75, height=2, units="in", dpi=300, filename=args$output_plots[1])
  cat(paste("-plot saved at", args$output_plots[1], "\n"))
  ggsave(plot=p, width=3.75, height=2, units="in", dpi=300, filename=args$output_plots_paper[1])
  cat(paste("-plot saved at", args$output_plots_paper[1], "\n"))

  # draw inverted boxplots =============================================================================================
  p <- draw_inverted_boxplots_non_wgd(cohorts=args$cohorts, dfs_cna=dfs_cna)
  ggsave(args$output_plots[2], p, width=3.5, height=2, units="in")
  cat(paste("-plot saved at", args$output_plots[2], "\n"))
  ggsave(args$output_plots_paper[2], p, width=3.5, height=2, units="in")
  cat(paste("-plot saved at", args$output_plots_paper[2], "\n"))
  
  # genome fractions llg / mlg vs log / cn-loh =========================================================================

  # draw inverted boxplots - prism only, wgd vs non-wgd, genome fractions llg / mlg vs log / cn-loh
  tt_keep <- c("BLCA", "BRCA", "LUAD", "PAAD", "PRAD")
  df_cna <- dfs_cna$prism %>% filter(Tumor_Type %in% tt_keep)
  # df_cna <- bind_rows(df_cna, df_cna %>% mutate(Tumor_Type="All"))

  p  <- draw_inverted_boxplots_wgd(df_cna=df_cna, col_y1="GAIN:LL_ML_amplification", col_y2="LOSS:LOH_cnLOH",
                                   ylab1="Genome fraction covered by\nLLG or MLG",
                                   ylab2="Genome fraction covered by\nLOH or cn-LOH", 
                                   scale_y1_breaks=c(0.2, 0.4, 0.6), scale_y1_labels=c("20%", "40%", "60%"),
                                   scale_y2_breaks=c(0.2, 0.4, 0.6), scale_y2_labels=c("20%", "40%", "60%"))

  ggsave(args$output_plots[3], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots[3], "\n"))
  ggsave(args$output_plots_paper[3], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots_paper[3], "\n"))

  # genome fractions HLG vs HD =========================================================================================
  p  <- draw_inverted_boxplots_wgd(df_cna=df_cna, col_y1="GAIN:HL_amplification", col_y2="LOSS:Deletion",
                                   ylab1="Genome fraction covered by\nHLG",
                                   ylab2="Genome fraction covered by\nHD", 
                                   scale_y1_breaks=c(0.2e-2, 0.4e-2, 0.6e-2),
                                   scale_y1_labels=c("0.20%", "0.40%", "0.60%"),
                                   scale_y2_breaks=c(0.2e-2, 0.4e-2, 0.6e-2),
                                   scale_y2_labels=c("0.20%", "0.40%", "0.60%"))
  ggsave(args$output_plots[4], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots[4], "\n"))
  ggsave(args$output_plots_paper[4], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots_paper[4], "\n"))

  # exact medians
  df_med <- df_cna %>% mutate(Subcohort=ifelse(WGD==0, paste(Cohort, "diploid"), paste(Cohort, "polyploid")))
    group_by(Tumor_Type, Subcohort) %>%
    summarize(Median_HLG=median(`GAIN:HL_amplification`), Median_HD=median(`LOSS:Deletion`)) %>% 
    mutate(Median_HLG=paste0(round(Median_HLG*100, 3), "%"), Median_HD=paste0(round(Median_HD*100,3), "%"))

  # number of chr arm gains/losses =====================================================================================
  df_arm_gain <- dfs_arm$prism %>% filter(copy_number %in% c(1,2)) %>%
    group_by(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% summarize(N_Chr_Arm_Gain=n())
  df_arm_loss <- dfs_arm$prism %>% filter(copy_number %in% c(-1)) %>%
    group_by(Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode) %>% summarize(N_Chr_Arm_Loss=n())
  df_cna <- left_join(df_cna, df_arm_gain)
  df_cna <- left_join(df_cna, df_arm_loss)
  df_cna <- df_cna %>% mutate(N_Chr_Arm_Gain=ifelse(is.na(N_Chr_Arm_Gain), 0, N_Chr_Arm_Gain)) %>%
    mutate(N_Chr_Arm_Loss=ifelse(is.na(N_Chr_Arm_Loss), 0, N_Chr_Arm_Loss))

  p  <- draw_inverted_boxplots_wgd(df_cna=df_cna, col_y1="N_Chr_Arm_Loss", col_y2="N_Chr_Arm_Gain",
                                   ylab1="Number of chromosome\narm losses",
                                   ylab2="Number of chromosome\narm gains",
                                   scale_y1_breaks=c(10,20,30),
                                   scale_y1_labels=c("10","20","30"),
                                   scale_y2_breaks=c(10,20,30),
                                   scale_y2_labels=c("10","20","30"))
  ggsave(args$output_plots[5], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots[5], "\n"))
  ggsave(args$output_plots_paper[5], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots_paper[5], "\n"))

  # select samples for combining mutations and cnas ====================================================================
  dfs_sam_com <- lapply(args$cohorts, function(x) inner_join(dfs_sam_cna[[x]], dfs_sam_mut[[x]]) %>%
                        mutate(Use_combination=ifelse(Use_selections==1 & Use_mutpanning==1, 1, 0)))
  dfs_sam_com <- setNames(dfs_sam_com, args$cohorts)

  # apply selection on cnas
  dfs_cna_com <- lapply(args$cohorts, function(x) select_samples_and_add_cln_data(df_dat=dfs_cna[[x]], df_cln=dfs_cln[[x]],
                                                                                  df_sam=dfs_sam_com[[x]],
                                                                                  name="combination", cols_cln=c(),
                                                                                  cohort=x))
  dfs_cna_com <- setNames(dfs_cna_com, args$cohorts)

  # apply selection on muts
  dfs_mut_com <- lapply(args$cohorts, function(x) select_samples_and_add_cln_data(df_dat=dfs_mut[[x]], df_cln=dfs_cln[[x]],
                                                                                  df_sam=dfs_sam_com[[x]],
                                                                                  name="combination",
                                                                                  cohort=x))
  dfs_mut_com <- setNames(dfs_mut_com, args$cohorts)

  # add number of mutated oncogenes and tumorsuppressors to df_cna_com tables
  df_gen <- load_table(args$gen_table)

  # oncogenes
  dfs_cna_com <- lapply(args$cohorts, function(x) add_number_genes_mut(df_cna=dfs_cna_com[[x]], df_mut=dfs_mut_com[[x]],
                                                                       df_gen=df_gen, type="og"))
  dfs_cna_com <- setNames(dfs_cna_com, args$cohorts)

  # tumorsuppressors
  dfs_cna_com <- lapply(args$cohorts, function(x) add_number_genes_mut(df_cna=dfs_cna_com[[x]], df_mut=dfs_mut_com[[x]],
                                                                       df_gen=df_gen, type="tsg"))
  dfs_cna_com <- setNames(dfs_cna_com, args$cohorts)

  # draw inverted boxplots - prism only, wgd vs non-wgd, og vs tsg
  tt_keep <- c("BLCA", "BRCA", "LUAD", "PAAD", "PRAD")
  df_cna <- dfs_cna_com$prism %>% filter(Tumor_Type %in% tt_keep)
  p  <- draw_inverted_boxplots_wgd(df_cna=df_cna, col_y1="Count_og", col_y2="Count_tsg",
                                   ylab1="Number of oncogenes\naltered ",
                                   ylab2="Number of tumorsuppressors\naltered",
                                   scale_y1_breaks=c(1,2,4), scale_y1_labels=c("1", "2", "4"),
                                   scale_y2_breaks=c(1,2,4), scale_y2_labels=c("1", "2", "4"))
  ggsave(args$output_plots[6], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots[6], "\n"))
  ggsave(args$output_plots_paper[6], p, width=3.5, height=1.75, units="in")
  cat(paste("-plot saved at", args$output_plots_paper[6], "\n"))
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw plot of ploidy densities.')
  parser$add_argument("--cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--sam_cna_tables", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/somatic_cnas/selection/selection_samples_prism.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_met500.tsv",
                                "../../../results/somatic_cnas/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--sam_mut_tables", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/somatic_mutations/selection/selection_samples_prism.tsv",
                                "../../../results/somatic_mutations/selection/selection_samples_met500.tsv",
                                "../../../results/somatic_mutations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--cln_tables", nargs="+", help="Paths to tables of clinical attributes.",
                      default=c("../../../data/prism/clinical/curated/cln_prism_in_design_curated.tsv",
                                "../../../data/met500/clinical/curated/cln_met500_in_design_curated.tsv",
                                "../../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv"))
  parser$add_argument("--cna_tables", nargs="+", help="Paths to tables of CNA calls.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls_summary_statistics.tsv.gz"))
  parser$add_argument("--arm_tables", nargs="+", help="Paths to tables of CNA chr arm calls.",
                      default=c("../../../data/prism/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz",
                                "../../../data/met500/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz",
                                "../../../data/tcga/wes/somatic_cna/somatic_calls_per_chr_arm.tsv.gz"))
  parser$add_argument("--mut_tables", nargs="+", help="Paths to tables of mutation calls.",
                      default=c("../../../data/prism/wes/somatic_maf/somatic_calls_union_ann.maf.gz",
                                "../../../data/met500/wes/somatic_maf/somatic_calls_union_ann.maf.gz",
                                "../../../data/tcga/wes/somatic_maf/somatic_calls_union_ann.maf.gz"))
  parser$add_argument("--gen_table", help="Path to table of genes.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/somatic_cnas/selection/selection_tumor_types.tsv")
  parser$add_argument("--output_plots", nargs="+", type="character", help="Path to output plots.",
                      default=c("../../../results/somatic_cnas/other_plots/density_ploidies.pdf",
                                "../../../results/somatic_cnas/other_plots/boxplot_losses_and_gains_all.pdf",
                                "../../../results/somatic_cnas/other_plots/boxplot_llg_mlg_and_loh_cnloh_prism_wgd.pdf",
                                "../../../results/somatic_cnas/other_plots/boxplot_hlg_and_hd_prism_wgd.pdf",
                                "../../../results/somatic_cnas/other_plots/boxplot_chr_gains_and_chr_losses_prism_wgd.pdf",
                                "../../../results/somatic_cnas/other_plots/boxplot_ogs_and_tsgs_prism_wgd.pdf"))
  parser$add_argument("--output_plots_paper", nargs="+", type="character",
                      default=c("../../../results/figures_paper/F2c.eps",
                                "../../../results/figures_paper/F2d.eps",
                                "../../../results/figures_paper/FS11a.pdf",
                                "../../../results/figures_paper/FS11b.pdf",
                                "../../../results/figures_paper/FS11c.pdf",
                                "../../../results/figures_paper/FS11d.pdf"))
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
