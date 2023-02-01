# @created: 17 Mar 22
# @modified: 01 Jun 22
# @authors: Yoann Pradat
#
# Draw violin plots showing age and drug count distribution per tumor type. 

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# functions ============================================================================================================

main <- function(args){
  # load data
  df_cln <- load_cln(args$cohort)

  col_sid <- "Subject_Id"
  col_drug <- "Drugs_Before_Biopsy"
  col_type <- "Project_TCGA_More"
  col_type_original <- "Project_TCGA_More_Original"
  col_age <- "Age_At_Biopsy"
  col_time <- "Survival_Time"
  col_status <- "Survival_Status"
  col_met <- "N_Metastatic_Sites"
  cols_keep <- c(col_sid, col_drug, col_type, col_type_original, col_age, col_time, col_status, col_met)

  # add additional tumor type correspond to full cohort
  df_cln <- bind_rows(df_cln, df_cln %>% mutate(!!col_type:=paste("All", toupper(args$cohort))))

  # group Not_TCGA tumor types into "Rare subtypes - Not_TCGA"
  rare_type <- "Rare subtypes"
  unknown_type <- "Unknown primary"
  df_cln <- df_cln %>% mutate(!!col_type:=ifelse(grepl("Not_TCGA", .data[[col_type]]), rare_type, .data[[col_type]]))
  df_cln <- df_cln %>% mutate(!!col_type:=ifelse(grepl("Unknown", .data[[col_type]]), unknown_type, .data[[col_type]]))

  # set Censored to NA
  df_cln <- df_cln %>% mutate(!!col_time:=ifelse(.data[["Survival_Status"]]=="Censored", NA, .data[[col_time]]))

  # transform type to add number of patients
  df_cln <- df_cln %>% group_by_at(all_of(col_type)) %>% mutate(Count=n(), .groups="drop")
  df_cln[[col_type_original]] <- df_cln[[col_type]]
  df_cln[[col_type]] <- paste0(df_cln[[col_type]], " (n = ", df_cln[["Count"]], ")")
  df_cln <- df_cln %>% mutate(!!col_time:=.data[[col_time]]/30)

  # drop patients with no drugs data except for All
  df_cln <- df_cln %>% select(all_of(cols_keep)) %>% filter(!is.na(.data[[col_drug]]) || .data[[col_type]]=="All")

  # drop tumor types with insufficiently many patient
  df_cln <- df_cln %>% group_by_at(all_of(col_type)) %>% filter(n() >= args$min_count) %>% ungroup()

  # set tumor type as factor which chosen order
  # put Rare subtypes and Unknown primary at the end
  ordered_types <- df_cln %>% group_by_at(all_of(col_type)) %>% summarize(Count=n()) %>% arrange(desc(Count)) %>%
    pull(var=col_type)
  types_end <- c(ordered_types[grepl(rare_type, ordered_types)],ordered_types[grepl(unknown_type, ordered_types)])
  ordered_types <- c(ordered_types[!ordered_types %in% types_end], types_end)

  df_cln[[col_type]] <- factor(df_cln[[col_type]], levels=rev(ordered_types))

  # split aggregated drugs into multiple rows and distinguish between None (no drug) and other drug names
  df_drug <- df_cln %>% separate_rows({{col_drug}}, sep="\\|")
  df_drug_ez <- df_drug %>% filter(.data[[col_drug]]=="None")
  df_drug_nz <- df_drug %>% filter(.data[[col_drug]]!="None")

  # build table containg drug counts per patient
  col_count <- paste0("Count_", col_drug)
  df_count_ez <- df_drug_ez %>% select(all_of(c(col_sid, col_type))) %>%
    mutate(!!col_count:=0)
  df_count_nz <- df_drug_nz %>% group_by_at(all_of(c(col_sid, col_type))) %>%
    summarize(!!col_count:=n(), .groups="drop")
  df_count <- bind_rows(df_count_nz, df_count_ez)
  df_count[[col_type]] <- factor(df_count[[col_type]], levels=rev(ordered_types))

  # set theme
  theme_set(theme_bw() +
            theme(axis.text = element_text(size=6, family="Helvetica", face="plain", color="black"),
                  axis.title = element_text(size=6, family="Helvetica", face="plain", color="black"),
                  legend.title = element_text(size=6, family="Helvetica", face="plain", color="black"),
                  legend.text = element_text(size=6, family="Helvetica", face="plain", color="black"),
                  legend.key.size = unit(0.2, "cm")))

  # colors
  colors_type <- load_colors(col_type)
  colors_type[[paste("All", toupper(args$cohort))]] <- "#D81E5B"
  type_original_to_new <- df_cln %>% select(all_of(c(col_type_original, col_type))) %>% distinct()

  colors_type <- colors_type[type_original_to_new[[col_type_original]]]
  names(colors_type) <- type_original_to_new[[col_type]]
  colors_fill <- list(all_tumor_type="#D81E5B", per_tumor_type="#6BAED6")

  # draw plot
  p1 <- ggplot(df_cln, aes_string(col_type, col_age, fill=col_type)) + 
    geom_violin(scale="area", size=0.1) +
    stat_summary(fun.y=median, geom="point", shape=20, size=1) +
    theme(panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(color="black", size=0.1, linetype="dashed"),
          panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank(), legend.position="bottom") + 
    ggtitle("") + xlab("") + ylab("Age at biopsy\n(Years)") + 
    scale_y_continuous(expand=c(0.01, 0.01), breaks=c(40,60,80)) +
    scale_fill_manual(values=colors_type) +
    theme(legend.position="none") +
    coord_flip() + 
    theme(axis.text=element_text(size=6), axis.title=element_text(size=6),
          plot.margin=unit(c(0, 0.1, 0.1, 0), "cm"))


  # draw plot
  p2 <- ggplot(df_cln, aes_string(col_type, col_met, fill=col_type)) + 
    geom_boxplot(outlier.size=0.5, outlier.shape=1, outlier.stroke=0.1, fatten=NULL, lwd=0.1) +
    stat_summary(fun.y=median, geom="point", shape=20, size=1) +
    theme(panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(color="black", size=0.1, linetype="dashed"),
          panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank(), legend.position="bottom") + 
    ggtitle("") + xlab("") + ylab("# Metastases\n") + 
    scale_y_continuous(expand=c(0.01, 0.01), limits=c(1,8), breaks=c(2,4,6)) +
    scale_fill_manual(values=colors_type) +
    theme(legend.position="none") +
    coord_flip() + 
    labs(x=NULL) +
    theme(axis.text.y=element_blank(), axis.text.x=element_text(size=6), axis.title=element_text(size=6),
          plot.margin=unit(c(0, 0.1, 0.1, 0), "cm"))

  # draw plot
  p3 <- ggplot(df_count, aes_string(col_type, col_count, fill=col_type)) + 
    geom_violin(scale="area", size=0.1) +
    stat_summary(fun.y=median, geom="point", shape=20, size=1) +
    theme(panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(color="black", size=0.1, linetype="dashed"),
          panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank(), legend.position="bottom") + 
    ggtitle("") + xlab("") + ylab("Drugs before\nbiopsy") + 
    scale_y_continuous(expand=c(0.01, 0.01), breaks=c(4,8,12)) +
    scale_fill_manual(values=colors_type) +
    theme(legend.position="none") +
    coord_flip() + 
    labs(x=NULL) +
    theme(axis.text.y=element_blank(), axis.text.x=element_text(size=6), axis.title=element_text(size=6),
          plot.margin=unit(c(0, 0.1, 0.1, 0), "cm"))


  # draw plot
  p4 <- ggplot(df_cln, aes_string(col_type, col_time, fill=col_type)) + 
    geom_violin(scale="area", size=0.1) +
    stat_summary(fun.y=median, geom="point", shape=20, size=1) +
    theme(panel.grid.major.y=element_blank(),
          panel.grid.major.x=element_line(color="black", size=0.1, linetype="dashed"),
          panel.grid.minor=element_blank()) +
    theme(axis.ticks=element_blank(), legend.position="bottom") + 
    ggtitle("") + xlab("") + ylab("Survival from\nbiopsy (Months)") + 
    scale_y_continuous(expand=c(0.01, 0.01)) +
    scale_y_log10(limits=c(0.4,80), breaks=c(0.5,6,90), labels=c(0.5, 6, 90)) +
    scale_fill_manual(values=colors_type) +
    theme(legend.position="none") +
    coord_flip() + 
    labs(x=NULL) +
    theme(axis.text.y=element_blank(), axis.text.x=element_text(size=6), axis.title=element_text(size=6),
          plot.margin=unit(c(0, 0.1, 0.1, 0), "cm"))


  # save figures
  p <- plot_grid(p1,p2,p3,p4, ncol=4, rel_widths=c(3.5, 1.25, 1.25, 1.25), labels=NULL)
  ggplot2::ggsave(args$outputs[1], p, width=90, height=60, units="mm", dpi=300)
  cat(paste("-plot saved at", args$outputs[1], "\n"))

  if (!is.null(args$outputs_paper)){
    ggplot2::ggsave(args$outputs_paper[1], p, width=90, height=60, units="mm", dpi=300)
    cat(paste("-plot saved at", args$outputs_paper[1], "\n"))
  }

  # save table of medians
  df_med <- df_cln %>% group_by_at(all_of(col_type)) %>%
    summarize(Median_Survival_Time=median(Survival_Time, na.rm=T), Median_Age=median(Age_At_Biopsy, na.rm=T)) %>%
    arrange(desc(Project_TCGA_More))

  df_med <- left_join(df_med, df_count %>% group_by_at(all_of(col_type)) %>%
    summarize(Median_Count_Drugs=median(Count_Drugs_Before_Biopsy, na.rm=T)))

  write.table(df_med, args$output_table, sep="\t", row.names=F, quote=F)
  cat(paste("-table saved at", args$output_table, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw violin plots showing age and drug count distribution per tumor type.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--min_count", type="integer", default=10,
                      help="Minimal number of patients in a given tumor type for it be showin in the plot")
  parser$add_argument('--outputs', type="character", nargs="+", help='Path to output plot.',
                  default=c("../../../results/data_overview/violins/violins_age_drug_surv_met.pdf"))
  parser$add_argument('--outputs_paper', type="character", nargs="+", help='Path to output plot for the paper.',
                  default=NULL)
  parser$add_argument('--output_table', type="character", help='Path to output plot.',
                  default="../../../results/data_overview/violins/violins_medians.tsv")
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
