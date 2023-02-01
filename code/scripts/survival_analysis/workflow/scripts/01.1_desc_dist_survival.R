# @created: 20 Oct 21
# @modified: 21 Oct 21
# @authors: Yoann Pradat
#
# Describe the distribution of overall survival in PRISM and per tumor type.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggridges))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

source("workflow/functions/utils.R")

# functions ============================================================================================================

draw_plot_dist <- function(df, output, width, height, col_x="Survival_Time", col_y="Tumor_Type", 
                           col_z="Survival_Status", col="All", fill="All", xlab="Survival time (months)"){

  ggplot(df, aes_string(
      x=col_x, y=col_y, 
      col=col, fill=fill)) + 
      facet_wrap(~df[[col_z]], scales = "free_x", nrow = 1) + 
      stat_density_ridges(quantile_lines = TRUE, quantiles = 2, alpha=0.2,
                          jittered_points = TRUE,
                          position = position_points_jitter(width = 0.05, height = 0),
                          point_shape = '|', point_size = 1, point_alpha = 1) + 
      theme_ridges() + theme(
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
      xlab(xlab) + ylab(gsub("_", " ", col_y))
}


draw_plot_hist <- function(df, col_x="Count_Drugs_Before_Biopsy", col_y="Tumor_Type", col_z="Variable", fill="All",
                           xlab="# drugs received"){

  ggplot(df, aes_string(x=col_x)) + 
      facet_grid(df[[col_y]]~df[[col_z]], scales="free") + 
      geom_bar(aes(fill=All, color=All), alpha=0.2) + 
      scale_y_continuous(breaks=c(0, NULL), minor_breaks=c(NULL), position="right") +
      theme_ridges() + theme(
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          strip.text.y.right = element_text(angle=0),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
      xlab(xlab) + ylab("")
}


main <- function(args){
  # load tables
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  df_bio <- load_bio(args$cohort)

  # transform days to months
  df_cln$Survival_Time <- df_cln$Survival_Time/30

  # build table for plot 1
  cols_plot <- c("Subject_Id", "Survival_Time", "Tumor_Type", "Survival_Status")
  df_plot <- df_cln[,cols_plot] %>% filter(!is.na(Survival_Time))
  df_plot_1 <- df_plot %>% mutate(All="No") %>% group_by(Tumor_Type) %>% filter(n()>args$min_count)
  df_plot_2 <- df_plot %>% mutate(All="Yes", Tumor_Type="All tumors")
  df_plot <- bind_rows(df_plot_1, df_plot_2)

  # choose row order
  tumor_types_ord <- df_plot_1 %>% filter(Survival_Status=="Deceased") %>%
    group_by(Tumor_Type) %>% summarize(Median=quantile(Survival_Time, 0.5)) %>% arrange(Median) %>% 
    pull(Tumor_Type)
  levels_ord <- rev(c("All tumors", tumor_types_ord))
  df_plot["Tumor_Type"] <- factor(df_plot[["Tumor_Type"]], levels = levels_ord)

  # draw
  p1 <- draw_plot_dist(df_plot)
  p1 <- p1 + scale_x_continuous(limits = c(0.5, NA), trans="log", breaks=c(0.5, 6, 120))
  ggsave(filename=args$output[[1]], plot=p1, width=args$width, height=args$height, units="in")

  cols_x <- c("Drugs_Before_Biopsy", "Metastatic_Sites")
  xlabs <- c("# Drugs received", "# Metastatic sites")
  names <- c("Drugs", "Metastatic sites")
  outputs <- args$outputs[2:3]
  
  for (i in 1:2){
    # build table for hist plot
    col_x <- cols_x[[i]]
    col_count_x <- paste("Count", col_x, sep="_")
    df_cln <- add_col_count(df_cln, col_x)

    cols_plot <- c("Subject_Id", col_count_x, "Tumor_Type")
    df_plot <- df_cln[,cols_plot] %>% filter(!is.na(.data[[col_count_x]]))
    df_plot_1 <- df_plot %>% mutate(All="No") %>% group_by(Tumor_Type) %>% filter(n()>args$min_count)
    df_plot_2 <- df_plot %>% mutate(All="Yes", Tumor_Type="All tumors")
    df_plot <- bind_rows(df_plot_1, df_plot_2)
    df_plot[["Variable"]] <- names[[i]]
    levels_ord <- c("All tumors", tumor_types_ord)
    df_plot["Tumor_Type"] <- factor(df_plot[["Tumor_Type"]], levels = levels_ord)

    p2 <- draw_plot_hist(df_plot, col_x=col_count_x, xlab=xlabs[[i]])
    p2 <- p2 + scale_x_continuous(limits=c(NA,12), breaks=seq(1,12,2))

    ggsave(filename=outputs[[i]], plot=p2, width=args$width, height=args$height, units="in")
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe distribution of survival times.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument('--min_count', type="integer", default=10,
                      help="Minimum number of tumors in a given class.")
  parser$add_argument("--width", type="double", default=8, help="Width of output plots in inches.")
  parser$add_argument("--height", type="double", default=10, help="Height of output plots in inches.")
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plot.",
                      default=c("../../../results/survival_analysis/description/dist_survival_times.pdf",
                                "../../../results/survival_analysis/description/dist_counts_drug.pdf",
                                "../../../results/survival_analysis/description/dist_counts_mets.pdf"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  for (name in names(args)){
    if (!is.null(args[[name]]) && length(args[[name]])==1){
      if (args[[name]]=="None"){
        args[[name]] <- NULL
      }
    }
  }

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
