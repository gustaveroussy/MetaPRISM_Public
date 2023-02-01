# @modified: 25 Jun 21
# @modified: 17 Aug 21
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))

# functions ============================================================================================================

draw_burden_plot_fusion <- function(dfs_burden, cohorts2colors){
  cohorts <- names(dfs_burden)
  offset <- 0.1
  tt_order <- compute_groups_order_burden_plot(df=bind_rows(dfs_burden),
                                               col_burden="Fusion_Burden",
                                               col_groups="Tumor_Type")

  dfs_plot <- list()
  for (stack in cohorts){
    dfs_plot[[stack]] <- compute_coordinates_burden_plot(dfs_burden[[stack]], groups_order=tt_order,
                                                         col_burden="Fusion_Burden", col_groups="Tumor_Type",
                                                         offset=offset)
  }

  fig <- draw_burden_plot(dfs=dfs_plot, groups=tt_order, stacks=cohorts, yrange=c(10^-1, 10^3), yname="Fusion burden",
                          ytitle="Fusions prevalence<br>(number of fusions per tumor)", marker_opacity=1,
                          stacks2colors=cohorts2colors, offset=offset, min_median_size=10)
  fig
}

main <- function(args){
  # load
  dfs_burden <- setNames(lapply(args$input_tables, read_tsv, col_types=cols()), args$input_cohorts)
  dfs_burden <- lapply(dfs_burden, function(df) df %>% replace(.==0, 0.5))

  # plot
  cohorts2colors <- load_colors(sheet="Global")[args$input_cohorts]
  fig <- draw_burden_plot_fusion(dfs_burden, cohorts2colors)

  # save
  orca(fig, args$output, width=args$output_width, height=args$output_height)
  cat(paste("-plot saved at", args$output, "\n"))
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Make plot of fusions burden.')
  parser$add_argument("--input_cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--input_tables", nargs="+", help="Paths to input tables.",
                        default=c("workflow/results/tables/prism/burden_table_sel.tsv",
                                  "workflow/results/tables/met500/burden_table_sel.tsv",
                                  "workflow/results/tables/tcga/burden_table_sel.tsv"))
  parser$add_argument("--output_width", type="integer", default=1200, help="Width of output plot in pixels.")
  parser$add_argument("--output_height", type="integer", default=800, help="Width of output plot in pixels.")
  parser$add_argument("--output", type="character", help="Path where output plot will be saved.",
                      default="workflow/results/plots/plot_burden_tumor_types_sel.pdf")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
