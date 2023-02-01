# @created: 07 Feb 22
# @modified: 12 May 22
# @authors: Yoann Pradat
#
# Forest-like tables showing the characteristics of the distributions of coefficients estimates of a given model.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(forester))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tools))

# functions ============================================================================================================

format_table <- function(df_cov, tab_size=2, ci_type="Percentile"){
  tab_1 <- paste0(rep(" ", tab_size), collapse="")
  tab_2 <- paste0(rep(" ", tab_size*2), collapse="")

  # add class for plot
  # for clinical covariates, only show the first level
  df_cov <- df_cov %>%
    mutate(Class_Plot=ifelse(Class_Lvl_1=="Clinical", Class_Lvl_1, paste(Class_Lvl_1, Class_Lvl_2, sep="-"))) %>%
    mutate(Class_1_Plot=gsub("_", " ", Class_Lvl_1)) %>%
    mutate(Class_2_Plot=ifelse(Class_Lvl_1=="Clinical", NA, paste0(tab_1, gsub("_", " ", Class_Lvl_2)))) %>%
    arrange(Class_Plot) %>%
    mutate(Covariate_Plot=gsub("_", " ", Plot_Name)) %>%
    mutate(Covariate_Plot=ifelse(is.na(Class_2_Plot), paste0(tab_1, Covariate_Plot), paste0(tab_2, Covariate_Plot)))

  cols_plot <- colnames(df_cov)[grepl("_Plot$", colnames(df_cov))]
  cols_plot <- c(cols_plot, "N_Error", "N_Warning", "N_Selected", "N_Not_Selected")
  df_cov <- df_cov %>% mutate(N_Not_Selected = N_Repeats - N_Error - N_Warning - N_Selected)
  col_ci_low <- paste0("Coefficient_CI_Low_", ci_type)
  col_ci_high <- paste0("Coefficient_CI_High_", ci_type)
  cols_plot <- c(cols_plot, "Coefficient", col_ci_low, col_ci_high)
  df_plot <- df_cov %>% select(all_of(cols_plot)) %>%
    rename(Coefficient_CI_Low=.data[[col_ci_low]], Coefficient_CI_High=.data[[col_ci_high]])

  # add lines with coefficient estimates for delimiting classes
  # add level 1 classes
  dfs_plot <- list()
  for (class_plot in unique(df_plot$Class_Plot)){
    df_plot_a <- df_plot %>% filter(Class_Plot==class_plot)
    dfs_plot_b <- list()
    for (i in 1:length(unlist(strsplit(class_plot, "-")))){
      dfs_plot_b[[i]] <- data.frame(Covariate_Plot=unique(df_plot_a[[paste0("Class_", i, "_Plot")]]))
    }
    df_plot_b <- Reduce(bind_rows, dfs_plot_b)
    dfs_plot[[class_plot]] <- bind_rows(df_plot_b, df_plot_a)
  }

  df_plot <- Reduce(bind_rows, dfs_plot) %>% rename(Covariate=Covariate_Plot)
  df_plot <- df_plot %>% distinct()

  df_plot
}


draw_plot <- function(df_plot, output, ci_conf=0.95, outcome="survival", ...){
  # draw forest plot
  xmin <- max(-10, floor(1.1*min(df_plot$Coefficient_CI_Low, na.rm=T)))
  xmax <- min(10, ceiling(1.1*max(df_plot$Coefficient_CI_High, na.rm=T)))

  if (outcome=="survival"){
    arrow_labels <- c("Better survival", "Worse survival")
    coeff_name <- "Log hazard ratio"
  } else {
    arrow_labels <- c(paste("Not", outcome, "resistant"), paste(outcome, "resistant"))
    coeff_name <- "Coefficient"
  }

  forester(left_side_data=df_plot[,1,drop=F],
           estimate=df_plot$Coefficient,
           ci_low=df_plot$Coefficient_CI_Low,
           ci_high=df_plot$Coefficient_CI_High,
           display=FALSE,
           # font_family="sans", # for reason that I ignore, using 'sans' provokes a mispositioning
           null_line_at=0,
           xlim=c(xmin, xmax),
           x_scale_linear=T,
           estimate_precision=2,
           arrows=T,
           arrow_labels=arrow_labels,
           estimate_col_name=paste0(coeff_name, " (", round(ci_conf*100, 1), "% CI)"),
           render_as=file_ext(output),
           file_path=output, ...)
}


categorize_covariates <- function(df_cov_final, df_cov_remov, class_lvls=c("Class_Lvl_1", "Class_Lvl_2")){
  # align column types before binding
  cov_final_types <- sapply(colnames(df_cov_final), function(x) class(df_cov_final[[x]]))
  cov_remov_types <- sapply(colnames(df_cov_remov), function(x) class(df_cov_remov[[x]]))
  for (col in names(cov_remov_types)){
    if (col %in% colnames(df_cov_final)){
      col_final_type <- cov_final_types[[col]]
      col_remov_type <- cov_remov_types[[col]]
      if (col_final_type!=col_remov_type){
        df_cov_remov[[col]] <- eval(parse(text=(paste0("as.",col_final_type,"(df_cov_remov[[col]])"))))
      }
    }
  }

  # concat tables
  if (nrow(df_cov_remov)>0){
    df_cov_all <- bind_rows(df_cov_final, df_cov_remov)
  } else {
    df_cov_all <- df_cov_final %>% mutate(Reason=NA)
  }
  df_cov_all <- df_cov_all %>% unite("Class_Cat", all_of(class_lvls), sep="-", remove=F)
  
  # manual grouping of some categories
  df_cov_all <- df_cov_all %>% mutate(Class_Cat=ifelse(Class_Lvl_1=="Clinical", "Clinical", Class_Cat)) %>%
    mutate(Class_Cat=ifelse(grepl("Total", Class_Lvl_2) & grepl("DNA", Class_Lvl_1), "DNA-Total_Counts", Class_Cat)) %>%
    mutate(Class_Cat=ifelse(grepl("Total", Class_Lvl_2) & grepl("RNA", Class_Lvl_1), "RNA-Total_Counts", Class_Cat)) %>%
    mutate(Class_Cat=ifelse(grepl("CNA_Summary_Statistics", Class_Lvl_2), "CNA_Summary_Stats", Class_Cat))

  # add classes and colors for the pie plots
  df_cov_all <- df_cov_all %>% mutate(Class_Pie=ifelse(!is.na(Reason), Reason, NA))
  mask_final <- is.na(df_cov_all$Reason)

  reasons_all <- unique(df_cov_remov$Reason)
  reason_redundant <- reasons_all[grepl("^Redundant R", reasons_all)]
  reason_nonzero <- reasons_all[grepl("non-zero", reasons_all)]

  class2color_pie <- list("Exactly redundant"="#9E4770")
  if (length(reason_redundant)>0) class2color_pie[[reason_redundant]] <- "#FFC8DD"

  class2color_pie[["Uniquely valued"]] <- "#D7B377"
  if (length(reason_nonzero)>0) class2color_pie[[reason_nonzero]] <- "#E5E5E5"

  selection_threshold <- seq(0,1,by=0.25)
  colors_selection <- colorRampPalette(c("#d8f3dc", "#2d6a4f"))(length(selection_threshold)+1)

  for (i in 1:(length(selection_threshold)+1)){
    if (i==1){
      df_cov_all <- df_cov_all %>% mutate(Class_Pie=ifelse(mask_final & Selected==0, "0%", Class_Pie))
      class2color_pie[["0%"]] <- colors_selection[i]
    } else if (i==length(selection_threshold)+1) {
      df_cov_all <- df_cov_all %>% mutate(Class_Pie=ifelse(mask_final & Selected==1, "100%", Class_Pie))
      class2color_pie[["100%"]] <- colors_selection[i]
    } else {
      mask_left <- df_cov_all$Selected > selection_threshold[i-1]
      if (i==length(selection_threshold)){
        mask_right <- df_cov_all$Selected < selection_threshold[i]
        name <- paste0("]", selection_threshold[i-1], ",", selection_threshold[i], "[")
      } else {
        mask_right <- df_cov_all$Selected <= selection_threshold[i]
        name <- paste0("]", selection_threshold[i-1], ",", selection_threshold[i], "]")
      }
      df_cov_all <- df_cov_all %>% mutate(Class_Pie=ifelse(mask_final & mask_left & mask_right, name, Class_Pie))
      class2color_pie[[name]] <- colors_selection[i]
    }
  }

  stopifnot(sum(is.na(df_cov_all$Class_Pie))==0)
  list(df_cov_all=df_cov_all, class2color_pie=class2color_pie)
} 


make_count <- function(df, field){
  df_count <- df %>%
    rename(Label={{ field }}) %>%
    replace_na(list(Label="N/A")) %>%
    group_by(Label) %>%
    summarize(Count=n())  %>% 
    mutate(Percentage=Count/sum(Count)) %>%
    ungroup()

  df_count
}


pie_count_plotly <- function(fig, df, field, row, column, pal=NULL, values2colors, textinfo='label'){
  df_count <- make_count(df, field)

  if(!is.null(values2colors)){
    df_count <- df_count %>%
      mutate(Color=recode(Label, !!!values2colors))
  } else {
    df_count <- df_count %>%
      mutate(Color=get_label_colors(Label, pal=pal))
  }

  # reorder to match order of values2colors
  df_count$Label <- factor(df_count$Label, levels=names(values2colors))
  df_count <- df_count %>% arrange(Label)

  fig <- fig %>% add_pie(labels=df_count$Label, values=df_count$Count, textposition="none",
                         textinfo=textinfo,
                         text=df_count$Percentage,
                         sort=F,
                         hovertemplate = paste0(paste0("<b>%{label}</b>"),
                                                '<br><br>Count: %{value}<br>',
                                                'Percent: %{text:.3f}',
                                                "<extra></extra>"),
                         hole=0.6,
                         marker = list(line=list(color='#FFFFFF', width=1), colors=df_count$Color),
                         domain = list(row=row, column=column))

  fig
}


pie_counts_plotly <- function(dfs, main, anns, field, values2colors, fontsize=10, xgap=0.25, max_per_row=4){
  annotations <- list()
  n_plot <- length(dfs)

  font <- list(size=fontsize)
  font_anns <- list(size=fontsize*1.25, family="Arial")
  font_title <- list(size=fontsize*2, family="Arial")

  fig <- plot_ly()

  # get grid layout
  if (n_plot < max_per_row){
    columns <- n_plot
  } else {
    columns <- max_per_row
  }
  rows <- ceiling(n_plot/max_per_row)


  for (i in 1:n_plot){
    df <- dfs[[i]]
    ann <- anns[[i]]
    row <- floor((i-1)/max_per_row)
    column <- (i-1) %% max_per_row

    fig <- pie_count_plotly(fig, df=df, field=field, values2colors=values2colors, row=row, column=column)
    annotation <- list(text=ann, showarrow=F, font=font_anns, xanchor="center", xref=paste0("x",column+1),
                       yref=paste0("y",row+1))
    annotations[[i]] <- annotation
  }


  fig <- fig %>% layout(title = list(text=main, font=font_title, xref="paper", yref="paper", x=0.5, y=1),
                        font = font, showlegend = T, annotations=annotations,
                        grid = list(rows=rows, columns=columns, xgap=xgap))

  emptyaxis <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)
  for (row in 1:rows){
    eval(parse(text=paste0("fig <- fig %>% layout(yaxis", row, "=emptyaxis)")))
  }
  for (column in 1:columns){
    eval(parse(text=paste0("fig <- fig %>% layout(xaxis", column, "=emptyaxis)")))
  }

  fig
}


main <- function(args){
  # load tables
  df_met <- read_tsv(args$input_met, progress=F, show_col_types=F)
  df_cov <- read_tsv(args$input_cov, progress=F, show_col_types=F)
  df_cov_final <- read_tsv(args$cov_final, progress=F, show_col_types=F)
  df_cov_remov <- read_tsv(args$cov_remov, progress=F, show_col_types=F)

  # map codes back to names and add additional info
  df_cov <- left_join(df_cov %>% rename(Code=Covariate), df_cov_final, by="Code")

  # create new Selected frequency that includes models with errors or warnings
  df_cov <- df_cov %>% mutate(Selected_WE=(N_Selected + N_Warning + N_Error)/N_Repeats)

  # split covariates from the model based on the selection frequency
  df_cov_in_plot <- df_cov %>% filter(Selected_WE >= args$threshold)
  df_cov_out_plot <- df_cov %>% filter(Selected_WE < args$threshold)

  # get table ready for plot
  df_plot <- format_table(df_cov_in_plot)
  df_plot$Index <- rev(0:(nrow(df_plot)-1))

  # draw side plot showing percentage of selection
  if (args$threshold!=1){
    df_bar <- df_plot %>% select(Index, N_Selected, N_Warning, N_Error, N_Not_Selected) %>% 
      gather("Category", "# Models", -Index)
    df_bar$Category <- factor(df_bar$Category, levels=c("N_Not_Selected", "N_Error","N_Warning", "N_Selected"))
    levels(df_bar$Category) <- c("Not selected", "Selected & model error", "Selected & model warning", "Selected")
    add_plot <- ggplot(df_bar, aes(fill=Category, x=`# Models`, y=Index)) + geom_col(orientation="y") +
      theme_classic() + # base theme
      scale_fill_manual(values = c("#FE5F55", "#F6AE2D", "#C7EFCF", "#157A6E"),
                        breaks = c("Selected & model error", "Selected & model warning", "Not selected", "Selected")) +
        theme(text = ggplot2::element_text(family="mono", size=12),
              axis.line.x = ggplot2::element_line(colour="black"),
              axis.ticks.length.x = grid::unit(.07, "in"),
              axis.title.y = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_blank(),
              axis.ticks.y = ggplot2::element_blank(),
              axis.line.y = ggplot2::element_blank(),
              panel.background = ggplot2::element_rect(fill = "transparent"),
              plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              legend.background = ggplot2::element_rect(fill = "transparent"),
              legend.box.background = ggplot2::element_blank(),
              legend.title = element_blank()) +
        scale_x_continuous(n.breaks=3)

  } else {
    add_plot <- NULL
  }

  # draw forest plot
  dir.create(file.path(dirname(args$output_forester)), showWarnings=F, recursive=T)
  draw_plot(df_plot, output=args$output_forester, ci_conf=args$ci_conf, outcome=args$outcome, add_plot=add_plot,
            add_plot_width=35)

  # draw donut plots showing proportion of feature removed (with reason for removal) and proportion features
  # selected at different thresholds. each donut represent
  df_cov_final <- left_join(df_cov_final, df_cov[,c("Covariate", "Selected")], by="Covariate")
  out <- categorize_covariates(df_cov_final=df_cov_final, df_cov_remov=df_cov_remov)

  df_cov_all <- out$df_cov_all
  class2color_pie <- out$class2color_pie
  classes_cat <- unique(df_cov_all$Class_Cat)
  dfs_cov_all <- lapply(classes_cat, function(cl) df_cov_all %>% filter(Class_Cat==cl))
  dfs_cov_all <- setNames(dfs_cov_all, classes_cat)

  field <- "Class_Pie"
  main <- NULL
  anns <- sapply(classes_cat, function(x) paste0(paste(unlist(str_split(x, "-")), collapse="\n"),
                                                 "\n", paste(nrow(dfs_cov_all[[x]]), "features")))
  max_per_row <- 3

  fig <- pie_counts_plotly(dfs=dfs_cov_all, main=main, anns=anns, field=field, values2colors=class2color_pie,
                           fontsize=15, xgap=0.25, max_per_row=max_per_row)
  fig <- fig %>% layout(legend=list(font=list(family="Helvetica", color="black", size=24), x=1, y=0.5))
  
  width <- 200 + 500*min(length(dfs_cov_all), max_per_row)
  height <- 400*ceiling(length(dfs_cov_all)/max_per_row)

  dir.create(file.path(dirname(args$output_pieplots)), showWarnings=F, recursive=T)
  save_image(fig, args$output_pieplots, width=width, height=height)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  results_folder <- "../../../results/survival_analysis"
  data_folder <- file.path(results_folder, "data_prism/sub_features/BRCA_dna_and_rna/cln_biol_1_dna_1_rna_1")
  models_folder <- file.path(results_folder, "models_prism/BRCA_dna_and_rna/cln_biol_1_dna_1_rna_1/none_coxph_standard")
  plots_folder <- file.path(results_folder, "plots_prism/BRCA_dna_and_rna/cln_biol_1_dna_1_rna_1")

  parser <- ArgumentParser(description='Pool model quality metrics and coefficients estimates across bootstraps.')
  parser$add_argument("--input_met", type="character", help="Path to input table of pooled model quality metrics.",
                      default=file.path(models_folder, "mets.pooled_ax_rep.tsv.gz"))
  parser$add_argument("--input_cov", type="character", help="Path to input table of pooled model coefficients estimates.",
                      default=file.path(models_folder, "covs.pooled_ax_rep.tsv.gz"))
  parser$add_argument("--cov_final", type="character", help="Path to table mapping covariate codes to names.",
                      default=file.path(data_folder, "processed/covs.final.tsv.gz"))
  parser$add_argument("--cov_remov", type="character", help="Path to table mapping covariate codes to names.",
                      default=file.path(data_folder, "processed/covs.removed.tsv.gz"))
  parser$add_argument("--threshold", type="double", default=0.2,
                      help="Minimum percentage of selection of a covariate to be included inthe plot.")
  parser$add_argument("--ci_conf", type="double", default=0.95,
                      help="Confidence level of intervals (for plot labelling).")
  parser$add_argument("--outcome", type="character", default="survival",
                      help="Outcome name (for plot labelling).")
  parser$add_argument("--output_forester", type="character", help="Path to output plot.",
                default=file.path(plots_folder, "forester_coeffs_none_coxph_standard.pdf"))
  parser$add_argument("--output_pieplots", type="character", help="Path to output plot.",
                default=file.path(plots_folder, "pieplots_coeffs_none_coxph_standard.pdf"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args = parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
