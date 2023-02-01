suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

get_mutations_classifications <- function(){
  c(Upstream="5'UTR",
    Upstream="5'Flank",
    Downstream="3'UTR",
    Downstream="3'Flank",
    Missense="De_novo_Start_InFrame",
    Missense="De_novo_Start_OutOfFrame",
    FSIndel="Frame_Shift_Del",
    FSIndel="Frame_Shift_Ins",
    FSIndel="frameshift deletion",
    FSIndel="frameshift insertion",
    FSIndel="frameshift_variant",
    IFIndel="In_Frame_Del",
    IFIndel="In_Frame_Ins",
    IFIndel="nonframeshift insertion",
    IFIndel="nonframeshift deletion",
    StartSite="Start_Codon_Del",
    Indel="Start_Codon_Ins",
    Nonstop="Stop_Codon_Del",
    Indel="Stop_Codon_Ins",
    StartSite="Start_Codon_SNP",
    StartSite="Translation_Start_Site",
    StartSite="start_lost",
    Missense="Missense_Mutation",
    Missense="nonsynonymous SNV",
    Missense="missense_variant",
    Nonsense="Nonsense_Mutation",
    Nonsense="stopgain",
    Nonsense="stop_gained",
    Nonstop="Nonstop_Mutation",
    Silent="Silent",
    Silent="Intron",
    Splice="Splice_Site",
    Splice="Splice_Region",
    Splice="splice_acceptor_variant",
    Splice="splice_donor_variant",
    Splice="splicing")
}


get_mut_count <- function(df_mut, col_evt, col_evt_classes, evt_classes, col_gby="Tumor_Type"){
  variants_filtered <- c("IGR", "lincRNA", "RNA")
  df_mut_filt <- df_mut %>%
    filter(!.data[[col_evt]] %in% c("IGR", "lincRNA", "RNA"))

  evt_classes_filt <- evt_classes[evt_classes %in% unique(df_mut_filt[[col_evt_classes]])]

  df <- df_mut_filt %>%
    mutate(!!col_evt_classes:=fct_recode(.data[[col_evt_classes]], !!!evt_classes_filt)) %>%
    mutate(!!col_evt_classes:=as.character(.data[[col_evt_classes]])) %>%
    group_by(.data[[col_gby]], Subject_Id, Hugo_Symbol, .data[[col_evt_classes]]) %>%
    summarize(Count=n(), .groups="keep")

  df <- df  %>%
    group_by(.data[[col_gby]], Subject_Id, Hugo_Symbol) %>%
    mutate(Count_Subject=sum(Count)) %>%
    ungroup() %>%
    mutate(!!col_evt_classes:=ifelse(Count_Subject > 1, "Multihit", .data[[col_evt_classes]]), Count=1) %>%
    select(-c(Count_Subject)) %>%
    distinct() %>%
    group_by(.data[[col_gby]], Hugo_Symbol, .data[[col_evt_classes]]) %>%
    summarize(Count=sum(Count), .groups="keep")

  df
}


check_min_counts <- function(min_counts, cohorts){
  if (length(min_counts)==1){
    min_counts <- rep(min_counts, length.out=length(cohorts))
  } else {
    stopifnot(length(min_counts)==length(cohorts))
  }

  setNames(rep(min_counts, length.out=length(cohorts)), cohorts)
}


order_dfs <- function(dfs, col_var){
  order_var <- dfs$count_row %>% arrange(Count) %>% pull(col_var) 

  for (name in names(dfs)){
    if (col_var %in% colnames(dfs[[name]])){
      dfs[[name]] <- dfs[[name]] %>% arrange(match(.data[[col_var]], order_var))
    }
  }

  dfs
}


draw_facetted_heatmap <- function(dfs, col_var, names2colors, colors_palette_heatmap, alpha_left, alpha_heatmap, 
                                  names_left_pvals=NULL, names_heatmap_pvals=NULL, col_stack=NULL, stacks2colors=NULL,
                                  ...){
  names2plots <- list(left_pvals=names_left_pvals,
                      heatmap="percent",
                      heatmap_hover="count",
                      heatmap_pvals=names_heatmap_pvals,
                      middle_right_bar="count_row")
  if (!is.null(col_stack)){
    names2plots[["extreme_right_bar"]] <- "fraction_row"
  }

  # order rows for the heatmap
  dfs <- order_dfs(dfs, col_var=col_var)

  figs <- facetted_heatmap_barplots(dfs=dfs, names2plots=names2plots, col_var=col_var, names2colors=names2colors,
                                    colors_palette_heatmap=colors_palette_heatmap,
                                    alpha_left=alpha_left, alpha_heatmap=alpha_heatmap, col_stack=col_stack,
                                    width_one=50, height_one=20, stacks2colors=stacks2colors, ...)

  # Assemble figures
  # fig_left <- subplot(figs[["left_pvals"]], figs[["heatmap"]], margin=0.01, widths=c(0.15, 0.85))
  fig_left <- figs[["heatmap"]]

  if (is.null(col_stack)){
    fig_right <- figs[["middle_right_bar"]]
    fig <- subplot(fig_left, fig_right, margin=0.01, widths=c(0.8, 0.2))
  } else {
    fig_right <- subplot(figs[["middle_right_bar"]], figs[["extreme_right_bar"]], margin=0.01, widths=c(0.5,0.5))
    fig <- subplot(fig_left, fig_right, margin=0.05, widths=c(0.6, 0.4))
  }
  fig <- fig %>% layout(plot_bgcolor='rgba(0,0,0,0)')

  fig
}
