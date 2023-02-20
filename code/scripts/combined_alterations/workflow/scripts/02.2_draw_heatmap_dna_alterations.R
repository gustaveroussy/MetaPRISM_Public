# @created: 28 Sep 21
# @modified: 14 Dec 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xlsx))

# functions ============================================================================================================

load_oncogenes_or_tumorsuppressor <- function(gene_type){
  df_genes <- load_from_data("./resources/curated/cancer_genes_curated.tsv", progress=F)

  if (gene_type=="oncogenes"){
    genes <- df_genes %>%
      filter(Is_Oncogene=="Yes" & Inclusion_Level==1) %>%
      pull(Hugo_Symbol)
  } else if (gene_type=="tumorsuppressors"){
    genes <- df_genes %>%
      filter(Is_Tumor_Suppressor_Gene=="Yes" & Inclusion_Level==1 & Is_Oncogene=="No") %>%
      pull(Hugo_Symbol)
  } else {
    stop("- Unsupported value of --gene_type. Choose 'oncogenes' or 'tumorsuppressors.")
  }


  genes
}


select_multiple_breakpoints_fusions <- function(df_fus){
  # add a Coverage column that will be the max of all *_Total_Reads columns
  cols_total <- colnames(df_fus)[grepl("_Total_Reads$", colnames(df_fus))]
  df_fus <- df_fus %>% mutate(Coverage=pmax(!!!syms(cols_total), na.rm=T))

  # select 1 breakpoint per fusion with the following rules
  # - prioritize breakpoints seen by multiple callers
  # - if ties left, prioritize breakpoints with the highest Coverage

  # add number of callers
  df_fus$N_Algos_Wt_Breakpoint <- as.vector(sapply(df_fus$Algo_Wt_Breakpoint, function(x) 
                                                   length(unlist(str_split(x, "\\|")))))
  df_fus <- df_fus %>% unite(Row_Id_Wo_Breakpoint, Sample_Id, Fusion_Id, sep="_", remove=F) %>%
    arrange(Row_Id_Wo_Breakpoint, desc(N_Algos_Wt_Breakpoint), desc(Coverage))
  df_fus <- df_fus[!duplicated(df_fus$Row_Id_Wo_Breakpoint),]

  df_fus
}


select_bidirectional_fusions <- function(df_fus){
  df_fus$Meta_Id <- apply(df_fus[,c("Gene_1", "Gene_2")], 1, function(x) paste(sort(x), collapse="/"))
  df_fus$Fusion_Id_Reverse <- apply(df_fus[,c("Gene_1", "Gene_2")], 1, function(x) paste(rev(x), collapse="--")) 

  df_counts_fwd <- df_fus %>%
    group_by(Fusion_Id, Meta_Id, Gene_1, Gene_2) %>% summarize(n_fwd=n(), .groups="keep") %>% arrange(-n_fwd)

  df_counts_rev <- df_counts_fwd %>%
    rename(n_rev=n_fwd) %>%
    unite(Fusion_Id, c("Gene_2", "Gene_1"), sep="--")

  df_counts <- left_join(df_counts_fwd, df_counts_rev, by=c("Meta_Id", "Fusion_Id")) %>% replace(is.na(.), 0)

  df_translate <- df_counts %>%
    group_by(Meta_Id) %>% filter(n() > 1) %>% arrange(Meta_Id, n_fwd) %>% group_by(Meta_Id) %>%
    slice(1) %>% unite(Fusion_Id_Good, c("Gene_2", "Gene_1"), sep="--")

  # add direction
  df_fus <- df_fus %>% mutate(Direction="fwd")

  df_fus_good <- df_fus %>% filter(!Fusion_Id %in% df_translate$Fusion_Id)
  df_fus_revr <- df_fus %>% filter(Fusion_Id %in% df_translate$Fusion_Id) %>%
    mutate(Tmp=Fusion_Id, Fusion_Id=Fusion_Id_Reverse, Fusion_Id_Reverse=Tmp) %>%
    mutate(Tmp=Gene_1, Gene_1=Gene_2, Gene_2=Tmp) %>%
    mutate(Direction="rev") %>%
    select(-c(Tmp))

  # assemble and remove duplicated fusions in same sample. Record the information that the fusion is seen in 
  # both directions in the sample
  df_fus <- bind_rows(df_fus_good, df_fus_revr)
  df_fus_dir <- df_fus %>% group_by(Sample_Id, Fusion_Id) %>%
    summarize(Direction=paste0(Direction, collapse="|"), Coverage=paste0(Coverage, collapse="|"), .groups="keep")
  df_fus_sel <- df_fus %>% group_by(Sample_Id, Fusion_Id) %>% slice(which.max(Coverage)) %>%
    select(-c(Direction, Coverage))

  left_join(df_fus_sel, df_fus_dir, by=c("Sample_Id", "Fusion_Id"))
}


simplify_annotations_custom <- function(df_fus){
  annots_remove <- c("Chitars_All_Human_v5.0", "ONE_PARTNER_IS_DRIVER")
  df_fus$Annotations_Custom_Simple <- as.vector(sapply(df_fus$Annotations_Custom, function(x){
                                                         annots <- unlist(str_split(x, "\\|"));
                                                         annots <- annots[!annots %in% annots_remove];
                                                         if (length(annots)>2){
                                                           return("3 databases or more")
                                                         } else {
                                                           if (length(annots)==0){
                                                             return("Novel")
                                                           } else {
                                                             return(paste0(annots, collapse="|"))
                                                           }
                                                         }}))
  df_fus
}


collapse_to_gene <- function(df_fus){
  df_fus <- df_fus %>% mutate(Hugo_Symbol=Fusion_Id) %>% separate_rows(Hugo_Symbol, sep="--")
  levels_annots <- unique(df_fus$Annotations_Custom_Simple)
  levels_annots_mid <- levels_annots[!levels_annots %in% c("Novel", "3 databases or more")]
  levels_annots_mid <- levels_annots_mid[order(nchar(levels_annots_mid), levels_annots_mid)]
  levels_annots <- c("Novel", levels_annots_mid, "3 databases or more")
  df_fus$Annotations_Custom_Simple <- factor(df_fus$Annotations_Custom_Simple, levels=levels_annots)

  df_fus <- df_fus %>% arrange(desc(Annotations_Custom_Simple), desc(Coverage))
  df_fus <- df_fus %>% unite(Row_Id_Hugo_Symbol, Sample_Id, Hugo_Symbol, sep="_", remove=F)
  df_fus <- df_fus[!duplicated(df_fus$Row_Id_Hugo_Symbol),] %>% select(-Row_Id_Hugo_Symbol)

  df_fus
}


load_alterations_and_select_samples <- function(args, select_samples=T){
  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_pid <- "Sample_Id_DNA_P"
  col_tid <- "Sample_Id_DNA_T"
  col_nid <- "Sample_Id_DNA_N"
  col_rid <- "Sample_Id_RNA_T"
  col_sub <- "Subject_Id"
  col_sam <- "Sample_Id"

  if (args$sample_type=="DNA_T__DNA_N"){
    col_sam_cln <- col_pid
  } else if (args$sample_type=="DNA_T"){
    col_sam_cln <- col_tid
  } else if (args$sample_type=="DNA_N"){
    col_sam_cln <- col_nid
  } else if (args$sample_type=="RNA_T"){
    col_sam_cln <- col_rid
  }

  n_cohorts <- length(args$cohorts)
  dfs_alt <- list()

  for (i in 1:n_cohorts){
    cohort <- args$cohorts[[i]]
    input_table <- args$alterations[[i]]
    input_sample <- args$samples[[i]]
    cat(paste("-processing",  cohort, "...\n"))

    # load alterations
    cat("-loading alterations \n")
    df_alt <- load_table(input_table, guess_max=1e4)
    
    # add Sample_Id if it was in rownmaes
    if (colnames(df_alt)[1]=="...1"){
      df_alt <- df_alt %>% rename(!!col_sam:=`...1`)
    }

    # select in-design samples
    df_cln <- load_cln(study=cohort)
    df_cln <- df_cln %>% rename(Tumor_Type=Project_TCGA_More)

    if (args$sample_type=="DNA_T__DNA_N"){
      df_cln <- df_cln %>% unite(!!col_pid, all_of(c(col_tid, col_nid)), sep="_vs_", remove=F)
    }

    mask_in <- df_alt[[col_sam]] %in% df_cln[[col_sam_cln]][!is.na(df_cln[[col_sam_cln]])] 
    df_alt <- df_alt[mask_in,]
    cat("-selected", paste0(sum(mask_in), "/", length(mask_in)) , "alterations from in-design", args$sample_type,
        "samples ...\n")

    # add clinical data if missing
    cols_cln <- c(col_sub, "Tumor_Type")
    cols_cln <- setdiff(cols_cln, colnames(df_alt))

    if (length(cols_cln) > 0) {
      df_alt <- left_join(df_alt, df_cln[,c(cols_cln, col_sam_cln)], by=setNames(col_sam_cln, col_sam))
    }

    # select samples using table of samples
    if (select_samples){
      # select samples
      df_sam <- load_table(input_sample)
      df_sam <- df_sam %>% filter(.data[[paste0("Use_", args$selection_name)]]==1)
      mask_in <- df_alt[[col_sam]] %in% df_sam[[col_sam]]
      df_alt <- df_alt[mask_in,]
      cat("-selected", paste0(sum(mask_in), "/", length(mask_in)) , "alterations from selected", args$sample_type,
          "samples ...\n")
      cat("-number of unique subjects", df_alt %>% distinct(.data[[col_sub]]) %>% nrow(), "\n")
    }

    # if RNA_T, select 1 breakpoint per fusion and deal with the fusions seen in 2 directions in one sample
    if (args$sample_type=="RNA_T" & (args$gene_type=="genes" | args$gene_type=="fusions")){
      df_alt <- select_multiple_breakpoints_fusions(df_fus=df_alt)
      df_alt <- select_bidirectional_fusions(df_fus=df_alt)
      df_alt <- simplify_annotations_custom(df_fus=df_alt)
    }

    dfs_alt[[cohort]] <- df_alt
  }


  # When aggregating by Hugo_Symbol, a sample may have multiple fusions involving the same gene but with different
  # values for Annotations_Custom_Simple. In that case, known fusions have priority over Novel. For instance, 
  # if a sample has PTEN involved in a known COSMIC fusion and a novel fusion, annotation for this sample and this
  # gene will be "COSMIC_curated_v95" in the plot. 
  if (args$sample_type=="RNA_T" & args$gene_type=="genes"){
    known_fusions <- unique(dfs_alt$prism %>% group_by(Fusion_Id) %>% filter(n()>=4) %>% pull(var="Fusion_Id"))

    for (cohort in args$cohorts){
      df_alt <- dfs_alt[[cohort]]
      df_alt <- df_alt %>% filter(!Fusion_Id %in% known_fusions)
      df_alt <- collapse_to_gene(df_fus=df_alt)
      dfs_alt[[cohort]] <- df_alt
    }
  }

  dfs_alt
}


get_alterations_classifications <- function(){
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
    Splice="splicing",
    Amplification="Amplification",
    Deletion="Deletion",
    "Novel"="Novel",
    "ChimerKB_v4.0"="ChimerKB_v4.0",
    "Chitars_Cancer_v5.0"="Chitars_Cancer_v5.0",
    "COSMIC_curated_v95"="COSMIC_curated_v95",
    "ChimerKB_v4.0|Chitars_Cancer_v5.0"="ChimerKB_v4.0|Chitars_Cancer_v5.0",
    "ChimerKB_v4.0|COSMIC_curated_v95"="ChimerKB_v4.0|COSMIC_curated_v95",
    "Chitars_Cancer_v5.0|COSMIC_curated_v95"="Chitars_Cancer_v5.0|COSMIC_curated_v95",
    "3 databases or more"="3 databases or more")
}


get_alt_count <- function(df_alt, col_evt, col_evt_classes, evt_classes, col_gby="Tumor_Type"){
  variants_filtered <- c("IGR", "lincRNA", "RNA")
  df_alt_filt <- df_alt %>%
    filter(!.data[[col_evt]] %in% c("IGR", "lincRNA", "RNA"))

  if (!is.null(col_evt_classes)){
    evt_classes_filt <- evt_classes[evt_classes %in% unique(df_alt_filt[[col_evt_classes]])]

    df <- df_alt_filt %>%
      mutate(!!col_evt_classes:=fct_recode(.data[[col_evt_classes]], !!!evt_classes_filt)) %>%
      mutate(!!col_evt_classes:=as.character(.data[[col_evt_classes]])) %>%
      group_by(.data[[col_gby]], Subject_Id, .data[[col_evt]], .data[[col_evt_classes]]) %>%
      summarize(Count=n(), .groups="keep")

    df <- df  %>%
      group_by(.data[[col_gby]], Subject_Id, .data[[col_evt]]) %>%
      mutate(Count_Subject=sum(Count)) %>%
      ungroup() %>%
      mutate(!!col_evt_classes:=ifelse(Count_Subject > 1, "Multihit", .data[[col_evt_classes]]), Count=1) %>%
      select(-c(Count_Subject)) %>%
      distinct() %>%
      group_by(.data[[col_gby]], .data[[col_evt]], .data[[col_evt_classes]]) %>%
      summarize(Count=sum(Count), .groups="keep")
  } else {
    df <- df_alt_filt %>%
      group_by(.data[[col_gby]], Subject_Id, .data[[col_evt]]) %>%
      summarize(Count=n(), .groups="keep")
  }

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


get_tables_for_plot <- function(dfs_alt, df_tt_count, cohorts_keep, evt_keep, tt_keep, col_evt, col_tt, col_stack,
                                cohorts_a=c("prism", "met500"), cohorts_b=c("tcga", "tcga"), min_counts_evt=NULL,
                                max_evt=NULL, n_cores=4){
  dfs_plot <- list()
  
  # event classes for right bar plot
  if (!is.null(col_stack)){
    evt_classes <- get_alterations_classifications()
  } else {
    evt_classes <- NULL
  }

  # For each cohort, compute tables of event counts/percentages aggregated by Tumor_Type
  for (cohort in cohorts_keep){
    df_alt <- dfs_alt[[cohort]]
    df_alt <- df_alt %>% filter(.data[[col_evt]] %in% evt_keep)

    # fill Variant classifications for Amplification/Deletion
    if ("Alteration_Category" %in% colnames(df_alt)){
      df_alt[df_alt$Alteration_Category=="Amplification", "Variant_Classification"] <- "Amplification"
      df_alt[df_alt$Alteration_Category=="Deletion", "Variant_Classification"] <- "Deletion"
    }

    df_alt_count <- get_alt_count(df_alt=df_alt, col_evt=col_evt, col_evt_classes=col_stack, evt_classes=evt_classes)
    df_tt_count_cohort <- df_tt_count %>% rename(Count=.data[[paste0(toupper(cohort), "_", args$selection_name)]])
    dfs_tabs <- get_tables_for_facetted_heatmap_barplots(df_alt_count, df_tt_count_cohort, tt_keep, col_evt, 
                                                         col_tt, col_stack)
    if ("percent" %in% names(dfs_tabs)){
      dfs_tabs$percent <- dfs_tabs$percent  %>% replace(is.na(.), 0)
    }
    dfs_plot[[cohort]] <- dfs_tabs
  }

  # Select list of genes for the plot. This step is important as it conditions how many tests will be performed
  in_plot_evt <- select_in_plot_evt(dfs_plot, min_counts_evt, col_evt, max_evt=max_evt, max_evt_cohort="prism")
  dfs_plot <- lapply(dfs_plot, function(dfs) select_vals_in_dfs(dfs, col_evt, in_plot_evt))

  # differential testing per tumor type
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b, test="fisher", col_evt=col_evt,
                               in_plot_evt=in_plot_evt, tt_keep=tt_keep, n_cores=n_cores, suffix="")

  # differential testing across tumor types 
  dfs_plot <- add_pvals_tables(dfs_plot, cohorts_a=cohorts_a, cohorts_b=cohorts_b, test="cmh", col_evt=col_evt,
                               in_plot_evt=in_plot_evt, tt_keep=tt_keep, n_cores=n_cores, suffix="_strats")

  dfs_plot
}
 

save_tables_xlsx <- function(dfs_plot, col_tt, filepath){
  df_desc <- data.frame(Supplementary=c("Tables for facetted heatmap"))
  write.xlsx(df_desc, file=filepath, sheetName="description", row.names=FALSE)

  for (name_a in names(dfs_plot)){
    for (name_b in names(dfs_plot[[name_a]])){
      sheet_name <- paste(name_a, name_b, sep="_")
      df_name <- dfs_plot[[name_a]][[name_b]]
      
      if(!is.null(df_name)){
        if (name_b=="count_col") df_name <- df_name %>% rownames_to_column(var=col_tt)
        write.xlsx(as.data.frame(df_name), file=filepath, sheetName=sheet_name, row.names=FALSE, append=T)
      }
    }
  }

  cat(paste("-INFO: excel workbook saved at", filepath), "\n")
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


main <- function(args){
  # load alterations and tt counts
  dfs_alt <- load_alterations_and_select_samples(args)
  df_tt_count <- load_table(args$counts)

  # define util variable names
  col_tt <- "Tumor_Type"

  if (args$sample_type=="RNA_T"){
    if (args$gene_type=="fusions"){
      col_evt <- "Fusion_Id"
      col_stack <- "Annotations_Custom_Simple"
    } else if (args$gene_type=="genes") {
      col_evt <- "Hugo_Symbol"
      col_stack <- "Annotations_Custom_Simple"
    } else if (args$gene_type=="tmes"){
      col_evt <- "Label"
      col_stack <- NULL
    }
  } else {
    col_evt <- "Hugo_Symbol"
    col_stack <- "Variant_Classification"
  }

  # read input min counts
  min_counts_evt <- check_min_counts(args$min_counts_evt, args$cohorts)

  # select data for columns
  cohorts_keep <- c("prism", "met500", "tcga")

  if (args$sample_type=="DNA_T__DNA_N"){
    evt_keep <- load_oncogenes_or_tumorsuppressor(args$gene_type)
  } else {
    evt_keep <- Reduce(union, lapply(dfs_alt, function(df) unique(df[[col_evt]])))
  }
  df_tt_count <- df_tt_count %>% filter(.data[[paste0("Use_", args$selection_name)]]==1)
  tt_keep <- df_tt_count %>% pull(.data[[col_tt]])

  dfs_plot <- get_tables_for_plot(dfs_alt=dfs_alt, df_tt_count=df_tt_count, cohorts_keep=cohorts_keep,
                                  evt_keep=evt_keep, tt_keep=tt_keep, col_evt=col_evt, col_tt=col_tt,
                                  col_stack=col_stack, cohorts_a=c("prism", "met500"), cohorts_b=c("tcga", "tcga"), 
                                  min_counts_evt=min_counts_evt, max_evt=args$max_evt, n_cores=args$n_cores)


  # save all tables in an excel workbook
  save_tables_xlsx(dfs_plot=dfs_plot, col_tt=col_tt, filepath=args$output_tables)

  # aggregate tables across cohorts to prepare for the plot
  cohorts_keep_order <- c("prism", "met500", "tcga")
  tt_keep_order <- df_tt_count %>% arrange(desc(.data[[paste0("PRISM_", args$selection_name)]])) %>%
    pull(.data[[col_tt]])
  dfs_plot_agg <- dfs_plot$prism

  for (name_comp in names(dfs_plot)[grepl("_vs_", names(dfs_plot))]){
    for (subname in names(dfs_plot[[name_comp]])){
      new_name <- paste(name_comp, subname, sep="_")
      dfs_plot_agg[[new_name]] <- dfs_plot[[name_comp]][[subname]]
    }
  }

  for (name in names(dfs_plot_agg)){
    df_for_plot <- dfs_plot_agg[[name]]
    colnames_df <- colnames(df_for_plot)
    if (length(intersect(colnames_df, tt_keep_order))>0){
      colnames_df_other <- setdiff(colnames_df, tt_keep_order)
      df_for_plot <- df_for_plot[c(colnames_df_other, tt_keep_order)]
    }
    dfs_plot_agg[[name]] <- df_for_plot
  }

  # order rows for the heatmap
  dfs_plot_agg <- order_dfs(dfs_plot_agg, col_var=col_evt)

  # define arguments for the function draw_facetted_heatmap_barplots_2
  cohorts2colors <- load_colors(sheet="Global")[c("prism", "met500")]
  tumors2colors <- load_colors(sheet="Project_TCGA_More")

  if (args$sample_type=="DNA_T__DNA_N" | args$sample_type=="DNA_N"){
    stacks2colors <- load_colors(sheet="Variant_Classification_Custom")
  } else if (args$sample_type=="RNA_T" & (args$gene_type=="genes" | args$gene_type=="fusions")){
    stacks2colors <- load_colors(sheet="Fusion_Annotations_Custom")
  } else {
    stacks2colors <- NULL
  }

  if (args$sample_type=="DNA_T__DNA_N"){
    if (args$gene_type=="oncogenes"){
      x_tick_heatmap_side <- "top"
      colors_palette_heatmap <- "Reds"
    } else {
      x_tick_heatmap_side <- "bottom"
      colors_palette_heatmap <- "Blues"
    }
  } else if (args$sample_type=="DNA_N") {
    x_tick_heatmap_side <- "bottom"
    colors_palette_heatmap <- "Greens"
  } else if (args$sample_type=="RNA_T"){
    if (args$gene_type=="fusions"){
      colors_palette_heatmap <- "Reds"
      x_tick_heatmap_side <- "top"
    } else if (args$gene_type=="genes") {
      colors_palette_heatmap <- "Reds"
      x_tick_heatmap_side <- "bottom"
    } else {
      x_tick_heatmap_side <- "top"
      colors_palette_heatmap <- "Blues"
    }
  }

  names2colors <- list(met500_vs_tcga_qvals_strats=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals_strats=cohorts2colors[["prism"]],
                       met500_vs_tcga_qvals=cohorts2colors[["met500"]],
                       prism_vs_tcga_qvals=cohorts2colors[["prism"]],
                       colnames=tumors2colors)

  names2plots <- list(heatmap="percent",
                      heatmap_hover="count",
                      heatmap_pvals=c("prism_vs_tcga_qvals", "met500_vs_tcga_qvals"),
                      left_pvals=c("prism_vs_tcga_qvals_strats", "met500_vs_tcga_qvals_strats"),
                      middle_right_bar="count_row",
                      count_col="count_col")
  if (!is.null(col_stack)) names2plots[["extreme_right_bar"]] <- "fraction_row"


  if (args$sample_type=="RNA_T" & args$gene_type %in% c("fusions", "genes")){
    fonts <- list(x_tick_heatmap=list(size=12, color="black", family="Helvetica"),
                  y_tick_heatmap=list(size=10, color="black", family="Helvetica"),
                  x_tick_row_bar=list(size=14, color="black", family="Helvetica"),
                  z_heatmap=list(size=14, color="black", family="Helvetica"),
                  cohort_size=list(size=14, color="black", family="Helvetica"),
                  legend=list(size=12, color="black", family="Helvetica"))
  } else if (args$sample_type=="RNA_T" & args$gene_type %in% c("tmes")){
    fonts <- list(x_tick_heatmap=list(size=10, color="black", family="Helvetica"),
                  y_tick_heatmap=list(size=8, color="black", family="Helvetica"),
                  x_tick_row_bar=list(size=8, color="black", family="Helvetica"),
                  z_heatmap=list(size=10, color="black", family="Helvetica"),
                  cohort_size=list(size=10, color="black", family="Helvetica"),
                  legend=list(size=10, color="black", family="Helvetica"))
  } else {
    fonts <- list(x_tick_heatmap=list(size=12, color="black", family="Helvetica"),
                  y_tick_heatmap=list(size=9, color="black", family="Helvetica"),
                  x_tick_row_bar=list(size=9, color="black", family="Helvetica"),
                  z_heatmap=list(size=10, color="black", family="Helvetica"),
                  cohort_size=list(size=12, color="black", family="Helvetica"),
                  legend=list(size=10, color="black", family="Helvetica"))
  }

  # draw figure
  figs <- draw_facetted_heatmap_barplots_1(dfs=dfs_plot_agg, 
                                           col_var=col_evt,
                                           names2plots=names2plots,
                                           names2colors=names2colors,
                                           colors_palette_heatmap=colors_palette_heatmap,
                                           alpha_left=0.1, 
                                           alpha_heatmap=0.1,
                                           width_edges_heatmap=0.3,
                                           col_stack=col_stack,
                                           width_one=args$output_plot_width_one,
                                           height_one=args$output_plot_height_one,
                                           fonts=fonts,
                                           stacks2colors=stacks2colors,
                                           add_colors_names=F,
                                           x_tick_heatmap_side=x_tick_heatmap_side,
                                           add_cohorts_sizes=T)

  # assemble figures
  if (args$sample_type=="DNA_T__DNA_N" | args$sample_type=="DNA_N"){
    margin <- 0.09
  } else if (args$sample_type=="RNA_T" & args$gene_type=="tmes"){
    margin <- 0.04
  } else {
    margin <- 0.11
  }
  
  if (args$sample_type=="RNA_T" & args$gene_type=="tmes"){
    fig_left <- figs[["heatmap"]]
    fig_right <- figs[["middle_right_bar"]] 
  } else {
    fig_left <- figs[["heatmap"]]
    fig_right <- subplot(figs[["middle_right_bar"]], figs[["extreme_right_bar"]], margin=0.01, widths=c(0.5,0.5))
  }

  if (args$sample_type=="DNA_T__DNA_N" | args$sample_type=="DNA_N"){
    fig <- subplot(fig_left, fig_right, margin=margin, widths=c(0.65, 0.35))
    fig <- fig %>% layout(plot_bgcolor='rgba(0,0,0,0)')
  } else if (args$sample_type=="RNA_T" & args$gene_type=="tmes") {
    fig <- subplot(fig_left, fig_right, margin=margin, widths=c(0.85, 0.15))
    fig <- fig %>% layout(plot_bgcolor='rgba(0,0,0,0)')
  } else {
    fig <- subplot(fig_left, fig_right, margin=margin, widths=c(0.75, 0.25))
    fig <- fig %>% layout(plot_bgcolor='rgba(0,0,0,0)')
  }

  # save
  width <- args$output_plot_width_one*ncol(dfs_plot_agg$percent)
  height <- args$output_plot_height_one*nrow(dfs_plot_agg$percent)
  cat(paste("figure width:", width, "and height:", height,"\n"))

  if (args$sample_type=="DNA_N"){
    height <- height + 60
  } else if (args$sample_type=="DNA_T__DNA_N") {
    if (args$gene_type=="tumorsuppressors"){
      height <- height + 40
    }
  } else if (args$sample_type=="RNA_T"){
    if (args$gene_type=="fusions" | args$gene_type=="genes"){
      height <- height + 200
    } else {
      height <- height + 120
    }
  }

  save_image(fig, args$output_plot, width=width, height=height)
  cat(paste("-plot saved at", args$output_plot, "\n"))

  if (!is.null(args$output_plot_paper)){
    save_image(fig, args$output_plot_paper, width=width, height=height)
    cat(paste("-plot saved at", args$output_plot_paper, "\n"))
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw heatmap of alterations.')
  parser$add_argument("--cohorts", nargs="+", help="Names of input cohorts.",
                      default=c("prism", "met500", "tcga"))
  parser$add_argument("--samples", nargs="+", help="Paths to input sample tables.",
                      default=c("../../../results/combined_alterations/selection/selection_samples_prism.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_met500.tsv",
                                "../../../results/combined_alterations/selection/selection_samples_tcga.tsv"))
  parser$add_argument("--alterations", nargs="+", help="Paths to input alterations tables.",
                default=c("../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv",
                          "../../../results/combined_alterations/alterations/aggregated_alterations_met500.tsv",
                          "../../../results/combined_alterations/alterations/aggregated_alterations_tcga.tsv"))
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/combined_alterations/selection/selection_tumor_types.tsv")
  parser$add_argument("--min_counts_evt", type="integer", nargs="+",
                      default=c(10,7,90),
                      help="Criteria for building the list of genes in the plot.")
  parser$add_argument("--max_evt", type="integer", default=NULL,
                      help="If not NULL, this option sets the maximum number of rows displayed in heatmap.")
  parser$add_argument("--gene_type", type="character", help="Type of events considered in the plot.",
                      default="oncogenes")
  parser$add_argument("--sample_type", type="character", help="Choose DNA_T, RNA_T or DNA_N.",
                      default="DNA_T__DNA_N")
  parser$add_argument("--selection_name", type="character", help="Name of the selection in selection tables.",
                      default="heatmap_dna")
  parser$add_argument("--n_cores", type="integer", default=4, help="Number of cores for parallel computations.")
  parser$add_argument("--output_tables", type="character",
                      help="Path to where output tables are saved.",
                      default="../../../results/combined_alterations/heatmap_dna/tables_oncogenes.xlsx")
  parser$add_argument("--output_plot_width_one", type="integer", default=45, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot_height_one", type="integer", default=16, help="Width of output plot in pixels.")
  parser$add_argument("--output_plot", type="character",
                      help="Path to where output facetted heatmap plot is saved.",
                      default="../../../results/combined_alterations/heatmap_dna/heatmap_oncogenes.svg")
  parser$add_argument("--output_plot_paper", type="character",
                      help="Path to where output facetted heatmap plot is saved.",
                      default=NULL)
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
