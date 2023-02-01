# @created: 18 Oct 21
# @modified: 20 Oct 21
# @authors: Yoann Pradat
#
# Describe the distribution of metastatic sites according to tumor types.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# functions ============================================================================================================

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

  fig <- fig %>% add_pie(labels=df_count$Label, values=df_count$Count, textposition="outside",
                         textinfo=textinfo,
                         text=df_count$Percentage,
                         hovertemplate = paste0(paste0("<b>%{label}</b>"),
                                                '<br><br>Count: %{value}<br>',
                                                'Percent: %{text:.3f}',
                                                "<extra></extra>"),
                         hole=0.4,
                         marker = list(line=list(color='#FFFFFF', width=1), colors=df_count$Color),
                         domain = list(row=row, column=column))

  fig
}

add_drug_classes <- function(df, df_drug, col_drug_name, col_drug_class){
  df_drug <- df_drug %>%
    rename(!!col_drug_class:=Class, !!col_drug_name:=Drug)

  # add column col_drug_class to df
  cols_id <- colnames(df)[!colnames(df) %in% col_drug_name]
  df_sep <- df %>% separate_rows(.data[[col_drug_name]], sep="\\|")
  df_sep <- left_join(df_sep, df_drug[,c(col_drug_name, col_drug_class)], by=col_drug_name)

  df_sep %>% group_by_at(all_of(cols_id)) %>%
    summarize(!!col_drug_class:=paste0(sort(unique(.data[[col_drug_class]])), collapse="|"),
              !!col_drug_name:=paste0(sort(unique(.data[[col_drug_name]])), collapse="|"), .groups="keep") %>% 
    replace(.=="", NA) %>%
    ungroup()
}


get_matching_drug_classes <- function(df_drug, drug_class){
  drug_class_lvls <- unlist(strsplit(drug_class, split=" - "))
  n_lvls <- length(drug_class_lvls)
  if (n_lvls==1){
    drug_classes <- df_drug %>% filter(Class_Lvl_1==drug_class_lvls[[1]]) %>% pull(Class)
  } else {
    intermediate <- paste(drug_class_lvls[1:(n_lvls-1)], collapse=" - ")
    last_levels <- unlist(strsplit(drug_class_lvls[n_lvls], split="/"))
    df_drug_sub <- df_drug %>% filter(grepl(paste0("^", intermediate), Class))
    drug_classes <- c()
    for (last_level in last_levels){
      regex <- paste0(" - ", last_level, "|/", last_level)
      drug_classes <- c(drug_classes, df_drug_sub %>% filter(grepl(regex, Class)) %>% pull(Class) %>% unique())
    }
  }

  unique(drug_classes)
}


main <- function(args){
  # load tables
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  df_alts <- load_table(args$alts_table) %>% filter(Cohort==args$cohort)

  df_drug <- load_table(args$drug_table)
  cols_levels <- sort(colnames(df_drug)[grepl("Class_Lvl", colnames(df_drug))])
  df_drug <- df_drug %>% unite("Class", all_of(cols_levels), sep=" - ", remove=F)
  df_drug <- df_drug %>% mutate(Class=gsub(" - NA", "", Class)) %>%
    mutate(Class=ifelse(Class %in% c("", "NA"), NA, Class))

  # add Class column
  for (mode in c("Sen", "Res")){
    col_drug_name <- paste(mode, "Drug_Simple", sep="_")
    col_drug_class <- paste(mode, "Class_Simple", sep="_")
    df_alts <- add_drug_classes(df_alts, df_drug, col_drug_name, col_drug_class)
  }
  df_cln <- add_drug_classes(df_cln, df_drug, "Drugs_Before_Biopsy", "Classes_Before_Biopsy")

  # select tumor types
  if (tolower(args$tumor_types)!="all"){
    df_cln <- df_cln %>% filter(Tumor_Type %in% toupper(args$tumor_types))
  }

  # filter out patients with no treatments data
  df_cln <- df_cln %>% filter(!is.na(Systematic_Treatment_Before_Biopsy))

  # select alterations
  samples_id <- union(df_cln$Sample_Id_DNA_T, df_cln$Sample_Id_RNA_T)
  df_alts <- df_alts %>% filter(Sample_Id %in% samples_id)

  # plot data
  subjects_id <- df_cln %>% filter(Sample_Type=="DNA_N|DNA_T|RNA_T") %>%
    pull(Subject_Id)

  drug_classes <- c("Hormonotherapy", "Platinum_salts", "Targeted_Therapy - BRAF", "Targeted_Therapy - EGFR")
  mains <- paste("Resistance events to", drug_classes)
  mains <- gsub("_", " ", mains)

  for (i in 1:length(drug_classes)){
    drug_class <- drug_classes[[i]]
    main <- mains[[i]]
    output <- args$outputs[[i]]

    drug_class_match <- get_matching_drug_classes(df_drug, drug_class)
    regex <- paste(drug_class_match, collapse="|")
    df_sub <- df_alts %>% filter(grepl(regex, Res_Class_Simple)) %>%
      filter(Subject_Id %in% subjects_id)
    df_sub <- df_sub %>% unite(Mutation_Detail, c(Hugo_Symbol, Alteration), sep=" - ") %>%
      mutate(Tumor_Type = paste0("(", Tumor_Type, ")")) %>%
      unite(Mutation_Tumor, c(Mutation_Detail, Tumor_Type), sep=" ")

    # plot settings
    fontsize <- 10
    font <- list(size=fontsize)
    font_anns <- list(size=fontsize*1.25, family="Arial")
    font_title <- list(size=fontsize*2, family="Arial")

    # draw
    field <- "Mutation_Tumor"
    ann <- paste(length(unique(df_sub$Sample_Id)), "samples")
    annotation <- list(text=ann, showarrow=F, font=font_anns, xanchor="center", xref="x", yref="y")
    emptyaxis <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)

    fig <- plot_ly()
    fig <- pie_count_plotly(fig, df=df_sub, field=field, values2colors=NULL, row=0, column=0, pal="Dark2")
    fig <- fig %>% layout(title = list(text=main, font=font_title, xref="paper", yref="container", x=0.5, y=0.99),
                          font = font, showlegend = F, annotations=list(annotation),
                          grid = list(rows=1, columns=1),  xaxis = emptyaxis, yaxis = emptyaxis)

    fig <- fig %>% layout(margin=list(l=50, b=75, top=75, right=50))
    orca(fig, output, width=520, height=350, verbose=F)
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe mutational landscape of APOBEC-rich and APOBEC-poor tumors.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument("--drug_table", type="character",
                      default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx",
                      help="Path to table linking drug names to drug classes. Required if --class_names is used.")
  parser$add_argument("--tumor_types", type="character", nargs="+", default="All", 
                      help="List of tumor types to be considered.")
  parser$add_argument('--alts_table', type="character", help='Path to the table of annotated genomic alterations.',
                      default="../../../results/oncokb_civic/alterations/aggregated_alterations.tsv")
  parser$add_argument("--outputs", type="character", help="Path to output plot.",
    default=c("../../../results/data_overview/description/resistance/hormonotherapy_resistance_mutation.pdf",
              "../../../results/data_overview/description/resistance/platinum_salts_resistance_mutation.pdf",
              "../../../results/data_overview/description/resistance/targeted_egfr_resistance_mutation.pdf",
              "../../../results/data_overview/description/resistance/targeted_braf_resistance_mutation.pdf"))
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
