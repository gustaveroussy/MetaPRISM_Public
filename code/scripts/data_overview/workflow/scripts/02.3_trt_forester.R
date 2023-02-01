# @created: 18 Oct 21
# @modified: 17 Mar 22
# @authors: Yoann Pradat
#
# Describe the distribution of metastatic sites according to tumor types.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DescTools))
suppressPackageStartupMessages(library(extrafont))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(Exact))
suppressPackageStartupMessages(library(forester))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))

# functions ============================================================================================================

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


add_statistics <- function(table, subgroups, c_as, m_a, c_bs, m_b, estimate=c("odds_ratio", "diff_prop"), conf_int=F,
                           conf_lvl=0.95, progress=F, n_cores=1){
  stopifnot(length(subgroups)==length(c_as))
  stopifnot(length(subgroups)==length(c_bs))
  estimate <- match.arg(estimate)

  df_count_a <- data.frame(Count=c_as, row.names=subgroups)
  df_count_b <- data.frame(Count=c_bs, row.names=subgroups)
  margins_a <- t(data.frame(Count=m_a))
  margins_b <- t(data.frame(Count=m_b))

  if (estimate=="diff_prop"){
    # estimate difference in proportion of the subgroup in populations A and B using a Fischer-Booschloo test
    statistics <- get_pvals_fisher(df_count_a, df_count_b, margins_a, margins_b, conf_int=conf_int, conf_lvl=conf_lvl,
                                   progress=progress, n_cores=n_cores)

    for (subgroup in subgroups){
      c_a <- df_count_a[subgroup, Count]
      c_b <- df_count_b[subgroup, Count]
      est <- statistics$est[subgroup,"Count"]
      pval <- statistics$pval[subgroup,"Count"]

      if (conf_int){
        cil <- statistics$cil[subgroup,"Count"]
        cih <- statistics$cih[subgroup,"Count"]
        table <- bind_rows(table, tibble(Subgroup=subgroup, Present_A=c_a, Present_B=c_b,
                                         P_val=pval, Estimate=est, CI_Low=cil, CI_High=cih))
      } else {
        table <- bind_rows(table, tibble(Subgroup=subgroup, Present_A=c_a, Present_B=c_b,
                                         P_val=pval, Estimate=est))
      }
    }

  } else {
    for (subgroup in subgroups){
      c_a <- df_count_a[subgroup, "Count"]
      c_b <- df_count_b[subgroup, "Count"]
      data <- matrix(c(c_a, m_a - c_a, c_b, m_b - c_b), byrow=T, nrow=2)

      if (conf_int){
        out <- OddsRatio(x=data, digits=3, method="wald", conf.level=conf_lvl)
        est <- out[["odds ratio"]]
        cil <- out[["lwr.ci"]]
        cih <- out[["upr.ci"]]
        table <- bind_rows(table, tibble(Subgroup=subgroup, Present_A=c_a, Present_B=c_b,
                                         Estimate=est, CI_Low=cil, CI_High=cih))
      } else {
        out <- OddsRatio(x=data, digits=3, method="wald", conf.level=NULL)
        est <- out[["odds ratio"]]
        table <- bind_rows(table, tibble(Subgroup=subgroup, Present_A=c_a, Present_B=c_b, Estimate=est))
      }
    }
  }

  table
}

add_values_with_zero_count <- function(df_count, col, vals){
  groups_exp <- c("A", "B")
  for (val in vals){
    groups_val <- df_count %>% filter(.data[[col]]==val) %>% pull(Group)
    groups_mis <- setdiff(groups_exp, groups_val)
    for (group in groups_mis){
      df_count <- bind_rows(df_count, tibble(Group=group, Present=0, Col=val) %>% rename(!!col:=Col))
    }
  }
  df_count
}

add_statistics_mutational_signatures <- function(table, df_cln, args, show_main_count=T){
  df_sigs <- load_table(args$sigs_table)
  samples <- intersect(colnames(df_sigs[2:ncol(df_sigs)]), df_cln$Sample_Id_DNA_T)
  df_sigs <- df_sigs[,c("Signature", samples)]
  df_sigs <- df_sigs %>% pivot_longer(!Signature, names_to="Sample_Id_DNA_T", values_to="Count")
  cols_cln <- c("Sample_Id_DNA_T", "Group", "Subject_Id")
  df_sigs <- left_join(df_sigs, df_cln[,cols_cln], by="Sample_Id_DNA_T")

  df_u <- df_sigs %>% distinct(Subject_Id, Group)
  m_a <- sum(df_u$Group=="A")
  m_b <- sum(df_u$Group=="B")
  if (show_main_count){
    table <- bind_rows(table, tibble(Subgroup="Mutational signatures", A=m_a, B=m_b))
  } else {
    table <- bind_rows(table, tibble(Subgroup="Mutational signatures", A=NA, B=NA))
  }

  col_sel <- "Signature"
  vals_sel <- args$sigs_names

  df_count <- df_sigs %>% filter(.data[[col_sel]] %in% vals_sel) %>% 
    group_by(Group, .data[[col_sel]]) %>% summarize(Present=sum(Count>0), .groups="keep")
  df_count <- add_values_with_zero_count(df_count, col_sel, vals_sel)
  df_count <- df_count %>% arrange(.data[[col_sel]])
  c_as <- df_count %>% filter(Group=="A") %>% pull(Present)
  c_bs <- df_count %>% filter(Group=="B") %>% pull(Present)
  add_statistics(table, sort(vals_sel), c_as, m_a, c_bs, m_b, conf_int=T, conf_lvl=0.95)
}


add_statistics_clinical <- function(table, subgroup, df_cln, col_sel, min_count, show_main_count=T){
  df_cln <- df_cln %>% filter(!is.na(.data[[col_sel]])) %>% mutate(!!col_sel:=as.character(.data[[col_sel]]))
  vals_sel <- df_cln %>% group_by(.data[[col_sel]]) %>% summarize(n=n()) %>% filter(n>=min_count) %>%
    pull(.data[[col_sel]])

  df_u <- df_cln %>% distinct(Subject_Id, Group)
  m_a <- sum(df_u$Group=="A")
  m_b <- sum(df_u$Group=="B")
  if (show_main_count){
    table <- bind_rows(table, tibble(Subgroup=subgroup, A=m_a, B=m_b))
  } else {
    table <- bind_rows(table, tibble(Subgroup=subgroup, A=NA, B=NA))
  }

  df_count <- df_cln %>% filter(.data[[col_sel]] %in% vals_sel) %>% 
    group_by(Group, .data[[col_sel]]) %>% summarize(Present=n(), .groups="keep")
  df_count <- add_values_with_zero_count(df_count, col_sel, vals_sel)
  df_count <- df_count %>% arrange(.data[[col_sel]])
  c_as <- df_count %>% filter(Group=="A") %>% pull(Present)
  c_bs <- df_count %>% filter(Group=="B") %>% pull(Present)
  add_statistics(table, sort(vals_sel), c_as, m_a, c_bs, m_b, conf_int=T, conf_lvl=0.95)
}


add_statistics_alterations <- function(table, subgroup, df_cln, df_alts, col_row, mode=c("Res", "Sen"),
                                       show_main_count=T, sample_types="DNA_N|DNA_T|RNA_T", drug_names=NULL,
                                       drug_classes=NULL, genes=NULL, tiers=NULL, alterations=NULL, pathways=NULL,
                                       df_pathways=NULL, min_count=5){
  mode <- match.arg(mode)
  col_drug_name <- paste(mode, "Drug_Simple", sep="_")
  col_drug_class <- paste(mode, "Class_Simple", sep="_")
  col_tier <- paste(mode, "Level_Simple", sep="_")

  # select samples
  df_cln <- df_cln %>% filter(Sample_Type %in% sample_types)
  samples_id <- union(df_cln$Sample_Id_DNA_T, df_cln$Sample_Id_RNA_T)
  df_alts <- df_alts %>% filter(Sample_Id %in% samples_id)

  # add Group column
  df_alts <- left_join(df_alts, df_cln[,c("Subject_Id", "Group")], by="Subject_Id")

  # count events per drug class
  if (!is.null(drug_classes)){
    col_sel <- col_drug_class
    if (drug_classes=="All"){
      vals_sel <- df_alts %>% separate_rows(.data[[col_drug_class]], sep="\\|") %>%
        filter(!is.na(.data[[col_drug_class]])) %>%
        group_by(.data[[col_drug_class]]) %>%
        summarize(Count=n()) %>%
        filter(!.data[[col_sel]] %in% c("Exclude", "TO BE FILLED"), Count>=min_count) %>%
        pull(.data[[col_sel]])
    } else {
      vals_sel <- drug_classes
    }

    df_alts <- df_alts %>% filter(grepl(paste(vals_sel, collapse="|"), .data[[col_sel]]))
    if (col_row==col_sel) vals_row <- vals_sel
  }

  if (!is.null(drug_names)){
    col_sel <- col_drug_name
    vals_sel <- drug_names
    df_alts <- df_alts %>% filter(grepl(paste(vals_sel, collapse="|"), .data[[col_sel]]))
    if (col_row==col_sel) vals_row <- vals_sel
  }

  if (!is.null(genes)){
    col_sel <- "Hugo_Symbol"
    vals_sel <- genes
    df_alts <- df_alts %>% filter(grepl(paste(vals_sel, collapse="|"), .data[[col_sel]]))
    if (col_row==col_sel) vals_row <- vals_sel
  }


  if (!is.null(tiers)){
    col_sel <- col_tier
    vals_sel <- tiers
    df_alts <- df_alts %>% filter(grepl(paste(vals_sel, collapse="|"), .data[[col_sel]]))
    if (col_row==col_sel) vals_row <- vals_sel
  }


  df_u <- df_cln %>% distinct(Subject_Id, Group)
  m_a <- sum(df_u$Group=="A")
  m_b <- sum(df_u$Group=="B")
  if (show_main_count){
    table <- bind_rows(table, tibble(Subgroup=subgroup, A=m_a, B=m_b))
  } else {
    table <- bind_rows(table, tibble(Subgroup=subgroup, A=NA, B=NA))
  }

  df_count <- tibble()
  for (val_row in vals_row){
    df_count_row <- df_alts %>% filter(grepl(val_row, .data[[col_row]])) %>% 
      group_by(Group) %>% summarize(Present=n(), .groups="keep")
    df_count_row[[col_row]] <- val_row
    df_count <- bind_rows(df_count, df_count_row)
  }
  df_count <- add_values_with_zero_count(df_count, col_row, vals_row)
  df_count <- df_count %>% arrange(.data[[col_row]])
  c_as <- df_count %>% filter(Group=="A") %>% pull(Present)
  c_bs <- df_count %>% filter(Group=="B") %>% pull(Present)
  add_statistics(table, sort(vals_row), c_as, m_a, c_bs, m_b, conf_int=T, conf_lvl=0.95)
}


draw_plot <- function(table, A_name, B_name, output){
  table <- table %>% mutate(Estimate=ifelse(Present_A==0&Present_B==0, NA, Estimate)) %>%
    mutate(CI_Low=ifelse(Present_A==0&Present_B==0, NA, CI_Low)) %>%
    mutate(CI_High=ifelse(Present_A==0&Present_B==0, NA, CI_High))

  table$Subgroup <- gsub("_", " ", table$Subgroup)
  table$Subgroup <- ifelse(is.na(table$Present_A), 
                        table$Subgroup,
                        paste0("   ", table$Subgroup))
  table <- table %>% mutate(A=ifelse(is.na(Present_A), A, Present_A), B=ifelse(is.na(Present_B), B, Present_B))
  table <- table %>% rename(!!A_name:=A, !!B_name:=B)

  # draw forest plot
  xmin <- min(0.05, 0.5*min(table$Estimate, na.rm=T))
  xmax <- max(20, 2*max(table$Estimate, na.rm=T))

  forester(left_side_data = table[,1:3],
           estimate = table$Estimate,
           ci_low = table$CI_Low,
           ci_high = table$CI_High,
           display = FALSE,
           font_family = "sans",
           null_line_at=1,
           xlim = c(xmin,xmax),
           x_scale_linear=F,
           estimate_precision = 2,
           arrows=F,
           estimate_col_name = "Estimate (95% CI)",
           file_path = output)
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
  df_alts <- load_table(args$alts_table)

  if (!is.null(args$drug_table)){
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
  }

  if (!is.null(args$pathways_table)){
    df_pathways <- load_table(args$pathways_table)
  }

  # select tumor types
  if (tolower(args$tumor_types)!="all"){
    df_cln <- df_cln %>% filter(Tumor_Type %in% toupper(args$tumor_types))
  }

  # filter out patients with no treatments data
  df_cln <- df_cln %>% filter(!is.na(Systematic_Treatment_Before_Biopsy))

  # select alterations
  mask_dna_and_rna <- grepl("DNA_T\\|RNA_T", df_cln$Sample_Type)
  samples_id <- union(df_cln[mask_dna_and_rna,]$Sample_Id_DNA_T, df_cln[mask_dna_and_rna,]$Sample_Id_RNA_T)
  df_alts <- df_alts %>% filter(Sample_Id %in% samples_id)

  # if drug selection based on classes, add Classes_Before_Biopsy
  if (!is.null(args$drug_classes) & !("Classes_Before_Biopsy" %in% colnames(df_cln))){
    df_cln <- add_drug_classes(df_cln, df_drug, col_drug_name="Drugs_Before_Biopsy", col_drug_class="Classes_Before_Biopsy")
  }

  if (args$plot_mode==1){
    # Mode 1 ===========================================================================================================

    # cohort separation based on drug names or classes
    if (!is.null(args$drug_classes)){
      col_group <- "Classes_Before_Biopsy"
      vals_group <- unique(Reduce(c, sapply(args$drug_classes, function(d) get_matching_drug_classes(df_drug, d))))
    } else if (!is.null(args$drug_names)){
      col_group <- "Drugs_Before_Biopsy"
      vals_group <- args$drug_names
    } else {
        stop("One of --drug_names or --drug_classes must be not NULL")
    }

    # separate the cohort into 2 groups.
    # - A: those that received one ("union") or all ("intersection") of the drugs (names or classes)
    # - B: those that are not A

    if (args$combination=="union"){
      regex <- paste0(vals_group, collapse="|")
      df_cln <- df_cln %>% mutate(Group=ifelse(grepl(regex, .data[[col_group]]), "A", "B"))
    } else if (args$combination=="intersection") {
      mask <- Reduce("&", lapply(vals_group, function(c) grepl(c, df_cln[[col_group]])))
      df_cln <- df_cln %>% mutate(Group=ifelse(mask, "A", "B"))
    } else {
      stop("Choose 'union' or 'intersection' for --combination")
    }

    # Table of clinical characteristics
    table <- tibble(Subgroup=character(0), A=numeric(0), B=numeric(0), Present_A=numeric(0), Present_B=numeric(0))
    table <- bind_rows(table, tibble(Subgroup="All patients", A=sum(df_cln$Group=="A"), B=sum(df_cln$Group=="B")))
    if (args$tumor_types=="All" | length(args$tumor_types)!=1)
      table <- add_statistics_clinical(table, "Tumor type", df_cln, "Tumor_Type", min_count=10, show_main_count=F)
    if (!is.null(args$sigs_names))
      table <- add_statistics_mutational_signatures(table, df_cln, args, show_main_count=T)
    draw_plot(table, args$A_name, args$B_name, args$output[[1]])

    # Table of molecular characteristics
    table <- tibble(Subgroup=character(0), A=numeric(0), B=numeric(0), Present_A=numeric(0), Present_B=numeric(0))
    table <- bind_rows(table, tibble(Subgroup="All patients", A=sum(df_cln$Group=="A"), B=sum(df_cln$Group=="B")))
    table <- add_statistics_alterations(table, "Resistance alterations", df_cln, df_alts, col_row="Res_Class_Simple",
                                        mode="Res", drug_classes="All", min_count=5)
    draw_plot(table, args$A_name, args$B_name, args$output[[2]])

  } else if (args$plot_mode==2){
    # Mode 2 ===========================================================================================================

    table <- tibble(Subgroup=character(0), A=numeric(0), B=numeric(0), Present_A=numeric(0), Present_B=numeric(0))

    if (!is.null(args$drug_classes)){
      col_group <- "Classes_Before_Biopsy"
      tiers <- c("Tier1", "Tier2", "Tier3")

      for (drug_class in args$drug_classes){
        vals_group <- get_matching_drug_classes(df_drug, drug_class)
        regex <- paste0(vals_group, collapse="|")
        df_cln <- df_cln %>% mutate(Group=ifelse(grepl(regex, .data[[col_group]]), "A", "B"))
        table <- add_statistics_alterations(table, drug_class, df_cln, df_alts, col_row="Res_Class_Simple", mode="Res", 
                                            drug_classes=regex, min_count=1)
        table$Subgroup[nrow(table)] <- "Resistance"

        # genes <- df_alts %>% filter(grepl(regex, Res_Class_Simple)) %>%  pull(Hugo_Symbol) %>% unique()
        # table <- add_statistics_alterations(table, drug_class, df_cln, df_alts, mode="Res", col_row="Hugo_Symbol",
        #                                     drug_classes=regex, genes=genes, min_count=1)
        table <- add_statistics_alterations(table, drug_class, df_cln, df_alts, mode="Res", col_row="Res_Level_Simple",
                                            drug_classes=regex, tiers=tiers, min_count=1)
        table <- table %>% distinct()
      }

      table$Subgroup <- ifelse(grepl(paste(tiers, collapse="|"), table$Subgroup),
                               paste0("   ", table$Subgroup),
                               table$Subgroup)
      old2new <- list("Immune checkpoint inhibitor - PD1/PD-L1"="Immunotherapy - Checkpoint_inhibitor - PD1/PD-L1")
      table <- table %>% mutate(Subgroup=fct_recode(Subgroup, !!!old2new))
      draw_plot(table, args$A_name, args$B_name, args$output[[1]])

    } else if (!is.null(args$drug_names)){
      col_group <- "Drugs_Before_Biopsy"
      vals_group <- args$drug_names
    } else {
        stop("One of --drug_names or --drug_classes must be not NULL")
    }
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Describe odds ratios for different events in 2 populations A and B.')
  parser$add_argument("--cohort", type="character", default="prism", help="Cohort name.")
  parser$add_argument("--plot_mode", type="integer", default=1,
                      help=paste("Choose 1 for fixed A and B groups for odds ratios. Choose 2 for changing A and B",
                                "populations for each subset of odds ratio."))
  parser$add_argument("--combination", type="character", default="union", help="Choose 'union' or 'intersection'.")
  parser$add_argument("--drug_classes", type="character", nargs="*",
                      default=c("Platinum_salts"),
                      help=paste("List of drug classes to be considered for separating the patients into 2 groups: those",
                                 "that received none of these drugs and those that received 1 (if 'union') or all",
                                 "(if 'intersection') of them. You can use --drug_names instead if you want to be",
                                 "specific"))
  parser$add_argument("--drug_names", type="character", nargs="*", default=NULL, 
                      help=paste("List of drugs to be considered for separating the patients into 2 groups: those",
                                 "that received none of these drugs and those that received 1 (if 'union') or all",
                                 "(if 'intersection') of them. You can use --drug_classes instead if all drugs fall",
                                 "within the same therapeutic class."))
  parser$add_argument("--drug_table", type="character",
                      default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx",
                      help="Path to table linking drug names to drug classes. Required if --class_names is used.")
  parser$add_argument("--tumor_types", type="character", nargs="+", default="LUAD", 
                      help="List of tumor types to be considered.")
  parser$add_argument("--sigs_table", type="character", help="Path to table of mutational signatures.",
    default=paste0("../../../results/mutational_signatures/projection_known_signatures/",
      "MutationalPatterns/counts_signatures_cosmic_sbs_96_v3.2_sbs_96_min_mut_sparse_sigprofilerjulia_prism.tsv"))
  parser$add_argument("--sigs_names", type="character", nargs="*", help="Signatures to be considered.",
                      default=c("SBS31", "SBS35"))
  parser$add_argument('--alts_table', type="character", help='Path to the table of annotated genomic alterations.',
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv")
  parser$add_argument("--A_name", type="character", default="Platinum Tx", 
                      help="Name displayed on the graph for group A.")
  parser$add_argument("--B_name", type="character", default="No Platinum Tx", 
                      help="Name displayed on the graph for group B.")
  parser$add_argument("--outputs", type="character", nargs="+", help="Path to output plots.",
    default=c("../../../results/data_overview/treatments/forester/LUAD/platinum_salts_LUAD_clinical.pdf",
              "../../../results/data_overview/treatments/forester/LUAD/platinum_salts_LUAD_molecular.pdf"))
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
