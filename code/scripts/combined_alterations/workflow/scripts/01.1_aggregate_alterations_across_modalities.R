# @created: 17 Sep 21
# @modified: 14 Dec 22
# @authors: Yoann Pradat
#
# Aggregate alterations (copy-number, mutations and fusions) annotated with OncoKb and CiVIC annotations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

select_in_design_samples <- function(df, df_cln, cat="mut"){
  if (is.null(df)) return(NULL)

  if (cat %in% c("cna", "msi", "mut", "tmb")){
    df_cln <- df_cln %>% unite("Sample_Id", Sample_Id_DNA_T, Sample_Id_DNA_N, sep="_vs_", remove=F)
    mask <- df[["Sample_Id"]] %in% df_cln[["Sample_Id"]]
    n_pair_bef <- length(unique(df$Sample_Id))
    df <- df[mask,]
    n_pair_aft <- length(unique(df$Sample_Id))
    cat(paste("-INFO: selected", paste0(n_pair_aft,"/",n_pair_bef), cat, "lines from pairs in cln table\n"))
    df_cln <- df_cln %>% select(-Sample_Id)

  } else if (cat %in% c("exp", "fus")){
    mask <- df[["Sample_Id"]] %in% df_cln[["Sample_Id_RNA_T"]]
    n_sam_bef <- length(unique(df$Sample_Id))
    df <- df[mask,]
    n_sam_aft <- length(unique(df$Sample_Id))
    cat(paste("-INFO: selected", paste0(n_sam_aft,"/",n_sam_bef), cat, "lines from samples in cln table\n"))
  }

  df
}


process_cna <- function(df_cna){
  if (is.null(df_cna)) return(NULL)

  # process cnas
  df_cna <- df_cna %>% rename(ALTERATION_CATEGORY=Alteration, SAMPLE_ID=Sample_Id, HUGO_SYMBOL=Hugo_Symbol,
                              SUBJECT_ID=Subject_Id)

  df_cna
}


process_fus <- function(df_fus){
  if (is.null(df_fus)) return(NULL)

  # process fusions: disregard breakpoints data and keep only data about partners
  cols_fus_gen <- c("Sample_Id", "Subject_Id", "Fusion_Id", "Gene_1", "Gene_Id_1", "Gene_2", "Gene_Id_2",
                    "Gene_Type_1", "Gene_Type_2", "Annotations", "Annotations_Custom", "Algo_Wo_Breakpoint")
  cols_fus_okb <- intersect(colnames(df_fus), get_oncokb_columns())
  cols_fus_civ <- intersect(colnames(df_fus), get_civic_columns())
  cols_fus_keep <- c(cols_fus_gen, cols_fus_okb, cols_fus_civ)
  df_fus <- df_fus %>% select(all_of(cols_fus_keep)) %>% distinct()
  stopifnot(nrow(df_fus %>% select(Sample_Id, Fusion_Id) %>% distinct())==nrow(df_fus))
  df_fus <- df_fus %>% mutate(ALTERATION_CATEGORY="Fusion") %>%
    rename(SAMPLE_ID=Sample_Id, FUSION=Fusion_Id, SUBJECT_ID=Subject_Id)

  df_fus
}


process_msi <- function(df_msi){
  if (is.null(df_msi)) return(NULL)

  # process msi
  # - add OncoKB annotations manually https://www.oncokb.org/gene/Other%20Biomarkers/MSI-H
  coadread_okb_type <- c("COADREAD", "COAD", "CAIS", "MACR", "READ", "SRCCR")
  coadread_okb_drug <- "Pembrolizumab,Ipilimumab+Nivolumab,Nivolumab"
  solid_okb_drug <- "Pembrolizumab"

  df_msi <- df_msi %>% 
    mutate(LEVEL_1=ifelse(MSKCC_Oncotree %in% coadread_okb_type & Status=="Unstable", coadread_okb_drug,
                          ifelse(Status=="Unstable", solid_okb_drug, NA))) %>%
    mutate(ALTERATION_CATEGORY="MSI", SAMPLE_ID=Sample_Id, HUGO_SYMBOL="MSI High",
           ONCOGENIC="Oncogenic", GENE_IN_ONCOKB=TRUE, VARIANT_IN_ONCOKB=TRUE, MUTATION_EFFECT="Unknown",
           SUBJECT_ID=Subject_Id)
  df_msi <- df_msi %>% filter(!is.na(LEVEL_1))

  df_msi
}


process_mut <- function(df_mut){
  if (is.null(df_mut)) return(NULL)

  # process mutations:
  # - where Variant_Classification is Splice_Site or Splice_Region, set PROTEIN_CHANGE to Splice_Site
  # - where TP53, replace PROTEIN_CHANGE by DNA Binding Domain where applicable
  df_mut <- df_mut %>% rowwise %>% mutate(ALTERATION_CATEGORY=get_alteration_category_mut(Variant_Classification)) %>%
    rename(SAMPLE_ID=Sample_Id, SUBJECT_ID=Subject_Id, HUGO_SYMBOL=Hugo_Symbol, PROTEIN_CHANGE=HGVSp_Short)
  df_mut[["PROTEIN_CHANGE"]] <- gsub("%3D", "=", df_mut[["PROTEIN_CHANGE"]])
  mask_tp53 <- df_mut[["HUGO_SYMBOL"]]=="TP53"
  mask_dnabd <- (df_mut$Start_Position >= 7577149) & (df_mut$End_Position <- 7578443)
  df_mut[mask_tp53 & mask_dnabd, "PROTEIN_CHANGE_MORE"] <- "DNA Binding Domain"

  mask_splice <- grepl("splice", df_mut[["PROTEIN_CHANGE"]])
  mask_splice <- mask_splice | grepl("Splice|Silent", df_mut[["Variant_Classification"]])
  df_mut[mask_splice, "PROTEIN_CHANGE"] <- "Splice_Site"
  df_mut[mask_splice, "Variant_Classification"] <- "Splice_Site"

  df_mut
}


process_tmb <- function(df_tmb){
  if (is.null(df_tmb)) return(NULL)

  # process tmb
  # - add OncoKB annotations manually https://www.oncokb.org/gene/Other%20Biomarkers/TMB-H
  df_bed <- read.table(args$target_bed)
  df_bed <- df_bed %>% mutate(Size=V3-V2+1)
  target_size <- sum(df_bed$Size)/1e6 

  okb_drug <- "Pembrolizumab"

  df_tmb <- df_tmb %>% 
    mutate(LEVEL_1=ifelse(N_Mutations/target_size > 10, okb_drug, NA)) %>%
    mutate(ALTERATION_CATEGORY="TMB", SAMPLE_ID=Sample_Id, HUGO_SYMBOL="TMB High", SUBJECT_ID=Subject_Id,
           ONCOGENIC="Likely oncogenic", GENE_IN_ONCOKB=TRUE, VARIANT_IN_ONCOKB=TRUE, MUTATION_EFFECT="Unknown")
  df_tmb <- df_tmb %>% filter(!is.na(LEVEL_1))

  df_tmb
}


process_exp_arv7 <- function(df_exp_arv7){
  if (is.null(df_exp_arv7)) return(NULL)

  # process exp ARV7
  civ_subtype <- "Prostate Cancer"
  civ_drug <- "Enzalutamide,Abiraterone"

  df_exp_arv7 <- df_exp_arv7 %>% 
    mutate(`Predictive:N:B`=ifelse(grepl(civ_subtype, Civic_Disease) & ARV7ratio >= 0.05, civ_drug, NA)) %>%
    mutate(ALTERATION_CATEGORY="Exp", SAMPLE_ID=Sample_Id, HUGO_SYMBOL="AR", ALTERATION="AR-V7 High",
           ONCOGENIC="Oncogenic", GENE_IN_ONCOKB=TRUE, VARIANT_IN_ONCOKB=TRUE, MUTATION_EFFECT="Gain-of-function",
           SUBJECT_ID=Subject_Id)
  df_exp_arv7 <- df_exp_arv7 %>% filter(!is.na(`Predictive:N:B`))

  df_exp_arv7
}


correct_civic_annotations <- function(df_agg, df_cln){
  # correct some civic annotations that are unspecific
  # FGFR3 mutation is Tier1 in BLCA only for the Fusion and R248C, S249C, G370C, Y373C
  # Other alterations
  mask_blca <- grepl("Bladder Carcinoma", df_agg$CIViC_Matching_Disease)
  mask_fgfr3 <- df_agg$HUGO_SYMBOL=="FGFR3" 
  mask_fgfr3_tier1_mut <- df_agg$HUGO_SYMBOL=="FGFR3" & grepl("R248C|S249C|G370C|G372C|Y373C|Y375C", df_agg$PROTEIN_CHANGE)
  mask_fgfr3_tier1_fus <- grepl("FGFR3--TACC3", df_agg$FUSION)

  mask_fgfr3_tier1 <- mask_blca & mask_fgfr3 & (mask_fgfr3_tier1_mut | mask_fgfr3_tier1_fus) 
  mask_fgfr3_tier_other <- mask_blca & mask_fgfr3 & (!mask_fgfr3_tier1_mut & !mask_fgfr3_tier1_fus) 
  df_agg[mask_fgfr3_tier_other, "Predictive:P:A"] <- NA

  # KRAS Exon 2 mutation is Tier2 in COAD only for codon 2 or 13. Other Exon 2 mutation are Tier3 only 
  mask_coad <- grepl("Colorectal", df_agg$CIViC_Matching_Disease)
  mask_mut <- grepl("Mut", df_agg$ALTERATION_CATEGORY)
  mask_kras <- grepl("KRAS", df_agg$HUGO_SYMBOL)
  mask_exon2 <- grepl("2/", df_agg$EXON)
  mask_tier2_mut <- grepl("[A-Z*]12[A-Z*]|[A-Z*]13[A-Z*]", df_agg$PROTEIN_CHANGE)

  mask_kras_tier2 <- mask_coad & mask_mut & mask_kras & mask_exon2 & mask_tier2_mut 
  mask_kras_tier_other <- mask_coad & mask_mut & mask_kras & mask_exon2 & !mask_tier2_mut 
  df_agg[mask_kras_tier_other,"Predictive:P:B"] <- NA

  # BRAF V600E mutation in Lung Non-small Cell Carcinoma does not cause resistance to dabrafenib.
  # It is rather KRAS mutations. For sake of clarity, only patients for which dabrafenib is indicated, i.e patients
  # with a detected BRAF V600 mutation, will have an annotation of dabrafenib resistance.
  subjects_nsclc <- df_cln %>% 
    filter(MSKCC_Oncotree %in% c("LUAD", "LUSC", "CCLC", "SGTTL", "NSCLCPD", "LUACC", "LUAS")) %>%
    pull(var="Subject_Id")
  mask_nsclc <- df_agg$SUBJECT_ID %in% subjects_nsclc
  mask_kras <- grepl("KRAS", df_agg$HUGO_SYMBOL) & grepl("G12D", df_agg$PROTEIN_CHANGE)
  mask_braf <- grepl("BRAF", df_agg$HUGO_SYMBOL) & grepl("V600", df_agg$PROTEIN_CHANGE)

  mask_kras_add <- mask_nsclc & mask_kras
  mask_braf_rmv <- mask_nsclc & mask_braf
  if (sum(mask_kras_add)>0){
    df_agg[mask_kras_add,"Predictive:N:C"] <- unlist(lapply(df_agg[mask_kras_add,"Predictive:N:C", drop=T], function(x) {
                                        if (is.na(x)) {"Dabrafenib"} else {paste(x, "Dabrafenib", sep=";")}}))
  }
  if (sum(mask_braf_rmv)){
    df_agg[mask_braf_rmv,"Predictive:N:C"] <- gsub("Dabrafenib", "", df_agg[mask_braf_rmv,"Predictive:N:C", drop=T])
  }

  df_agg
}


get_oncokb_columns <- function(){
  c("GENE_IN_ONCOKB", "VARIANT_IN_ONCOKB", "MUTATION_EFFECT", "ONCOGENIC", "TX_CITATIONS",
    "LEVEL_1","LEVEL_2","LEVEL_3A","LEVEL_3B","LEVEL_4","LEVEL_R1","LEVEL_R2","LEVEL_R3", "HIGHEST_LEVEL",
    "LEVEL_Dx1","LEVEL_Dx2","LEVEL_Dx3", "HIGHEST_DX_LEVEL",
    "LEVEL_Px1","LEVEL_Px2","LEVEL_Px3", "HIGHEST_PX_LEVEL")
}


get_civic_columns <- function(){
  c("CIViC_Matching_Disease", "CIViC_Matching_Type", "CIViC_Matching_Gene_Variant", "CIViC_Matcing_Evidence_Id",
    "CIViC_Matching_Citation",
    "Predictive:N:A","Predictive:N:B","Predictive:N:C","Predictive:N:D","Predictive:N:E",
    "Predictive:P:A","Predictive:P:B","Predictive:P:C","Predictive:P:D","Predictive:P:E",
    "Diagnostic:N:A","Diagnostic:N:B","Diagnostic:N:C","Diagnostic:N:D","Diagnostic:N:E",
    "Diagnostic:P:A","Diagnostic:P:B","Diagnostic:P:C","Diagnostic:P:D","Diagnostic:P:E",
    "Prognostic:N:A","Prognostic:N:B","Prognostic:N:C","Prognostic:N:D","Prognostic:N:E",
    "Prognostic:P:A","Prognostic:P:B","Prognostic:P:C","Prognostic:P:D","Prognostic:P:E")
}


get_alteration_category_mut <- function(x){
  if (grepl("Ins", x)){
    return("Ins")
  } else if (grepl("Del", x)) {
    return("Del")
  } else {
    return("Mut")
  }
}

aggregate_alterations <- function(df_cna=NULL, df_fus=NULL, df_msi=NULL, df_mut=NULL, df_tmb=NULL, df_exp_arv7=NULL){
  cat("-aggregating alterations across modalities...")
  cols_civic <- get_civic_columns()
  cols_oncokb <- get_oncokb_columns()
  cols_keep <- c("SUBJECT_ID", "SAMPLE_ID", "ALTERATION_CATEGORY", cols_oncokb, cols_civic)

  if (!is.null(df_cna)){
    cols_keep_cna <- intersect(colnames(df_cna), c(cols_keep, "HUGO_SYMBOL"))
    df_cna <- df_cna[cols_keep_cna]
  } else {
    df_cna <- data.frame()
  }

  if (!is.null(df_fus)){
    cols_keep_fus <- intersect(colnames(df_fus), c(cols_keep, "FUSION"))
    df_fus <- df_fus[cols_keep_fus]
  } else {
    df_fus <- data.frame()
  }

  if (!is.null(df_msi)){
    cols_keep_msi <- intersect(colnames(df_msi),
                               c(cols_keep, "HUGO_SYMBOL"))
    df_msi <- df_msi[cols_keep_msi]
  } else {
    df_msi <- data.frame()
  }

  if (!is.null(df_mut)){
    cols_keep_mut <- intersect(colnames(df_mut),
                               c(cols_keep, "HUGO_SYMBOL", "PROTEIN_CHANGE", "PROTEIN_CHANGE_MORE",
                                 "Variant_Classification", "EXON"))
    df_mut <- df_mut[cols_keep_mut]
  } else {
    df_mut <- data.frame()
  }

  if (!is.null(df_tmb)){
    cols_keep_tmb <- intersect(colnames(df_tmb),
                               c(cols_keep, "HUGO_SYMBOL"))
    df_tmb <- df_tmb[cols_keep_tmb]
  } else {
    df_tmb <- data.frame()
  }

  if (!is.null(df_exp_arv7)){
    cols_keep_exp_arv7 <- intersect(colnames(df_exp_arv7),
                               c(cols_keep, "HUGO_SYMBOL", "ALTERATION"))
    df_exp_arv7 <- df_exp_arv7[cols_keep_exp_arv7]
  } else {
    df_exp_arv7 <- data.frame()
  }

  cat("done!\n")
  bind_rows(df_cna, df_fus, df_msi, df_mut, df_tmb, df_exp_arv7)
}


make_identifier <- function(df, cols){
  df[cols] <- df[cols] %>% replace(is.na(.), "N/A")
  df <- df %>% unite("Row_Identifier", all_of(cols), remove=F, sep="/")
  df
}


collapse_oncokb_annotations <- function(df_agg){
  cat("-collapsing oncokb annotations...")
  col_cite <- "TX_CITATIONS"

  #### sensitivity
  levels_sen <- sort(colnames(df_agg)[grepl("^(?=LEVEL_[0-9]).*", colnames(df_agg), perl=T)])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Oncokb_Sen_Level=NA, Oncokb_Sen_Drug=NA, Oncokb_Sen_Citations=NA)
  for (level in levels_sen){
    df_agg_best <- df_agg_best %>%
      mutate(Oncokb_Sen_Level=ifelse(is.na(Oncokb_Sen_Level)&!is.na(.data[[level]]), level, Oncokb_Sen_Level)) %>%
      mutate(Oncokb_Sen_Drug=ifelse(Oncokb_Sen_Level==level, .data[[level]], Oncokb_Sen_Drug)) %>%
      mutate(Oncokb_Sen_Citations=ifelse(Oncokb_Sen_Level==level, .data[[col_cite]], Oncokb_Sen_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Oncokb_Sen_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Oncokb_Sen_Level=NA, Oncokb_Sen_Drug=NA, Oncokb_Sen_Citations=NA)
  dfs_agg_all <- list()
  for (level in levels_sen){
    df_agg_level <- df_agg_all %>%
      mutate(Oncokb_Sen_Level=level) %>%
      mutate(Oncokb_Sen_Drug=.data[[level]]) %>%
      mutate(Oncokb_Sen_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Oncokb_Sen_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Oncokb_Sen_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Oncokb_Sen_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Oncokb_Sen_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  sen_level_simple <- list(LEVEL_1="Tier1", LEVEL_2="Tier1", LEVEL_3A="Tier2", LEVEL_3B="Tier3", LEVEL_4="Tier3")
  df_agg$Oncokb_Sen_Level_Simple <- recode(df_agg$Oncokb_Sen_Level, !!!sen_level_simple)

  #### resistance
  levels_res <- sort(colnames(df_agg)[grepl("^LEVEL_R", colnames(df_agg))])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Oncokb_Res_Level=NA, Oncokb_Res_Drug=NA, Oncokb_Res_Citations=NA)
  for (level in levels_res){
    df_agg_best <- df_agg_best %>%
      mutate(Oncokb_Res_Level=ifelse(is.na(Oncokb_Res_Level)&!is.na(.data[[level]]), level, Oncokb_Res_Level)) %>%
      mutate(Oncokb_Res_Drug=ifelse(Oncokb_Res_Level==level, .data[[level]], Oncokb_Res_Drug)) %>%
      mutate(Oncokb_Res_Citations=ifelse(Oncokb_Res_Level==level, .data[[col_cite]], Oncokb_Res_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Oncokb_Res_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Oncokb_Res_Level=NA, Oncokb_Res_Drug=NA, Oncokb_Res_Citations=NA)
  dfs_agg_all <- list()
  for (level in levels_res){
    df_agg_level <- df_agg_all %>%
      mutate(Oncokb_Res_Level=level) %>%
      mutate(Oncokb_Res_Drug=.data[[level]]) %>%
      mutate(Oncokb_Res_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Oncokb_Res_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Oncokb_Res_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Oncokb_Res_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Oncokb_Res_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  res_level_simple <- list(LEVEL_R1="Tier1", LEVEL_R2="Tier2", LEVEL_R3="Tier3")
  df_agg$Oncokb_Res_Level_Simple <- recode(df_agg$Oncokb_Res_Level, !!!res_level_simple)

  cat("done!\n")
  cols <- c("Oncokb_Annotated",
            "Oncokb_Sen_Level", "Oncokb_Sen_Level_Simple", "Oncokb_Sen_Drug", "Oncokb_Sen_Drug_Best_Level", 
            "Oncokb_Sen_Citations", "Oncokb_Res_Level", "Oncokb_Res_Level_Simple", "Oncokb_Res_Drug",
            "Oncokb_Res_Drug_Best_Level", "Oncokb_Res_Citations", "Row_Id")
  df_agg %>% select(all_of(cols))
}


collapse_civic_annotations <- function(df_agg){
  cat("-collapsing civic annotations...")
  col_cite <- "CIViC_Matching_Citation"

  #### sensitivity
  levels_sen <- sort(colnames(df_agg)[grepl("^Predictive:P", colnames(df_agg), perl=T)])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Civic_Sen_Level=NA, Civic_Sen_Drug=NA, Civic_Sen_Citations=NA)
  for (level in levels_sen){
    df_agg_best <- df_agg_best %>%
      mutate(Civic_Sen_Level=ifelse(is.na(Civic_Sen_Level)&!is.na(.data[[level]]), level, Civic_Sen_Level)) %>%
      mutate(Civic_Sen_Drug=ifelse(Civic_Sen_Level==level, .data[[level]], Civic_Sen_Drug)) %>%
      mutate(Civic_Sen_Citations=ifelse(Civic_Sen_Level==level, .data[[col_cite]], Civic_Sen_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Civic_Sen_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Civic_Sen_Level=NA, Civic_Sen_Drug=NA, Civic_Sen_Citations=NA)
  dfs_agg_all <- list()
  for (level in levels_sen){
    df_agg_level <- df_agg_all %>%
      mutate(Civic_Sen_Level=level) %>%
      mutate(Civic_Sen_Drug=.data[[level]]) %>%
      mutate(Civic_Sen_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Civic_Sen_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Civic_Sen_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Civic_Sen_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Civic_Sen_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  sen_level_simple <- list(`Predictive:P:A`="Tier1", `Predictive:P:B`="Tier2", `Predictive:P:C`="Tier3",
                           `Predictive:P:D`="Tier3", `Predictive:P:E`="Tier3")
  df_agg$Civic_Sen_Level_Simple <-recode(df_agg$Civic_Sen_Level, !!!sen_level_simple)

  #### resistance
  levels_res <- sort(colnames(df_agg)[grepl("^Predictive:N", colnames(df_agg))])

  # for each alteration, select the annotation with the best level
  df_agg_best <- df_agg %>% mutate(Civic_Res_Level=NA, Civic_Res_Drug=NA, Civic_Res_Citations=NA)
  for (level in levels_res){
    df_agg_best <- df_agg_best %>%
      mutate(Civic_Res_Level=ifelse(is.na(Civic_Res_Level)&!is.na(.data[[level]]), level, Civic_Res_Level)) %>%
      mutate(Civic_Res_Drug=ifelse(Civic_Res_Level==level, .data[[level]], Civic_Res_Drug)) %>%
      mutate(Civic_Res_Citations=ifelse(Civic_Res_Level==level, .data[[col_cite]], Civic_Res_Citations))
  }
  df_agg_best <- df_agg_best %>% mutate(Civic_Res_Drug_Best_Level=1)

  # for each alteration, combine annotations from all levels
  df_agg_all <- df_agg %>% mutate(Civic_Res_Level=NA, Civic_Res_Drug=NA, Civic_Res_Drug=NA)
  dfs_agg_all <- list()
  for (level in levels_res){
    df_agg_level <- df_agg_all %>%
      mutate(Civic_Res_Level=level) %>%
      mutate(Civic_Res_Drug=.data[[level]]) %>%
      mutate(Civic_Res_Citations=.data[[col_cite]]) %>%
      filter(!is.na(Civic_Res_Drug))
    if (nrow(df_agg_level)>0){
      dfs_agg_all[[level]] <- df_agg_level
    }
  }
  df_agg_all <- bind_rows(dfs_agg_all)
  df_agg_all <- df_agg_all %>% mutate(Civic_Res_Drug_Best_Level=0)

  # combine and remove duplicates on all columns except Civic_Res_Drug_Best_Level
  cols_distinct <- colnames(df_agg_best %>% select(-Civic_Res_Drug_Best_Level))
  df_agg <- bind_rows(df_agg_best, df_agg_all)
  df_agg <- df_agg %>% distinct(across(all_of(cols_distinct)), .keep_all=T)

  # recode oncokb to escat levels
  res_level_simple <- list(`Predictive:N:A`="Tier1", `Predictive:N:B`="Tier2", `Predictive:N:C`="Tier3",
                           `Predictive:N:D`="Tier3", `Predictive:N:E`="Tier3")
  df_agg$Civic_Res_Level_Simple <- recode(df_agg$Civic_Res_Level, !!!res_level_simple)

  cat("done!\n")
  cols <- c("Civic_Annotated",
            "Civic_Sen_Level", "Civic_Sen_Level_Simple", "Civic_Sen_Drug", "Civic_Sen_Drug_Best_Level",
            "Civic_Sen_Citations", "Civic_Res_Level", "Civic_Res_Level_Simple", "Civic_Res_Drug",
            "Civic_Res_Drug_Best_Level", "Civic_Res_Citations", "Row_Id")
  df_agg %>% select(all_of(cols))
}


get_alteration <- function(df_agg){
  regex <- "(?<=\\/)[0-9]*"
  cols <- c("Hugo_Symbol", "Alteration_Category", "Alteration", "Alteration_Detail")
  df_agg %>% 
    mutate(Alteration=NA) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Deletion", "Del", Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Amplification", "Amp", Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Fusion", "Fus", Alteration)) %>% 
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="MSI", "MSI High", Alteration)) %>% 
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="TMB", "TMB High", Alteration)) %>% 
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Exp", ALTERATION, Alteration)) %>% 
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Mut",gsub("^p.","",PROTEIN_CHANGE),Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Ins",
                             paste("Exon", gsub(regex,"",EXON,perl=T), "Ins"), Alteration)) %>%
    mutate(Alteration=ifelse(is.na(Alteration)&ALTERATION_CATEGORY=="Del",
                             paste("Exon", gsub(regex,"",EXON,perl=T), "Del"), Alteration)) %>%
    mutate(Alteration_Detail=ifelse(ALTERATION_CATEGORY=="Ins",gsub("^p.","",PROTEIN_CHANGE),NA)) %>%
    mutate(Alteration_Detail=ifelse(ALTERATION_CATEGORY=="Del",gsub("^p.","",PROTEIN_CHANGE),Alteration_Detail)) %>%
    mutate(Alteration=ifelse(!is.na(PROTEIN_CHANGE_MORE), PROTEIN_CHANGE_MORE, Alteration)) %>%
    mutate(Alteration_Detail=ifelse(!is.na(PROTEIN_CHANGE_MORE), PROTEIN_CHANGE, Alteration_Detail)) %>%
    mutate(Alteration=ifelse(!is.na(Alteration), gsub("\\/","",Alteration,perl=T), NA)) %>%
    mutate(Alteration=ifelse(is.na(Alteration), Variant_Classification, Alteration)) %>%
    mutate(HUGO_SYMBOL=ifelse(ALTERATION_CATEGORY=="Fusion", FUSION, HUGO_SYMBOL)) %>% 
    rename(Alteration_Category=ALTERATION_CATEGORY, Hugo_Symbol=HUGO_SYMBOL) %>% 
    select(all_of(cols))
}


harmonize_drugs <- function(x, df_drug){
  x <- x[!is.na(x)]
  if (length(x)==0){
    return(NA)
  } else {
    vals <- unlist(strsplit(x, ",|;|\\+|\\|"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- toupper(vals)
    if (!all(vals %in% df_drug$Drug)){
      vals_mis <- setdiff(vals, df_drug$Drug)
      warning(paste("-the following drugs are not in the drug table:", paste(vals_mis, collapse=",")))
    } else {
      vals_mis <- c()
    }

    x <- toupper(x)
    df_recode <- df_drug %>% filter(Drug %in% vals) %>% select(Drug, DCI) %>% distinct()
    df_recode <- bind_rows(df_recode, tibble(Drug=vals_mis, DCI=vals_mis))
    for (i in 1:nrow(df_recode)){
      x <- gsub(df_recode[i,"Drug"], df_recode[i,"DCI"], x)
    }

    vals <- unlist(strsplit(x, ",|;"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- sort(unique(vals))
    return(paste0(vals, collapse="|"))
  }
}


union_drugs <- function(vals){
  vals <- vals[!is.na(vals)]
  if (length(vals)==0){
    return(NA)
  } else {
    vals <- paste0(vals, collapse=",")
    vals <- unlist(strsplit(vals, ",|;"))
    vals <- sapply(vals, function(v) toupper(trimws(v)))
    vals <- sort(unique(vals))
    return(paste0(vals, collapse="|"))
  }
}

save_table <- function(table, output){
  if (grepl(".gz$", output)){
    write.table(table, file=gsub(".gz$", "", output), sep="\t", quote=F, row.names=F, na="")
    system(paste("gzip", gsub(".gz$", "", output)))
  } else {
    write.table(table, file=output, sep="\t", quote=F, row.names=F, na="")
  }
  cat(paste("-file saved at", output, "\n"))
}


main <- function(args){
  df_bio <- load_table(args$bio) 
  df_cln <- load_table(args$cln) %>% rename(Tumor_Type=Project_TCGA_More)
  df_cna <- load_table(args$cna, guess_max=1e5)
  df_fus <- load_table(args$fus, guess_max=1e5)
  df_msi <- load_table(args$msi, guess_max=1e5)
  df_mut <- load_table(args$mut, guess_max=1e5)
  df_tmb <- load_table(args$tmb, guess_max=1e5)
  df_exp_arv7 <- load_table(args$exp_arv7, guess_max=1e5)
  df_gen <- load_table(args$gen, guess_max=1e5)

  # select driver genes
  df_gen <- df_gen %>% filter(Inclusion_Level==1)
  cat(paste("-INFO: selected", nrow(df_gen), "driver genes\n"))

  # select on driver genes
  if (!is.null(df_cna)){
    mask <- df_cna$Hugo_Symbol %in% df_gen$Hugo_Symbol
    cat(paste("-INFO: selected", paste0(sum(mask), "/", length(mask)), "lines from CNA table using driver genes list\n"))
    df_cna <- df_cna %>% filter(mask)
  }

  if (!is.null(df_mut)){
    mask <- df_mut$Hugo_Symbol %in% df_gen$Hugo_Symbol
    cat(paste("-INFO: selected", paste0(sum(mask), "/", length(mask)), "lines from MUT table using driver genes list\n"))
    df_mut <- df_mut %>% filter(mask)
  }

  # add Sample_Id and/or Subject_Id where missing
  cols_bio <- c("Subject_Id", "Sample_Id")
  cols_cln <- c("Subject_Id", "Tumor_Type", "MSKCC_Oncotree", "Civic_Disease")
  df_cln_dna <- df_cln %>% filter(grepl("DNA_N\\|DNA_T", Sample_Type)) %>% 
    unite("Sample_Id", Sample_Id_DNA_T, Sample_Id_DNA_N, sep="_vs_", remove=F)
  df_cln_rna <- df_cln %>% filter(grepl("RNA_T", Sample_Type)) %>% 
    mutate(Sample_Id=Sample_Id_RNA_T)

  if (!is.null(df_cna)){
    df_cna <- df_cna %>% unite("Sample_Id", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_")
    df_cna <- left_join(df_cna, df_cln_dna[,c("Sample_Id", "Subject_Id")], by="Sample_Id")
  }

  if (!is.null(df_msi)){
    df_msi <- left_join(df_msi, df_bio[,cols_bio] %>% distinct(), by=c("Tumor_Sample_Id"="Sample_Id"))
    df_msi <- df_msi %>% unite("Sample_Id", Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_")
    df_msi <- left_join(df_msi, df_cln[,cols_cln] %>% distinct(), by=c("Subject_Id"))
  }

  if (!is.null(df_mut)){
    df_mut <- df_mut %>% unite("Sample_Id", Tumor_Sample_Barcode, Matched_Norm_Sample_Barcode, sep="_vs_")
    df_mut <- left_join(df_mut, df_cln_dna[,c("Sample_Id", "Subject_Id")], by="Sample_Id")
  }

  if (!is.null(df_tmb)){
    df_tmb <- df_tmb %>% unite("Sample_Id", Tumor_Sample_Id, Normal_Sample_Id, sep="_vs_")
    df_tmb <- left_join(df_tmb, df_cln[,cols_cln] %>% distinct(), by=c("Subject_Id"))
  }

  if (!is.null(df_fus)){
    df_fus <- left_join(df_fus, df_cln_rna[,c("Sample_Id", "Subject_Id")], by="Sample_Id")
  }
  if (!is.null(df_exp_arv7)) df_exp_arv7 <- left_join(df_exp_arv7, df_cln[,cols_cln] %>% distinct(), by=c("Subject_Id"))

  # select in-design samples
  df_cna <- select_in_design_samples(df_cna, df_cln, cat="cna")
  df_fus <- select_in_design_samples(df_fus, df_cln, cat="fus")
  df_msi <- select_in_design_samples(df_msi, df_cln, cat="msi")
  df_mut <- select_in_design_samples(df_mut, df_cln, cat="mut")
  df_tmb <- select_in_design_samples(df_tmb, df_cln, cat="tmb")
  df_exp_arv7 <- select_in_design_samples(df_exp_arv7, df_cln, cat="exp")

  # process each modality
  df_cna <- process_cna(df_cna)
  df_fus <- process_fus(df_fus)
  df_msi <- process_msi(df_msi)
  df_mut <- process_mut(df_mut)
  df_tmb <- process_tmb(df_tmb)
  df_exp_arv7 <- process_exp_arv7(df_exp_arv7)

  # aggregate
  df_agg <- aggregate_alterations(df_cna=df_cna,
                                  df_fus=df_fus,
                                  df_msi=df_msi,
                                  df_mut=df_mut,
                                  df_tmb=df_tmb,
                                  df_exp_arv7=df_exp_arv7)
  df_agg <- df_agg %>% distinct()

  # get table of alterations
  df_agg_alteration <- get_alteration(df_agg)

  # add indicators of oncokb and civic annotation
  mask_okb <- !is.na(df_agg["ONCOGENIC"])
  mask_civ <- !is.na(df_agg["CIViC_Matching_Gene_Variant"])
  df_agg[mask_okb, "Oncokb_Annotated"] <- "Yes"
  df_agg[mask_civ, "Civic_Annotated"] <- "Yes"

  # correct some civic annotations that are unspecific
  df_agg <- correct_civic_annotations(df_agg, df_cln)

  # add row id because collapse of oncokb and civic annotations may produce multiple
  # rows per row id in case multiple levels of annotation coexist for the same alteration
  df_agg <- df_agg %>% mutate(Row_Id=1:nrow(df_agg))

  # collapse annotations
  df_agg_oncokb <- collapse_oncokb_annotations(df_agg)
  df_agg_civic <- collapse_civic_annotations(df_agg)

  # add Tumor_Type columns
  df_agg <- df_agg %>% rename(Sample_Id=SAMPLE_ID, Subject_Id=SUBJECT_ID)
  cols_cln <- c("Subject_Id", "Tumor_Type", "MSKCC_Oncotree", "Civic_Disease")
  df_agg <- left_join(df_agg, df_cln[,cols_cln] %>% distinct(), by="Subject_Id")
  cols_agg <- c("Variant_Classification")

  # column bind and joins
  df_fin <- bind_cols(df_agg[c("Sample_Id", cols_cln, cols_agg, "Row_Id")], df_agg_alteration)
  cols_fin <- colnames(df_fin)
  df_fin <- make_identifier(df_fin, setdiff(colnames(df_fin), "Row_Id"))

  df_agg_oncokb_best <- df_agg_oncokb %>% filter(Oncokb_Sen_Drug_Best_Level==1, Oncokb_Res_Drug_Best_Level==1)
  df_agg_oncokb_oth <- df_agg_oncokb %>% filter(Oncokb_Sen_Drug_Best_Level!=1 | Oncokb_Res_Drug_Best_Level!=1)
  df_agg_civic_best <- df_agg_civic %>% filter(Civic_Sen_Drug_Best_Level==1, Civic_Res_Drug_Best_Level==1)
  df_agg_civic_oth <- df_agg_civic %>% filter(Civic_Sen_Drug_Best_Level!=1 | Civic_Res_Drug_Best_Level!=1)

  df_fin_best <- full_join(df_fin, df_agg_oncokb_best, by="Row_Id")
  df_fin_best <- full_join(df_fin_best, df_agg_civic_best, by="Row_Id")

  df_fin_oth_oncokb <- right_join(df_fin, df_agg_oncokb_oth, by="Row_Id")
  df_fin_oth_civic <- right_join(df_fin, df_agg_civic_oth, by="Row_Id")
  df_fin_oth <- full_join(df_fin_oth_oncokb, df_fin_oth_civic) 
  df_fin <- bind_rows(df_fin_best, df_fin_oth) %>% arrange(Row_Identifier)

  df_fin$Row_Id <- NULL
  df_fin <- df_fin %>% distinct()

  # for CNA, set Variant_Classification to Deletion or Amplification
  mask_null <- is.na(df_fin$Variant_Classification)
  df_fin[mask_null, "Variant_Classification"] <- df_fin[mask_null, "Alteration_Category"]

  # # restrict table to sample from subjects in df_cln
  # samples_keep <- c(df_cln$Sample_Id_DNA_T, df_cln$Sample_Id_RNA_T)
  # samples_keep <- samples_keep[which(!is.na(samples_keep))]
  # df_fin <- df_fin %>% filter(Subject_Id %in% df_cln$Subject_Id) %>% filter(Sample_Id %in% samples_keep)

  # consensus level by taking the best of oncokb and civic
  cols <- c("Oncokb_Sen_Level_Simple", "Civic_Sen_Level_Simple")
  df_fin$Sen_Level_Simple <- pmin(df_fin[[cols[1]]], df_fin[[cols[2]]], na.rm=T)
  cols <- c("Oncokb_Res_Level_Simple", "Civic_Res_Level_Simple")
  df_fin$Res_Level_Simple <- pmin(df_fin[[cols[1]]], df_fin[[cols[2]]], na.rm=T)

  # harmonize drug names
  df_drug <- load_table(args$drug)
  cols_drug <- c("Oncokb_Sen_Drug", "Oncokb_Res_Drug", "Civic_Sen_Drug", "Civic_Res_Drug")
  for (col_drug in cols_drug){
    df_uni <- tibble(Old=unique(df_fin[[col_drug]]))
    df_uni$New <- unlist(lapply(df_uni$Old, function(x) harmonize_drugs(x, df_drug)))
    df_fin <- left_join(df_fin, df_uni, by=setNames(c("Old"), col_drug))
    df_fin <- df_fin %>% select(-all_of(col_drug)) %>% rename(!!col_drug:=New)
  }

  # # consensus drug by taking union of oncokb and civic
  # df_fin <- df_fin %>% rowwise() %>%
  #   mutate(Res_Drug_Simple=union_drugs(c(Oncokb_Res_Drug, Civic_Res_Drug))) %>%
  #   mutate(Sen_Drug_Simple=union_drugs(c(Oncokb_Sen_Drug, Civic_Sen_Drug)))

  # select only annotations with best level
  df_fin_best <- df_fin %>% filter_at(vars(ends_with("_Best_Level")), all_vars(. == 1))

  # save tables
  save_table(df_fin, args$output_all)
  save_table(df_fin_best, args$output_best)
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregate alterations.')
  parser$add_argument("--bio", type="character", help="Path to input curated biospecimen table.",
                      default="../../../data/tcga/clinical/curated/bio_tcga_in_design_curated.tsv")
  parser$add_argument("--cln", type="character", help="Path to input curated clincal table.",
                      default="../../../data/tcga/clinical/curated/cln_tcga_in_design_curated.tsv")
  parser$add_argument("--cna", type="character", help="Path to input annotated CNAs table.",
                      default="../../../data/tcga/wes/somatic_cna/somatic_calls_union_ann.tsv.gz")
  parser$add_argument("--fus", type="character", help="Path to input annotated fusions table.",
                      default="../../../data/tcga/rna/fusions/tcga_annotated_filtered_union_ann.tsv.gz")
  parser$add_argument("--msi", type="character", help="Path to input table of MSI.",
                      default="../../../data/tcga/wes/somatic_msi/somatic_msi.tsv")
  parser$add_argument("--mut", type="character", help="Path to input annotated mutations table.",
                      default="../../../data/tcga/wes/somatic_maf/somatic_calls_union_ann.maf.gz")
  parser$add_argument("--tmb", type="character", help="Path to input table of TMB.",
                      default="../../../data/tcga/wes/summary/somatic_maf.tsv")
  parser$add_argument("--exp_arv7", type="character", nargs="*", help="Path to table of AR-V7 expression.",
                      default=NULL)
  parser$add_argument("--gen", type="character", help="Path to input list of genes to be considered.",
                      default="../../../data/resources/curated/cancer_genes_curated.tsv")
  parser$add_argument("--target_bed", type="character", help="Path to target bed file.",
                      default="../../../data/resources/target_files/all_targets_intersect_padded_10n.bed")
  parser$add_argument('--drug', type="character", help='Path to table of drugs.',
                      default="../../../data/resources/drug_tables/Table_Drugs_v7.xlsx")
  parser$add_argument("--output_best", type="character", help="Paths to output aggregated table",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_tcga.tsv")
  parser$add_argument("--output_all", type="character", help="Paths to output aggregated table",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations_tcga_all.tsv")
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  for (arg_name in names(args)){
    if (!is.null(args[[arg_name]])){
      if (is.list(args[[arg_name]]) & length(args[[arg_name]])==0){
        args[[arg_name]] <- NULL
      } else if (args[[arg_name]]==""){
        args[[arg_name]] <- NULL
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
