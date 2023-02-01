# @modified: 06 Jul 21
# @modified: 19 Oct 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

remove_multi_polymorphisms <- function(df_maf){
  if ("Variant_Type" %in% colnames(df_maf)){
    types_allowed <- c("DEL", "INS", "SNP")
    mask_in <- df_maf$Variant_Type %in% types_allowed
    cat(paste0("-INFO: excluded ", sum(!mask_in), "/",  length(mask_in), " lines not in ", 
               paste(types_allowed, collapse=";"), ".\n"))
    df_maf <- df_maf[mask_in,]
  }

  df_maf
}


select_vc <- function(df_maf, vc_list){
  if ("Variant_Classification" %in% colnames(df_maf)){
    mask_in <- df_maf$Variant_Classification %in% vc_list
    cat(paste0("-INFO: excluded ", sum(!mask_in), "/",  length(mask_in), " lines not in ", 
               paste(vc_list, collapse=";"), ".\n"))
    df_maf <- df_maf[mask_in,]
  }

  df_maf
}


get_tcga_tt_to_mutpanning_tt <- function(){
  list("ANUS - Not_TCGA"=NULL,
       "BLCA"="Bladder",
       "BLCA - Not_TCGA"="Bladder",
       "BRCA"="Breast",
       "CESC"="Cervix",
       "CHOL"="Cholangio",
       "COAD"="Colorectal",
       "ESCA"="Gastroesophageal",
       "GBM"="Brain",
       "HNAC - Not_TCGA"="AdenoidCystic",
       "HNSC"="HeadNeck",
       "KICH"="KidneyNonClear",
       "KIRC"="KidneyClear",
       "LIHC"="Liver",
       "LUAD"="LungAD",
       "LUSC"="LungSC",
       "LUNEC - Not_TCGA"=NULL,
       "MESO"="Pleura",
       "MISC - Not_TCGA"=NULL,
       "N/A"="NULL",
       "OV"="Ovarian",
       "PAAD"="Pancreas",
       "PRAD"="Prostate",
       "SARC"="Sarcoma",
       "SARC - Not_TCGA"="Sarcoma",
       "SI - Not_TCGA"=NULL,
       "SKCM"="Skin",
       "STAD"="Gastroesophageal",
       "TGCT"="TesticularGermCell",
       "THCA"="Thyroid",
       "THYM"="Thymus",
       "UCEC"="Endometrium",
       "Unknown_Primary"=NULL,
       "UVM"="UvealMelanoma")
}


get_genes_list <- function(df_maf, args){
  genes_list <- c()
  if ("Tumor_Sample_Barcode" %in% colnames(df_maf)){
    col_sam <- "Tumor_Sample_Barcode"
  } else {
    col_sam <- "Sample_Id"
  }

  if ("mutpanning_pairs" %in% args$genes_selection_mode){
    df_mutpanning_pairs <- load_table(args$genes_selection_mutpanning_pairs)
    tcga_tt_to_mutpanning_tt <- get_tcga_tt_to_mutpanning_tt()
    mutpanning_tt <- sapply(unique(df_maf$Tumor_Type),
                            function(tt) tcga_tt_to_mutpanning_tt[[tt]], USE.NAMES=F)
    if (is.null(mutpanning_tt)){
      stop("-the tumor type specified has no corresponding cancer type in the mutpanning pairs table.")
    } else {
      genes_list_mutpanning <- df_mutpanning_pairs %>%
        filter(`Cancer Type` %in% mutpanning_tt) %>%
        pull(Gene)
      genes_list <- union(genes_list, genes_list_mutpanning)
    }
  }

  if ("mutpanning_run" %in% args$genes_selection_mode){
    df_mutpanning_run <- load_table(args$genes_selection_mutpanning_run)
    genes_list_mutpanning <- df_mutpanning_run %>%
      filter(FDR <= 0.25) %>%
      pull(Name)
    genes_list <- union(genes_list, genes_list_mutpanning)
  }

  if ("min_freq" %in% args$genes_selection_mode){
    n_tumors <- length(unique(df_maf[[col_sam]]))
    n_tumors_min <- ceiling(args$genes_selection_min_freq * n_tumors)
    genes_list_freq <- df_maf %>%
      filter(Hugo_Symbol != "Unknown") %>%
      distinct(.data[[col_sam]], Hugo_Symbol) %>%
      group_by(Hugo_Symbol) %>%
      summarize(n_tumors=n()) %>%
      filter(n_tumors >= n_tumors_min) %>%
      pull(Hugo_Symbol)
    genes_list <- union(genes_list, genes_list_freq)
  }

  if ("min_count" %in% args$genes_selection_mode){
    n_tumors <- length(unique(df_maf[[col_sam]]))
    n_tumors_min <- args$genes_selection_min_count
    genes_list_freq <- df_maf %>%
      filter(Hugo_Symbol != "Unknown") %>%
      distinct(.data[[col_sam]], Hugo_Symbol) %>%
      group_by(Hugo_Symbol) %>%
      summarize(n_tumors=n()) %>%
      filter(n_tumors >= n_tumors_min) %>%
      pull(Hugo_Symbol)
    genes_list <- union(genes_list, genes_list_freq)
  }

  if ("list" %in% args$genes_selection_mode){
    genes_list <- union(genes_list, args$genes_selection_list)
  }

  genes_list
}


select_genes <- function(df_maf, genes_list){
  df_maf_sel <- df_maf %>% filter(Hugo_Symbol %in% genes_list)
  cat(paste0("-removed ", nrow(df_maf)-nrow(df_maf_sel), "/", nrow(df_maf), " variants using Hugo_Symbol\n"))
  df_maf_sel
}


get_vc_to_simple_vc <- function(){
  list(Frame_Shift_Del="Frameshift_Indel",
       Frame_Shift_Ins="Frameshift_Indel",
       Splice_Site="Splice_Site",
       Translation_Start_Site="Other",
       Nonsense_Mutation="Nonsense",
       Nonstop_Mutation="Nonsense",
       In_Frame_Del="Inframe_Indel",
       In_Frame_Ins="Inframe_Indel",
       Missense_Mutation="Missense",
       Start_Codon_Del="Start_Stop_Codon",
       Start_Codon_SNP="Start_Stop_Codon",
       Stop_Codon_Del="Start_Stop_Codon",
       Stop_Codon_Ins="Start_Stop_Codon",
       Amp="Amplification",
       Del="Deletion",
       `MSI High`="MSI",
       `TMB High`="TMB",
       Fusion="Fusion")
}


simplify_variant_classification <- function(df_maf, vc_to_simple_vc){
  for (vc_value in unique(df_maf$Variant_Classification)){
    if (vc_value %in% names(vc_to_simple_vc)){
      simple_vc_value <- vc_to_simple_vc[[vc_value]]
    } else {
      simple_vc_value <- "Other"
    }
    df_maf[df_maf$Variant_Classification==vc_value, "Variant_Classification"] <- simple_vc_value
  }

  df_maf
}


get_left_bar_data <- function(df_maf, genes_list){
  if ("Tumor_Sample_Barcode" %in% colnames(df_maf)){
    col_sam <- "Tumor_Sample_Barcode"
  } else {
    col_sam <- "Sample_Id"
  }

  df_maf %>%
    filter(Hugo_Symbol %in% genes_list) %>%
    distinct(.data[[col_sam]], Hugo_Symbol) %>%
    group_by(Hugo_Symbol) %>%
    summarize(Count=n())
}


get_right_bar_data <- function(df_maf, genes_list, genes_selection_mutpanning_run){
  if (is.null(genes_selection_mutpanning_run)){
    return(NULL)
  }
  df_mutpanning_run <- load_table(genes_selection_mutpanning_run)
  df_mutpanning_run %>%
    rename(Hugo_Symbol=Name) %>%
    filter(Hugo_Symbol %in% genes_list) %>%
    arrange(Significance, FDR) %>%
    distinct(Hugo_Symbol, .keep_all=T) %>%
    mutate(`-log10 (fdr)`=ifelse(FDR==0, 100, -log10(FDR))) %>%
    select(c(Hugo_Symbol, `-log10 (fdr)`))
}


save_table <- function(df, name, output){
  write.table(df, output, quote=F, row.names=F, sep="\t")
  cat(paste("-input", name, "file saved at", output, "\n"))
}


main <- function(args){
  # load mutations
  df_maf <- load_table(args$mutations)

  # select samples
  col_sam_cln <- paste0("Sample_Id_", args$sample_type)
  df_cln <- load_cln(study=args$cohort)
  df_cln <- df_cln %>% filter(!is.na(.data[[col_sam_cln]]))

  # if possible select sample pairs
  if (args$sample_type=="DNA_T"){
    col_pair <- "DNA_P"
    cols_psb <- c("Tumor_Sample_Barcode", "Matched_Norm_Sample_Barcode")
    cols_psi <- c("Tumor_Sample_Id", "Normal_Sample_Id")
    cols_pair <- NULL
    if (all(cols_psb %in% colnames(df_maf))){
      cols_pair <- cols_psb
    } else if (all(cols_psi %in% colnames(df_maf))){
      cols_pair <- cols_psi
    }

    if (!is.null(cols_pair)){
      df_maf <- df_maf %>% unite(!!col_pair, all_of(cols_pair), sep="_vs_", remove=F)
      df_cln <- df_cln %>% unite(!!col_pair, Sample_Id_DNA_T, Sample_Id_DNA_N, sep="_vs_", remove=F)
      mask_in <- df_maf[[col_pair]] %in% df_cln[[col_pair]]
      cat(paste0("-INFO: selected ", sum(mask_in), "/",  length(mask_in), " lines from DNA_T/DNA_N in-design pairs.\n"))
      df_maf <- df_maf[mask_in,]
      df_maf[[col_pair]] <- NULL
    }
  }

  # select by tumor type and algorithm
  df_cln <- df_cln %>% filter(Project_TCGA_More==args$tumor_type)
  samples_cln <- df_cln[[col_sam_cln]]

  df_bio <- load_bio(study=args$cohort)
  df_bio <- df_bio %>% filter(if_all(all_of(args$algorithms), ~ .x == 1))
  samples_bio <- df_bio$Sample_Id
  
  samples <- intersect(samples_cln, samples_bio)
  cat(paste0("-INFO: selected ", length(samples), " samples from ", args$tumor_type, " and ", 
             paste0(args$algorithms, collapse=";"), ".\n"))

  df_cln <- df_cln[df_cln[[col_sam_cln]] %in% samples,]
  df_cln[["Tumor_Sample_Barcode"]] <- df_cln[[col_sam_cln]]

  col_tsb <- "Tumor_Sample_Barcode"
  col_sid <- "Sample_Id"

  if (col_tsb %in% colnames(df_maf)){
    col_sam <- col_tsb
  } else if (col_sid %in% colnames(df_maf)){
    col_sam <- col_sid
  }
  mask_in <- df_maf[[col_sam]] %in% samples
  cat(paste0("-INFO: selected ", sum(mask_in), "/",  length(mask_in), " lines from ", length(samples), " samples.\n"))
  df_maf <- df_maf[mask_in,]
  
  # select variants
  df_maf <- remove_multi_polymorphisms(df_maf)

  # fill NA Variant_Classification if Any
  if (all(c("Variant_Classification", "Alteration", "Alteration_Category") %in% colnames(df_maf))){
    mask_null <- is.na(df_maf$Variant_Classification) |  df_maf$Variant_Classification=="N/A"
    if (args$sample_type=="DNA_T"){
      df_maf[mask_null, "Variant_Classification"] <- df_maf[mask_null, "Alteration"]
    } else if (args$sample_type=="RNA_T"){
      df_maf[mask_null, "Variant_Classification"] <- df_maf[mask_null, "Alteration_Category"]
    }
  }
  df_maf <- select_vc(df_maf, args$vc_selection_list)

  # select genes
  genes_list <- get_genes_list(df_maf, args)
  df_genes <- data.frame(Hugo_Symbol=genes_list)

  # remap values of Variant_Classification
  vc_to_simple_vc <- get_vc_to_simple_vc()
  df_maf <- simplify_variant_classification(df_maf, vc_to_simple_vc)

  # side bar plots
  df_left_bar_data <- get_left_bar_data(df_maf, genes_list)

  if (!is.null(args$output_rbar)){
    df_right_bar_data <- get_right_bar_data(df_maf, genes_list, args$genes_selection_mutpanning_run)
    save_table(df_right_bar_data, "right bar data", args$output_rbar)
  }
  
  # add columns required for read.maf function
  if (c("Alteration_Category") %in% colnames(df_maf)){
    old2new <- list("Amplification"="CNV", "Del"="DEL", "Deletion"="CNV",  "Exp"="EXP", "Fusion"="FUS",
                    "Ins"="INS", "MSI"="MSI", "Mut"="SNP", "TMB"="TMB")
    df_maf[["Variant_Type"]] <- unlist(plyr::revalue(df_maf[["Alteration_Category"]], old2new))
  }

  cols_req <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2",
                "Variant_Classification", "Variant_Type")
  for (col_req in cols_req){
    if (!col_req %in% colnames(df_maf))
      df_maf[[col_req]] <- NA
  }

  if (!"Tumor_Sample_Barcode" %in% colnames(df_maf)){
    df_maf$Tumor_Sample_Barcode <- df_maf$Sample_Id
  }

  # add fake lines for samples without mutations
  samples_not_maf <- setdiff(samples, df_maf$Tumor_Sample_Barcode)
  df_maf_fake <- tibble(Tumor_Sample_Barcode=samples_not_maf)
  df_maf <- bind_rows(df_maf, df_maf_fake)

  # save
  save_table(df_genes, "genes list", args$output_gen)
  save_table(df_cln, "cln", args$output_cln)
  save_table(df_maf, "maf", args$output_maf)
  save_table(df_left_bar_data, "left bar data", args$output_lbar)
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Build tables for oncoplot via maftools.')
  parser$add_argument("--mutations", type="character", help="Path to input mutations table.",
                      default="../../../results/somatic_mutations/mutpanning/prism/LUAD_inputs/mutations_LUAD.maf")
  parser$add_argument("--cohort", type="character", help="Name of the cohort.",
                      default="prism")
  parser$add_argument("--tumor_type", type="character", help="Name of the tumor type analyzed.",
                      default="LUAD")
  parser$add_argument("--sample_type", type="character", help="Choose 'DNA_T' or 'RNA_T'.",
                      default="DNA_T")
  parser$add_argument("--algorithms", type="character", nargs="+", help="Algorithms for samples selection.",
                      default=c("WES_somatic_maf"))
  parser$add_argument("--vc_selection_list", type="character", nargs="*",
                      help="Only variants with a Variant_Classification in this list will be included in the plot.",
                      default=c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Translation_Start_Site",
                                "Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                                "Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins", "Amp", "Del",
                                "MSI High", "TMB High"))
  parser$add_argument("--genes_selection_mode", type="character", nargs="*",
                      help=paste("How to select genes for the plot. Possible values are 'min_freq', 'min_count'",
                                 ", 'mutpanning_pairs' or 'list'. If you specfify multiple modes, the union of the",
                                 "genes from each mode is selected."),
                      default=c("mutpanning_pairs", "mutpanning_run"))
  parser$add_argument("--genes_selection_min_freq", type="double", default=0.01,
                      help="Used only if --genes_selection_mode contains 'min_freq'.")
  parser$add_argument("--genes_selection_min_count", type="integer", default=5,
                      help="Used only if --genes_selection_mode contains 'min_count'.")
  parser$add_argument("--genes_selection_list", type="character", default=NULL, nargs="*",
                      help=paste("Used only if --genes_selection_mode contains 'list'.",
                                 "Variants located in genes in this list will be included in the plot."))
  parser$add_argument("--genes_selection_mutpanning_pairs", type="character",
                      help=paste("Used only if --genes_selection_mode contains 'mutpanning_pairs'.",
                                 "Path to results table from mutpanning paper with gene-tumor type associations."),
                      default=NULL)
  parser$add_argument("--genes_selection_mutpanning_run", type="character", nargs="*",
                      help=paste("Used only if --genes_selection_mode contains 'mutpanning_run'.",
                                 "Variants located in genes in this list will be included in the plot."),
              default="../../../results/somatic_mutations/mutpanning/prism/LUAD/SignificanceFiltered/SignificanceLUAD.txt")
  parser$add_argument("--output_gen", type="character", help="Path to output genes list file.",
                      default="../../../results/somatic_mutations/oncoplots/genes_prism_LUAD.tsv")
  parser$add_argument("--output_cln", type="character", help="Path to output clinical file.",
                      default="../../../results/somatic_mutations/oncoplots/clinical_prism_LUAD.tsv")
  parser$add_argument("--output_maf", type="character", help="Path to output mutations file.",
                      default="../../../results/somatic_mutations/oncoplots/mutations_prism_LUAD.maf")
  parser$add_argument("--output_lbar", type="character", help="Path to output right bar data file.",
                      default="../../../results/somatic_mutations/oncoplots/left_bar_data_prism_LUAD.tsv")
  parser$add_argument("--output_rbar", type="character", help="Path to output left bar data file.",
                      default="../../../results/somatic_mutations/oncoplots/right_bar_data_prism_LUAD.tsv")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
