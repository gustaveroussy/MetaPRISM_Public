# @created: 09 Jun 21
# @modified: 24 Feb 22
# @authors: Yoann Pradat

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

## Arriba

preprocess_fusions_arriba <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Chr_1", "Breakpoint_1", "Strand_1", "AR_Strand_Fusion_1",
            "Gene_2", "Chr_2",  "Breakpoint_2", "Strand_2", "AR_Strand_Fusion_2",
            "AR_Split_Reads", "AR_Spanning_Reads", "AR_Total_Reads", "AR_Direction_1", "AR_Direction_2",
            "AR_Site_1", "AR_Site_2", "AR_Type", "AR_Confidence", "AR_Reading_Frame", "AR_Call")

  df_fus <- df_fus %>% rename(Gene_1=`#gene1`, AR_Site_1=site1, AR_Direction_1=direction1,
                              Gene_2=gene2, AR_Site_2=site2, AR_Direction_2=direction2, AR_Type=type,
                              AR_Confidence=confidence, AR_Spanning_Reads=discordant_mates,
                              AR_Reading_Frame=reading_frame)

  df_fus <- df_fus %>%
    separate(breakpoint1, c("Chr_1", "Breakpoint_1"), sep=":") %>%
    separate(breakpoint2, c("Chr_2", "Breakpoint_2"), sep=":") %>%
    separate(`strand1(gene/fusion)`, c("Strand_1", "AR_Strand_Fusion_1"), sep="/") %>%
    separate(`strand2(gene/fusion)`, c("Strand_2", "AR_Strand_Fusion_2"), sep="/") %>%
    mutate(AR_Split_Reads=split_reads1+split_reads2) %>%
    mutate(AR_Total_Reads=AR_Split_Reads+AR_Spanning_Reads)

  df_fus <- df_fus %>%
    separate_rows(Gene_1, sep=",") %>%
    mutate(Gene_1=gsub("\\([^\\)]*\\)", "", Gene_1))

  df_fus <- df_fus %>%
    separate_rows(Gene_2, sep=",") %>%
    mutate(Gene_2=gsub("\\([^\\)]*\\)", "", Gene_2))

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus <- df_fus %>% filter(AR_Confidence %in% c("medium", "high")) %>%
    filter(AR_Site_1!="intergenic" & AR_Site_2!="intergenic")
  df_fus["AR_Call"] <- 1
  df_fus[cols]
}

## DEEPEST PNAS 2019

preprocess_fusions_deepest_pnas_2019 <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Gene_Id_1", "Chr_1", "Breakpoint_1", "Strand_1",
            "Gene_2", "Chr_2", "Gene_Id_2", "Breakpoint_2", "Strand_2",
            "Deepest_Pnas_2019_Type", "Deepest_Pnas_2019_Total_Reads", "Deepest_Pnas_2019_Score",
            "Deepest_Pnas_2019_Call")

  df_sam <- read_tsv("../../../data/tcga/rna/deepest_pnas_2019/sample_list.tsv", show_col_types=F, progress=F)
  df_sam <- df_sam %>% filter(Sample %in% df_fus$Sample)

  df_fus <- left_join(df_fus, df_sam[,c("Sample", "Sample_Id")], by="Sample") %>% distinct()
  df_fus <- df_fus %>% rename(Gene_Id_1=gene_id_1, Gene_Id_2=gene_id_2,
                              Deepest_Pnas_2019_Total_Reads=numReads, Deepest_Pnas_2019_Score=junction_cdf)

  junction_format <- c("Chr", "Gene", "Breakpoint", "Strand")
  df_fus <- df_fus %>%
    separate(junction, c(paste0(junction_format, "_1"), paste0(junction_format, "_2"), "Deepest_Pnas_2019_Type"), sep=":")

  # extract sample ids analyzed
  df_ids <- load_from_data("tcga/rna/deepest_pnas_2019/sample_list.tsv")

  if (!"Sample_Id" %in% colnames(df_fus)){
    df_ids <- df_ids %>% rename(File_Name_Legacy=`RNA-seq file`)
    df_bio <- load_bio(study="tcga", mode="all")
    df_bio <- df_bio %>% separate_rows(File_Id_Legacy, File_Name_Legacy, File_Size_Legacy, sep="\\|") %>%
      mutate(File_Name_Legacy=gsub(".gz", "", File_Name_Legacy))
    cols_bio <- c("File_Name_Legacy", "Sample_Id")
    df_ids <- left_join(df_ids, df_bio[cols_bio], by="File_Name_Legacy")

    # some sample ids are 16 characters long, others are 28 characters long
    # make them all 16 characters long and match back to 28 characters long using df_ids
    df_fus <- left_join(df_fus, df_ids, by="Sample")
  }

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["Deepest_Pnas_2019_Call"] <- 1
  df_fus[cols]
}

## Ericscript

preprocess_fusions_ericscript <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Chr_1", "Breakpoint_1", "Strand_1",
            "Gene_2", "Chr_2",  "Breakpoint_2", "Strand_2",
            "ES_Split_Reads", "ES_Spanning_Reads", "ES_Total_Reads", "ES_Score", "ES_Fusion_Type",
            "ES_Gene_Expr_1", "ES_Gene_Expr_2", "ES_Gene_Expr_Fusion", "ES_Call")

  df_fus <- df_fus %>% rename(Gene_1=GeneName1, Chr_1=chr1, Breakpoint_1=Breakpoint1, Strand_1=strand1,
                              Gene_2=GeneName2, Chr_2=chr2, Breakpoint_2=Breakpoint2, Strand_2=strand2,
                              ES_Split_Reads=crossingreads, ES_Spanning_Reads=spanningreads, 
                              ES_Score=EricScore, ES_Fusion_Type=fusiontype, ES_Gene_Expr_1=GeneExpr1,
                              ES_Gene_Expr_2=GeneExpr2, ES_Gene_Expr_Fusion=GeneExpr_Fused)

  df_fus <- df_fus %>% mutate(ES_Total_Reads=ES_Split_Reads+ES_Spanning_Reads)

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["ES_Call"] <- 1
  df_fus[cols]
}

## Fusioncatcher

preprocess_fusions_fusioncatcher <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Gene_Id_1", "Chr_1", "Breakpoint_1", "Strand_1",
            "Gene_2", "Gene_Id_2", "Chr_2", "Breakpoint_2", "Strand_2",
            "FC_Spanning_Reads", "FC_Call")

  df_fus <- df_fus %>% 
    rename(Gene_1=`Gene_1_symbol(5end_fusion_partner)`, Gene_2=`Gene_2_symbol(3end_fusion_partner)`) %>%
    rename(Gene_Id_1=`Gene_1_id(5end_fusion_partner)`, Gene_Id_2=`Gene_2_id(3end_fusion_partner)`) %>%
    rename(FC_Spanning_Reads=Spanning_unique_reads)

  df_fus <- df_fus %>%
    separate(`Fusion_point_for_gene_1(5end_fusion_partner)`, c("Chr_1", "Breakpoint_1", "Strand_1"), sep=":") %>%
    separate(`Fusion_point_for_gene_2(3end_fusion_partner)`, c("Chr_2", "Breakpoint_2", "Strand_2"), sep=":")

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["FC_Call"] <- 1
  df_fus[cols]
}

## Pizzly

extract_pizzly_gene_id <- function(gene_id){
  unlist(strsplit(gene_id, "\\."))[[1]]
}


preprocess_fusions_pizzly <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Gene_Id_1", "Gene_2", "Gene_Id_2",
            "PZ_Split_Reads", "PZ_Spanning_Reads", "PZ_Total_Reads", "PZ_Call")

  df_fus <- df_fus %>% rename(Gene_1=geneA.name, Gene_Id_1=geneA.id, Gene_2=geneB.name, Gene_Id_2=geneB.id,
                              PZ_Split_Reads=splitcount, PZ_Spanning_Reads=paircount)

  df_fus <- df_fus %>% rowwise() %>% mutate(Gene_Id_1=extract_pizzly_gene_id(Gene_Id_1),
                                            Gene_Id_2=extract_pizzly_gene_id(Gene_Id_2))

  df_fus$PZ_Total_Reads <- df_fus$PZ_Split_Reads + df_fus$PZ_Spanning_Reads
  df_fus <- df_fus %>% filter(PZ_Split_Reads > 0, PZ_Spanning_Reads > 0)

  # remove duplicates
  cols_id <- c("Sample_Id", "Gene_1", "Gene_2")
  col_ord <- "PZ_Total_Reads"
  df_fus <- df_fus %>% group_by_at(vars(all_of(cols_id))) %>% slice_max(col_ord)

  df_fus["PZ_Call"] <- 1
  df_fus[cols]
}

## PRADA NAR 2018

preprocess_fusions_prada_nar_2018 <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Chr_1", "Breakpoint_1", "Prada_Nar_2018_Strand_Fusion_1",
            "Gene_2", "Chr_2", "Breakpoint_2", "Prada_Nar_2018_Strand_Fusion_2",
            "Prada_Nar_2018_Split_Reads", "Prada_Nar_2018_Spanning_Reads", "Prada_Nar_2018_Total_Reads",
            "Prada_Nar_2018_Reading_Frame", "Prada_Nar_2018_Call")

  # add ids
  df_sam <- read_tsv("../../../data/tcga/rna/prada_nar_2018/sample_list.tsv", show_col_types=F, progress=F)
  df_fus_a <- df_fus %>% filter(str_length(Sample)>=20)
  df_fus_b <- df_fus %>% filter(str_length(Sample)<20)

  df_fus_a <- df_fus_a %>% mutate(Sample_Id=str_sub(Sample, start=1, end=20))
  df_sam <- df_sam %>% mutate(Sample=str_sub(Sample_Id, start=1, end=16)) %>%
    filter(Sample %in% df_fus_b$Sample)
  df_fus_b <- left_join(df_fus_b, df_sam[,c("Sample", "Sample_Id")], by="Sample") %>% distinct()
  df_fus <- bind_rows(df_fus_a, df_fus_b)

  df_fus <- df_fus %>% rename(Gene_1=Gene_A, Chr_1=`Gene A Chr`, Prada_Nar_2018_Strand_Fusion_1=`Gene A Strand`,
                              Gene_2=Gene_B, Chr_2=`Gene B Chr`, Prada_Nar_2018_Strand_Fusion_2=`Gene B Strand`,
                              Prada_Nar_2018_Split_Reads=`Discordant Read Pairs`,
                              Prada_Nar_2018_Spanning_Reads=`Junction Spanning Reads`,
                              Prada_Nar_2018_Reading_Frame=`Frame Prediction`)

  df_fus <- df_fus %>%
    separate_rows(Junction, `Transcriptional Alleleic Fraction`, sep="\\|") %>%
    separate(Junction, c("Junction_Clean", "Prada_Nar_2018_Spanning_Reads_Bis"), sep=",") %>%
    separate(Junction_Clean, c("Junction_1", "Junction_2"), sep="_") %>%
    separate(Junction_1, c("Gene_1_Bis", "Chr_1_Bis", "Breakpoint_1"), sep=":") %>%
    separate(Junction_2, c("Gene_2_Bis", "Chr_2_Bis", "Breakpoint_2"), sep=":")

  df_fus <- df_fus %>% mutate(Prada_Nar_2018_Total_Reads=Prada_Nar_2018_Split_Reads+Prada_Nar_2018_Spanning_Reads)

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["Prada_Nar_2018_Call"] <- 1
  df_fus[cols]
}

## Squid

preprocess_fusions_squid <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Chr_1", "Breakpoint_1", "SQ_Strand_Fusion_1",
            "Gene_2", "Chr_2", "Breakpoint_2", "SQ_Strand_Fusion_2", "SQ_Call")

  df_fus <- df_fus %>% filter(Type=="fusion-gene") %>%
    rename(Chr_1=`# chrom1`, Chr_2=chrom2, SQ_Strand_Fusion_1=strand1, SQ_Strand_Fusion_2=strand2)

  df_fus <- df_fus %>% 
    mutate(Breakpoint_1=ifelse(SQ_Strand_Fusion_1=="-", start1, end1)) %>%
    mutate(Breakpoint_2=ifelse(SQ_Strand_Fusion_2=="-", start2, end2)) %>%
    mutate(Multiple_Gene_Pairs=ifelse(grepl(",", FusedGenes), "Yes", "No")) %>%
    separate_rows(FusedGenes, sep=",") %>%
    separate(FusedGenes, c("Gene_2", "Gene_1"), sep=":")

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)

  df_fus["SQ_Call"] <- 1
  df_fus[cols]
}

## Star-Fusion

mutate_starfusion_annots <- function(annots){
  regex_annots <- "[0-9\\+\\-\\:[A-Z][a-z]\\_\\[\\]\\.]+"
  splits <- unlist(lapply(strsplit(annots, "^\\[|,|\\]$"), function(s) str_extract(s, regex_annots)))
  splits <- splits[!is.na(splits)]
  paste(splits, collapse="|")
}


extract_starfusion_fusion_type <- function(annots){
  regex_fusiontype <- "(?=INTRACHROMOSOMAL)[\\w\\[\\:\\.\\]\\|\\-\\+]+$|(?=INTERCHROMOSOMAL)[\\w\\[\\:\\.\\]\\|\\-\\+]+$"
  fusiontype <- str_extract(annots, regex_fusiontype)

  regex_fusiontype <- "(?=\\[).*?(?<=\\])"
  fusiontype <- str_remove_all(fusiontype, regex_fusiontype)

  regex_fusiontype <- "(?=\\:).*?(?<=\\:)"
  fusiontype <- str_remove_all(fusiontype, regex_fusiontype)

  regex_fusiontype <- "\\+|\\-"
  str_remove_all(fusiontype, regex_fusiontype)
}


extract_starfusion_annots <- function(annots){
  regex_fusiontype <- "^[\\w\\[\\:\\.\\]\\|\\-\\+]+(?=INTRACHROMOSOMAL)|^[\\w\\[\\:\\.\\]\\|\\-\\+]+(?=INTERCHROMOSOMAL)"
  annots <- str_extract(annots, regex_fusiontype)
  str_remove(annots, "\\|$")
}


preprocess_fusions_starfusion <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Gene_Id_1", "Chr_1", "Breakpoint_1", "Strand_1",
            "Gene_2", "Gene_Id_2", "Chr_2", "Breakpoint_2", "Strand_2",
            "SF_Split_Reads", "SF_Spanning_Reads", "SF_Total_Reads", "SF_Splice_Type", "SF_Large_Anchor_Support",
            "SF_Fusion_Type", "SF_Annots", "SF_FFPM", "SF_Call")

  df_fus <- df_fus %>% rename(SF_Splice_Type=SpliceType, SF_Large_Anchor_Support=LargeAnchorSupport,
                              SF_Annots=annots, SF_FFPM=FFPM)

  df_fus <- df_fus %>%
    separate(LeftGene, c("Gene_1", "Gene_Id_1"), sep="\\^") %>%
    separate(RightGene, c("Gene_2", "Gene_Id_2"), sep="\\^") %>%
    separate(LeftBreakpoint, c("Chr_1", "Breakpoint_1", "Strand_1"), sep=":") %>%
    separate(RightBreakpoint, c("Chr_2", "Breakpoint_2", "Strand_2"), sep=":") %>%
    mutate(SF_Split_Reads=JunctionReadCount, SF_Spanning_Reads=SpanningFragCount) %>%
    mutate(SF_Splice_Type=ifelse(SF_Splice_Type=="YES_LDAS", "Yes", "No"))


  df_fus <- df_fus %>% rowwise() %>%
    mutate(SF_Annots=mutate_starfusion_annots(SF_Annots)) %>%
    mutate(SF_Fusion_Type=extract_starfusion_fusion_type(SF_Annots)) %>%
    mutate(SF_Annots=extract_starfusion_annots(SF_Annots))
                            
  df_fus <- df_fus %>% mutate(SF_Total_Reads=SF_Split_Reads+SF_Spanning_Reads)

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["SF_Call"] <- 1
  df_fus[cols]
}


## STARFUSION CELL 2018

preprocess_fusions_starfusion_cell_2018 <- function(df_fus){
  cols <- c("Sample_Id", "Gene_1", "Chr_1", "Breakpoint_1", "Strand_1",
            "Gene_2", "Chr_2", "Breakpoint_2", "Strand_2",
            "Split_Reads", "Spanning_Reads", "Total_Reads", "Starfusion_Cell_2018_Call")

  df_fus <- df_fus %>% mutate(Sample=str_sub(Sample, start=1, end=20))
  df_fus <- df_fus %>% rename(Sample_Id=Sample, Split_Reads=Junction, Spanning_Reads=Spanning)

  df_fus <- df_fus %>%
    separate(Fusion, c("Gene_1", "Gene_2"), sep="--") %>%
    separate(Breakpoint1, c("Chr_1", "Breakpoint_1", "Strand_1"), sep=":") %>%
    separate(Breakpoint2, c("Chr_2", "Breakpoint_2", "Strand_2"), sep=":") %>%
    mutate(Chr_1=gsub("chr", "", Chr_1), Chr_2=gsub("chr", "", Chr_2)) %>%
    mutate(Total_Reads=Split_Reads+Spanning_Reads)

  cols_num <- c("Breakpoint_1", "Breakpoint_2")
  df_fus <- df_fus %>% mutate_at(all_of(cols_num), as.numeric)
  df_fus <- df_fus %>% mutate_at(vars(ends_with("Reads")), as.numeric)

  df_fus["Starfusion_Cell_2018_Call"] <- 1
  df_fus[cols]
}

# Merge

preprocess_fusions <- function(df_fus, algo){
  cat(paste("-preprocessing fusions called by", algo, "...")) 
  if (algo=="arriba"){
    df_out <- preprocess_fusions_arriba(df_fus)
  } else if (algo=="deepest_pnas_2019"){
    df_out <- preprocess_fusions_deepest_pnas_2019(df_fus)
  } else if (algo=="ericscript"){
    df_out <- preprocess_fusions_ericscript(df_fus)
  } else if (algo=="fusioncatcher"){
    df_out <- preprocess_fusions_fusioncatcher(df_fus)
  } else if (algo=="pizzly"){
    df_out <- preprocess_fusions_pizzly(df_fus)
  } else if (algo=="prada_nar_2018"){
    df_out <- preprocess_fusions_prada_nar_2018(df_fus)
  } else if (algo=="squid"){
    df_out <- preprocess_fusions_squid(df_fus)
  } else if (algo=="starfusion"){
    df_out <- preprocess_fusions_starfusion(df_fus)
  } else if (algo=="starfusion_cell_2018"){
    df_out <- preprocess_fusions_starfusion_cell_2018(df_fus)
  } else {
    stop(paste("-unsupported value of algo:", algo))
  }
  cat("done!\n")
  df_out
}


order_fusions_columns <- function(df_fus){
  cols_first <- c("Sample_Id", "Fusion_Id",
                  "Gene_1", "Gene_Id_1", "Chr_1", "Breakpoint_1", "Strand_1", 
                  "Gene_2", "Gene_Id_2", "Chr_2", "Breakpoint_2", "Strand_2")
  cols_first <- intersect(cols_first, colnames(df_fus))

  cols_algos <- colnames(df_fus)[grepl("_Call$", colnames(df_fus))]
  cols_algos <- intersect(cols_algos, colnames(df_fus))
  df_fus[cols_algos] <- df_fus[cols_algos] %>% replace(is.na(.), 0)

  cols_other <- setdiff(colnames(df_fus), c(cols_first, cols_algos))

  df_fus[c(cols_first, cols_algos, cols_other)]
}


getmode <- function(v) {
   uniqv <- unique(v[!is.na(v)])
   uniqv[which.max(tabulate(match(v, uniqv)))]
}


extract_harmonized_column <- function(df, col){
  prefixes <- c("AR", "ES", "FC", "PZ", "SQ", "SF")
  cols_w_prefix <- paste(prefixes, col, sep="_")
  cols_w_prefix <- intersect(colnames(df), cols_w_prefix)

  df[[col]] <- apply(df[cols_w_prefix], 1, getmode)
  df[cols_w_prefix] <- NULL

  df
}


harmonize_column <- function(df, col_a, col_b){
  cols_order <- colnames(df)
  df_ab <- df %>% select(all_of(c(col_a, col_b))) %>% distinct() %>% filter(!is.na(.data[[col_b]]))

  if (length(unique(df_ab[[col_a]]))!=nrow(df_ab)){
    col_a_dup <- df_ab %>% group_by(.data[[col_a]]) %>% filter(n()>=2) %>% pull(.data[[col_a]])
    col_a_dup <- unique(col_a_dup)
    cat(paste("-WARNING: the following values of", col_a, "have more >=2 different values of", col_b, "\n"))
    cat(paste0("\t", paste0(col_a_dup, collapse="\n\t")), "\n")

    df <- df %>% mutate(Row_Order=1:nrow(df))
    df_wt_dup <- df %>% filter(.data[[col_a]] %in% col_a_dup)
    df_wo_dup <- df %>% filter(!.data[[col_a]] %in% col_a_dup)
    df_ab_wo_dup <- df_ab %>% filter(!.data[[col_a]] %in% col_a_dup)
    df_wo_dup[[col_b]] <- NULL
    df_wo_dup <- left_join(df_wo_dup, df_ab_wo_dup, by=col_a)

    df <- bind_rows(df_wt_dup, df_wo_dup) %>% arrange(Row_Order) %>% select(-Row_Order)
  } else {
    df[[col_b]] <- NULL
    df <- left_join(df, df_ab, by=col_a)
  }

  df[cols_order]
}


update_gene_symbol <- function(study, df_fus){
  filepath_tmp <- paste0("fusions_", study, ".tsv")
  filepath_upd <- paste0("fusions_", study, "_updated.tsv")
  cwd <- getwd()

  for (i in c(1,2)){
    write.table(df_fus, filepath_tmp, quote=F, row.names=F, sep="\t")
    col_gen <- paste0("Gene_", i)
    col_gid <- paste0("Gene_Id_", i)
    col_chr <- paste0("Chr_", i)

    cat(paste("-updating gene symbol", col_gen, "...\n"))

    if (!col_gid %in% colnames(df_fus)){
      col_gid <- "None"
    }
    if (!col_chr %in% colnames(df_fus)){
      col_chr <- "None"
    }

    script <- paste0("update_", study, ".sh")
    sink(script)
    cat("#!/bin/bash\n")
    cat("source ~/.bash_profile\n")
    cat("source activate pyMetaPrism\n")
    cat(paste("python",  file.path(cwd, "workflow/scripts/update_gene_symbols.py"), "\\\n",
              paste("--input", filepath_tmp, "\\\n"),
              paste("--mode", "None",  "\\\n"),
              paste("--gene_name", col_gen, "\\\n"),
              paste("--gene_ensembl_id", col_gid,  "\\\n"),
              paste("--chr_name", col_chr,  "\\\n"),
              paste("--hgnc", args$hgnc, "\\\n"),
              paste("--gencode", args$gencode, "\\\n"),
              paste("--output", filepath_upd, "\n")))
    cat("conda deactivate\n")
    sink()
    oe <- system(paste("/bin/bash", script), intern=T)
    cat(oe, "\n")
    df_fus <- read_tsv(filepath_upd, show_col_types=F, progress=F)[colnames(df_fus)]
  }

  file.remove(filepath_tmp)
  file.remove(filepath_upd)
  file.remove(script)

  df_fus
}


max_combiner <- function(x){
  if (all(is.na(x))){
    return(NA)
  } else {
    return(max(x, na.rm=T))
  }
}


check_unique_row_id <- function(df, cols_row_id, msg){
  df <- df %>% unite(Row_Id, cols_row_id, sep="--", remove=F)
  if (length(unique(df$Row_Id)) < nrow(df)){
    print(msg)
    print("-duplicates are left, they will be removed by combining duplicated rows")

    df_dup <- df %>% group_by(Row_Id) %>% filter(n()>1) %>% arrange(Row_Id) %>% ungroup()
    df_no_dup <- df %>% group_by(Row_Id) %>% filter(n()==1) %>% ungroup()
    print(paste("-number of duplicated rows:", length(unique(df_dup$Row_Id))))

    df_dup <- df_dup %>% group_by(Row_Id) %>% summarise_all(~max_combiner(.x)) %>% ungroup()
    df <- bind_rows(df_no_dup, df_dup)
  }

  df
}


join_iteratively <- function(dfs, names, cols_row_id){
  if (length(names)>0){
    df <- dfs[[names[1]]]
    if (length(names)>1){
      for (name in names[2:length(names)]){
        cols_com <- intersect(colnames(df), colnames(dfs[[name]]))
        df <- full_join(df, dfs[[name]], by=cols_com)

        # check that no duplicates are left
        msg <- paste("-WARNING: merging of", paste0(names[1:(which(name==names)-1)], collapse=";"), "and", name,
                     "introduced duplicates")
        df <- check_unique_row_id(df, cols_row_id, msg)
      }
    }
  } else {
    df <- data.frame()
  }
  
  df
}

ids_1_not_2 <- c('TCGA-CH-5744-01A-11R--NT5DC1--RDH11--116101023--67691244',
 'TCGA-C5-A7X5-01A-11R--CYRIB--ASAP1--129939608--130152805',
 'TCGA-D5-6929-01A-31R--DNM1L--TSPAN11--32679465--30953981',
 'TCGA-05-4420-01A-01R--FOXK2--NARF--82520307--82472564',
 'TCGA-4G-AAZO-01A-12R--CLASP1--FBXO11--121605701--47820456',
 'TCGA-3X-AAVE-01A-11R--PLCXD2--TMPRSS7--111732703--112038072',
 'TCGA-EJ-5518-01A-01R--SUCO--RASAL2--172533497--178390100',
 'TCGA-EJ-5518-01A-01R--SUCO--RABGAP1L--172551626--174231145',
 'TCGA-EJ-5518-01A-01R--CHRM3--SYT14--239545749--210094322',
 'TCGA-GV-A3QK-01B-11R--KDM2A--GSTP1--67121358--67583580',
 'TCGA-GV-A3QK-01B-11R--DHRS2--GSTM4--23630235--109671287',
 'TCGA-C8-A134-01A-11R--RHOBTB1--ATP6V1C2--60943971--10768719',
 'TCGA-2F-A9KP-01A-11R--ABHD17C--PRAMEF6--80696019--12942540',
 'TCGA-DC-6157-01A-11R--NAPB--CST9L--23402993--23566087',
 'TCGA-C5-A7X5-01A-11R--MADD--DDB2--47296055--47234573',
 'TCGA-W5-AA2U-01A-11R--OLA1--CAB39--174123595--230798729',
 'TCGA-CH-5744-01A-11R--ECSIT--RGL3--11513056--11411444',
 'TCGA-GV-A3QK-01B-11R--KDM2A--GSTP1--67121358--67584134',
 'TCGA-AG-3902-01A-01R--ACACA--MED1--37358999--39432016',
 'TCGA-B6-A0I1-01A-11R--TIMM50--IFNL1--39483156--39296805',
 'TCGA-W5-AA39-01A-11R--KIAA1671--NF2--25112462--29639090',
 'TCGA-3X-AAVC-01A-21R--ZNF311--TTBK1--28998734--43282727',
 'TCGA-EJ-5518-01A-01R--TOR1AIP2--IFI16--179877239--159049432',
 'TCGA-86-8585-01A-11R--EXOC4--SUPT3H--133480138--44829876',
 'TCGA-IQ-7630-01A-11R--PSAP--ASCC1--71851182--72213331',
 'TCGA-DM-A28F-01A-11R--HECTD4--TMEM116--112264084--112005303',
 'TCGA-W5-AA2U-01A-11R--NCOA4--RSU1--46023330--16593496',
 'TCGA-86-8585-01A-11R--FRS2--OS9--69470530--57715760',
 'TCGA-EJ-5518-01A-01R--CEP350--GPR161--179955142--168114549',
 'TCGA-86-8585-01A-11R--FRS2--OS9--69470530--57695780',
 'TCGA-EJ-5518-01A-01R--RASAL2--RABGAP1L--178390206--174241483',
 'TCGA-C5-A7X5-01A-11R--KIZ--SLPI--21136552--45253733',
 'TCGA-AG-3902-01A-01R--ACACA--MED1--37358999--39431188',
 'TCGA-4G-AAZO-01A-12R--RASAL2--COPA--178283691--160340294',
 'TCGA-MP-A4T9-01A-11R--TNRC18--TGS1--5370375--55824581',
 'TCGA-AF-5654-01A-01R--GATB--SH3D19--151758772--151226086',
 'TCGA-K4-A5RH-01A-11R--MPRIP--TOM1L2--17042971--17848859',
 'TCGA-W5-AA2U-01A-11R--RSU1--LARP4B--16817315--817889',
 'TCGA-BI-A0VS-01A-11R--FRMD6--TMX1--51489420--51243856',
 'TCGA-KQ-A41O-01A-12R--ARAP1--SLC22A8--72693325--63014962')


aggregate_tables_callers <- function(study, algos=c("arriba", "ericscript", "starfusion")){
  # preprocess per algo
  dfs_fus <- list()
  for (algo in algos){
    cat(paste("-processing algo...", algo, "\n"))
    df_fus <- load_rna_fus(study=study, mode=algo, progress=F)
    df_fus <- preprocess_fusions(df_fus, algo)
    df_fus <- update_gene_symbol(study, df_fus)
    dfs_fus[[algo]] <- df_fus %>% distinct()
  }
  
  # merge all algo calls
  cols_row_id <- c("Sample_Id", "Gene_1", "Gene_2", "Breakpoint_1", "Breakpoint_2")
  cols_gid <- c("Gene_Id_1", "Gene_Id_2")
  cols_brk <- c("Breakpoint_1", "Breakpoint_2")

  # split algos with breakpoints from algos withouth breakpoints
  # algos without breakpoints must be merged in the last step
  algos_wt_breakpoint <- c()
  algos_wo_breakpoint <- c()
  for (algo in algos){
    if (all(cols_brk %in% colnames(dfs_fus[[algo]]))){
      algos_wt_breakpoint <- c(algos_wt_breakpoint, algo)
    } else {
      algos_wo_breakpoint <- c(algos_wo_breakpoint, algo)
    }
  }

  # split algos with gene ids from algos withouth gene ids 
  # algos without gene ids must be merged in the last step
  algos_wt_gene_id <- c()
  algos_wo_gene_id <- c()

  for (algo in algos){
    if (all(cols_gid %in% colnames(dfs_fus[[algo]]))){
      algos_wt_gene_id <- c(algos_wt_gene_id, algo)
    } else {
      algos_wo_gene_id <- c(algos_wo_gene_id, algo)
    }
  }

  # reorder algos to put algos wo breakpoint at the end
  algos_wt_gene_id <- c(setdiff(algos_wt_gene_id, algos_wo_breakpoint), intersect(algos_wt_gene_id, algos_wo_breakpoint))
  algos_wo_gene_id <- c(setdiff(algos_wo_gene_id, algos_wo_breakpoint), intersect(algos_wo_gene_id, algos_wo_breakpoint))
  
  df_fus_wt_gene_id <- join_iteratively(dfs_fus, algos_wt_gene_id, cols_row_id)
  df_fus_wo_gene_id <- join_iteratively(dfs_fus, algos_wo_gene_id, cols_row_id)
  dfs_fus_gene_id <- list(wt_gene_id=df_fus_wt_gene_id, wo_gene_id=df_fus_wo_gene_id)
  df_fus <- join_iteratively(dfs_fus_gene_id, names(dfs_fus_gene_id), cols_row_id)

  # add Fusion_Id
  df_fus$Row_Id <- NULL
  df_fus <- df_fus %>% unite(Fusion_Id, c("Gene_1", "Gene_2"), sep="--", remove=F)

  # Harmonize
  if ("Gene_Id_1" %in% colnames(df_fus)) df_fus <- harmonize_column(df_fus, col_a="Gene_1", col_b="Gene_Id_1")
  if ("Gene_Id_2" %in% colnames(df_fus)) df_fus <- harmonize_column(df_fus, col_a="Gene_2", col_b="Gene_Id_2")

  # order columns
  df_fus <- order_fusions_columns(df_fus)

  df_fus
}


main <- function(args){
  df_fus <- aggregate_tables_callers(study=args$cohort, algos=args$algos)
  if (grepl(".gz$", args$output)){
    write.table(df_fus, gsub(".gz$", "", args$output), quote=F, row.names=F, sep="\t")
    system(paste("gzip",gsub(".gz$", "", args$output)))
  } else {
    write.table(df_fus, args$output, quote=F, row.names=F, sep="\t")
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregated fusion calls for one cohort across different callers.')
  parser$add_argument("-c", "--cohort", type="character", default="tcga_validation", help="Name of the cohort.")
  parser$add_argument("-a", "--algos", nargs="*", help="Name of the cohort.",
                      default=c("arriba", "ericscript", "fusioncatcher", "pizzly", "squid", "starfusion"))
  parser$add_argument("--hgnc", type="character", help="Path to hgnc table",
                      default="../../../data/resources/hgnc/hgnc_all_symbols_03012022.tsv")
  parser$add_argument('--gencode', type="character", help='Path to gencode table.',
                      default="../../../data/resources/gencode/gencode.v27.annotation_genes.tsv")
  parser$add_argument("-o", "--output", type="character", help="Path where output fusion table will be saved.")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")
  print(args)

  main(args)
}
