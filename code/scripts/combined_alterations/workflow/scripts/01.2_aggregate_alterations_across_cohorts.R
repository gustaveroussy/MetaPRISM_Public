# @created: 20 Sep 21
# @modified: 17 Jun 22
# @authors: Yoann Pradat
#
# Aggregate alterations(copy-number, mutations and fusions) annotated with OncoKb and CiVIC annotations across all
# cohorts and group alterations with compatible ESCAT levels across tumor types. Grouped alterations will be used
# as rows in the heatmaps.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(httr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(rprismtools))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(tidyr))

# functions ============================================================================================================

combine_alteration_profiles <- function(df_alt_1, df_alt_2, col_level){
  df_alt_1 <- df_alt_1 %>% rename(Alt_1=.data[[col_level]])
  df_alt_2 <- df_alt_2 %>% rename(Alt_2=.data[[col_level]])
  df_alt <- full_join(df_alt_1, df_alt_2, by="Tumor_Type")
  df_alt %>% mutate(!!col_level:=ifelse(is.na(Alt_1), Alt_2, Alt_1)) %>%
    select(Tumor_Type, .data[[col_level]])
}


get_alteration_profile <- function(df_agg_cat, alteration, col_level, tumor_types){
  df_alt <- df_agg_cat %>% filter(Alteration==alteration) %>%
    distinct(Tumor_Type, .data[[col_level]]) %>% select(Tumor_Type, .data[[col_level]])
  df_alt_mis <- data.frame(Tumor_Type=setdiff(tumor_types, df_alt$Tumor_Type))
  bind_rows(df_alt, df_alt_mis) %>% arrange(Tumor_Type)
}


are_compatible_profiles <- function(df_alt_1, df_alt_2, col_level){
  df_alt_1 <- df_alt_1 %>% rename(Alt_1=.data[[col_level]])
  df_alt_2 <- df_alt_2 %>% rename(Alt_2=.data[[col_level]])
  df_alt <- full_join(df_alt_1, df_alt_2, by="Tumor_Type")
  n_row_dis <- df_alt %>% filter(!is.na(Alt_1), !is.na(Alt_2), Alt_1!=Alt_2) %>% nrow()
  n_row_dis == 0
}

get_profile_counts <- function(profile, col_level, tiers){
  counts <- c()
  for (tier in tiers){
    counts <- c(counts, sum(profile[[col_level]]==tier, na.rm=T))
  }
  counts
}

is_lower_lexicographically <- function(a, b){
  stopifnot(length(a)==length(b))
  n <- length(a)
  for (i in 1:n){
    if (a[i]<b[i]) return(T)
    else if (a[i]>b[i]) return(F)
  }

  return(F)
}


sort_lexicographically <- function(vecs, decreasing=T){
  if (length(vecs)==1){
    return(vecs)
  } else {
    names_vecs <- names(vecs)
    has_changed <- F
    for (i in 2:length(vecs)){
      if (decreasing){
        if (is_lower_lexicographically(vecs[[names_vecs[i-1]]], vecs[[names_vecs[i]]])){
          has_changed <- T
          name_temp <- names_vecs[i]
          names_vecs[i] <- names_vecs[i-1]
          names_vecs[i-1] <- name_temp
        }
      } else {
        if (is_lower_lexicographically(vecs[[names_vecs[i]]], vecs[[names_vecs[i-1]]])){
          has_changed <- T
          name_temp <- names_vecs[i]
          names_vecs[i] <- names_vecs[i-1]
          names_vecs[i-1] <- name_temp
        }
      }
    }
    if (!has_changed){
      return(vecs[names_vecs])
    } else {
      sort_lexicographically(vecs[names_vecs], decreasing)
    }
  }
}


check_level_uniformity <- function(df_agg_sub, col_level){
  df_count <- df_agg_sub %>% distinct(Tumor_Type, .data[[col_level]]) %>%
    group_by(Tumor_Type) %>%
    summarize(Count=n(), .groups="keep")

  mean(df_count$Count==1)==1
}


choose_one_direction_per_fusion <- function(df_agg){
  df_fus <- df_agg %>% filter(Alteration_Category=="Fusion")
  df_oth <- df_agg %>% filter(Alteration_Category!="Fusion")
  df_fus[["Fusion_Id_Ord"]] <- sapply(df_fus[["Hugo_Symbol"]],
                                     function(x) paste(sort(unlist(str_split(x, "--"))), collapse="--"))
  n_row_bef <- nrow(df_fus)
  df_fus <- df_fus %>% unite("ESCAT", Sen_Level_Simple, Res_Level_Simple, sep=",", remove=F)
  df_fus <- df_fus %>% arrange(Fusion_Id_Ord) %>% group_by(Sample_Id, Fusion_Id_Ord) %>%
    slice_min(ESCAT, with_ties=F) %>% ungroup()
  df_fus <- df_fus %>% select(-Fusion_Id_Ord, -ESCAT)
  n_row_aft <- nrow(df_fus)
  cat(paste("-INFO: selected", paste0(n_row_aft, "/", n_row_bef), "fusions by selecting one direction per sample",
            "and per fusion\n"))
  bind_rows(df_oth, df_fus)
}


choose_one_gene_per_fusion <- function(df_agg, col_level){
  df_agg_nfus <- df_agg %>% filter(Alteration_Category!="Fusion") %>%
    mutate(Hugo_Symbol_Original=Hugo_Symbol)
  df_agg_fus <- df_agg %>% filter(Alteration_Category=="Fusion")
  df_agg_fus <- df_agg_fus %>% separate(Hugo_Symbol, c("Gene_1", "Gene_2"), sep="--", remove=F)

  regex <- "^AS(?<=-))|(?=-)AS$"
  genes_all <- union(df_agg_fus$Gene_1, df_agg_fus$Gene_2)
  genes_occ <- sapply(genes_all, function(g) df_agg_fus %>% filter(Gene_1==g|Gene_2==g) %>% nrow())
  df_genes_occ <- data.frame(Occurrence=genes_occ) %>% rownames_to_column("Gene")
  df_genes_occ <- as_tibble(df_genes_occ) %>% arrange(desc(Occurrence))

  df_agg_fus <- df_agg_fus %>% mutate(Hugo_Symbol_Final=NA)

  # special cases
  # BRCA1--NF1 in PRAD is Tier1 sensitivity assign it to BRCA1 and not NF1
  mask <- df_agg_fus$Hugo_Symbol=="BRCA1--NF1" & df_agg_fus$Tumor_Type=="PRAD"
  df_agg_fus[mask, "Hugo_Symbol_Final"] <- "BRCA1"

  for (gene in df_genes_occ$Gene){
    df_agg_sub <- df_agg_fus %>% filter(grepl(gene, Hugo_Symbol), is.na(Hugo_Symbol_Final))

    if (nrow(df_agg_sub) > 0){
      df_agg_fus <- df_agg_fus %>%
        mutate(Hugo_Symbol_Final=ifelse(is.na(Hugo_Symbol_Final)&(Gene_1==gene|Gene_2==gene),gene,Hugo_Symbol_Final))
      
      # If the fusions involving the gene `gene` are not all at the same level of annotation within at least one 
      # tumor type, it is necessary to consider it as a separate event from other fusions involving this gene.
      # Fusions involving this gene will be scattered on different lines in the heatmap plot according to the level of
      # annotation. For the fusions with the best level annotation, the alteration recorded is the full name of
      # the fusion with both parnters.
      if (!check_level_uniformity(df_agg_sub, col_level)){
        # only preserve fusion name in fusions having different annotation level within the same tumor type

        # df_dis <- df_agg_sub %>% group_by(Tumor_Type, .data[[col_level]]) %>% summarize(n=n(), .groups="keep") %>%
        #   group_by(Tumor_Type) %>% filter(n()>1)
        # tumor_types_dis <- df_dis %>% distinct(Tumor_Type) %>% pull()
        # for (tumor_type in tumor_types_dis){
        #   levels_all <- df_dis %>% filter(Tumor_Type==tumor_type) %>% distinct(.data[[col_level]]) %>% pull()
        #   levels_but_last <- sort(levels_all)[1:length(levels_all)-1]

        #   for (level in levels_but_last){
        #     fusions_tumor_type_level <- df_agg_fus %>%
        #       filter(Tumor_Type==tumor_type, .data[[col_level]]==level, (Gene_1==gene|Gene_2==gene)) %>%
        #       pull(Hugo_Symbol)
        #     df_agg_fus <- df_agg_fus %>%
        #       mutate(Alteration=ifelse(Hugo_Symbol %in% fusions_tumor_type_level, paste("Fus", Hugo_Symbol), Alteration))
        #   }
        # }

        # or preserve fusion name of all fusions involving this gene if at least one tumor type has different annotation
        # level
        df_agg_fus <- df_agg_fus %>% 
          mutate(Alteration=ifelse(Gene_1==gene|Gene_2==gene, paste("Fus", Hugo_Symbol), Alteration))
      }
    }
  }

  # fill na if left
  df_agg_fus <- df_agg_fus %>% mutate(Hugo_Symbol_Final=ifelse(is.na(Hugo_Symbol_Final), Gene_1, Hugo_Symbol_Final))
  df_agg_fus <- df_agg_fus %>% select(-Gene_1,-Gene_2) %>% 
    rename(Hugo_Symbol_Original=Hugo_Symbol, Hugo_Symbol=Hugo_Symbol_Final)

  bind_rows(df_agg_nfus, df_agg_fus)
}


collapse_prefix_names <- function(alteration_group_name, prefix){
  name_subparts <- unlist(str_split(alteration_group_name, "\\|"))
  mask_prefix <- grepl(paste0("^", prefix), name_subparts)
  name_subparts_mask <- paste0(prefix, paste0(gsub(prefix, "", name_subparts[mask_prefix]), collapse="/"))
  name_subparts_not_mask <- paste0(name_subparts[!mask_prefix], collapse="|")
  alteration_group_name <- paste0(c(name_subparts_mask, name_subparts_not_mask), collapse="|")

  alteration_group_name
}


add_alteration_group_name <- function(df_agg_sub_nna, alteration_groups){
  if (length(alteration_groups)==1){
    # 1 group = 1 line for this gene
    df_agg_sub_nna <- df_agg_sub_nna %>% mutate(Alteration_Group_Name=Hugo_Symbol)
  } else {
    df_agg_sub_nna <- df_agg_sub_nna %>% mutate(Alteration_Simple=NA, Alteration_Group=NA)

    # add column stating to which ordered group each alteration belongs
    for (i in 1:length(alteration_groups)){
      df_agg_sub_nna <- df_agg_sub_nna %>% 
        mutate(Alteration_Group=ifelse(Alteration %in% alteration_groups[[i]], i, Alteration_Group))
    }

    # simplify categories. Short insertions and deletions will be categorized as mut
    df_agg_sub_nna <- df_agg_sub_nna %>% 
      mutate(Alteration_Category_Original=Alteration_Category) %>%
      mutate(Alteration_Category=ifelse(Alteration_Category %in% c("Del", "Ins"), "Mut", Alteration_Category)) %>%
      mutate(Alteration_Category=ifelse(Alteration_Category=="Fusion", "Fus", Alteration_Category)) %>%
      group_by(Alteration_Group) %>% mutate(Group_Size=n()) %>% ungroup()

    for (alteration_category in df_agg_sub_nna$Alteration_Category){
      # apply specific rules
      # - if any of the group contains all "Mut", mut alterations become "Mut"
      # - if more than 2 groups contains "Mut", mut alterations from last one become "Mut Other" unless it contains only
      # 1 alteration
      # - apply same rule for Other
      df_cat_counts <- df_agg_sub_nna %>% distinct(Alteration, Alteration_Category, Alteration_Group) %>% 
        filter(Alteration_Category==alteration_category) %>% group_by(Alteration_Group) %>% summarize(n=n())

      if (nrow(df_cat_counts)==1){
        cat_group <- df_cat_counts %>% pull(Alteration_Group)
        df_agg_sub_nna <- df_agg_sub_nna %>% 
          mutate(Alteration_Simple=ifelse(Alteration_Category==alteration_category & Alteration_Group==cat_group,
                                          alteration_category,Alteration_Simple))
      } else if (nrow(df_cat_counts)>=2) {
        # consider last group
        cat_group <- max(df_cat_counts %>% pull(Alteration_Group))
        df_agg_sub_nna <- df_agg_sub_nna %>% 
          mutate(Alteration_Simple=ifelse(Alteration_Category==alteration_category & Alteration_Group==cat_group &
                                          Group_Size > 1, paste(alteration_category, "Other"),Alteration_Simple))
      }
    }

    # build Alteration_Group_Name
    df_agg_sub_nna <- df_agg_sub_nna %>% 
      mutate(Alteration_Simple=ifelse(is.na(Alteration_Simple),Alteration,Alteration_Simple)) %>%
      mutate(Alteration_Simple=ifelse(Alteration_Simple=="Fusion","Fus",Alteration_Simple)) %>%
      mutate(Alteration_Simple=ifelse(Alteration_Simple=="Amplification","Amp",Alteration_Simple)) %>%
      mutate(Alteration_Simple=ifelse(Alteration_Simple=="Deletion","Del",Alteration_Simple)) %>%
      mutate(Alteration_Simple=gsub("Exon", "Ex.", Alteration_Simple)) %>%
      select(-Group_Size)

    # try to sort alterations by increasing amino acid number and then alphabetically if possible
    df_agg_sub_nna["AA_Pos"] <- as.numeric(str_extract(df_agg_sub_nna$Alteration_Simple, "(?<=^[A-Z])([0-9]+)"))
    df_group_name <- df_agg_sub_nna %>% arrange(Alteration_Simple) %>% arrange(AA_Pos) %>%
      distinct(Alteration_Simple, Alteration_Group) %>%  group_by(Alteration_Group) %>%
      summarize(Alteration_Group_Name=paste(Alteration_Simple, collapse="|"))
    df_agg_sub_nna <- df_agg_sub_nna %>% select(-AA_Pos)

    # special groups names for some genes
    hugo_symbol <- unique(df_agg_sub_nna$Hugo_Symbol)
    alteration_group_names <- df_group_name$Alteration_Group_Name

    if (hugo_symbol=="KRAS"){
      alteration_group_names_new <- c()
      for (alteration_group_name in alteration_group_names){
        for (prefix in c("G13", "G12")){
          if (grepl(prefix, alteration_group_name) & alteration_group_name!="G12C"){
            alteration_group_name <- collapse_prefix_names(alteration_group_name, prefix)
          }
        }
        alteration_group_names_new <- c(alteration_group_names_new, alteration_group_name)
      }
    } else if (hugo_symbol=="EGFR"){
      alteration_group_names_new <- c()
      for (alteration_group_name in alteration_group_names){
        for (prefix in c("G719")){
          if (grepl(prefix, alteration_group_name)){
            alteration_group_name <- collapse_prefix_names(alteration_group_name, prefix)
          }
        }
        alteration_group_names_new <- c(alteration_group_names_new, alteration_group_name)
      }
    } else {
      alteration_group_names_new <- alteration_group_names
    }
    
    df_group_name$Alteration_Group_Name <- alteration_group_names_new

    # if multiple cat Other in last group, rename to hugo_symbol Other
    group <- max(df_agg_sub_nna %>% pull(Alteration_Group))
    group_name <- df_group_name %>% filter(Alteration_Group==cat_group) %>% pull(Alteration_Group_Name)
    if (str_count(group_name, "Other") > 1){
      df_group_name <- df_group_name %>% 
        mutate(Alteration_Group_Name=ifelse(Alteration_Group==group, "Other", Alteration_Group_Name))
    }

    df_agg_sub_nna <- left_join(df_agg_sub_nna, df_group_name, by="Alteration_Group")
    df_agg_sub_nna <- df_agg_sub_nna %>% select(-Alteration_Simple, -Alteration_Group, -Alteration_Category)
    df_agg_sub_nna <- df_agg_sub_nna %>% rename(Alteration_Category=Alteration_Category_Original)
  }

  df_agg_sub_nna
}


get_alterations_groups <- function(alterations, alteration_profiles, col_level){
  if (length(alterations)==0){
    alteration_groups <- list()
    alteration_groups_profiles <- list()
  } else {
    alteration_groups <- setNames(list(alterations[1]), alterations[1])
    alteration_groups_profiles <- setNames(list(alteration_profiles[[alterations[1]]]), alterations[1])

    if (length(alterations)>1){
      for (alt_2 in alterations[2:length(alterations)]){
        alt_2_grouped <- F
        for (alt_1 in names(alteration_groups)){
          if (are_compatible_profiles(alteration_groups_profiles[[alt_1]], alteration_profiles[[alt_2]], col_level)){
            alteration_groups[[alt_1]] <- c(alteration_groups[[alt_1]], alt_2)
            alteration_groups_profiles[[alt_1]] <- combine_alteration_profiles(alteration_groups_profiles[[alt_1]],
                                                                               alteration_profiles[[alt_2]],
                                                                               col_level)
            alt_2_grouped <- T
            break
          }
        }
        if (!alt_2_grouped){
          alteration_groups[[alt_2]] <- alt_2
          alteration_groups_profiles[[alt_2]] <- alteration_profiles[[alt_2]]
        }
      }
    }
  }

  list(groups=alteration_groups, groups_profiles=alteration_groups_profiles)
}


get_alterations_groups_special <- function(alterations, alteration_profiles, col_level, hugo_symbol){
  if (col_level=="Sen_Level_Simple" & hugo_symbol=="KRAS"){
    mask_tier1 <- sapply(alterations, function(x) any(grepl("Tier1", alteration_profiles[[x]][[col_level]])))
    mask_tier2 <- sapply(alterations, function(x) any(grepl("Tier2", alteration_profiles[[x]][[col_level]])))
    mask_tier3_coad <- sapply(alterations, function(x) {
                           lvl <- alteration_profiles[[x]][alteration_profiles[[x]]$Tumor_Type=="COAD",][[col_level]];
                           if(is.na(lvl)){ return(FALSE)} else {return(lvl=="Tier3")}})

    alterations_tier1 <- alterations[mask_tier1 & !mask_tier3_coad]
    alterations_tier2 <- alterations[mask_tier2 & !mask_tier1 & !mask_tier3_coad]
    alterations_other <- alterations[(!mask_tier2 & !mask_tier1) | mask_tier3_coad]

    groups_tier1 <- get_alterations_groups(alterations_tier1, alteration_profiles, col_level)
    groups_tier2 <- get_alterations_groups(alterations_tier2, alteration_profiles, col_level)
    groups_other <- get_alterations_groups(alterations_other, alteration_profiles, col_level)

    groups <- list(groups_tier1, groups_tier2, groups_other)
  } else if (col_level=="Res_Level_Simple" & hugo_symbol == "TP53"){
    alterations_spe <- c("DNA Binding Domain")
    alterations_not_spe <- setdiff(alterations, alterations_spe)

    groups_spe <- get_alterations_groups(alterations_spe, alteration_profiles, col_level)
    groups_not_spe <- get_alterations_groups(alterations_not_spe, alteration_profiles, col_level)

    groups <- list(groups_spe, groups_not_spe)
  } else if (col_level=="Sen_Level_Simple" & hugo_symbol == "EGFR"){
    alterations_a <- c("G719A", "G719C", "S768I", "L861Q")
    alterations_b <- c("Exon 19 Del", "L858R")
    alterations_c <- c("Exon 20 Ins")
    alterations_d <- c("T790M")
    alterations_oth <- setdiff(alterations, c(alterations_a, alterations_b, alterations_c, alterations_d))

    groups_a <- get_alterations_groups(alterations_a, alteration_profiles, col_level)
    groups_b <- get_alterations_groups(alterations_b, alteration_profiles, col_level)
    groups_c <- get_alterations_groups(alterations_c, alteration_profiles, col_level)
    groups_d <- get_alterations_groups(alterations_d, alteration_profiles, col_level)
    groups_oth <- get_alterations_groups(alterations_oth, alteration_profiles, col_level)

    groups <- list(groups_a, groups_b, groups_c, groups_d, groups_oth)
  } else {
    mask_tier1 <- sapply(alterations, function(x) any(grepl("Tier1", alteration_profiles[[x]][[col_level]])))
    alterations_tier1 <- alterations[mask_tier1]
    alterations_other <- alterations[!mask_tier1]

    groups_tier1 <- get_alterations_groups(alterations_tier1, alteration_profiles, col_level)
    groups_other <- get_alterations_groups(alterations_other, alteration_profiles, col_level)

    groups <- list(groups_tier1, groups_other)
  }

  groups
}


group_alterations_according_to_level <- function(df_agg_sub_nna, col_level, tumor_types, hugo_symbol){
  # get list of alterations
  alterations <- df_agg_sub_nna %>%  distinct(Alteration, Alteration_Category) %>% 
    arrange(desc(Alteration_Category)) %>% pull(Alteration)

  # get profile of each alteration
  alteration_profiles <- lapply(alterations, function(x) get_alteration_profile(df_agg=df_agg_sub_nna, 
                                                                                alteration=x,
                                                                                col_level=col_level,
                                                                                tumor_types=tumor_types))
  alteration_profiles <- setNames(alteration_profiles, alterations)

  # sort alterations by Tier
  mask_tier1 <- sapply(alterations, function(x) any(grepl("Tier1", alteration_profiles[[x]][[col_level]])))
  mask_tier2 <- sapply(alterations, function(x) any(grepl("Tier2", alteration_profiles[[x]][[col_level]])))
  mask_tier3 <- sapply(alterations, function(x) any(grepl("Tier3", alteration_profiles[[x]][[col_level]])))
  alterations_tier1 <- alterations[mask_tier1]
  alterations_tier2 <- alterations[mask_tier2 & !mask_tier1]
  alterations_tier3 <- alterations[mask_tier3 & !mask_tier2 & !mask_tier1]
  alterations_notier <- alterations[!mask_tier3 & !mask_tier2 & !mask_tier1]
  alterations <- c(alterations_tier1, alterations_tier2, alterations_tier3, alterations_notier)

  # for some genes, perform some special grouping
  hugo_symbol_special_sen <- c("FGFR3", "KRAS", "EGFR")
  hugo_symbol_special_res <- c("TP53")
  special_sen <- col_level=="Sen_Level_Simple" & hugo_symbol %in% hugo_symbol_special_sen
  special_res <- col_level=="Res_Level_Simple" & hugo_symbol %in% hugo_symbol_special_res

  if (special_sen | special_res){
    groups <- get_alterations_groups_special(alterations, alteration_profiles, col_level, hugo_symbol)
  } else {
    groups <- list(get_alterations_groups(alterations, alteration_profiles, col_level))
  }

  alteration_groups <- do.call(c, lapply(groups, function(x) x[["groups"]]))
  alteration_groups_profiles <- do.call(c, lapply(groups, function(x) x[["groups_profiles"]]))

  # order alteration groups
  tiers <- df_agg_sub_nna[col_level] %>% distinct() %>% pull()
  tiers <- sort(tiers)
  alteration_groups_counts <- lapply(alteration_groups_profiles, get_profile_counts, col_level=col_level, tiers=tiers)
  alteration_groups_counts <- sort_lexicographically(alteration_groups_counts)
  alteration_groups <- alteration_groups[names(alteration_groups_counts)]

  # add alteration groups names
  df_agg_sub_nna <- add_alteration_group_name(df_agg_sub_nna, alteration_groups)
  col_alteration_level <- paste0("Alteration_", col_level)
  df_agg_sub_nna <- df_agg_sub_nna %>% rename(!!col_alteration_level:=Alteration_Group_Name)

  df_agg_sub_nna
}


categorize_alterations_according_to_level_hugo_symbol <- function(df_agg_level, hugo_symbol, col_level){
  tumor_types <- df_agg_level %>% distinct(Tumor_Type) %>% arrange(Tumor_Type) %>% pull()
  df_agg_sub <- df_agg_level %>% filter(Hugo_Symbol==hugo_symbol) %>% arrange(desc(Alteration_Category), .data[[col_level]])
  df_agg_sub <- df_agg_sub %>% mutate(!!col_level:=ifelse(is.na(.data[[col_level]]), "N/A", .data[[col_level]]))
  df_agg_sub <- df_agg_sub %>% distinct(Hugo_Symbol, Alteration, Alteration_Category, Tumor_Type, .data[[col_level]])
  
  df_agg_sub_na <- df_agg_sub %>% filter(.data[[col_level]]=="N/A") %>% 
    mutate(!!paste0("Alteration_", col_level):="Oncogenic, no Level")
  df_agg_sub_nna <- df_agg_sub %>% filter(.data[[col_level]]!="N/A")

  if (nrow(df_agg_sub_nna)>0){
    df_agg_sub_nna <- group_alterations_according_to_level(df_agg_sub_nna, col_level, tumor_types, hugo_symbol)
  }

  df_agg_sub <- bind_rows(df_agg_sub_nna, df_agg_sub_na)
  df_agg_sub <- df_agg_sub %>% mutate(!!col_level:=ifelse(.data[[col_level]]=="N/A", NA, .data[[col_level]]))
  df_agg_sub
}


categorize_alterations_according_to_level <- function(df_agg_level, col_level){
  hugo_symbols <- df_agg_level %>% distinct(Hugo_Symbol) %>% pull()

  df_agg_level_cat <- tibble()
  i <- 0
  cat("-grouping alterations for each gene...")
  for (hugo_symbol in hugo_symbols){
    i <- i+1
    df_agg_level_cat_hs <- categorize_alterations_according_to_level_hugo_symbol(df_agg_level, hugo_symbol, col_level)
    df_agg_level_cat <- bind_rows(df_agg_level_cat, df_agg_level_cat_hs)
    progress(i, length(hugo_symbols))
  }
  cat("done!\n")

  df_agg_level_cat 
}


check_and_clean_incoherences <- function(df_agg, raise_error=F){
  n_unique_combs <- df_agg %>% distinct(Tumor_Type, Hugo_Symbol, Alteration) %>% nrow()
  n_unique_combs_level <- df_agg %>%
    distinct(Tumor_Type, Hugo_Symbol, Alteration, Sen_Level_Simple, Res_Level_Simple)  %>% nrow()

  # df_agg %>%
  #   unite("Combination_Lvl", Tumor_Type, Hugo_Symbol, Alteration, Sen_Level_Simple, Res_Level_Simple, remove=F) %>%
  #   unite("Combination", Tumor_Type, Hugo_Symbol, Alteration, remove=F) %>%
  #   select(Combination, Combination_Lvl) %>% distinct() %>%
  #   group_by(Combination) %>% filter(n()>1) %>% arrange(Combination)

  if (n_unique_combs!=n_unique_combs_level){
    if (!raise_error){
      warning(paste("-there are discrepancies in the annotations!", n_unique_combs, "unique combs vs",
                    n_unique_combs_level, " unique combs with level."))
    } else {
      stop(paste("-there are discrepancies in the annotations!", n_unique_combs, "unique combs vs",
                    n_unique_combs_level, " unique combs with level."))
    }

    df_agg_cat <- df_agg %>%
      select(Tumor_Type, Hugo_Symbol, Alteration, Sen_Level_Simple, Res_Level_Simple)  %>% distinct()
    df_agg_cat <- df_agg_cat %>% mutate(Sen_Level_Simple=ifelse(is.na(Sen_Level_Simple), "Tier4", Sen_Level_Simple)) %>%
      mutate(Res_Level_Simple=ifelse(is.na(Res_Level_Simple), "Tier4", Res_Level_Simple))
    df_agg_cat <- df_agg_cat %>% group_by(Tumor_Type, Hugo_Symbol, Alteration) %>%
      filter(Sen_Level_Simple==min(Sen_Level_Simple), Res_Level_Simple==min(Res_Level_Simple)) %>% ungroup()
    df_agg_cat <- df_agg_cat %>% mutate(Sen_Level_Simple=ifelse(Sen_Level_Simple=="Tier4", NA, Sen_Level_Simple)) %>%
      mutate(Res_Level_Simple=ifelse(Res_Level_Simple=="Tier4", NA, Res_Level_Simple))
    stopifnot(nrow(df_agg_cat)==n_unique_combs)

    df_agg <- df_agg %>% select(-Sen_Level_Simple, -Res_Level_Simple)
    df_agg <- left_join(df_agg , df_agg_cat, by=c("Tumor_Type", "Hugo_Symbol", "Alteration"))
  }
  
  df_agg
}


main <- function(args){
  # load alterations
  dfs_agg <- setNames(lapply(args$tables, function(x) load_table(x)), args$cohorts)
  dfs_agg <- setNames(lapply(args$cohorts, function(x) dfs_agg[[x]] %>% mutate(Cohort=x)), args$cohorts)

  # aggregate cohorts
  df_agg <- bind_rows(dfs_agg)
  df_agg <- df_agg %>% distinct()

  # restrict the table to samples from selected tumor types
  # reason for this is that aggregation results may vary with additional tumor types
  df_tt_count <- load_table(args$counts)
  df_tt_count <- df_tt_count %>% filter(Use_heatmap_all==1)
  tt_keep <- df_tt_count %>% pull(var="Tumor_Type")
  df_agg <- df_agg %>% filter(Tumor_Type %in% tt_keep)

  # for representation purposes, some levels are manually modified.
  # - sometimes, different events are annotated at different ESCAT levels because of differences
  #   in intervals associated to some events in CIViC database from the same gene (ex: FBXW7). This
  #   has no biological meaning.
  # - sometimes, a single event has different ESCAT level. For NF1 for instance, all fusions involving
  #   NF1 are tier 3 except for BRCA1--NF1 which is tier 1 in PRAD. In this case, the fusions should have
  #   been attributed to BRCA1 and not NF1. 

  # TSC1
  mask <- df_agg$Hugo_Symbol=="TSC1" & df_agg$Tumor_Type=="BLCA" & grepl("Del", df_agg$Alteration_Category)
  df_agg[mask, "Sen_Level_Simple"] <- "Tier2"
  # df_agg[mask, "Sen_Drug_Simple"] <- "EVEROLIMUS"
  
  # FBXW7
  mask <- df_agg$Hugo_Symbol=="FBXW7" & grepl("\\*|Exon 12", df_agg$Alteration)
  df_agg[mask, "Sen_Level_Simple"] <- "Tier2"
  # df_agg[mask, "Sen_Drug_Simple"] <- "MTOR INHIBITOR"

  mask <- df_agg$Hugo_Symbol=="FBXW7" & df_agg$Tumor_Type=="BLCA" & grepl("Splice_Site", df_agg$Alteration)
  df_agg[mask, "Sen_Level_Simple"] <- "Tier2"
  # df_agg[mask, "Sen_Drug_Simple"] <- "SIROLIMUS"

  mask <- df_agg$Hugo_Symbol=="FBXW7" & df_agg$Tumor_Type=="COAD" & grepl("\\*|Exon 12", df_agg$Alteration) 
  df_agg[mask, "Res_Level_Simple"] <- "Tier2"
  # df_agg[mask, "Res_Drug_Simple"] <- "PANITUMUMAB|CETUXIMAB"

  # KRAS
  mask <- df_agg$Hugo_Symbol=="KRAS" & df_agg$Alteration=="3'UTR"
  df_agg[mask, "Alteration"] <- "F156L"

  # clean incoherences
  df_agg <- check_and_clean_incoherences(df_agg)

  # for fusions seen in both directions, select only one direction using ESCAT levels
  df_agg <- choose_one_direction_per_fusion(df_agg)

  # add columns Hugo_Symbol_Sen_Level_Simple, Alteration_Sen_Level_Simple and
  # and Hugo_Symbol_Res_Level_Simple, Alteration_Res_Level_Simple giving final
  # names of alterations for the heatmap
  for (col_level in c("Sen_Level_Simple", "Res_Level_Simple")){
    # for fusions, keep only 1 of the partner genes
    cat(paste("-choosing one gene per fusion for annotation", col_level, "..."))
    df_agg_level <- choose_one_gene_per_fusion(df_agg, col_level)
    stopifnot(nrow(df_agg_level)==nrow(df_agg))
    cat("done!\n")

    # categorize alterations according to annotation level across tumor types
    df_agg_level_cat <- categorize_alterations_according_to_level(df_agg_level, col_level)
    df_agg_level_cat <- df_agg_level_cat %>% replace(.=="N/A", NA)

    # sanity checks
    stopifnot(df_agg_level_cat %>% group_by(Tumor_Type, Hugo_Symbol, Alteration_Category, Alteration) %>%
              filter(n()>1) %>% nrow()==0)

    cols_common <- intersect(colnames(df_agg_level), colnames(df_agg_level_cat))
    df_agg_level <- left_join(df_agg_level, df_agg_level_cat, by=cols_common)
    stopifnot(nrow(df_agg_level)==nrow(df_agg))

    df_agg_level <- df_agg_level %>%
      rename(!!paste0("Hugo_Symbol_", col_level):=Hugo_Symbol, Hugo_Symbol=Hugo_Symbol_Original)
    df_agg_level <- df_agg_level %>% 
      mutate(Alteration=ifelse(Alteration_Category=="Fusion", "Fus", Alteration))
    cols_common <- intersect(colnames(df_agg), colnames(df_agg_level))
    df_agg <- left_join(df_agg, df_agg_level, by=cols_common)
  }

  # save
  df_agg <- df_agg %>% arrange(Subject_Id, Alteration_Category)
  write.table(df_agg, file=args$output, sep="\t", quote=F,row.names=F)
  cat(paste("-file saved at", args$output, "\n"))
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Aggregate alterations.')
  parser$add_argument("--cohorts", nargs="+", help="Names of inputs cohorts.",
                      default=c("met500", "prism", "tcga"))
  parser$add_argument("--counts", type="character", help="Path to count table of tumor types.",
                      default="../../../results/combined_alterations/selection/selection_tumor_types.tsv")
  parser$add_argument("--tables", nargs="+", help="Path to tables of aggregated alterations of each cohort.",
                      default=c("../../../results/combined_alterations/alterations/aggregated_alterations_met500.tsv",
                                "../../../results/combined_alterations/alterations/aggregated_alterations_prism.tsv",
                                "../../../results/combined_alterations/alterations/aggregated_alterations_tcga.tsv"))
  parser$add_argument("--output", type="character", help="Path to output aggregated table",
                      default="../../../results/combined_alterations/alterations/aggregated_alterations.tsv")
  parser$add_argument('--log', type="character", help='Path to log file.')
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)

  main(args)
}
