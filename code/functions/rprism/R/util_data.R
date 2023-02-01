#' Function for selecting tumor types.
#'
#' Using the column \code{col_tt}, this function subsets each input dataframe to the values of \code{col_tt} with at
#' least \code{tt_min_size} entries.
#'
#' @param dfs A named list of \code{data.frame}. They must contain columns \code{col_tt}.
#' @param tt_min_size Numeric value.
#' @param tt_drop Values that should be discarded regarless of the size.
#' @return a list containing a named list of dataframes and a named list of character vectors, each with the sames names
#' as \code{dfs}.
#'
#' @importFrom dplyr filter n_distinct
#' @importFrom magrittr %>%
#' @importFrom tidyr unite replace_na
#' @importFrom rlang .data
#'
#' @author Yoann Pradat
#' @export
select_tumor_types <- function(dfs, col_tt, tt_min_size, tt_drop){
  dfs_tt_count <- list()
  tts_keep <- list()

  for (name in names(dfs)){
    df_tt_count <- dfs[[name]] %>%
      mutate(!!col_tt:=replace_na(.data[[col_tt]], "N/A")) %>%
      group_by(.data[[col_tt]]) %>%
      summarize(Count=n_distinct(Subject_Id))

    tts_keep[[name]] <- df_tt_count %>%
      filter(!.data[[col_tt]] %in% tt_drop) %>%
      filter(Count >= tt_min_size) %>%
      pull(.data[[col_tt]])

    dfs_tt_count[[name]] <- df_tt_count %>%
      rename(!!toupper(name) := Count)
  }

  list(dfs_tt_count=dfs_tt_count, tts_keep=tts_keep)
}


#' Add columns
#' 
#' Add annotations from bio and cln files.
#' @param df the dataframe to be annotated
#' @param df_bio the dataframe of biospecimen annotations
#' @param df_cln the dataframe of subjects annotations
#' @param col_on_bio column in the dataframe to be used to connect to df_bio
#' @param col_on_cln column in the dataframe to be used to connect to df_cln
#' @param cols_bio columns from the df_bio dataframe to be added to df
#' @param cols_cln columns from the df_cln dataframe to be added to df
#'
#' @importFrom dplyr distinct
#' @author Yoann Pradat
#' @export
add_annotations_bio_cln <- function(df, df_bio=NULL, df_cln=NULL, col_on_bio=NULL, col_on_cln=NULL, cols_bio=NULL,
                                    cols_cln=NULL){
  if(!is.null(cols_bio)){
    df <- left_join(df, df_bio[union(cols_bio, "Sample_Id")] %>% distinct(), by=setNames("Sample_Id", col_on_bio))
  }

  if(!is.null(cols_cln)){
    df <- left_join(df, df_cln[union(cols_cln, "Subject_Id")] %>% distinct(), by=setNames("Subject_Id", col_on_cln))
  }

  df
}


#' Preprocess mutations table
#'
#' Add bio and/or cln attributes to the mutations table and perform a selection of one pair tumor/normal per subject
#' if requested.
#' 
#' @param df_mut data.frame of mutations
#' @param cohort name of the cohort
#' @param cols_bio (optional) list ist of bio attributes
#' @param cols_cln (optional) list of cln attributes.
#' @param select_pairs (optional) If set to TRUE, 1 pair tumor/normal is selected for each subject.
#' @param selection_mut (optional)
#'   Here are the available modes
#'  \enumerate{
#'    \item{'all'}{no selection of variants}
#'    \item{'annotated'}{no selection of variants}
#'    \item{'non_synonymous'}{select variants for which the value of 'Variant_Classification' is one of
#'             * Frame_Shift_Del
#'             * Frame_Shift_Ins,
#'             * Splice_Site
#'             * Translation_Start_Site
#'             * Nonsense_Mutation
#'             * Nonstop_Mutation
#'             * In_Frame_Del
#'             * In_Frame_Ins
#'             * Missense_Mutation
#'             * Start_Codon_Del
#'             * Start_Codon_SNP
#'             * Stop_Codon_Del
#'             * Stop_Codon_Ins}
#'    \item{'truncating'}{select variants for which the value of 'Variant_Classification' is one of
#'             * Frame_Shift_Del
#'             * Frame_Shift_Ins,
#'             * Splice_Site
#'             * Nonsense_Mutation}
#' }
#' @param verbose (optional) Should info messages be printed?
#' @return A dataframe of mutations with possibly additional attributes and possibly a reduced number of lines.
#'
#' @import dplyr
#' @importFrom tidyr unite
#'
#' @export
preprocess_wes_mut <- function(df_mut, cohort, cols_bio=c(), cols_cln=c(), select_pairs=FALSE, selection_mut="all",
                               verbose=TRUE){

  col_tsb <- "Tumor_Sample_Barcode"
  col_nsb <- "Matched_Norm_Sample_Barcode"
  col_tid <- "Tumor_Sample_Id"
  col_nid <- "Normal_Sample_Id"
  col_id <- "Sample_Id"

  if (all(c(col_tid, col_nid) %in% colnames(df_mut))){
    mut_mode <- "somatic"
    cols_rmv_id <- c()
  } else  if (all(c(col_tsb, col_nsb) %in% colnames(df_mut))){
    mut_mode <- "somatic"
    df_mut[[col_tid]] <- df_mut[[col_tsb]]
    df_mut[[col_nid]] <- df_mut[[col_nsb]]
    cols_rmv_id <- c(col_tid, col_nid)
  } else if (col_nid %in% colnames(df_mut)){
    mut_mode <- "germline"
    cols_rmv_id <- c()
  } else if (col_nsb %in% colnames(df_mut)){
    mut_mode <- "germline"
    df_mut[[col_nid]] <- df_mut[[col_nsb]]
    cols_rmv_id <- c(col_nid)
  } else if (col_id %in% colnames(df_mut)){
    mut_mode <- "germline"
    df_mut[[col_nid]] <- df_mut[[col_id]]
    cols_rmv_id <- c(col_nid)
  } else {
    stop("-unrecognized names for sample identifiers")
  }

  # bio attributes
  cols_bio_min <- c("Sample_Id", "Subject_Id")
  cols_bio_use <- union(cols_bio, cols_bio_min)
  cols_bio_rmv <- setdiff(cols_bio_use, cols_bio)
  df_mut <- df_mut %>% select(-one_of(intersect(colnames(df_mut), cols_bio_use)))
  df_bio <- load_bio(cohort)[cols_bio_use] %>% distinct()

  if (mut_mode=="somatic"){
    by = setNames(col_id, col_tid)
  } else {
    by = setNames(col_id, col_nid) 
  }
  df_mut <- left_join(df_mut, df_bio, by=by, keep=T)

  # cln attributes
  cols_cln_min <- c("Subject_Id", "Sample_Id_DNA_T", "Sample_Id_DNA_N")
  cols_cln_use <- union(cols_cln, cols_cln_min)
  cols_cln_rmv <- setdiff(cols_cln_use, cols_cln)
  df_cln <- load_cln(cohort)[cols_cln_use] %>% distinct()
  df_cln <- df_cln[rowSums(is.na(df_cln[cols_cln_min]))==0,]
  df_mut <- df_mut %>% select(-one_of(setdiff(intersect(colnames(df_mut), cols_cln_use), "Subject_Id")))
  df_mut <- left_join(df_mut, df_cln, by="Subject_Id")

  # normal or tumor/normal pairs selection
  if (select_pairs){
    col_tsid_m <- "Tumor_Sample_Id"
    col_tsid_c <- "Sample_Id_DNA_T"
    col_nsid_m <- "Normal_Sample_Id"
    col_nsid_c <- "Sample_Id_DNA_N"
    col_pid <- "Tumor_vs_Normal"

    if (mut_mode == "somatic"){

      df_mut <- df_mut %>% unite({{col_pid}}, {{col_tsid_m}}, {{col_nsid_m}}, remove=F, sep="_vs_")
      df_cln <- df_cln %>% unite({{col_pid}}, {{col_tsid_c}}, {{col_nsid_c}}, remove=F, sep="_vs_")
      df_mut <- df_mut %>% filter(.data[[col_pid]] %in% df_cln[[col_pid]])

      # info messages
      if (verbose){
        cat(paste("-number of unique pairs: ", length(unique(df_mut[[col_pid]]))), "\n")
        cat(paste("-number of unique subjects: ", length(unique(df_mut[["Subject_Id"]]))), "\n")
      }

      df_mut <- df_mut %>% select(-{{col_pid}})

    } else if (mut_mode == "germline") {
      df_mut <- df_mut %>% filter(.data[[col_nsid_m]] %in% df_cln[[col_nsid_c]])

      if (verbose){
        cat(paste("-number of unique samples: ", length(unique(df_mut[[col_nsid_m]]))), "\n")
        cat(paste("-number of unique subjects: ", length(unique(df_mut[["Subject_Id"]]))), "\n")
      }
    }
  }

  # remove unwanted attributes
  cols_rmv <- setdiff(union(cols_bio_rmv, cols_cln_rmv), union(cols_bio, cols_cln))
  cols_rmv <- union(cols_rmv, cols_rmv_id)
  cols_rmv <- intersect(colnames(df_mut), cols_rmv)
  df_mut <- df_mut %>% select(-one_of(cols_rmv))

  # variants selection
  if (!selection_mut %in% c("all", "annotated")){
    if (selection_mut=="non_synonymous"){
      col_keep <- "Variant_Classification"
      val_keep <- c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site","Splice_Region","Translation_Start_Site",
            "Nonsense_Mutation","Nonstop_Mutation","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
            "Start_Codon_Del","Start_Codon_SNP","Stop_Codon_Del","Stop_Codon_Ins","Translation_Start_Site")
    } else if (selection_mut=="truncating"){
      col_keep <- "Variant_Classification"
      val_keep <- c("Frame_Shift_Del","Frame_Shift_Ins","Splice_Site", "Nonsense_Mutation", "Splice_Site",
                    "Splice_Region", "Nonstop_Mutation", "Stop_Codon_Del", "Stop_Codon_Ins")
    } else {
      stop(paste("-unsupported value",  selection_mut, "for argument 'selection_mut'."))
    }

    df_mut <- df_mut %>% filter(.data[[col_keep]] %in% val_keep)
    if (verbose) cat(paste("-number of selected variants:", nrow(df_mut), "\n"))
  }


  df_mut
}


#' Split targeted therapy classes with multiple targets
#'
#' @description Split multiple-targets therapies into multiple single-target therapies.
#' @param x a string containing entries (e.g drug names) separated by sep
#' @param sep the separator, default is "|"
#' @return a string
#'
#' @importFrom stringr str_split
#' @export
split_targeted_therapy_targets <- function(x, sep="|"){
  if (sep %in% c("\\", "^", "$", ".", "?", "*", "|", "+", "(", ")", "[")){
    regex <- paste0("\\", sep)
  } else {
    regex <- sep
  }
  classes <- unlist(strsplit(x, regex))
  prefix <- "Targeted_Therapy - "
  classes_new <- c()
  classes_old <- c()
  for (c in classes){
    if (grepl(paste0("^", prefix), c)){
      targets <- unlist(strsplit(gsub(prefix, "", c), "/"))
      if (length(targets) > 1){
        classes_new <- c(classes_new, paste0(prefix, targets))
        classes_old <- c(classes_old, c)
      }
    }
  }
  classes <- union(setdiff(classes, classes_old), classes_new)
  paste(classes, collapse=sep)
}
