# @created: 05 Nov 21
# @modified: 16 May 22
# @authors: Yoann Pradat
#
# Considering a sub table with selected features, prepare it for being fed to a survival model. The preparation steps
# include
#  - setting types of covariates according to their category
#  - analyzing patterns of missing values
#  - dropping covariates with too few non-zero values
#  - grouping small-size categories of qualitative features into "Other". This is required to avert
#    numerical issues in cross-validation where some splits would have columns of 0s.
#  - checking for redundant features before imputation
#  - imputing missing values
#  - encoding qualitative features by spreading them into into dummy columns
#  - adding interaction terms if any
#  - checking for redundant features after spreading and interactions

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Hmisc))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(mice))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(yaml))

source("../common/functions/model_utils.R")

# functions ============================================================================================================

remove_too_many_na <- function(df_dat, df_cov, df_cov_rmv, thresh=0.5){
  means_na <- colMeans(is.na(df_dat))
  covs_rmv <- names(means_na[means_na >= thresh])

  df_cov_rmv_na <- df_cov %>% filter(Covariate %in% covs_rmv) %>% 
    mutate(Reason=paste("More than", paste0(thresh*100, "%"), "NA"))
  df_cov_rmv <- bind_rows(df_cov_rmv, df_cov_rmv_na)
  df_cov <- df_cov %>% filter(!Covariate %in% covs_rmv)
  df_dat <- df_dat %>% select(-all_of(covs_rmv))
  if (length(covs_rmv) > 0) {
    cat(paste0("-INFO: removed the following columns with too many NAs:\n\t", paste0(covs_rmv, collapse="\n\t"), "\n"))
  }

  list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv)
}


remove_exactly_redundant <- function(df_dat, df_cov, df_cov_rmv){
  covs_names <- df_cov$Covariate
  covs_redun <- df_cov$Redundancy
  n_covs <- df_cov %>% nrow()
  
  covs_rmv <- c()
  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_redun <- covs_redun[[j]]
    if (!is.null(cov_redun)){
      if (grepl("&", cov_redun)){
        cov_redun_names <- unlist(str_split(cov_redun, "\\&"))
        if (all(cov_redun_names %in% covs_names)) covs_rmv <- c(covs_rmv, cov_name)
      } else if (grepl("|", cov_redun)){
        cov_redun_names <- unlist(str_split(cov_redun, "\\|"))
        if (any(cov_redun_names %in% covs_names)) covs_rmv <- c(covs_rmv, cov_name)
      } else {
        if (cov_redun %in% covs_names) covs_rmv <- c(covs_rmv, cov_name)
      }
    }
  }


  df_cov_rmv_exc <- df_cov %>% filter(Covariate %in% covs_rmv) %>% mutate(Reason="Exactly redundant")
  df_cov_rmv <- bind_rows(df_cov_rmv, df_cov_rmv_exc)
  df_cov <- df_cov %>% filter(!Covariate %in% covs_rmv)
  df_dat <- df_dat %>% select(-all_of(covs_rmv))
  if (length(covs_rmv) > 0) {
    cat(paste0("-INFO: removed the following exactly redundant columns:\n\t", paste0(covs_rmv, collapse="\n\t"), "\n"))
  }

  list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv)
}


get_covs_reprs <- function(df_cov){
  cols_levels <- sort(colnames(df_cov)[grepl("Class_Lvl", colnames(df_cov))])
  df_cov <- df_cov %>% unite("Class", all_of(cols_levels), sep=" - ", remove=F)
  df_cov <- df_cov %>% mutate(Class=gsub(" - NA", "", Class)) %>%
    mutate(Class=ifelse(Class %in% c("", "NA"), NA, Class))
  df_cov$Prefix <- sapply(df_cov$Covariate, function(x) unlist(str_split(x, "_"))[1], USE.NAMES=F)
  df_cov <- df_cov %>% mutate(Prefix=ifelse(grepl("^DNA|^RNA", Class), Class, Prefix))
  covs_reprs <- df_cov %>% distinct(Prefix, Class, .keep_all=T) %>% pull(Covariate)

  covs_reprs
}


draw_plot_missing_values <- function(df_dat, df_cov, output, point_size=10, text_size=4){
  covs_names <- get_covs_reprs(df_cov)
  count_na <- colSums(is.na(df_dat[,covs_names]))
  count_na <- count_na[count_na > 0]
  df_count_na <- data.frame(x=names(count_na), y=count_na)
  names(df_count_na) <- c("Covariate", 'Count_NA')
  df_count_na$Covariate <- factor(df_count_na$Covariate, levels=df_count_na$Covariate[order(df_count_na$Count_NA)])

  if (nrow(df_count_na)>0){
    theme_set(theme_bw())
    plot <- ggplot(df_count_na, aes(x=Covariate, y=Count_NA, label=Count_NA)) + 
      geom_point(stat='identity', fill="black", size=point_size)  +
      geom_segment(aes(y=0, x=Covariate, yend=Count_NA, xend=Covariate),  color = "black") +
      geom_text(color="white", size=text_size) +
      labs(title="Number of missing values per covariate") +
      coord_flip()
    
    ggsave(filename=output, plot=plot)
    cat(paste("-INFO: plot of missing values saved at", output, "\n"))
  } else {
    pdf(output)
    dev.off()
    cat(paste("-INFO: no missing values so empty plot saved at", output, "\n"))
  }
}


drop_uniquely_valued_covs <- function(df_dat, df_cov, df_cov_rmv){
  covs_names <- df_cov$Covariate
  n_unique_per_cov <- df_dat %>% summarize_all(n_distinct, na.rm=T)

  covs_rmv <- c()
  for (cov_name in covs_names){
    if (n_unique_per_cov[[cov_name]]<=1){
      covs_rmv <- c(covs_rmv, cov_name)
      df_dat[[cov_name]] <- NULL
    }
  }

  df_cov_rmv_unq <- df_cov %>% filter(Covariate %in% covs_rmv) %>% mutate(Reason="Uniquely valued")
  df_cov_rmv <- bind_rows(df_cov_rmv, df_cov_rmv_unq)
  df_cov <- df_cov %>% filter(!Covariate %in% covs_rmv)
  if (length(covs_rmv)>0){
    cat(paste("-INFO: removed variable(s) with only 1 unique non-NA value:\n\t"))
    cat(paste(covs_rmv, collapse="\n\t"), "\n")
  }

  list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv)
}


drop_min_size_cats <- function(df_dat, df_cov, df_cov_rmv, min_size=args$min_size){
  if (min_size==1) return(list(df_dat=df_dat, df_cov=df_cov))

  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()

  covs_rmv <- c()
  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    if (cov_type=="Continuous" | cov_type=="Binary"){
      n_nonzero <- df_dat %>% filter(.data[[cov_name]]!=0) %>% nrow()
      if (n_nonzero < min_size){
        covs_rmv <- c(covs_rmv, cov_name)
        df_dat[[cov_name]] <- NULL
      }
    }
  }

  df_cov_rmv_min <- df_cov %>% filter(Covariate %in% covs_rmv) %>% 
    mutate(Reason=paste("<", min_size, "non-zero values"))
  df_cov_rmv <- bind_rows(df_cov_rmv, df_cov_rmv_min)
  df_cov <- df_cov %>% filter(!Covariate %in% covs_rmv)
  if (length(covs_rmv)>0){
    cat(paste("-INFO: removed binary or continuous variables with strictly less than", min_size, "non-zero values:\n\t"))
    cat(paste(covs_rmv, collapse="\n\t"), "\n")
  }

  list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv)
}


simplify_others_classes <- function(x, sep="|"){
  if (sep %in% c("\\", "^", "$", ".", "?", "*", "|", "+", "(", ")", "[")){
    regex <- paste0("\\", sep)
  } else {
    regex <- sep
  }
  classes <- unlist(strsplit(x, regex))
  prefix <- "Others"
  classes_new <- c()
  classes_old <- c()
  for (c in classes){
    if (grepl(paste0("^", prefix), c)){
      classes_new <- "Other"
      classes_old <- c(classes_old, c)
    }
  }
  classes <- union(setdiff(classes, classes_old), classes_new)
  paste(classes, collapse=sep)
}


set_type_covs <- function(df_dat, df_cov, discrete_ordered_as_factor=T){
  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()
  cat("-setting numerical types of all covariates ...")

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    if (cov_type=="Continuous"){
      df_dat[[cov_name]] <- as.numeric(df_dat[[cov_name]])
    } else if (cov_type=="Binary"){
      df_dat[[cov_name]] <- as.numeric(df_dat[[cov_name]])
    } else if (cov_type=="Discrete_Ordered"){
      if (discrete_ordered_as_factor){
        levels <- sort(unique(df_dat[[cov_name]]))
        df_dat[[cov_name]] <- factor(df_dat[[cov_name]], ordered=T, levels=levels)
      } else {
        df_dat[[cov_name]] <- as.numeric(df_dat[[cov_name]])
      }
    } else if (cov_type=="Discrete_Exclusive"|cov_type=="Discrete_Not_Exclusive"|cov_type=="Binary"){
      df_dat[[cov_name]] <- as.factor(df_dat[[cov_name]])
    } else if (cov_type=="Date"){
      df_dat[[cov_name]] <- as.Date(df_dat[[cov_name]], format="%d/%m/%Y")
    } else if (cov_type=="Text"){
      stop("-cannot transform text variables, please preprocess it and change its nature")
    }
  }
  cat("done!\n")

  df_dat
}


remove_predicted_redundant <- function(df_dat, df_cov, df_cov_rmv, indices=NULL, tlinear=F, nk=0, r2=0.9, pr=T,
                                       step="before"){
  if (nrow(df_cov)==1){
    return(list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv))
  } else {
    if (is.null(indices)){
      indices <- seq(1, nrow(df_dat))
    }
    # check redundancy
    # uses flexible parametric additive models (see 'areg' and its use of regression splines) to determine how well each
    # variable can be predicted from the remaining variables.  Variables are dropped in a stepwise fashion, removing the
    # most predictable variable at each   step. The remaining variables are used to predict.  The process continues
    # until no variable still in the list of predictors can be predicted with an R^2 or adjusted R^2 of at least 'r2' or
    # until dropping the variable with the highest R^2 (adjusted or ordinary) would cause a variable that was dropped
    cat(paste("-fitting independent supervised models to identify redundant features", step, "imputation ...\n"))
    redundant <- redun(~., data=df_dat[indices, df_cov$Code], tlinear=tlinear, nk=nk, r2=r2, pr=pr)

    codes_redun <- redundant$Out
    covs_redun <- df_cov %>% filter(Code %in% codes_redun) %>% pull(Covariate)

    if (length(covs_redun)>0)
      cat(paste0("\t->INFO: the following columns are redundant:\n\t",  paste(covs_redun, collapse="\n\t")), "\n")

    # align column types in case they are discrepant for unknown reasons
    for (name in colnames(df_cov_rmv)){
      df_cov_rmv[[name]] <- as.character(df_cov_rmv[[name]])
    }

    for (name in colnames(df_cov)){
      df_cov[[name]] <- as.character(df_cov[[name]])
    }

    df_cov_red <- df_cov %>% filter(Code %in% codes_redun) %>% 
      mutate(Reason=paste("Redundant R^2 >=", r2))
    df_cov_rmv <- bind_rows(df_cov_rmv, df_cov_red)
    df_dat <- df_dat %>% select(-all_of(codes_redun))
    df_cov <- df_cov %>% filter(!Code %in% codes_redun)

    return(list(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv))
  }
}


group_min_size_cats <- function(df_dat, df_cov, min_size=1, spread_discrete_ordered=F, sep_discrete_not_exclusive="\\|"){
  if (min_size==1) return(df_dat)
  cat(paste("-grouping categories of categorical variables with strictly less than", min_size, "occurrences ...\n"))

  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()
  covs_grouped <- c()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    if (cov_type=="Discrete_Exclusive"){
      if (length(levels(df_dat[[cov_name]]))>2){
        small_lvls <- df_dat %>% group_by_at(all_of(cov_name)) %>% summarise(n=n()) %>%
          filter(n<min_size) %>% pull(var=cov_name)
        if (length(small_lvls)>1){
          levels(df_dat[[cov_name]])[levels(df_dat[[cov_name]]) %in% small_lvls] <- "Other"
          cat(paste0("\t->INFO: grouped small-size categories of ", cov_name, ":\n\t",
                     paste(small_lvls, collapse="\n\t")), "\n")
          covs_grouped <- c(covs_grouped, cov_name)
        }
      }
    } else if (cov_type=="Discrete_Not_Exclusive"){
      df_dat[[cov_name]] <- as.character(df_dat[[cov_name]])
      small_lvls <- df_dat %>% separate_rows(all_of(cov_name), sep=sep_discrete_not_exclusive) %>%
        group_by_at(all_of(cov_name)) %>% summarise(n=n()) %>%
        filter(n<min_size) %>% pull(var=cov_name)
      if ("None" %in% small_lvls) small_lvls <- setdiff(small_lvls, "None")
      if (length(small_lvls)>1){
        regex <- paste0("(",paste0(small_lvls, collapse=")|("), ")")
        df_dat[[cov_name]] <- gsub(regex, "Other", df_dat[[cov_name]])
        cat(paste0("\t->INFO: grouped small-size categories of ", cov_name, ":\n\t",
                   paste(small_lvls, collapse="\n\t")), "\n")
        covs_grouped <- c(covs_grouped, cov_name)
      }
    } else if (cov_type=="Discrete_Ordered"){
      if (spread_discrete_ordered){
        lvls_bef <- levels(df_dat[[cov_name]])
        df_dat <- aggregate_small_strata(df_dat, cov_name, max_strata=Inf, min_size=min_size)
        lvls_aft <- levels(df_dat[[cov_name]])
        small_lvls <- setdiff(lvls_bef, lvls_aft)
        if (length(small_lvls)>0){
          cat(paste0("\t->INFO: grouped small-size categories of ", cov_name, ":\n\t",
                     paste(small_lvls, collapse="\n\t")), "\n")
          covs_grouped <- c(covs_grouped, cov_name)
        }
      }
    }
  }

  if (length(covs_grouped) == 0)
    cat(paste0("\t->INFO: there is no categorical covariate with strata of size smaller than ", min_size, "!\n"))
  
  df_dat
}


encode_discrete_exclusive <- function(df_dat, df_cov, cov_name, name_sep="_", remove_first=T, remove_col=T,
                                      ref_level=NULL){
  vec <- df_dat[[cov_name]]

  if (!is.factor(vec)){
    vec <- as.factor(vec)
  }

  if (is.null(ref_level)){
    if (remove_first){
      ref_level <- levels(vec)[1]
    }
  }

  if (!is.null(ref_level)){
    stopifnot(ref_level %in% levels(vec))
    dum_levels <- setdiff(levels(vec), ref_level)
  }

  df_dat_dum <- NULL
  for (dum_level in dum_levels){
    cov_dum <- paste0(cov_name, name_sep, dum_level)
    df_dat_dum <- tibble({{cov_dum}}:=as.numeric(vec==dum_level))
    df_cov_dum <- df_cov[df_cov$Covariate==cov_name,,drop=F] %>%
      mutate(Covariate=cov_dum, Plot_Name=paste(Plot_Name, dum_level),  Nature="Binary")

    df_dat <- bind_cols(df_dat, df_dat_dum)
    df_cov <- bind_rows(df_cov, df_cov_dum)
  }

  if (remove_col){
    df_dat[[cov_name]] <- NULL
    df_cov <- df_cov %>% filter(Covariate!=cov_name)
  }

  list(df_dat=df_dat, df_cov=df_cov)
}


encode_discrete_not_exclusive <- function(df_dat, df_cov, cov_name, sep="\\|", name_sep="_", remove_col=T){
  vec <- df_dat[[cov_name]]
  dum_levels <- df_dat %>% separate_rows(all_of(cov_name), sep=sep) %>% filter(!is.na(.data[[cov_name]])) %>%
    pull(.data[[cov_name]]) %>% unique()
  dum_levels <- setdiff(dum_levels, "None")

  for (dum_level in dum_levels){
    cov_dum <- paste0(cov_name, name_sep, dum_level)
    df_dat <- df_dat %>% mutate(!!cov_dum:=ifelse(grepl(dum_level, .data[[cov_name]]), 1, ifelse(is.na(.data[[cov_name]]), NA, 0)))
    df_dat <- df_dat %>% mutate(!!cov_dum:=ifelse(.data[[cov_name]]=="None", 0, .data[[cov_dum]]))

    df_cov_dum <- df_cov[df_cov$Covariate==cov_name,,drop=F] %>%
      mutate(Covariate=cov_dum, Plot_Name=paste(Plot_Name, dum_level), Nature="Binary")
    df_cov <- bind_rows(df_cov, df_cov_dum)
  }

  if (remove_col){
    df_dat[[cov_name]] <- NULL
    df_cov <- df_cov %>% filter(Covariate!=cov_name)
  }


  list(df_dat=df_dat, df_cov=df_cov)
}


encode_discrete_ordered <- function(df_dat, cov_name){
  if (is.ordered(df_dat[[cov_name]])){
    return(df_dat)
  } else {
    levels_ori <- df_dat %>% filter(!is.na(.data[[cov_name]])) %>% pull(.data[[cov_name]]) %>% unique()
    levels_ori <- sort(levels_ori)

    # try to convert to numeric if possible
    # if not possible, replace levels by levels positions
    n_na <- sum(is.na(df_dat[[cov_name]]))
    vec_num <- suppressWarnings(as.numeric(df_dat[[cov_name]]))
    levels_num <- sort(unique(vec_num))

    if (sum(is.na(df_dat[[cov_name]])) < sum(is.na(vec_num))){
      vec_num <- as.integer(factor(df_dat[[cov_name]], ordered=T, levels=levels_ori))
      levels_num <- seq(1, length(levels_ori))
    }

    df_dat[[cov_name]] <- factor(vec_num, ordered=T, levels=levels_num)

    return(df_dat)
  }
}


encode_discrete_covs <- function(df_dat, df_cov, spread_discrete_ordered=F, sep_discrete_not_exclusive="\\|",
                                 discrete_not_exclusive_only=F){
  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]

    if (cov_type=="Discrete_Not_Exclusive"){
      out <- encode_discrete_not_exclusive(df_dat, df_cov, cov_name, sep=sep_discrete_not_exclusive)
      df_dat <- out$df_dat
      df_cov <- out$df_cov
    } else if (cov_type=="Discrete_Exclusive" & !discrete_not_exclusive_only){
      out <- encode_discrete_exclusive(df_dat, df_cov, cov_name, remove_first=T)
      df_dat <- out$df_dat
      df_cov <- out$df_cov
    } else if (cov_type=="Discrete_Ordered" & !discrete_not_exclusive_only){
      if (spread_discrete_ordered){
        out <- encode_discrete_exclusive(df_dat, df_cov, cov_name, remove_first=T)
        df_dat <- out$df_dat
        df_cov <- out$df_cov
      } else {
        df_dat <- encode_discrete_ordered(df_dat, cov_name)
      }
    }
  }

  if (discrete_not_exclusive_only){
    cat("-encoded discrete not exclusive categorical variables into numerical covariables\n")
  } else {
    cat("-encoded all categorical variables into numerical covariables\n")
  }

  list(df_dat=df_dat, df_cov=df_cov)
}


choose_default_ref_levels <- function(df_dat, df_cov, spread_discrete_ordered){
  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  covs_ref_levels <- df_cov$Reference_Level
  n_covs <- df_cov %>% nrow()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    ref_level <- covs_ref_levels[[j]]

    if (is.na(ref_level) | !ref_level %in% df_dat[[cov_name]]){
      if (cov_type=="Discrete_Exclusive" | (cov_type=="Discrete_Ordered" & spread_discrete_ordered)){
        ref_level <- levels(df_dat[[cov_name]])[1]
        df_cov <- df_cov %>% mutate(Reference_Level=ifelse(Covariate==cov_name, ref_level, Reference_Level))
      }
    }
  }

  df_cov
}


get_covs_interactions <- function(covs, features, features_all, df_cov){
  for (cov in features){
    if (cov %in% names(features_all)){
      covs <- get_covs_interactions(covs, features_all[[cov]], features_all, df_cov)
    } else {
      if (grepl(":", cov)){
        covs <- c(covs, cov)
      }
    }
  }

  covs
}


add_covs_interactions <- function(covs_int, df_dat, df_cov, min_size){
  covs_int_new <- c()
  cat(paste0("-adding interactions if any ...\n"))

  for (cov_int in covs_int){
    # get list of covariate groups btw which to build interactions.
    cov_int_groups <- unlist(str_split(cov_int, ":"))

    # retrieve actual covariate name corresponding to covariate group (covariate class or single covariate)
    cov_int_lists <- lapply(cov_int_groups, function(cov){
      covs_more <- df_cov %>% filter(grepl(paste0("^", cov), df_cov$Class)) %>% pull(Covariate)
      if (length(covs_more) == 0){
        covs_more <- df_cov %>% filter(grepl(cov, df_cov$Covariate)) %>% pull(Covariate)
      }
      covs_more})

    # get all n-way interactions with n the number of covariate groups
    cov_combinations <- expand.grid(cov_int_lists, stringsAsFactors=F)

    # add only interactions with sufficient occurrences
    for (i in 1:nrow(cov_combinations)){
      cov_combination <- as.character(as.vector(cov_combinations[i,]))
      combination_name <- paste(cov_combination, collapse=":")
      combination_prod <- df_dat[,cov_combination] %>% mutate(Prod=Reduce(`*`, .)) %>% pull(Prod)
      if (sum(combination_prod!=0) >= min_size){
        covs_int_new <- c(covs_int_new, combination_name)
        df_dat[[combination_name]] <- combination_prod
        df_cov_combination <- data.frame(Covariate=combination_name) %>% 
          mutate(Nature="Continuous", Class_Lvl_1="Interaction", Class_Lvl_2=cov_int)
        df_cov <- bind_rows(df_cov, df_cov_combination)
      }
    }
  }

  if (length(covs_int_new) > 0) {
    cat(paste0("\t->INFO: added the following interactions:\n\t", paste0(covs_int_new, collapse="\n\t"), "\n"))
  }

  list(df_dat=df_dat, df_cov=df_cov)
}


special_grouping_covs <- function(df_dat, df_cov){
  # special case of Classes_Before_Biopsy. Split multiple targets therapies into multiple single-target therapies
  if ("Classes_Before_Biopsy" %in% df_cov$Covariate){
    cov_name <- "Classes_Before_Biopsy"
    df_dat[[cov_name]] <- sapply(df_dat[[cov_name]], function(x) simplify_others_classes(x), USE.NAMES=F) 
    df_dat[[cov_name]] <- sapply(df_dat[[cov_name]], function(x) split_targeted_therapy_targets(x), USE.NAMES=F)
    df_dat[[cov_name]][df_dat[[cov_name]]=="NA"] <- NA
    cat("-performed special grouping operation on values of Classes_Before_Biopsy\n")
  }

  df_dat
}


set_type_covs <- function(df_dat, df_cov){
  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    if (cov_type=="Continuous"){
      df_dat[[cov_name]] <- as.numeric(df_dat[[cov_name]])
    } else if (cov_type=="Binary"){
      df_dat[[cov_name]] <- as.numeric(df_dat[[cov_name]])
    } else if (cov_type=="Discrete_Ordered"){
      levels <- sort(unique(df_dat[[cov_name]]))
      df_dat[[cov_name]] <- factor(df_dat[[cov_name]], ordered=T, levels=levels)
    } else if (cov_type=="Discrete_Exclusive"|cov_type=="Discrete_Not_Exclusive"|cov_type=="Binary"){
      df_dat[[cov_name]] <- as.factor(df_dat[[cov_name]])
    } else if (cov_type=="Date"){
      df_dat[[cov_name]] <- as.Date(df_dat[[cov_name]], format="%d/%m/%Y")
    } else if (cov_type=="Text"){
      stop("-cannot transform text variables, please preprocess it and change its nature")
    }
  }

  cat("-set numerical types of all covariates\n")
  df_dat
}


perform_multiple_imputations <- function(df_dat, df_cov, covs_outcome, n_imputations, n_iterations, n_cores, seed){
  if (all(c("Survival_Time", "Survival_Status")  %in% covs_outcome)){
    # We use Nelson-Aalen estimates as a feature of the imputation models
    # Read "Imputing values for Cox model". Ian R. White and Patrick Royston, 2009.
    df_dat[,'Nelson_Aalen'] <- nelsonaalen(df_dat, Survival_Time, Survival_Status)
    df_dat_bimp <- df_dat[,c(df_cov$Code, covs_outcome, "Nelson_Aalen")]
    added_nelson_aalen <- T
  } else {
    df_dat_bimp <- df_dat[,c(df_cov$Code, covs_outcome)]
    added_nelson_aalen <- F
  }


  # Missing values imputation by Multiple Imputation Chained Equations
  # This is provided by mice package
  cat(paste("-INFO: table before imputation:", nrow(df_dat_bimp),"rows and", ncol(df_dat_bimp), "columns\n"))

  # imputation model
  # maybe collinear/redundant columns should be removed before imputation?
  set.seed(seed)
  ini <- mice(df_dat_bimp, maxit=0, remove.collinear=F)
  meth <- ini$meth
  pred <- ini$pred

  if (all(c("Survival_Time", "Survival_Status")  %in% covs_outcome)){
    pred[,"Survival_Time"] <- 0 # White and Royston suggestion (2009)
  }

  # prepare for parallel execution
  cl <- makePSOCKcluster(n_cores)
  registerDoParallel(cl)
  capt <- clusterEvalQ(cl, library(doParallel))
  capt <- clusterEvalQ(cl, library(mice))
  capt <- clusterEvalQ(cl, library(Hmisc))

  cat(paste("-running", n_imputations, "imputations in parallel on", n_cores, "cores..."))
  imp_pars <- 
    foreach(x = 1:n_imputations) %dopar%{
      mice(df_dat_bimp, m=1, pred=pred, meth=meth, maxit=n_iterations, print=F, seed=seed*x,
           donors=5, nnet.MaxNWts = 5000)
    }
  cat("done!\n")
  stopCluster(cl)

  # merge output from parallel execution
  imp_merged <- imp_pars[[1]]
  if (n_imputations >= 2){
    for (n in 2:length(imp_pars)){
      imp_merged <- ibind(imp_merged, imp_pars[[n]])
    }
  }

  # Recover all imputations in long dataframe
  # This will allow us to do some transformations of 
  # the dataset before working on models

  df_dat_aimp_long <- complete(imp_merged, action='long', include=T)
  nrow_one_iter <- nrow(df_dat_aimp_long)/(n_imputations+1)
  stopifnot(nrow_one_iter==nrow(df_dat_bimp))

  if (added_nelson_aalen){
    # Remove Nelson-Aalen estimate
    df_dat_aimp_long$Nelson_Aalen <- NULL
  }

  df_dat_aimp_long
}


main <- function(args){
  # load data and covariates tables
  df_dat <- load_table(args$input_dat)
  df_cov <- load_table(args$input_cov)
  sep_discrete_not_exclusive <- "\\|"
  spread_discrete_ordered <- tolower(args$spread_discrete_ordered)=="yes"
  df_cov_rmv <- tibble()

  cat(paste("-INFO: table after loading:", nrow(df_dat),"rows and", ncol(df_dat), "columns\n"))

  # process before imputation ==========================================================================================

  df_cov <- unite_class_levels(df_cov)

  # outcome covs
  covs_outcome <- df_cov %>% filter(grepl("Outcome", Class)) %>% pull(Covariate)
  df_cov_outcome <- df_cov %>% filter(Covariate %in% covs_outcome)
  df_cov <- df_cov %>% filter(!Covariate %in% covs_outcome)

  # drop cases with missing outcome
  df_dat <- df_dat[rowSums(is.na(df_dat[, covs_outcome]))==0,]

  # outcome covs transformation if any
  if ("Survival_Status" %in% covs_outcome){
    df_dat <- df_dat %>% mutate(Survival_Status=ifelse(Survival_Status=="Deceased", 1, 0))
  }

  # remove covariates with too many missing values
  out <- remove_too_many_na(df_dat, df_cov, df_cov_rmv)
  df_dat <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # remove covariates that are exactly redundant
  out <- remove_exactly_redundant(df_dat, df_cov, df_cov_rmv)
  df_dat <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # draw plot of missing values
  draw_plot_missing_values(df_dat, df_cov, args$output_pmc)

  # special grouping of some variables
  df_dat <- special_grouping_covs(df_dat, df_cov)

  # set types of covariates according to its nature (Continuous, Discrete_Exclusive/Ordered/Not_Exclusive)
  df_dat <- set_type_covs(df_dat, df_cov)

  # drop uniquely-valued covariates
  out <- drop_uniquely_valued_covs(df_dat, df_cov, df_cov_rmv)
  df_dat <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # drop small-size quantitative covariates
  out <- drop_min_size_cats(df_dat, df_cov, df_cov_rmv, min_size=args$min_size)
  df_dat <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # group small-size categorical covariates
  df_dat <- group_min_size_cats(df_dat, df_cov, min_size=args$min_size,
                                spread_discrete_ordered=spread_discrete_ordered, 
                                sep_discrete_not_exclusive=sep_discrete_not_exclusive)

  # encode discrete covariates before imputation
  # in this step only discrete not exclusive categorical variables are expanded
  # expanding discrete exclusive categorical variables may result in nonsense imputation for these variables
  # where a case would belong to more than 2 categories or to no category at all
  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  n_covs <- df_cov %>% nrow()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]

    if (cov_type=="Discrete_Not_Exclusive"){
      out <- encode_discrete_not_exclusive(df_dat, df_cov, cov_name, sep=sep_discrete_not_exclusive)
      df_dat <- out$df_dat
      df_cov <- out$df_cov
    }
  }
  cat("-INFO: encoded discrete not exclusive categorical variables into numerical covariables\n")

  # describe missing data patterns
  covs_reprs <- get_covs_reprs(df_cov)
  if (length(covs_reprs)>1){
    pdf(args$output_pmp)
    md_patterns <- md.pattern(df_dat[,covs_reprs], plot=T, rotate.names=T)
    dev.off()
    cat(paste0("-INFO: # patterns missingness: ", nrow(md_patterns), "\n"))
  } else {
    pdf(args$output_pmp)
    dev.off()
  }

  cat(paste0("-INFO: # samples if complete case: ", sum(cci(df_dat)), "\n"))
  cat(paste("-INFO: table before redundancy analysis:", nrow(df_dat),"rows and", ncol(df_dat), "columns\n"))

  # encode covariable names to algo-friendly codes
  out <- encode_covariable_names(df_dat, df_cov)
  df_dat <- out$df_dat
  df_cov <- out$df_cov

  # remove redundant features from data table
  out <- remove_predicted_redundant(df_dat=df_dat, df_cov=df_cov, df_cov_rmv=df_cov_rmv, r2=0.9, step="before")
  df_dat <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # impute multiple datasets ===========================================================================================

  if (nrow(df_cov)==1 | sum(colSums(is.na(df_dat)))==0){
    # let's not perform imputation in the case where only 1 coviarate is considered or when there is no missing data
    df_dat_aimp_long <- cbind(data.frame(.imp=rep(0, nrow(df_dat)), .id=seq(1, nrow(df_dat))), as.data.frame(df_dat))
  } else {
    df_dat_aimp_long <- perform_multiple_imputations(df_dat=df_dat, df_cov=df_cov,
                                                     covs_outcome=covs_outcome,
                                                     n_imputations=args$n_imputations,
                                                     n_iterations=args$n_iterations,
                                                     n_cores=args$n_cores,
                                                     seed=args$seed)
    
    df_dat_aimp_long["Subject_Id"] <- rep(df_dat[["Subject_Id"]], times=args$n_imputations+1)
  }

  # process after imputation ===========================================================================================

  # Expand categorical variables into dummies
  # We choose here what category is taken as reference (i.e what category we remove)
  # It is important to keep this in mind when interpreting hazard ratios from Cox model
  # 
  # As an example, if a factor has 3 levels a,b and c.
  # If we drop c, Cox solves betas such as log(lbda) = beta_a*x_a + beta_b*x_b
  # beta_a is the log hazard ratio (lhr) of a to c
  # beta_b is the lhr of b to c
  # 
  # If we drop a, Cox solves betas such as log(lbda) = beta_b*x_b + beta_c*x_c
  # beta_b is the lhr of b to a
  # beta_c is the lhr of c to a
  
  # Convert back from covariate codes to names
  out <- decode_covariable_names(df_dat_aimp_long, df_cov)
  df_dat_aimp_long <- out$df_dat
  
  # Choose reference levels
  df_cov <- choose_default_ref_levels(df_dat_aimp_long, df_cov, spread_discrete_ordered)

  covs_names <- df_cov$Covariate
  covs_types <- df_cov$Nature
  covs_ref_levels <- df_cov$Reference_Level
  n_covs <- df_cov %>% nrow()

  for (j in 1:n_covs){
    cov_name <- covs_names[[j]]
    cov_type <- covs_types[[j]]
    ref_level <- covs_ref_levels[[j]]

    if (cov_type=="Discrete_Exclusive"){
      out <- encode_discrete_exclusive(df_dat_aimp_long, df_cov, cov_name, remove_first=T, ref_level=ref_level)
      df_dat_aimp_long <- out$df_dat
      df_cov <- out$df_cov
    } else if (cov_type=="Discrete_Ordered"){
      if (spread_discrete_ordered){
        out <- encode_discrete_exclusive(df_dat_aimp_long, df_cov, cov_name, remove_first=T, ref_level=ref_level)
        df_dat_aimp_long <- out$df_dat
        df_cov <- out$df_cov
      } else {
        df_dat_aimp_long[[cov_name]] <- as.numeric(df_dat_aimp_long[[cov_name]])
        df_cov <- df_cov %>% mutate(Nature=ifelse(Covariate==cov_name, "Continuous", Nature))
      }
    }
  }
  cat("-INFO: encoded all categorical variables into numerical covariables\n")
  
  # set types
  df_dat_aimp_long <- set_type_covs(df_dat_aimp_long, df_cov)

  # add interaction terms if any
  features_all <- read_yaml(args$config_yaml)[[args$config_section]]$features
  df_cov <- unite_class_levels(df_cov)
  covs_int <- get_covs_interactions(covs=c(), features=args$features, features_all=features_all,  df_cov=df_cov)

  if (!is.null(covs_int)){
    out <- add_covs_interactions(covs_int=covs_int, df_dat=df_dat_aimp_long, df_cov=df_cov,
                                 min_size=args$min_size*(args$n_imputations+1))
    df_dat_aimp_long <- out$df_dat
    df_cov <- out$df_cov
  }
  
  # encode covariable names to algo-friendly codes once again
  df_cov$Code <- NULL
  out <- encode_covariable_names(df_dat_aimp_long, df_cov)
  df_dat_aimp_long <- out$df_dat
  df_cov <- out$df_cov

  # remove redundant features from data table
  out <- remove_predicted_redundant(df_dat=df_dat_aimp_long, df_cov=df_cov, df_cov_rmv=df_cov_rmv, r2=0.9,
                                    indices=seq(1, nrow(df_dat)), step="after")
  df_dat_aimp_long <- out$df_dat
  df_cov <- out$df_cov
  df_cov_rmv <- out$df_cov_rmv

  # min-max scaling
  n_cov <- nrow(df_cov)
  cov_codes <- df_cov$Code
  df_cov$Min_Value <- as.numeric(df_cov$Min_Value)
  df_cov$Max_Value <- as.numeric(df_cov$Max_Value)
  cov_mins <- df_cov$Min_Value
  cov_maxs <- df_cov$Max_Value
  df_cov <- df_cov %>% mutate(Min=NA, Max=NA)
  for (j in 1:n_cov){
    cov_code <- cov_codes[j]
    cov_max <- cov_maxs[j]
    cov_min <- cov_mins[j]
    if (is.na(cov_max))
      cov_max <- max(df_dat_aimp_long[[cov_code]], na.rm=T)
    if (is.na(cov_min))
      cov_min <- min(df_dat_aimp_long[[cov_code]], na.rm=T)
    df_dat_aimp_long[[cov_code]] <- (df_dat_aimp_long[[cov_code]] - cov_min)/(cov_max - cov_min)
    df_cov <- df_cov %>% mutate(Max=ifelse(Code==cov_code, cov_max, Max)) %>%
      mutate(Min=ifelse(Code==cov_code, cov_min, Min))
  }
  
  # save
  dir.create(args$output_dir, showWarnings=F)

  # imputed tables if imputation was performed
  if (max(df_dat_aimp_long[[".imp"]])==0){
    # no imputation was performed and using the "as.mids" function would raise the error
    # "Number of imputations (m) lower than 1"
    df_dat_cc <- df_dat_aimp_long[cci(df_dat_aimp_long),]
    df_dat_cc$.imp <- NULL
    df_dat_cc$.id <- NULL

    # complete case
    file <- file.path(args$output_dir, "data.complete_cases.tsv")
    write.table(df_dat_cc, file=file, row.names=F, sep="\t", quote=F)
    system(paste("gzip", file))
    cat(paste("-complete case data saved at", file, "\n"))
  } else {
    # create mice::mids object
    mids <- as.mids(df_dat_aimp_long)

    # save imputed
    for (l in 1:args$n_imputations){
      stopifnot(sum(is.na(complete(mids,l)))==0)
      file <- file.path(args$output_dir, paste0("data.imputed_",l,".tsv"))
      write.table(complete(mids,l), file=file, row.names=F, sep="\t", quote=F)
      system(paste("gzip", file))
      cat(paste("-imputed data", l, " saved at", file, "\n"))
    }
  }

  # covariates transformed
  file <- file.path(args$output_dir, "covs.final.tsv")
  write.table(df_cov, file=file, row.names=F, sep="\t", quote=F)
  system(paste("gzip", file))
  cat(paste("-updated covariates table saved at", file, "\n"))

  # covariates removed 
  file <- file.path(args$output_dir, "covs.removed.tsv")
  write.table(df_cov_rmv, file=file, row.names=F, sep="\t", quote=F)
  system(paste("gzip", file))
  cat(paste("-redundant covariates table saved at", file, "\n"))
}


# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  default_folder <- "../../../results/survival_analysis/data/sub_features/dna_cohort_dna_and_rna/cln_grim_1/original"

  parser <- ArgumentParser(description='Prepare tables with sub features selected.')
  parser$add_argument("--input_dat", type="character", help="Path to input data table.",
                    default=file.path(default_folder, "data.original.tsv.gz"))
  parser$add_argument("--input_cov", type="character", help="Path to input covariates table.",
                    default=file.path(default_folder, "covs.original.tsv.gz"))
  parser$add_argument("--config_yaml", type="character", default="config/config.yaml",
                      help="Path to the config.yaml file where subselections of features are defined.")
  parser$add_argument("--config_section", type="character", default="models",
                      help="Name of the section containing parameters for the models.")
  parser$add_argument("--features", type="character", help="Name of the selection of features.",
                      default="cln_grim_1")
  parser$add_argument("--spread_discrete_ordered", type="character", default="No",
                      help="Set to 'Yes' to spread ordered discrete covariates into dummy columns")
  parser$add_argument("--min_size", type="integer", default=5,
                      help="Categories of size smaller than this number will be grouped into 'Other'")
  parser$add_argument("--seed", type="integer", default=123, help="Seed for the multiple imputation process.")
  parser$add_argument("--n_imputations", type="integer", default=10, help="Number of imputed datasets.")
  parser$add_argument("--n_iterations", type="integer", default=15,
                      help="Number of iterations to generate each imputed dataset.")
  parser$add_argument("--n_cores", type="integer", default=4, help="Number of cores to be used.")
  parser$add_argument("--output_pmc", type="character", help="Path to output plot of missing data.",
                    default=file.path(default_folder, "plot_missing_data_counts.pdf"))
  parser$add_argument("--output_pmp", type="character", help="Path to output plot of missing data.",
                    default=file.path(default_folder, "plot_missing_data_patterns.pdf"))
  parser$add_argument("--output_dir", type="character", help="Path to output folder with imputed values.",
                    default=file.path(default_folder, "imputed"))
  parser$add_argument("--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
