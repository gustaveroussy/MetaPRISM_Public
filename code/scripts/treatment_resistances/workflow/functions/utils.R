suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(readxl))

aggregate_small_strata <- function(df_dat, col_strata, max_strata=9, min_size=10){
  mask_small <- TRUE
  has_changed <- TRUE
  if (!is.factor(df_dat[[col_strata]])) df_dat[[col_strata]] <- as.factor(df_dat[[col_strata]])

  while (sum(mask_small)!=0 & has_changed){
    df_cnt <- df_dat %>% group_by_at(all_of(col_strata)) %>% filter(!is.na(.data[[col_strata]])) %>%
      summarise(Count=n()) %>% arrange(.data[[col_strata]])
    stratas <- df_cnt[[col_strata]]
    stratas_init <- stratas
    mask_small <- df_cnt$Count < min_size
    
    k <- 1
    for (i in 1:nrow(df_cnt)){
      if (mask_small[[i]] & i >= k){
        mask_from_i <- mask_small[i:length(mask_small)]
        index_false <- which(!mask_from_i)
        if (length(index_false)==0){
          if (length(mask_from_i)==1 & i >= 2){
            stratas_agg <- stratas[(i-1):i]
            j <- i
          } else {
            j <- nrow(df_cnt)
            stratas_agg <- stratas[i:j]
          }
          break_loop <- T
        } else {
          j <- i + which(!mask_from_i)[[1]]-1
          stratas_agg <- stratas[i:j]
          break_loop <- F
        }
        name_agg_min <- min(as.numeric(unlist(strsplit(as.character(stratas_agg), "-"))))
        name_agg_max <- max(as.numeric(unlist(strsplit(as.character(stratas_agg), "-"))))
        name_agg <- paste0(unique(c(name_agg_min, name_agg_max)), collapse="-")
        levels(df_dat[[col_strata]])[levels(df_dat[[col_strata]]) %in% stratas_agg] <- name_agg
        levels(stratas)[levels(stratas) %in% stratas_agg] <- name_agg
        k <- j+1
        if (break_loop) break
      }
    }

    has_changed <- !setequal(stratas, stratas_init)
  }

  stratas <- df_cnt[[col_strata]] %>% unique()
  if (length(stratas)>max_strata){
    stratas_agg <- stratas[max_strata:length(stratas)]
    name_agg_min <- min(as.numeric(unlist(strsplit(as.character(stratas_agg), "-"))))
    name_agg_max <- max(as.numeric(unlist(strsplit(as.character(stratas_agg), "-"))))
    name_agg <- paste0(unique(c(name_agg_min, name_agg_max)), collapse="-")
    levels(df_dat[[col_strata]])[levels(df_dat[[col_strata]]) %in% stratas_agg] <- name_agg
  } 

  df_dat
}

load_table <- function(filepath, ...){
  if (grepl(".xlsx$", filepath)) {
    df <- read_xlsx(filepath, ...)
  } else if (grepl(".tsv$", filepath)) {
    df <- read_tsv(filepath, progress=F, show_col_types=F, ...)
  } else if (grepl(".csv$", filepath)) {
    df <- read_csv(filepath, progress=F, show_col_types=F, ...)
  } else {
    df <- read_delim(filepath, progress=F, show_col_types=F, ...)
  }

  df
}


encode_covariable_names <- function(df_dat, df_cov){
  n_covs <- nrow(df_cov)
  code_size <- ceil(log(n_covs)/log(26))
  covs_codes <- c()
  code_A <- utf8ToInt("A")
  for (j in 1:n_covs){
    cov_code <- paste0(sapply(code_size:1, function(p) intToUtf8(((j-1)/26^(p-1))%%26+code_A)), collapse="")
    covs_codes <- c(covs_codes, cov_code)
  }
  df_cov$Code <- covs_codes
  df_dat <- df_dat %>% rename_at(vars(df_cov$Covariate), ~ df_cov$Code)

  list(df_dat=df_dat, df_cov=df_cov)
}

decode_covariable_names <- function(df_dat, df_cov){
  stopifnot("Code" %in% colnames(df_cov))
  df_dat <- df_dat %>% rename_at(vars(df_cov$Code), ~ df_cov$Covariate)

  list(df_dat=df_dat, df_cov=df_cov)
}


unite_class_levels <- function(df){
  cols_levels <- sort(colnames(df)[grepl("Class_Lvl", colnames(df))])
  df <- df %>% unite("Class", all_of(cols_levels), sep="-", remove=F)
  df <- df %>% mutate(Class=gsub("-NA", "", Class)) %>%
    mutate(Class=ifelse(Class %in% c("", "NA"), NA, Class))

  df
}
