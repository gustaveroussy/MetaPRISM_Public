suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tableExtra))
suppressPackageStartupMessages(library(tibble))

get_plot_theme <- function(dscale, base_size=12, core_size=5, scale_breaks=seq(from=0, to=1, by=0.1),
                           color_breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 12)){
  # graphical parameters
  bg_colors <- c("#f8f9fa", "#e9ecef")
  color_palette <- c("#ffc651", "#ffa759", "#ff8962", "#ff6b6b", "#cc6999", "#9968c8", "#6767f8", "#4459ce",
                     "#224ba5", "#013d7c")
  
  # background colors
  bg_fill <- rep(NA, length.out=ncol(dscale))

  if (class(bg_colors)=="list"){
    for (name in names(bg_colors)){
      bg_fill[grep(name, colnames(dscale), ignore.case=T)] <- bg_colors[[name]]
    }
  } else {
    bg_fill <- bg_colors
  }
  

  theme <- ttheme_awesome(base_size=base_size,
                          rep_mode="col",
                          core_size=core_size, 
                          scale_breaks=scale_breaks,
                          color_palette=color_palette, 
                          color_breaks=color_breaks, 
                          legend_position="top_left",
                          core=list(bg_params=list(fill=bg_fill)))

  theme
}


get_plot_data <- function(df_count, rows_more=NULL, col_tt="Tumor_Type", col_sid="Subject_Id", dcolor_norm_factor=1,
                          dcolor_quantile=0.5, base_size=12, core_size=5){
  # dscale for plot
  if (!"n" %in% colnames(df_count))
    dscale <- df_count %>%
      group_by(.data[[col_tt]]) %>%
      mutate(n=n())
  else {
    dscale <- df_count
  }

  dscale <- dscale %>%
    group_by(.data[[col_tt]], n) %>%
    summarize_at(vars(-{{col_sid}}), ~sum(.x>0)) %>%
    mutate_at(vars(-{{col_tt}},-n), ~./n)

  # dcolor for plot
  dcolor <- data.frame(check.names=F)
  row_vals <- dscale[[col_tt]]
  col_vals <- colnames(df_count)[!colnames(df_count) %in% c(col_tt, col_sid, "n")]
  for (row_val in row_vals){
    dcolor_row <- setNames(list(row_val), col_tt)
    for (col_val in col_vals){
      df_sel <- df_count %>% filter(.data[[col_tt]]==row_val, .data[[col_val]]==1)
      dcolor_row <- c(dcolor_row, setNames(quantile(rowSums(df_sel[col_vals]), dcolor_quantile), col_val))
    }
    dcolor <- rbind(dcolor, data.frame(dcolor_row, check.names=F))
  }
  dcolor <- as_tibble(dcolor) %>% replace(is.na(.),0)

  # cols more
  cols_more <- list("n="=dscale$n)
  
  # dscale formatting
  dscale$n <- NULL
  dscale <- column_to_rownames(.data=dscale, var=col_tt)
  dscale <- t(as.matrix(dscale))

  # dcolor formatting
  dcolor$n <- NULL
  dcolor <- column_to_rownames(.data=dcolor, var=col_tt)
  dcolor <- t(as.matrix(dcolor))

  # theme
  theme <- get_plot_theme(dscale, base_size=base_size, core_size=core_size)

  list(dscale=dscale, dcolor=dcolor, cols_more=cols_more, rows_more=rows_more, theme=theme)
}

