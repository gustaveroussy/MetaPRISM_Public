
#' Render a table using DT
#' 
#' @param df data.frame to render
#' @param caption the table caption
#' @param full (optional) if T, render the full table
#' @param nrows (optional) number of rows to render. Ignored if full=T.
#' @param extensions (optional) passed to parameter 'extensions' of DT::datatable.
#' @param buttons (optional) passed to parameter 'options' of DT::datatable. Supported values are "copy", "csv" and
#'   "excel".
#'
#' @author Yoann Pradat
#'
#' @import DT
#' @import htmltools
#' @import htmlwidgets
#'
#' @export
render_table <- function(df, caption, full=F, nrows=5, extensions=c("Buttons", "Responsive"),
                         buttons=c("copy", "csv", "excel")){
  if (full){
    df_render <- df
  } else {
    if (nrows==-1){
      df_render <- df
    } else {
      df_render <- utils::head(df, nrows)
    }
  }

  DT::datatable(data=df_render, 
                caption=htmltools::tags$caption(style='caption-side: top; text-align: center; color:black; 
                                                font-size:150%;',
                                                caption),
                extensions=extensions,
                options=list(dom="Bfrtip", buttons=buttons,
                             initComplete = htmlwidgets::JS(
                                               "function(settings, json) {",
                                               "$(this.api().table().container()).css({'font-size': '10pt'});",
                                               "}")
                             )
  )
}
