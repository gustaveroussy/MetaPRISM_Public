#' Build a \code{ComplexHeatmap::Heatmap-class} object from a dataframe
#' 
#' Apply the \code{ComplexHeatmap::make_comb_mat} on a named list prepared from
#' a dataframe.
#'
#' @param df dataframe with one field for identifiers and one field to group identifiers
#' @param field_set Name of the field identifying groups.
#' @param field_identifier Name of the field identifier entries.
#' @return m a \code{ComplexHeatmap::Heatmap-class} object
#'
#' @importFrom ComplexHeatmap make_comb_mat
#' @author Yoann Pradat
#' @export
make_upset_m <- function(df, field_set, field_identifier){
  set_2_ids <- list()
  set_names <- unique(df[[field_set]])

  for (set_name in set_names){
    set_ids <- df[df[[field_set]]==set_name, field_identifier]
    if (length(set_ids) > 0){
      set_2_ids[[set_name]] <- set_ids
    }
  }

  make_comb_mat(set_2_ids)
}


#' Draw Upset plot with predefined them
#' 
#' Draw an upset plot using the function \code{ComplexHeatmap::Upset}. It takes as input 
#' an object produced by \code{ComplexHeatmap::make_comb_mat} and draws the plot 
#' with a predefined theme.
#'
#' @param m object produced by ComplexHeatmap::make_comb_mat
#' @param row_annot_fontsize (optional) fontsize of rows annotations
#' @param pt_size (optional) size for dots representing combination sets
#' @param lwd (optional) line width for the combination sets
#' @param height_top_annot (optional) height in inches
#' @param width_set_size (optional) height in inches
#' @param margin_row_text (optional) width in inches of the margin surrounding set names
#' @param ... Extra parameters passed to ComplexHeatmap::UpSet.
#'
#' @importFrom ComplexHeatmap set_size comb_size UpSet HeatmapAnnotation rowAnnotation anno_barplot anno_text draw
#' @importFrom ComplexHeatmap decorate_annotation column_order comb_degree set_name max_text_width
#' @author Yoann Pradat
#' @export
draw_upset_plot <- function(m, row_annot_fontsize=8, pt_size=grid::unit(3,"mm"), lwd=2, height_top_annot=3,
                            width_set_size=3, margin_row_text=4, ...){
  ss <- set_size(m)
  cs <- comb_size(m)
  ht <- UpSet(m, set_order=order(ss), comb_order=order(comb_degree(m), -cs), pt_size=pt_size, lwd=lwd,
              top_annotation=HeatmapAnnotation("Intersection"=anno_barplot(cs,
                                                                           ylim=c(0, max(cs)*1.1),
                                                                           border=FALSE,
                                                                           gp=grid::gpar(fill="#F5AA32"),
                                                                           height=grid::unit(height_top_annot,
                                                                                             "inches")),
                                               annotation_name_side="left",
                                               annotation_name_rot=90),
              left_annotation=rowAnnotation("Set size"=anno_barplot(set_size(m),
                                                                    axis_param=list(direction="reverse"),
                                                                    border=FALSE,
                                                                    gp=grid::gpar(fill="#F5AA32"),
                                                                    width=grid::unit(width_set_size, "inches")),
                                            set_name=anno_text(set_name(m),
                                                               location=0.5,
                                                               just="center",
                                                               gp=grid::gpar(fontsize=row_annot_fontsize),
                                                               width=max_text_width(set_name(m)) + 
                                                                 grid::unit(margin_row_text,"mm"))),
              show_row_names=F,
              right_annotation=NULL,
              ...)
  ht <- draw(ht)
  decorate_annotation("Intersection", {
                        grid::grid.text(cs[column_order(ht)],
                                        x=seq_along(column_order(ht)),
                                        y=grid::unit(cs[column_order(ht)], "native") + grid::unit(2, "pt"),
                                        default.units="native", just=c("left", "bottom"),
                                        gp=grid::gpar(fontsize=8, col="#404040"), rot=45)
              })
}
