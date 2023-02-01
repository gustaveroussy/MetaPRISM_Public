# ======================================================================================================================
#
# BURDEN PLOT 
#
# ======================================================================================================================

#' Draw multiple scatter plots showing cumulative distributions
#'
#' Scatter plots are organized by groups and stacks. The groups define the numbers of subplots while stacks define
#' the number of points clouds in each subplot.
#'
#' @param dfs A 2-level named list of dataframes. 1st level of names are stacks, 2nd level of names are groups.
#' @param groups A character vector.
#' @param stacks A character vector.
#' @param offset For group i, the scatter plot is drawn between i+offset and i+1-offset.
#' @param ytitle Title y-axis.
#' @param ytitlefont Font title y-axis.
#' @param yname Name of the y variable. Used for the hovering template.
#' @param yrange Range y-axis.
#' @param ylabelsfont Font of y-axis labels.
#' @param min_median_size Minimum number of points required in the points cloud to draw the median.
#' @param colors_bkgd Vector of colors for the background. If less than the number of grounds, values are recycled using 
#'   `rep` function.
#' @param markers_scatter List of markers characteristics.
#' @param lines_median List of median line characteristics.
#' @param lwd_axes Linewidth of line delimiting y anx x axes.
#' @param col_id Name of column id for hovering.
#'
#' @import plotly
#'
#' @author Yoann Pradat
#' @export
cumulative_scatter_plots <- function(dfs, groups, stacks, offset, ytitle, ytitlefont=list(size=14),
                                     yname="Burden", yrange=c(10^-3,10^4), ylabelsfont=list(size=12),
                                     min_median_size=5,
                                     colors_bkgd=c("#f3f3f3", "#fffff"),
                                     markers_scatter=list(list(size=5, color="black")),
                                     lines_median=list(list(color="red")), lwd_axes=2, col_id="Subject_Id"){

  n_groups <- length(groups)
  n_stacks <- length(stacks)

  # repeat parameters if required
  colors_bkgd <- rep(colors_bkgd, length.out=n_groups)
  markers_scatter <- rep(markers_scatter, length.out=n_stacks)
  lines_median <- rep(lines_median, length.out=n_stacks)

  # init empty plot
  fig <- plotly_empty()
  shapes <- list()
  showlegends <- rep(T, length.out=n_stacks)

  # scatter plot and background shape
  for (i in 1:n_groups){
    for (j in 1:n_stacks){
      df_group <- dfs[[stacks[j]]][[groups[i]]]
      marker_scatter <- markers_scatter[[j]]
      line_median <- lines_median[[j]]

      if (!is.null(df_group)){
        # scatter
        fig <- fig %>% add_trace(data=df_group, 
                                 name=toupper(stacks[j]),
                                 x=~X,
                                 y=~Y, 
                                 type="scatter",
                                 marker=marker_scatter,
                                 text=df_group[[col_id]],
                                 hovertemplate = paste0(paste0("<b>", toupper(stacks[j]), " - ", groups[i], "</b>"),
                                                        paste0('<br><br>', yname, ': %{y}<br>'),
                                                        paste0(col_id, ': %{text}'),
                                                        "<extra></extra>"),
                                 showlegend=showlegends[j],
                                 mode="markers")

        if (showlegends[j]){
          showlegends[j] <- F
        }

        # median line shape
        if (nrow(df_group) >= min_median_size){
          med <- stats::median(df_group[["Y"]])
          shape <- list(type="line",
                        x0=i+offset,
                        x1=i+1-offset,
                        y0=med,
                        y1=med,
                        line=line_median)

          shapes <- c(shapes, list(shape))
        }
      }
    }

    # bgkd shape
    shape <- list(type=rect, 
                  x0=i, 
                  x1=i+1, 
                  y0=yrange[1], 
                  y1=yrange[2], 
                  line=list(color='rgba(0,0,0,0)'),
                  fillcolor=colors_bkgd[i], 
                  layer='below')

    shapes <- c(shapes, list(shape))
  }

  # add grid lines
  for (j in seq(log10(yrange)[1]+1, log10(yrange)[2]-1, 1)){
    shape <- list(type="line",
                  x0=1,
                  x1=n_groups+1,
                  y0=10^j,
                  y1=10^j,
                  line=list(color="black", dash="dot", width=0.5))

    shapes <- c(shapes, list(shape))
  }

  # yaxis configuration
  yaxis <- list(title=list(text=ytitle, font=ytitlefont, standoff=10),
                tickmode="log",
                tick0=log10(yrange)[1],
                dtick=1,
                showticklabels=T,
                tickfont=ylabelsfont,
                zeroline=F,
                showgrid=F, 
                gridcolor="black",
                mirror=T, 
                linecolor="black",
                linewidth=lwd_axes,
                range=log10(yrange),
                type="log")


  # xaxis configuration
  xaxis <- list(title="",
                zeroline=F,
                showgrid=F,
                mirror=T,
                linecolor="black", 
                linewidth=lwd_axes,
                range=c(1,n_groups+1),
                showticklabels=F)

  # assemble pieces
  fig <- fig %>% layout(yaxis=yaxis, xaxis=xaxis, shapes=shapes, showlegend=T)

  fig
}

#' Compute the groups order for the burden plot.
#'
#' @param df A \code{data.frame} containing the columns \code{col_burden} and \code{col_groups}.
#' @param col_burden Name of the column holding burden values.
#' @param col_groups Name of the column holding groups values.
#' @param groups (optional) Character vector containing names of groups to be retained. 
#' @return a character vector containing ordered group names.
#'
#' @importFrom dplyr arrange bind_rows filter group_by n pull summarize
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median
#'
#' @export
compute_groups_order_burden_plot <- function(df, col_burden, col_groups, groups_keep=NULL){ 
  if (is.null(groups_keep)){
    groups_keep <- unique(df[[col_groups]])
  }

  groups_order <- df %>% 
    filter(.data[[col_groups]] %in% groups_keep) %>%
    group_by(.data[[col_groups]]) %>% 
    summarize(Median=median(.data[[col_burden]])) %>%
    arrange(Median) %>%
    pull(.data[[col_groups]])

  groups_order
}

#' Compute the coordinates of each point in the burden plot.
#'
#' @inheritParams compute_groups_order_burden_plot
#' @param groups_order A character vector containing ordered group names.
#' @param offset (optional) The j^th scatter plot is drawn between j+offset and j+1-offset.
#' @seealso [compute_groups_order_burden_plot()]
#'
#' @importFrom dplyr arrange filter
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats median
#'
#' @export
compute_coordinates_burden_plot <- function(df, groups_order, col_burden, col_groups, offset=0.1){
  dfs_groups <- list()
  for (i in 1:length(groups_order)){
    groups <- groups_order[i]
    if (groups %in% df[[col_groups]]){
      df_groups <- df %>% filter(.data[[col_groups]] == groups)
      df_groups <- df_groups %>% arrange(.data[[col_burden]])
      df_groups[,"X"] <- seq(i + offset, i + 1-offset, length=nrow(df_groups))
      df_groups[,"Y"] <- df_groups[col_burden]
      dfs_groups[[groups]] <- df_groups
    } else {
      dfs_groups[[groups]] <- NULL
    }
  }

  dfs_groups
}


#' Draw multiple scatter plots showing cumulative distributions
#'
#' Scatter plots are organized by groups and stacks. The groups define the numbers of subplots while stacks define
#' the number of points clouds in each subplot.
#'
#' @inheritParams cumulative_scatter_plots
#' @param stacks2colors A named list with names being values in \code{stacks} and values being color names.
#' @param marker_size Size of the scatter markers.
#' @param median_width Width of the median line.
#' @param margin Named list passed to the \code{margin} parameter of \code{plotly::layout}.
#' @param xlabelsfont Font of the x-axis labels.
#' @param legendfont Font of the legend.
#'
#' @import plotly
#'
#' @author Yoann Pradat
#' @export
draw_burden_plot <- function(dfs, groups, stacks, offset, stacks2colors, marker_size=4, marker_opacity=1,
                             median_width=2, yname="Burden", yrange=c(10^-3, 10^4), min_median_size=5, ytitle="Burden",
                             margin=list(t=250), lwd_axes=2, ytitlefont=list(size=18), xlabelsfont=list(size=20),
                             ylabelsfont=list(size=20), legendfont=list(size=18), col_id="Subject_Id"){
  # plot settings 
  markers_scatter <- list()
  lines_median <- list()

  for (stack in stacks){
    marker_scatter <- list(size=marker_size, color=stacks2colors[[stack]], opacity=marker_opacity)
    line_median <- list(color=toRGB(stacks2colors[[stack]]), width=2)

    markers_scatter <- c(markers_scatter, list(marker_scatter))
    lines_median <- c(lines_median, list(line_median))
  }

  # draw
  fig <- cumulative_scatter_plots(dfs=dfs, groups=groups, stacks=stacks, offset=offset, yrange=yrange, yname=yname,
                                  ytitle=ytitle, ytitlefont=ytitlefont, colors_bkgd=c("#f3f3f3", "#fffff"),
                                  min_median_size=min_median_size, ylabelsfont=ylabelsfont,
                                  markers_scatter=markers_scatter, lines_median=lines_median, lwd_axes=lwd_axes,
                                  col_id=col_id)

  # add names
  fig <- fig %>% layout(margin=margin)
  fig <- fig %>% add_annotations(x=seq(1,length(groups)+1,1)+0.5, y=1.01,
                                 text=groups, 
                                 textangle=90,
                                 font=xlabelsfont,
                                 xanchor="centered",
                                 yanchor="bottom",
                                 xref="x",
                                 yref="paper",
                                 showarrow=F)

  fig <- fig %>% layout(legend=list(orientation="h", xanchor="center", x=0.5, y=-0.02, font=legendfont,
                                    itemsizing="constant"))

  fig
}


# ======================================================================================================================
#
# FACETTED HEATMAP BARPLOTS FIGURE
#
# ======================================================================================================================

#' Colorscale vector
#'
#' Associate to each unique value in the numeric matrix \code{z} a color using a set of colors (a palette) and a set of 
#' limits defining intervals with the same color.
#'
#' @param z a numeric matrix or vector.
#' @param colors_limits (optional) A numeric vector, used to associate to each unique value in \code{z} an interval.
#' @param colors_palette (optional) A name ("Reds" or "Blues") or a character vector of color names/codes.
#' @return a \code{data.frame}.
#'
#' @importFrom scales rescale
#' @importFrom stats setNames
#'
#' @keywords internal
get_colorscale_htmp <- function(z, colors_limits=NULL, colors_palette=NULL){
  vals <- unique(scales::rescale(c(z)))
  o <- order(vals, decreasing = FALSE)

  if (is.null(colors_limits)){
    colors_limits <- c(0, 1e-9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20,
                       0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.65, 1+1e-9)
  } else {
    if (is.null(colors_palette)){
      stop("If you use colors_limits argument, you need to also specify the colors_palette argument as a list of colors.")
    }
  }

  if (is.null(colors_palette)){
    colors_palette <- "Reds"
  }

  if (length(colors_palette)==1 & colors_palette=="Reds"){
    colors_palette <- c("#FFFFFF", "#FFF3E7", "#FFECDB", "#FFE6CF", "#FEE3C9", "#FFDEC0","#FED6B1", "#FFD0A7",
                        "#FFCA9E", "#FEC592", "#FEC18A", "#FDBC82", "#FEAD6C", "#FDA560", "#FE9A51", "#FE9548",
                        "#FB8838", "#F67F31", "#F3792D", "#EE6F26", "#E8651D", "#E66019", "#E0550D", "#D94701")
  } else if (length(colors_palette)==1 & colors_palette=="Blues"){
    colors_palette <- c("#FFFFFF", "#F7FBFF", "#F2F8FC", "#ECF5FB", "#EAF3FA", "#E4EFF7", "#DEEBF5", "#D8E8F2",
                        "#D0E2EF", "#CDDFEC", "#C2DAE9", "#BBD6E7", "#B1D0E4", "#A7CBE2", "#97C3DF", "#81B7DA",
                        "#71B0D7", "#67AAD4", "#60A2D0", "#579ACB", "#4F92C7", "#478AC3", "#3C82BF", "#2573B6")
  } else if (length(colors_palette)==1 & colors_palette=="Greens"){
    colors_palette <- c("#FFFFFF", '#F7FCF5', '#F1F8EF', '#EAF5EA', '#E4F1E4', '#DDEEDE', '#D7EAD8', '#D1E6D3',
                        '#CAE3CD', '#C4DFC7', '#BDDCC1', '#B7D8BC', '#B1D4B6', '#AAD1B0','#A4CDAA', '#9EC9A5',
                        '#97C69F', '#91C299', '#8ABF93', '#84BB8E', '#7EB788', '#77B482', '#71B07C', '#64A971')
  } else {
    if (is.null(colors_palette)){
      stop("Specify a valid value for 'colors_palette'")
    }
  }

  cols <- sapply(vals, function(val) colors_palette[findInterval(val, colors_limits)])
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)

  colz
}

#' Colorscale vector
#'
#' Associate to each unique value in the numeric matrix \code{z} a color using a set of colors (a palette) and a set of 
#' limits defining intervals with the same color.
#'
#' @param z a numeric matrix or vector.
#' @param colors_limits (optional) A numeric vector, used to associate to each unique value in \code{z} an interval.
#' @param colors_palette (optional) A name ("Reds" or "Blues") or a character vector of color names/codes.
#' @return a \code{data.frame}.
#'
#' @importFrom scales rescale
#' @importFrom stats setNames
#'
#' @keywords internal
get_colorscale_htmp <- function(z, colors_limits=NULL, colors_palette=NULL){
  vals <- unique(scales::rescale(c(z)))
  o <- order(vals, decreasing = FALSE)

  if (is.null(colors_limits)){
    colors_limits <- c(0, 1e-9, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12, 0.14, 0.16, 0.18, 0.20,
                       0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.65, 1+1e-9)
  } else {
    if (is.null(colors_palette)){
      stop("If you use colors_limits argument, you need to also specify the colors_palette argument as a list of colors.")
    }
  }

  if (is.null(colors_palette)){
    colors_palette <- "Reds"
  }

  if (length(colors_palette)==1 & colors_palette=="Reds"){
    colors_palette <- c("#FFFFFF", "#FFF3E7", "#FFECDB", "#FFE6CF", "#FEE3C9", "#FFDEC0","#FED6B1", "#FFD0A7",
                        "#FFCA9E", "#FEC592", "#FEC18A", "#FDBC82", "#FEAD6C", "#FDA560", "#FE9A51", "#FE9548",
                        "#FB8838", "#F67F31", "#F3792D", "#EE6F26", "#E8651D", "#E66019", "#E0550D", "#D94701")
  } else if (length(colors_palette)==1 & colors_palette=="Blues"){
    colors_palette <- c("#FFFFFF", "#F7FBFF", "#F2F8FC", "#ECF5FB", "#EAF3FA", "#E4EFF7", "#DEEBF5", "#D8E8F2",
                        "#D0E2EF", "#CDDFEC", "#C2DAE9", "#BBD6E7", "#B1D0E4", "#A7CBE2", "#97C3DF", "#81B7DA",
                        "#71B0D7", "#67AAD4", "#60A2D0", "#579ACB", "#4F92C7", "#478AC3", "#3C82BF", "#2573B6")
  } else if (length(colors_palette)==1 & colors_palette=="Greens"){
    colors_palette <- c("#FFFFFF", '#F7FCF5', '#F1F8EF', '#EAF5EA', '#E4F1E4', '#DDEEDE', '#D7EAD8', '#D1E6D3',
                        '#CAE3CD', '#C4DFC7', '#BDDCC1', '#B7D8BC', '#B1D4B6', '#AAD1B0','#A4CDAA', '#9EC9A5',
                        '#97C69F', '#91C299', '#8ABF93', '#84BB8E', '#7EB788', '#77B482', '#71B07C', '#64A971')
  } else {
    if (is.null(colors_palette)){
      stop("Specify a valid value for 'colors_palette'")
    }
  }

  cols <- sapply(vals, function(val) colors_palette[findInterval(val, colors_limits)])
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)

  colz
}


#' Add text annotations on heatmap cells
#'
#' @param z a numeric matrix or vector.
#' @param fig plotly figure containing the heatmap.
#' @param black_white_cutoff (optional) Cutoff for deciding on text font color.
#' @param font (optional) Text font parameters.
#' @return the figure object \code{fig}.
#'
#' @keywords internal
add_z_annotations_htmp <- function(z, fig, black_white_cutoff=0.16, font=list(size=8)){
  text_anno <- apply(z, c(1,2), function(x) {if (x==0){return("")} else {return(paste0(round(x*100, 0)))}})
  text_anno <- apply(text_anno, c(1,2), function(x) {if (x==0){return("<1")} else {return(x)}})
  text_anno <- as.vector(t(text_anno))
  x_anno <- rep(seq(0, ncol(z)-1), nrow(z))
  y_anno <- rep(seq(0, nrow(z)-1), each=ncol(z))

  mask_black <- as.vector(t(z)) < black_white_cutoff
  mask_white <- as.vector(t(z)) >= black_white_cutoff

  if (sum(mask_black)>0){
    font_black <- font
    font_black$color <- "black"
    fig <- fig %>% add_annotations(x=x_anno[mask_black], y=y_anno[mask_black], text=text_anno[mask_black],
                                   showarrow=F, ax=0, ay=0, font=font_black)
  }
  if (sum(mask_white)>0){
    font_white <- font
    font_white$color <- "white"
    fig <- fig %>% add_annotations(x=x_anno[mask_white], y=y_anno[mask_white], text=text_anno[mask_white],
                                   showarrow=F, ax=0, ay=0, font=font_white)
  }

  fig
}


#' Add exponent star symbols on row names to show pvalues of some tests. 
#'
#' @param df_p a data.frame with first column containing original row names and second columns containing a pvalue.
#' @param color hexadecimal or string color.
#' @param alpha_left threshold for showing 1 star.
#' @param y_ann vector of current row names.
#' @param y vector of original row names.
#' @return a modified vector y_ann with additional symbols if any pvalue is below the threshold.
#'
#' @keywords internal
add_symbols_pvals <- function(df_p, color, alpha_left, y_ann, y){
  y_ann_new <- c()
  for (i in 1:length(y)){
    y_i <- y[i]
    y_ann_i <- y_ann[i]
    pval_i <- abs(df_p[df_p[colnames(df_p)[1]]==y_i,colnames(df_p)[2],drop=T])

    if (pval_i < 0.001){
      pattern_i <- "***"
    } else if (pval_i < 0.01) {
      pattern_i <- "**"
    } else if (pval_i < alpha_left) {
      pattern_i <- "*"
    } else {
      pattern_i <- ""
    }

    sup_i <- paste0("<span style='color:", color, "'>", pattern_i, "</span>")

    if (pattern_i==""){
      sup_i <- ""
    } else {
      if (y_i!=y_ann_i){
        sup_i <- paste0(",",sup_i)
      } else {
        sup_i <- paste0(sup_i)
      }
    }

    if (sup_i!=""){
      y_ann_new <- c(y_ann_new, paste0(y_ann_i,'<sup>', sup_i, '</sup>')) 
    } else {
      y_ann_new <- c(y_ann_new, y_ann_i) 
    }
  }
  
  y_ann_new
}


#' Set colors of labels
#'
#' @description Using a subpart of each label (split by "-" separator and take the first part), set the color using a
#' named list of colors.
#'
#' @param x vector of labels. 
#' @param colors named list of colors.
#' @param pos position of the substring to be considered.
#' @return a modified vector x with colors set using HTML syntax.
#'
#' @importFrom stringr str_split
#' @keywords internal
add_colors_names <- function(x, colors, pos=1){
  xnames <- sapply(x, function(e) trimws(unlist(str_split(e,"-"))[pos]), USE.NAMES=F)
  x_ann <- c()
  for (i in 1:length(x)){
    x_i <- x[i]
    xname_i <- xnames[i]
    x_ann_i <- paste0("<span style='color:", colors[[xname_i]], ";font-weight:bold'>", x_i, "</span>")
    x_ann <- c(x_ann,  x_ann_i)
  }

  x_ann
}


#' Configure x and y axes of the heatmap
#'
#' @param z a numeric matrix or vector.
#' @param fig plotly figure containing the heatmap.
#' @param black_white_cutoff (optional) Cutoff for deciding on text font color.
#' @param y_tick_font (optional) List with font parameters of y-axis ticks.
#' @param x_tick_font (optional) List with font parameters of x-axis ticks.
#' @param y_side (optional) Size of labels of y-axis.
#' @param x_side (optional) Size of labels of x-axis.
#' @param linewidth (optional) Linewidth of external borders (axes).
#' @return the figure object \code{fig}
#'
#' @keywords internal
configure_axes_htmp <- function(z, fig, y_tick_font=list(size=7), x_tick_font=list(size=12), y_side="left",
                                x_side="bottom", linewidth=0.01){
  # yaxis configuration
  yaxis <- list(title="",
                showticklabels=T,
                tickfont=y_tick_font,
                dtick=1,
                side=y_side,
                ticklen=0,
                zeroline=F,
                showline=T,
                linewidth=linewidth,
                linecolor="gray",
                showgrid=F, 
                visible=T,
                mirror=T, 
                range=c(-0.5, nrow(z)-0.5))


  # xaxis configuration
  xaxis <- list(title="",
                showticklabels=T,
                tickfont=x_tick_font,
                tickangle=-90,
                dtick=1,
                side=x_side,
                ticklen=0,
                zeroline=F,
                showline=T,
                linewidth=linewidth,
                linecolor="gray",
                showgrid=F,
                mirror=T,
                visible=T,
                range=c(-0.5,ncol(z)-0.5))


  fig %>% layout(yaxis=yaxis, xaxis=xaxis)
}


#' Configure x-axis for bar plots
#'
#' @param fig plotly figure containing the heatmap.
#' @param title (optional) Title of axis.
#' @param tickfont (optional) A list to specify tick font parameters. For instance use \code{list(size=12)} to specify
#' a tick font of size 12.
#' @param titlefont (optional) A list to specify title font parameters. For instance use \code{list(size=12)} to specify
#' a title font of size 12.
#' @param ... Extra parameters added to the xaxis configuration list.
#' @return the figure object \code{fig}
#'
#' @keywords internal
configure_xaxes_row_bar <- function(fig, title="", tickfont=list(size=8), titlefont=list(size=12), ...){
  # xaxis configuration
  xaxis <- list(title=title,
                titlefont=titlefont,
                showticklabels=T,
                tickfont=tickfont,
                tickcolor="gray",
                side="top",
                tickwidth=1,
                ticklen=5,
                zeroline=F,
                showline=T,
                linewidth=1,
                linecolor="gray",
                showgrid=F,
                mirror=F,
                visible=T)
  xaxis <- modifyList(xaxis, list(...))

  fig %>% layout(xaxis=xaxis)
}


#' Configure y-axis for bar plots
#'
#' @param fig plotly figure containing the heatmap.
#' @param title (optional) Title of axis.
#' @param tickfont (optional) A list to specify tick font parameters. For instance use \code{list(size=12)} to specify
#' a tick font of size 12.
#' @param titlefont (optional) A list to specify title font parameters. For instance use \code{list(size=12)} to specify
#' a title font of size 12.
#' @param ... Extra parameters added to the xaxis configuration list.
#' @return the figure object \code{fig}.
#'
#' @keywords internal
configure_yaxes_row_bar <- function(fig, title="", tickfont=list(size=8), titlefont=list(size=12), ...){
  # yaxis configuration
  yaxis <- list(title=title,
                titlefont=titlefont,
                showticklabels=F,
                tickfont=tickfont,
                ticklen=0,
                tickcolor="gray",
                zeroline=F,
                showline=F,
                linewidth=1,
                linecolor="gray",
                showgrid=F, 
                visible=F,
                mirror=F)
  yaxis <- modifyList(yaxis, list(...))

  fig %>% layout(yaxis=yaxis)
}


#' Build list of shapes
#'
#' Draw a rectangular shape with a specific color for highlighting cell edges.
#'
#' @param n_row Number of rows of the heatmap.
#' @param n_col Number of columns of the heatmap.
#' @param color Line color of the rectangle.
#' @param width Line width of the rectangle.
#' @return a list of shapes
#'
#' @keywords internal
shapes_cell_edges <- function(n_row, n_col, color, width=1){
  shapes <- list()

  for (i in 1:n_row){
    for (j in 1:n_col){
      x0 <- j-1-0.5
      x1 <- j-0.5
      y0 <- i-1-0.5
      y1 <- i-0.5

      shape <- list(type="rect", fillcolor=NULL, line=list(color=color, width=width),
                    x0=x0, x1=x1, y0=y0, y1=y1, xref="x", yref="y")
      shapes <- c(shapes, list(shape))
    }
  }

  shapes
}

#' Build list of shapes
#'
#' Draw a shape with a specific color according to the pvalue.
#'
#' @param dfs_pvals Names list of data.frames.
#' @param names Names of the list of data.frames.
#' @param names2colors Named list of colors.
#' @param alpha pvalue threshold for deciding to add a shape or not.
#' @param width Line width of the shape
#' @return a list of shapes
#'
#' @importFrom colorspace lighten
#'
#' @keywords internal
shapes_pvalues <- function(dfs_pval, names, names2colors, alpha=0.05, width=0.5){
  n_row <- nrow(dfs_pval[[names[1]]])
  n_col <- ncol(dfs_pval[[names[1]]])
  shapes <- list()

  for (i in 1:n_row){
    for (j in 1:n_col){
      offset <- 0.5
      margin <- 0.1
      for (name in names){
        pval <- dfs_pval[[name]][i,j,drop=T]
        color <- names2colors[[name]]
        if (abs(pval) < alpha){
          x0 <- j-1-0.5
          x1 <- j-0.5
          y0 <- i-1-0.5
          y1 <- i-0.5

          if (pval > 0){
            mx <- x0 + 0.75*(x1-x0)
            my <- y0 + 0.5*(y1-y0) - offset
            L1x <- mx + (x1-mx)/2
            L1y <- y1 - offset
            L2x <- x1
            L2y <- my
          } else {
            mx <- x0 + 0.75*(x1-x0)
            my <- y1 - offset
            L1x <- mx + (x1-mx)/2
            L1y <- y0 + 0.5*(y1-y0) - offset 
            L2x <- x1
            L2y <- my
          }
          shape <- list(type="path", path=paste(" M", mx, my, "L", L1x, L1y, "L", L2x, L2y,  "Z"),
                        fillcolor=lighten(color, amount=0.75), line=list(color=color, width=width))
          # shape <- list(type="rect", fillcolor=NULL, line=list(color=color, width=2),
          #               x0=j-1-0.5+offset, x1=j-0.5-offset, y0=i-1-0.5+offset, y1=i-0.5-offset,
          #               xref="x", yref="y")
          # offset <- offset + 0.1
          shapes <- c(shapes, list(shape))
        }
        offset <- offset - 0.5
      }
    }
  }

  shapes
}


#' Draw a facetted figure with template 1.
#'
#' @description
#' The main plot in the figure is the central heatmap. This heatmap serves to show the distribution of counts or
#' percentages of an event in a fixed number of variables and fixed number of samples.
#'
#' The figure may additionally contain a
#' \itemize{
#'   \item{a left double inverted barplot} {it is used to display -log10(pvalues) distribution from 
#' a per-variable test (e.g Fisher stratified by the samples).} 
#'   \item{a middle right barplot} {it is used to display total variable counts.} 
#'   \item{an extreme right barplot} {it is used to display additional information about a variable in the form of a
#'   stacked barplot. It is useful to show the breakdown of a total variable count into types.} 
#' }
#' 
#' @param dfs a named list of dataframes containing the values indicated in the list \code{names2plots}.
#' @param col_var A character vector defining the field used as variables.
#' @param names2plots a list with names among 'left_pvals', 'heatmap', 'heatmap_hover', 'heatmap_pvals',
#'   'middle_right_bar', 'extreme_right_bar'. The corresponding values must be in the list of names of \code{dfs}.
#' @param names2colors a list with names that must match values of the vectors names2plots$left_pvals and
#'    names2plots$heatmap_pvals.
#' @param colors_palette_heatmap (optional) The color palette to be used for the heatmap. You may use "Reds", "Blues"
#'   or a list of colors that matches the length of color_limits_heatmap.
#' @param colors_limits_heatmap (optional) The limit values for deciding on cell colors in the heatmap.
#' @param width_one (optional) The width of the plot is obtained as \code{width_one} multiplied by the number of columns
#'   in \code{z}.
#' @param height_one (optional) The height of the plot is obtained as \code{height_one} multiplied by the number of rows 
#'   in \code{z}.
#' @param alpha_left (optional) Threshold for drawing the bar for a significant variable pvalue.
#' @param alpha_heatmap (optional) Threshold for highlighting edges of a significant cell pvalue.
#' @param col_stack (optional) A character vector defining the field used for the extreme right stack barplot. Used
#'   only if 'extreme_right_bar' in the names of \code{names2plots}.
#' @param width_edges_heatmap (optional) Width of the line defining the edges of the heatmap.
#' @param fonts (optional) List of fonts. Required names are 'x_tick_heatmap', 'y_tick_heatmap', 'x_tick_rowbar'
#'   and 'legend'. Each font is a list with font parameters, e.g, 'size'.
#' @param bargap (optional) Spacing between bars of barplot. 
#' @param stacks2colors (optional) A list of colors with names corresponding to values for the extreme right stack
#'   barplot. Used only if 'extreme_right_bar' in the names of \code{names2plots}.
#' @param showlegend (optional) Set to FALSE to hide the legend.
#' @param add_colors_names (optional) Set to FALSE to not use colors in x-axis labels.
#' @param x_tick_heatmap_side (optional) Choose 'bottom' or 'top'.
#' @param x_title_row_bar (optional) Title for the x-axis of the side bar plot.
#' @param add_cohorts_sizes (optional) Set to TRUE to show cohort sizes above or below x tick labels. 
#' @return a list of plotly figure objects.
#'
#' @import plotly
#'
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @importFrom tidyr spread
#' @importFrom rlang .data
#'
#' @author Yoann Pradat
#' @export
draw_facetted_heatmap_barplots_1 <- function(dfs, col_var, names2plots, names2colors,
                                        colors_palette_heatmap="Reds", colors_limits_heatmap=NULL,
                                        col_pval="p.value", width_one=100, height_one=10, alpha_left=0.1, 
                                        alpha_heatmap=0.1, col_stack=NULL, width_edges_heatmap=0.05,
                                        fonts=list(z_heatmap=list(size=8),
                                                   x_tick_heatmap=list(size=10),
                                                   y_tick_heatmap=list(size=6),
                                                   x_tick_row_bar=list(size=8),
                                                   x_title_row_bar=list(size=10),
                                                   cohort_size=list(size=10),
                                                   legend=list(size=9)),
                                        bargap=0.1, stacks2colors=NULL, showlegend=TRUE, add_colors_names=TRUE,
                                        x_tick_heatmap_side="bottom", x_title_row_bar="", add_cohorts_sizes=F){
  # Shared x (samples) and y vectors (samples).
  y <- as.vector(dfs[[names2plots[["heatmap"]]]] %>% pull(col_var))
  x <- as.vector(dfs[[names2plots[["heatmap"]]]] %>% select(-c(all_of(col_var))) %>% names())
  if (add_colors_names){
    x_ann <- add_colors_names(x, names2colors[["colnames"]])
  } else {
    x_ann <- x
  }
  figs <- list()

  # # Left stratified pvalue bar plot
  # name <- "left_pvals"

  # if (name %in% names(names2plots) & !is.null(names2plots$left_pvals)){
  #   name_a <- names2plots$left_pvals[[1]]
  #   name_b <- names2plots$left_pvals[[2]]
  #   df_a <- dfs[[name_a]]
  #   df_b <- dfs[[name_b]]

  #   x_a <- df_a %>% mutate(Log_P=ifelse(.data[[col_pval]] <= alpha_left, -log10(.data[[col_pval]]), 0)) %>% pull(Log_P)
  #   x_b <- df_b %>% mutate(Log_P=ifelse(.data[[col_pval]] <= alpha_left, -log10(.data[[col_pval]]), 0)) %>% pull(Log_P)
  #   x_max <- max(max(x_a), max(x_b))
  #   hovertemplate <- paste0("<b>%{y}</b><br>",
  #                           '-log10(q-value): %{x:.3f}<br>',
  #                           "<extra></extra>")

  #   ## fig a
  #   fig_a <- plot_ly(x=x_a, y=y, type="bar", name=toupper(names2colors[[name_a]]),
  #                    orientation="h", marker=list(color=names2colors[[name_a]],
  #                                                 line=list(color="black",width=0)),
  #                    hovertemplate=hovertemplate, showlegend=F)

  #   if (max(x_a)==0){
  #     fig_a <- configure_xaxes_row_bar(fig_a, tickmode="array", tickvals=c(), ticktext=c(), range=c(x_max,0),
  #                                      tickfont=fonts$x_tick_row_bar)
  #   } else{
  #     fig_a <- configure_xaxes_row_bar(fig_a, tickmode="array", tickvals=c(max(x_a)), ticktext=c(round(max(x_a), 1)),
  #                                      range=c(x_max, 0), tickfont=fonts$x_tick_row_bar)
  #   }

  #   fig_a <- configure_yaxes_row_bar(fig_a, zeroline=F, visible=T, range=c(-0.5, nrow(y)-0.5), showticklabels=F,
  #                                    side="left", dtick=1, tickfont=fonts$y_tick_heatmap),
  #                                    linewidth=0, linecolor="gray", categoryorder="array", categoryarray=y)
  #   fig_a <- fig_a %>% layout(bargap=bargap)

  #   ## fig b
  #   fig_b <- plot_ly(x=x_b, y=y, type="bar", name=toupper(names2colors[[name_b]]),
  #                    orientation="h", marker=list(color=names2colors[[name_b]],
  #                                                 line=list(color="black",width=0)),
  #                    hovertemplate=hovertemplate,
  #                    showlegend=F)

  #   if (max(x_b)==0){
  #     fig_b <- configure_xaxes_row_bar(fig_b, tickmode="array", tickvals=c(), ticktext=c(), range=c(0, x_max),
  #                                      tickfont=fonts$x_tick_row_bar)
  #   } else{
  #     fig_b <- configure_xaxes_row_bar(fig_b, tickmode="array", tickvals=c(max(x_b)), ticktext=c(round(max(x_b), 1)),
  #                                      range=c(0,x_max), tickfont=fonts$x_tick_row_bar)
  #   }

  #   fig_b <- fig_b %>% layout(yaxis=list(zeroline=T, visible=T, range=c(-0.5, nrow(y)-0.5), linewidth=0.1,
  #                                        showticklabels=F, categoryorder="array", categoryarray=y), bargap=bargap)

  #   figs[[name]] <- subplot(fig_a, fig_b, margin=0.0)
  # }

  # Heatmap
  name <- "heatmap"
  name_hover <- "heatmap_hover"
  if (name %in% names(names2plots)){
    df <- dfs[[names2plots[[name]]]]
    df_hover <- dfs[[names2plots[[name_hover]]]]

    z <- as.matrix(df %>% select(-c(all_of(col_var))))
    text_hover <- as.matrix(df_hover %>% select(-c(all_of(col_var))))

    ## colors heatmap
    colz <- get_colorscale_htmp(z, colors_palette=colors_palette_heatmap, colors_limits=colors_limits_heatmap)

    ## if possible, add symbols to gene names for indicating significant stratified fisher tests
    y_ann <- c()
    if ("left_pvals" %in% names(names2plots) & !is.null(names2plots$left_pvals)){
      name_a <- names2plots$left_pvals[[1]]
      name_b <- names2plots$left_pvals[[2]]
      df_a <- dfs[[name_a]]
      df_b <- dfs[[name_b]]
      y_ann <- c()
      for (y_i in y){
        pval_a_i <- abs(df_a[df_a[colnames(df_a)[1]]==y_i,colnames(df_a)[2],drop=T])
        pval_b_i <- abs(df_b[df_b[colnames(df_b)[1]]==y_i,colnames(df_b)[2],drop=T])

        if (pval_a_i < 0.001){
          pattern_a_i <- "***"
        } else if (pval_a_i < 0.01) {
          pattern_a_i <- "**"
        } else if (pval_a_i < alpha_left) {
          pattern_a_i <- "*"
        } else {
          pattern_a_i <- ""
        }

        if (pval_b_i < 0.001){
          pattern_b_i <- "***"
        } else if (pval_b_i < 0.01) {
          pattern_b_i <- "**"
        } else if (pval_b_i < alpha_left) {
          pattern_b_i <- "*"
        } else {
          pattern_b_i <- ""
        }

        sup_a_i <- paste0("<span style='color:", names2colors[[name_a]], "'>", pattern_a_i, "</span>")
        sup_b_i <- paste0("<span style='color:", names2colors[[name_b]], "'>", pattern_b_i, "</span>")

        if (pattern_a_i=="" & pattern_b_i==""){
          sup_i <- ""
        } else if (pattern_a_i!="" & pattern_b_i!=""){
          sup_i <- paste0(sup_a_i,",",sup_b_i)
        } else if (pattern_a_i!=""){
          sup_i <- sup_a_i
        } else {
          sup_i <- sup_b_i
        }

        y_ann <- c(y_ann, paste0(y_i,'<sup>', sup_i, '</sup>')) 
      }
    } else {
      y_ann <- y
    }

    ## plot heatmap
    fig_htmp <- plot_ly(z=z, x=x_ann, y=y_ann, text=text_hover, type="heatmap",
                        colorscale=colz,
                        showscale=F,
                        hovertemplate=paste0(paste0('<b>%{x}', ' - ', "%{y}</b><br>"),
                                             'Percent: %{z}<br>',
                                             'Count: %{text}',
                                             "<extra></extra>"),
                        height=height_one*nrow(z), width=width_one*ncol(z))

    ## annotation text heatmap
    fig_htmp <- add_z_annotations_htmp(z, fig_htmp, black_white_cutoff=0.16, font=fonts$z_heatmap)

    ## axes configuration
    fig_htmp <- configure_axes_htmp(z, fig_htmp, x_tick_font=fonts$x_tick_heatmap,
                                    y_tick_font=fonts$y_tick_heatmap,
                                    x_side=x_tick_heatmap_side, y_side="right")
    
    ## additional shapes
    shapes <- c()

    ### shapes cell edges
    shapes_edges <- shapes_cell_edges(n_row=nrow(z), n_col=ncol(z), color="gray", width=width_edges_heatmap)
    shapes <- c(shapes, shapes_edges)

    ### shapes pvalues
    if ("heatmap_pvals" %in% names(names2plots)){
      dfs_pvals <- dfs[names2plots[["heatmap_pvals"]]]
      dfs_pvals <- lapply(dfs_pvals, function(df) {df %>% select(-c(all_of(col_var)))})
      shapes_pvals <- shapes_pvalues(dfs_pval=dfs_pvals, names=names2plots[["heatmap_pvals"]],
                                     names2colors=names2colors[names2plots[["heatmap_pvals"]]],
                                     alpha=alpha_heatmap)

      shapes <- c(shapes, shapes_pvals)
    }

    ### add shapes to figure
    fig_htmp <- fig_htmp %>% layout(shapes=shapes)

    figs[[name]] <- fig_htmp
  }

  # middle right bar
  name <- "middle_right_bar"

  if (name %in% names(names2plots)){
    percent <- dfs[[names2plots[[name]]]]$Percent
    text <- dfs[[names2plots[[name]]]]$Count

    fig_middle_right <- plot_ly(x=percent*100, y=y, type="bar", text=text,
                                orientation="h",
                                marker=list(color="#BC80BD", line=list(color="black",width=0)),
                                textposition="none",
                                hovertemplate=paste0("<b>%{y}</b><br>",
                                                     'Percent: %{x:.3f}<br>',
                                                     'Count: %{text}',
                                                     "<extra></extra>"),
                                showlegend=F)
    fig_middle_right <- configure_xaxes_row_bar(fig_middle_right, title=x_title_row_bar,
                                                titlefont=fonts$x_title_row_bar,
                                                tickfont=fonts$x_tick_row_bar, nticks=3)
    fig_middle_right <- configure_yaxes_row_bar(fig_middle_right, range=c(-0.5, nrow(y)-0.5),
                                                categoryorder="array", categoryarray=y)

    fig_middle_right <- fig_middle_right %>% layout(bargap=bargap)
    figs[[name]] <- fig_middle_right
  }

  # extreme right bar
  name <- "extreme_right_bar"

  if (name %in% names(names2plots) & !is.null(names2plots[[name]])){
    fraction_row <- dfs[[names2plots[[name]]]]
    fraction <- fraction_row %>%
      select(-c(Count)) %>%
      spread({{col_stack}}, Percent) %>%
      replace(is.na(.), 0) %>% arrange(match(.data[[col_var]], y))
    fraction$`0` <- NULL

    text <- fraction_row %>%
      select(-c(Percent)) %>%
      spread({{col_stack}}, Count) %>%
      replace(is.na(.), 0) %>% arrange(match(.data[[col_var]], y))
    text$`0` <- NULL

    stacks <- fraction %>% select(-c(all_of(col_var))) %>% names()

    fig_x_right <- plotly_empty()
    for (stack in stacks){
      fig_x_right <- fig_x_right %>% add_trace(x=fraction[[stack]], y=y, type="bar", orientation="h",
                                               name=stack, text=text[[stack]],
                                               textposition="none",
                                               hovertemplate=paste0("<b>%{y}</b><br>",
                                                                    paste0('Alteration: ', stack, '<br>'),
                                                                    'Percent: %{x:.3f}<br>',
                                                                    'Count: %{text}<br>',
                                                                    "<extra></extra>"),
                                               marker = list(color=stacks2colors[[stack]],
                                                             line=list(color="black", width=0)))
    }
    fig_x_right <- configure_xaxes_row_bar(fig_x_right, range=c(0,1), tickfont=fonts$x_tick_row_bar)
    fig_x_right <- configure_yaxes_row_bar(fig_x_right, range=c(-0.5, nrow(y)-0.5), categoryorder="array",
                                           categoryarray=y)
    fig_x_right <- fig_x_right %>% layout(barmode='stack', bargap=bargap, legend=list(font=fonts$legend),
                                          showlegend=showlegend)

    figs[[name]] <- fig_x_right
  }

  # add cohort sizes
  if (add_cohorts_sizes){
    name <- "count_col"
    rows <- rownames(dfs[[names2plots[[name]]]])
    vals <- as.vector(dfs[[names2plots[[name]]]]$Count)
    ords <- colnames(dfs[[names2plots[["heatmap"]]]])
    ords <- ords[2:length(ords)]
    mask <- match(ords, rows)
    rows <- rows[mask]
    vals <- vals[mask]
    if (add_colors_names){
      rows_ann <- add_colors_names(rows, names2colors[["colnames"]])
    } else {
      rows_ann <- rows
    }

    anns <- list()
    for (i in 1:length(rows)){
      if (x_tick_heatmap_side=="top"){
        y <- 1.01
        yanchor <- "bottom"
      } else if (x_tick_heatmap_side=="bottom") {
        y <- -0.01
        yanchor <- "top"
      }

      ann_i <- list(x=i-1,
                    y=y,
                    text=gsub(rows[i], paste0(vals[i]), rows_ann[i]),
                    textangle=-90,
                    xref="x",
                    yref="paper",
                    xanchor="middle",
                    yanchor=yanchor,
                    showarrow=F,
                    font=fonts$cohort_size)
      anns <- c(anns, list(ann_i))
    }
    
    figs[["heatmap"]] <- figs[["heatmap"]] %>% layout(annotations=anns)
    figs[["heatmap"]] <- figs[["heatmap"]] %>% layout(xaxis=list(ticklen=25, side=x_tick_heatmap_side,
                                                                 tickcolor="white"))
  }

  figs
}


#' Draw the legend manually
#'
#' @param labels2colors named list of legend colors and labels.
#' @param x_min x position at where legend should start.
#' @param x_max x position at where legend should end.
#' @param y_min y position at where legend should start.
#' @param y_max y position at where legend should end.
#' @param x_gap gap along x-axis between 2 legend items.
#' @param y_gap gap along y-axis between 2 legend items.
#' @param orientation choose "v" for vertical and "h" for horizontal.
#' @param x_to_y_ratio ratio for controlling the aspect of legend rectangles.
#' @param xref see shapes parameters in plotly.
#' @param yref see shapes parameters in plotly.
#' @param font_labels names list of font characteristics.
#' @return a list of shapes and annotations representing the legend.
#'
#' @author Yoann Pradat
#' @keywords internal
draw_legend_manual <- function(labels2colors, x_min, x_max, y_min, y_max, x_gap=NULL, y_gap=NULL, orientation="v",
                               x_to_y_ratio=1, xref="paper", yref="paper", font_labels=list(size=14)){
  shapes <- list()
  n_labs <- length(labels2colors)
  x_range <- x_max-x_min
  y_range <- y_max-y_min
  if (is.null(x_gap)) x_gap <- 0.05 * x_range
  if (is.null(y_gap)) y_gap <- 0.05 * y_range

  stopifnot(x_gap <= x_range)
  stopifnot(y_gap <= y_range)

  x0_labs <- c()
  x1_labs <- c()
  y0_labs <- c()
  y1_labs <- c()

  if (orientation=="v"){
    y_lab <- (y_range - (n_labs-1)*y_gap)/n_labs
    x_lab <- y_lab * x_to_y_ratio
    for (j in 1:n_labs){
      x0_labs <- c(x0_labs, x_min)
      x1_labs <- c(x1_labs, x_min + x_lab)
      y0_labs <- c(y0_labs, y_min + (y_lab + y_gap)*(j-1))
      y1_labs <- c(y1_labs, y_min + (y_lab + y_gap)*j - y_gap)
    }

    x_ann <- x0_labs + x_lab
    y_ann <- seq(y_min + y_lab/2, y_max - y_lab/2, length.out=n_labs)
    xanchor <- "left"
    yanchor <- "middle"
  } else if (orientation=="h"){
    x_lab <- (x_range - (n_labs-1)*x_gap)/n_labs
    y_lab <-  x_lab / x_to_y_ratio
    for (j in 1:n_labs){
      x0_labs <- c(x0_labs, x_min + (x_lab + x_gap)*(j-1))
      x1_labs <- c(x1_labs, x_min + (x_lab + x_gap)*j - x_gap)
      y0_labs <- c(y0_labs, y_min)
      y1_labs <- c(y1_labs, y_min + y_lab)
    }

    x_ann <- seq(x_min + x_lab/2, x_max - x_lab/2, length.out=n_labs)
    y_ann <- y0_labs - y_lab
    xanchor <- "middle"
    yanchor <- "top"
  }

  # legend symbols
  for (j in 1:n_labs){
    shape_label <- list(type="rect", fillcolor=labels2colors[[j]], line=list(color=labels2colors[[j]], width=0),
                        x0=x0_labs[j], x1=x1_labs[j], y0=y0_labs[j], y1=y1_labs[j],
                        xref=xref, yref=yref)

    shapes <- c(shapes, list(shape_label))
  }

  # legend annotations
  annotations <- list(x = x_ann,
                      y = y_ann,
                      text = names(labels2colors),
                      xref = xref,
                      yref = yref,
                      xanchor = xanchor,
                      yanchor = yanchor,
                      showarrow = F,
                      font = font_labels)


  list(shapes=shapes, annotations=annotations)
}


#' Draw a stacked barplot to be integrated above a heatmap.
#'
#' @inheritParams draw_facetted_heatmap_barplots_2
#' @param df_plot a data.frame with the stacks defined in the first column, the column labels defined in the column
#'    `col_x` and the values to be stacked in the `col_y` column.
#' @param stacks2colors a named list of colors.
#' @param x a vector of ordered x-axis labels matching values from `col_x` column.
#' @param x_ann a vector of ordered x-axis labels that do not need to match values from `col_x` column.
#' @param y a vector of ordered y-axis labels.
#' @param x_ann a vector of x-axis modified labels to be displayed as x label names.
#' @param legend_title title for the legend..
#' @param legend_x domain along the x-axis of the barplot in paper ref where the legend will be drawn.
#' @param legend_y domain along the y-axis of the barplot in paper ref where the legend will be drawn.
#' @param legend_x_to_y_ratio ratio controlling the aspect of legend items.
#' @param legend_y_gap gap between legend items along the y-axis.
#' @param legend_x_gap gap between legend items along the x-axis.
#' @param legend_orientation choose 'h' or 'v'.
#' @param xaxis_showticklabels set to T to display x label names (using `x_ann` if not NULL or `x` otherwise).
#' @return a plotly figure.
#'
#' @author Yoann Pradat
#' @keywords internal
draw_barplot_above_htmp <- function(df_plot, stacks2colors, fonts, linewidth, x=NULL, x_ann=NULL,
                                    xaxis_showticklabels=F, legend_title="", legend_x=c(-0.2, -0.2),
                                    legend_y=c(0.1, 0.9), legend_x_to_y_ratio=0.3, legend_y_gap=0.06, legend_x_gap=NULL,
                                    legend_orientation="v"){

  col_x <- "X"
  col_y <- "Y"

  if (is.null(x)) x <- colnames(df_plot)[2:ncol(df_plot)]
  if (is.null(x_ann)) x_ann <- x

  col_stack <- colnames(df_plot)[[1]]
  fraction <- df_plot %>% gather(key=!!col_x, value=!!col_y, -all_of(col_stack)) %>%
    spread({{col_stack}}, {{col_y}}) %>% arrange(match(.data[[col_x]], x))

  stacks <- fraction %>% select(-all_of(col_x)) %>% names()
  stacks_order <- names(stacks2colors)
  stacks_order <- stacks_order[which(stacks_order %in% stacks)]

  # init plot
  fig <- plotly_empty()

  # add bar trace for each stack
  for (stack in stacks_order){
    fig <- fig %>% add_trace(x=x_ann, y=fraction[[stack]]*100, type="bar", orientation="v",
                             name=stack, text=NULL,
                             textposition="none",
                             hovertemplate=paste0("<b>%{x}</b><br>",
                                                  paste0('Alteration: ', stack, '<br>'),
                                                  'Percent: %{y:.3f}<br>',
                                                  "<extra></extra>"),
                             marker = list(color=stacks2colors[[stack]],
                                           line=list(color="black", width=0)))
  }
  
  # x-axis configuration
  fig <- configure_xaxes_row_bar(fig, range=c(-0.5, length(x)-0.5), categoryorder="array", categoryarray=x_ann,
                                 showticklabels=xaxis_showticklabels, showline=T, visible=T, ticklen=0, mirror=T,
                                 linewidth=linewidth, tickfont=fonts$x_tick_row_bar, tickangle=-90)

  # y-axis configuration
  fig <- configure_yaxes_row_bar(fig, range=c(0,100), 
                                 showticklabels=T, showline=T, visible=T, title="Samples (%)", mirror=T,
                                 linewidth=linewidth, ticklen=5, tickwidth=linewidth, 
                                 tickfont=fonts$y_tick_row_bar, 
                                 titlefont=fonts$y_title_row_bar)

  # stack bar traces
  fig <- fig %>% layout(barmode='stack', bargap=0.1,showlegend=F)

  # draw legend manually
  legend_manual <- draw_legend_manual(labels2colors=rev(stacks2colors),
                                      y_min=legend_y[1], y_max=legend_y[2], y_gap=legend_y_gap,
                                      x_min=legend_x[1], x_max=legend_x[2], x_gap=legend_x_gap,
                                      x_to_y_ratio=legend_x_to_y_ratio,
                                      orientation=legend_orientation,
                                      font_labels=fonts$legend_label)

  annotation_title <- list(x=legend_x,
                           y=0.5,
                           text=legend_title,
                           textangle=-90,
                           xref="paper",
                           yref="paper",
                           xanchor="right",
                           yanchor="middle",
                           showarrow=F,
                           font=fonts$legend_title)

  fig <- fig %>% layout(shapes=legend_manual$shapes, annotations=legend_manual$annotations)
  fig <- fig %>% layout(annotations=annotation_title)

  fig
}


#' Draw a facetted figure with template 2.
#'
#' @description
#' Draw a facetted figure having a heatmap in its center. This heatmap serves to show the distribution of counts or
#' percentages of events (rows) per groups of observations (columns). The figure may additionally contain barplots 
#' positioned above or below the heatmap.
#' 
#' @param dfs a named list of dataframes containing the values indicated in the list \code{names2plots}.
#' @param col_var A character vector defining the field used as variables.
#' @param names2plots a list with names among 'heatmap', 'heatmap_color', 'heatmap_hover', 'left_pvals',
#'   'heatmap_pvals', 'barplot_top', 'barplot_bot'. The corresponding values must be in the list of names of \code{dfs}.
#' @param names2colors a list with names that must match values of the vectors names2plots$left_pvals and
#'   names2plots$heatmap_pvals.
#' @param linewidth (optional) Linewidth of external borders (axes).
#' @param width_one (optional) The width of the plot is obtained as \code{width_one} multiplied by the number of columns
#'   in \code{z}.
#' @param height_one (optional) The height of the plot is obtained as \code{height_one} multiplied by the number of rows 
#'   in \code{z}.
#' @param alpha_left (optional) Threshold for drawing the bar for a significant variable pvalue.
#' @param alpha_heatmap (optional) Threshold for highlighting edges of a significant cell pvalue.
#' @param width_edges_heatmap (optional) Width of the line defining the edges of the heatmap.
#' @param legend_titles (optional) Titles for legends left of barplots.
#' @param legend_x (optional) Position on the x-axis of the legend next to the heatmap barplots.
#' @param legend_x_to_y_ratio (optional) Ratio for controlling the aspect of legend items.
#' @param legend_y_gap (optional) For controlling vertical gap between legend items.
#' @param fonts (optional) List of fonts. Required names are 'x_tick_heatmap', 'y_tick_heatmap', 'x_tick_rowbar'
#'   and 'legend'. Each font is a list with font parameters, e.g, 'size'.
#' @param showlegend (optional) Set to FALSE to hide the legend.
#' @param add_colors_names (optional) Set to FALSE to not use colors in x-axis labels.
#' @return a list of plotly figure objects.
#'
#' @import plotly
#'
#' @importFrom dplyr pull
#' @importFrom tidyselect all_of
#' @importFrom tidyr spread
#' @importFrom rlang .data
#'
#' @author Yoann Pradat
#' @export
draw_facetted_heatmap_barplots_2 <- function(dfs, col_var, names2plots, names2colors, linewidth=2, width_one=100,
                                             height_one=10, alpha_left=0.1, alpha_heatmap=0.1, 
                                             width_edges_heatmap=0.05,
                                             legend_titles=c("ESCAT levels", "# of actionable<br>alterations"),
                                             legend_x=-0.25, legend_x_to_y_ratio=0.2, legend_y_gap=0.06,
                                             fonts=list(z_heatmap=list(size=8),
                                                        x_tick_heatmap=list(size=10),
                                                        y_tick_heatmap=list(size=10),
                                                        x_tick_row_bar=list(size=10),
                                                        y_title_row_bar=list(size=12),
                                                        cohort_size=list(size=10),
                                                        legend=list(size=12),
                                                        legend_label=list(size=12),
                                                        legend_title=list(size=16)),
                                             add_colors_names=TRUE){

  # Shared x and y vectors
  y <- as.vector(dfs[[names2plots[["heatmap"]]]] %>% pull(col_var))
  x <- as.vector(dfs[[names2plots[["heatmap"]]]] %>% select(-c(all_of(col_var))) %>% names())
  if (add_colors_names){
    x_ann <- add_colors_names(x, names2colors[["colnames"]])
  } else {
    x_ann <- x
  }
  figs <- list()

  # Heatmap
  name <- "heatmap"
  name_color <- "heatmap_color"
  name_hover <- "heatmap_hover"

  if (name %in% names(names2plots)){
    df <- dfs[[names2plots[[name]]]]
    df_color <- dfs[[names2plots[[name_color]]]]
    df_hover <- dfs[[names2plots[[name_hover]]]]

    z <- as.matrix(df %>% select(-c(all_of(col_var))))
    text_hover <- as.matrix(df_hover %>% select(-c(all_of(col_var))))

    z_color <- as.matrix(df_color %>% select(-c(all_of(col_var))))
    colors <- names2colors[["heatmap"]]
    names_in_z <- names(colors)[names(colors) %in% unique(as.vector(z_color))]
    colors_in_z <- colors[names_in_z]
    for (i in 1:length(colors_in_z)){
      name_i <- names(colors_in_z)[i]
      color_i <- colors_in_z[[name_i]]
      z_color[z_color==name_i] <- i/length(colors_in_z)
    }
    z_color <- apply(z_color, c(1,2), as.numeric)


    ## if possible, add symbols to gene names for indicating significant stratified fisher tests
    y_ann <- y
    if ("left_pvals" %in% names(names2plots)){
      for (name_p in names2plots[["left_pvals"]]){
        y_ann <- add_symbols_pvals(df_p=dfs[[name_p]], color=names2colors[[name_p]], alpha_left=alpha_left,
                                   y_ann=y_ann, y=y)
      }
    }

    ## init heatmap with colors
    col_1 <- c(rbind((1:length(colors_in_z)-1)/length(colors_in_z), (1:length(colors_in_z))/length(colors_in_z)))
    col_2 <- c(rbind(as.vector(unlist(colors_in_z)), as.vector(unlist(colors_in_z))))
    colorscale <- setNames(data.frame(col_1, col_2), NULL)
    fig_htmp <- plot_ly(z=z_color, x=x_ann, y=y_ann, text=text_hover, type="heatmap",
                        colorscale=colorscale,
                        showscale=F,
                        xgap=0.1,
                        ygap=0.1,
                        hovertemplate=paste0(paste0('<b>%{x}', ' - ', "%{y}</b><br>"),
                                             'Count: %{text}',
                                             "<extra></extra>"),
                        height=height_one*nrow(z), width=width_one*ncol(z)+500)

    ## annotation text heatmap
    fig_htmp <- add_z_annotations_htmp(z, fig_htmp, black_white_cutoff=0, font=fonts$z_heatmap)

    ## axes configuration
    fig_htmp <- configure_axes_htmp(z, fig_htmp, x_tick_font=fonts$x_tick_heatmap,
                                    y_tick_font=fonts$y_tick_heatmap, linewidth=linewidth)
    
    ## additional shapes
    shapes <- c()

    ### shapes cell edges
    shapes_edges <- shapes_cell_edges(n_row=nrow(z), n_col=ncol(z), color="gray", width=width_edges_heatmap)
    shapes <- c(shapes, shapes_edges)

    ### shapes pvalues
    if ("heatmap_pvals" %in% names(names2plots)){
      dfs_pvals <- dfs[names2plots[["heatmap_pvals"]]]
      dfs_pvals <- lapply(dfs_pvals, function(df) {df %>% select(-c(all_of(col_var)))})
      shapes_pvals <- shapes_pvalues(dfs_pval=dfs_pvals, names=names2plots[["heatmap_pvals"]],
                                     names2colors=names2colors[names2plots[["heatmap_pvals"]]],
                                     alpha=alpha_heatmap)

      shapes <- c(shapes, shapes_pvals)
    }

    ### add shapes to figure
    fig_htmp <- fig_htmp %>% layout(shapes=shapes, showlegend=F)

    figs[[name]] <- fig_htmp
  }

  # percent top and bot
  names_bars <- c("barplot_top", "barplot_bot")
  xaxis_showticklabels <- c(T, F)

  for (j in 1:length(names_bars)){
    figs[[names_bars[j]]] <- draw_barplot_above_htmp(df_plot=dfs[[names2plots[[names_bars[j]]]]], 
                                                     stacks2colors=names2colors[[names_bars[j]]], 
                                                     x=x, x_ann=x_ann, xaxis_showticklabels=xaxis_showticklabels[j],
                                                     fonts=fonts, linewidth=linewidth,
                                                     legend_title=legend_titles[j], legend_x=c(legend_x, legend_x),
                                                     legend_x_to_y_ratio=legend_x_to_y_ratio, legend_y_gap=legend_y_gap)
  }

  # add cohort sizes above top barplot
  name <- "count_col"
  rows <- rownames(dfs[[names2plots[[name]]]])
  vals <- as.vector(dfs[[names2plots[[name]]]])
  rows_ann <- add_colors_names(rows, names2colors[["colnames"]])

  anns <- list()
  for (i in 1:length(rows)){
    ann_i <- list(x=i-1,
                  y=1.05,
                  text=gsub(rows[i], paste0(vals[i]), rows_ann[i]),
                  textangle=-90,
                  xref="x",
                  yref="paper",
                  xanchor="middle",
                  yanchor="bottom",
                  showarrow=F,
                  font=fonts$cohort_size)
    anns <- c(anns, list(ann_i))
  }
  
  figs[["barplot_top"]] <- figs[["barplot_top"]] %>% layout(annotations=anns)
  figs[["barplot_top"]] <- figs[["barplot_top"]] %>% layout(xaxis=list(ticklen=30, side="top", tickcolor="white"))

  figs
}


#' Prepare tables for the facetted heatmap barplots
#'
#' Return a list of dataframes with events counts and percentages.
#'
#' @param df_evt_count The event counts aggreagted by tumor type.
#' @param df_tt_count Table of tumor type counts. It must contain the columns `Count` and col_tt.
#' @param tt_keep List of tumor types included in the plot.
#' @param col_evt Name of the column containing the events.
#' @param col_tt Name of the column containing tumor types.
#' @param col_evt_classes (optional) Name of the column containing subclasses of events.
#' @return a list of tables.
#'
#' @export
get_tables_for_facetted_heatmap_barplots <- function(df_evt_count, df_tt_count, tt_keep, col_evt, col_tt, 
                                                     col_evt_classes=NULL){
  # size of each tumor type
  df_count_tt <- df_tt_count %>%
    select(.data[[col_tt]], Count) %>%
    replace(is.na(.), 0) %>%
    column_to_rownames(var=col_tt)

  # counts evt
  df_count_evt <- df_evt_count %>%
    group_by(.data[[col_tt]], .data[[col_evt]]) %>%
    summarize(Count=sum(Count), .groups="keep") %>% 
    spread(.data[[col_tt]], Count) %>%
    replace(is.na(.), 0) %>%
    ungroup(.data[[col_evt]])

  # add columns if not events if any
  tt_keep_mis <- setdiff(tt_keep, colnames(df_count_evt))
  for (tt_mis in tt_keep_mis){
    df_count_evt <- df_count_evt %>% mutate(!!tt_mis:=0)
  }
  df_count_evt <- df_count_evt[,c(col_evt, tt_keep)]

  # percent evt
  tt_keep_inter <- intersect(tt_keep, colnames(df_count_evt))
  df_percent_evt <- df_count_evt %>% mutate(across(c(all_of(tt_keep_inter)), function(x)  x/df_count_tt[cur_column(),]))

  # total count evt
  df_tot_evt <- df_count_evt %>%
    mutate(Count=rowSums(across(where(is.numeric))), Percent=Count/sum(df_count_tt$Count)) %>%
    select(c(.data[[col_evt]], Count, Percent))

  # fraction evt classes
  if (!is.null(col_evt_classes)){
    df_fraction_evt_classes <- df_evt_count %>%
      group_by(.data[[col_evt]], .data[[col_evt_classes]]) %>%
      summarize(Count=sum(Count), .groups="keep") %>%
      group_by(.data[[col_evt]]) %>%
      mutate(Percent=Count/sum(Count)) %>%
      ungroup(.data[[col_evt]])
  } else {
    df_fraction_evt_classes <- NULL
  }

  list(count_col=df_count_tt,
       count_row=df_tot_evt,
       fraction_row=df_fraction_evt_classes,
       count=df_count_evt,
       percent=df_percent_evt)
}


#' Complete a column of a dataframe.
#'
#' @param df The dataframe.
#' @param col The name of the column.
#' @param vals Vals that should exist in the column. If values are missing, they are added as extra rows.
#' @return a data.frame
#'
#' @keywords internal
complete_col_df <- function(df, col, vals){
  vals_missing <- setdiff(vals, df[[col]])
  df_missing <- data.frame(matrix(0, nrow=length(vals_missing), ncol=ncol(df)))
  colnames(df_missing) <- colnames(df)
  df_missing[,col] <- vals_missing

  rbind(df, df_missing)
}


#' Subset a list of dataframes on a column.
#'
#' @param dfs List of dataframes.
#' @param col The name of the column.
#' @param vals Vals that should exist in the column. If values are missing, they are added as extra rows.
#' @return a data.frame
#'
#' @export
select_vals_in_dfs <- function(dfs, col, vals){
  for (name in names(dfs)){
    if (col %in% names(dfs[[name]])){
      dfs[[name]] <- complete_col_df(dfs[[name]], col, vals)
      dfs[[name]] <- dfs[[name]] %>% filter(.data[[col]] %in% vals) %>% arrange(.data[[col]])
    }
  }
  
  dfs
}

#' Get the union of values that meet threshold requirements.
#'
#' @param dfs_plot List of dataframes.
#' @param min_counts_evt List of minimum counts for each cohort.
#' @param col_evt Name of the column from which to extract a entries.
#' @param max_evt An integer giving the maximum number of events to be considered.
#' @param max_evt_cohort A name identifying the table in `dfs_plot` to be used for selecting the events to be kept
#'   in case there are more than `max_evt` selected after filtering by `min_counts_evt`.
#' @return a vector.
#'
#' @export
select_in_plot_evt <- function(dfs_plot, min_counts_evt, col_evt, max_evt=NULL, max_evt_cohort=NULL){
  cohorts <- names(dfs_plot)
  in_plot_evt <- c()
  for (cohort in cohorts){
    min_count <- min_counts_evt[[cohort]]
    in_plot_evt_cohort <- dfs_plot[[cohort]]$count_row %>% filter(Count >= min_count) %>% pull(.data[[col_evt]])
    in_plot_evt <- union(in_plot_evt, in_plot_evt_cohort)
  }

  if (!is.null(max_evt)){
    if (length(in_plot_evt) > max_evt){
      stopifnot(!is.null(max_evt_cohort))
      df_cnt <- dfs_plot[[max_evt_cohort]]$count_row %>% arrange(desc(Count))
      df_cnt <- df_cnt[1:max_evt,]
      in_plot_evt <- df_cnt %>% pull(.data[[col_evt]])
    }
  }

  sort(in_plot_evt)
}


# ======================================================================================================================
#
# SIMPLE HEATMAP
#
# ======================================================================================================================

#' Draw a simple heatmap with numbers
#'
#' @param df_z a data.frame with numeric values.
#' @param z_name Name of the values for the hovering template.
#'   'middle_right_bar', 'extreme_right_bar'. The corresponding values must be in the list of names of \code{dfs}.
#' @inheritParams get_colorscale_htmp
#' @inheritParams add_z_annotations_htmp
#' @param width_one (optional) The width of the plot is obtained as \code{width_one} multiplied by the number of columns
#'   in \code{z}.
#' @param height_one (optional) The height of the plot is obtained as \code{height_one} multiplied by the number of rows 
#'   in \code{z}.
#' @return a plotly figure object.
#'
#' @import plotly
#'
#' @author Yoann Pradat
#' @export
draw_numbered_heatmap <- function(df_z, z_name, colors_limits=NULL, colors_palette="Reds", width_one=50, height_one=15,
                                  black_white_cutoff=0.5, font=list(size=8)){
  y <- rownames(df_z)
  x <- colnames(df_z)

  # Heatmap
  z <- as.matrix(df_z)

  ## colors heatmap
  colz <- get_colorscale_htmp(z, colors_limits=colors_limits, colors_palette=colors_palette)

  ## plot heatmap
  fig_htmp <- plot_ly(z=z, x=x, y=y, type="heatmap",
                      colorscale=colz,
                      showscale=F,
                      hovertemplate=paste0(paste0('<b>%{x}', ' - ', "%{y}</b><br>"),
                                           paste0(z_name, ': %{z}<br>'),
                                           "<extra></extra>"),
                      height=height_one*nrow(z), width=width_one*ncol(z))

  ## annotation text heatmap
  fig_htmp <- add_z_annotations_htmp(z, fig_htmp, black_white_cutoff=black_white_cutoff, font=font)

  ## axes configuration
  fig_htmp <- configure_axes_htmp(z, fig_htmp)

  fig_htmp
}
