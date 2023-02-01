suppressPackageStartupMessages(library(survminer))

customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL){
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}


draw_km_curves <- function(df_sub, fit, output, width, height, title="Survival curves", xlim=88, break.time.by=12,
                      palette="lancet", tables.height=0.2, conf.int=F, legend.labs=c(seq(0,7), "8+"), annot=NULL,
                      annot_x=NULL, annot_y=NULL, ylab = "Survival probability", size_labels=12, ...){
  if (grepl(".pdf", output)){
    grDevices::pdf(file=output,
                   width=width,
                   height=height,
                   onefile=T)
  } else if (grepl(".svg", output)){
    grDevices::svg(file=output,
                   width=width,
                   height=height,
                   onefile=T)
  }

  theme_table <- theme(plot.title=element_blank(),
                       axis.line.x=element_blank(),
                       axis.line.y=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_text(size=size_labels),
                       axis.ticks.x=element_blank(),
                       axis.ticks.y=element_blank(),
                       axis.title.x=element_blank(),
                       axis.title.y=element_blank())

  ggsurv <- ggsurvplot(fit,                                   # survfit object with calculated statistics.
                       data = df_sub,                         # data used to fit survival curves.
                       risk.table = TRUE,                     # show risk table.
                       submain = title,
                       pval = TRUE,                           # show p-value of log-rank test.
                       conf.int = conf.int,                   # show confidence intervals for point estimates 
                       palette=palette,
                       xlim = c(0,xlim),                      # present narrower X axis, but not affect
                       # survival estimates.
                       xlab = "Time in months",               # customize X axis label.
                       ylab = ylab,                           # customize Y axis label.
                       break.time.by = break.time.by,         # break X axis in time intervals by 500.
                       tables.height = tables.height,         # the height of the risk table
                       legend.labs = legend.labs,             # change legend labels.
                       legend = "none",
                       risk.table.title = NULL,
                       risk.table.pos = "out",
                       tables.theme = theme_table,            # customize risk table with a theme.
                       ...
                       # ggtheme = theme_survminer(),         # customize plot and risk table with a theme.
                       # tables.theme = theme_cleantable(),   # customize plot and risk table with a theme.
                       # risk.table.y.text = FALSE,
                       # risk.table.x.text = FALSE,
                       # ncensor.plot = F,                      # plot the number of censored subjects at time t
                       # conf.int.style = "step",               # customize style of confidence intervals
  )

  # ggsurv$table <- ggsurv$table + theme(axis.text.y=element_blank())
  ggsurv$plot <- ggsurv$plot + theme(axis.title.y=element_text(size=size_labels),
                                     axis.title.x=element_text(size=size_labels),
                                     axis.text.y=element_text(size=size_labels),
                                     axis.text.x=element_text(size=size_labels)) +
                 scale_y_continuous(breaks=c(0,0.5,1))

  if (!is.null(annot)){
    ggsurv$plot <- ggsurv$plot + annotate(geom="text", x=annot_x, y=annot_y, label=annot, hjust=0)
  }

  ggsurv <- customize_labels(
    ggsurv,
    # font.title     = c(7, "plain", "black"),
    # font.subtitle  = c(7, "plain", "black"),
    # font.x         = c(10, "plain", "black"),
    # font.y         = c(10, "plain", "black"),
    # font.xtickslab = c(9),
    # font.ytickslab = c(9)
  )

  print(ggsurv, newpage=F)
  dev.off()
  
  cat(paste("-file saved at", output, "\n"))
}

