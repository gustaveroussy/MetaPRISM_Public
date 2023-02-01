# @created: 22 Nov 21
# @modified: 27 Dec 22
# @authors: Yoann Pradat
#
# Draw Kaplan-Meier curves where tumors are stratified according to the presence or absence of some genomic alterations.

suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(rprism))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(survival))

source("workflow/functions/km_curves.R")
source("workflow/functions/survival_models.R")

# functions ============================================================================================================

pie_count_plotly <- function(fig, df_count, pal=NULL, values2colors=NULL, textinfo='label', fontsize=8){
  font <- list(size=fontsize, family="Helvetica", color="black")
  font_anns <- list(size=fontsize*1.5, family="Helvetica", color="black")
  fig <- plot_ly()

  if(!is.null(values2colors)){
    df_count <- df_count %>%
      mutate(Color=recode(Label, !!!values2colors))
  } else {
    df_count <- df_count %>%
      mutate(Color=get_label_colors(Label, pal=pal))
  }

  # reorder to match order of values2colors
  df_count$Label <- factor(df_count$Label, levels=names(values2colors))
  df_count <- df_count %>% arrange(Label)

  fig <- fig %>% add_pie(labels=df_count$Label, values=df_count$Count,
                         textinfo=textinfo,
                         text=df_count$Percentage,
                         insidetextfont=font,
                         sort=F,
                         hovertemplate = paste0(paste0("<b>%{label}</b>"),
                                                '<br><br>Count: %{value}<br>',
                                                'Percent: %{text:.3f}',
                                                "<extra></extra>"),
                         hole=0.4,
                         marker = list(line=list(color='#FFFFFF', width=1), colors=df_count$Color))

  annotation <- list(text=sum(df_count$Count), showarrow=F, font=font_anns, xanchor="center", xref="x", yref="y")
  emptyaxis <- list(title = "", zeroline = FALSE, showline = FALSE, showticklabels = FALSE, showgrid = FALSE)

  fig <- fig %>% layout(font=font, showlegend=F, annotations=list(annotation), yaxis=emptyaxis, xaxis=emptyaxis)
  fig
}


encode_discrete_exclusive <- function(df_dat, cov_name, name_sep="_", remove_first=T, remove_col=T, ref_level=NULL){
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

    df_dat <- bind_cols(df_dat, df_dat_dum)
  }

  if (remove_col){
    df_dat[[cov_name]] <- NULL
  }

  df_dat
}


main <- function(args){
  # load tables
  df_cln <- load_cln(args$cohort) %>% rename(Tumor_Type=Project_TCGA_More)
  df_evt <- load_table(args$event_counts)

  # select tumor types
  if (tolower(args$tumor_types)=="dna_cohort"){
    tt_keep <- df_cln %>% filter(grepl("DNA_T", Sample_Type)) %>% 
      filter(!grepl("Not_TCGA", Tumor_Type)) %>% group_by(Tumor_Type) %>% filter(n() >= 10) %>% 
      select(Tumor_Type) %>% distinct() %>% pull(var="Tumor_Type")
    df_cln <- df_cln %>% filter(Tumor_Type %in% tt_keep)
  } else if (tolower(args$tumor_types)!="all"){
    df_cln <- df_cln %>% filter(Tumor_Type %in% toupper(args$tumor_types))
  }

  # select sample types
  mask_sts <- lapply(args$sample_types, function(st) grepl(st, df_cln$Sample_Type))
  mask_sts <- Reduce(function(m1, m2) m1 & m2, mask_sts)
  df_cln <- df_cln %>% filter(mask_sts)

  # transform survival columns
  df_cln <- df_cln %>% filter(!is.na(Survival_Time)) %>%
      mutate(Survival_Time=Survival_Time/30, Survival_Status=ifelse(Survival_Status=="Deceased", 1, 0))

  # if sample type is DNA_T only, remove Fus from alterations names
  if (setequal(args$sample_types, "DNA_T")){
    alterations <- setdiff(colnames(df_evt), "Subject_Id")
    alterations_fus_only <- alterations[grepl(" Fus$", alterations)]
    alterations_fus <- setdiff(alterations[grepl("Fus", alterations)], alterations_fus_only)
    for (alteration_fus in alterations_fus){
      alteration_fus_new <- gsub("Fus\\|", "", alteration_fus)
      alteration_fus_new <- gsub("\\|Fus", "", alteration_fus_new)
      df_evt <- df_evt %>% rename(!!alteration_fus_new:=.data[[alteration_fus]])
    }
    df_evt <- df_evt %>% select(-all_of(alterations_fus_only))
  }

  # granularity
  alterations <- setdiff(colnames(df_evt), "Subject_Id")
  alterations2genes <- sapply(alterations, function(x) unlist(str_split(x, " "))[1])
  genes <- unique(as.vector(alterations2genes))

  for (gene in genes){
    alterations_gene <- names(alterations2genes)[alterations2genes==gene]
    if (args$granularity=="gene"){
      df_evt <- df_evt %>% mutate(!!gene:=rowSums(across(all_of(alterations_gene)))) %>%
        select(-all_of(setdiff(alterations_gene, gene)))
      mask_wt <- df_evt[[gene]] == 0 
      df_evt[[gene]] <- as.character(df_evt[[gene]])
      df_evt[mask_wt, gene] <- "No oncogenic alt."
      df_evt[!mask_wt, gene] <- "Oncogenic alt."
    } else if (args$granularity=="alteration") {
      if (setequal(alterations_gene, gene)){
        mask_wt <- df_evt[[gene]] == 0 
        df_evt[[gene]] <- as.character(df_evt[[gene]])
        df_evt[mask_wt, gene] <- "No oncogenic alt."
        df_evt[!mask_wt, gene] <- "Oncogenic alt."
      } else {
        df_evt[[gene]] <- "No oncogenic alt."
        for (alteration in alterations_gene){
          mask_alt <- df_evt[[alteration]]==1
          df_evt[mask_alt, gene] <- alteration
          if (gene!=alteration){
            df_evt[alteration] <- NULL
          }
        }
      }
    } else {
      stop(paste("Unsupported value", args$granularity, "for --granularity option."))
    }
  }

  # merge
  df_cln <- left_join(df_cln, df_evt, by="Subject_Id")

  # remove NA
  df_cln <- df_cln %>% filter(!is.na(.data[[genes[1]]]))

  # select events with sufficient counts
  all_events <- setdiff(colnames(df_evt), "Subject_Id")
  df_cln_bin <- as_tibble(df_cln[all_events]!="No oncogenic alt.")
  df_cln_bin <- df_cln_bin %>% mutate(across(where(is.logical), as.numeric))
  sel_events <- all_events[colSums(df_cln_bin, na.rm=T)>=args$min_counts]
  df_cln <- df_cln %>% mutate_at(all_of(all_events), ~as.factor(.))
  tumortypes2colors <- load_colors("Project_TCGA_More")
  surv <- Surv(df_cln$Survival_Time, df_cln$Survival_Status)

  # draw km-curve and pie plot for each selected event
  for (sel_event in sel_events){
    # km-curve
    legend.labs <- levels(df_cln[[sel_event]])
    
    # change order of levels for colors
    if (setequal(legend.labs, c("No oncogenic alt.", "Oncogenic alt."))){
      legend.labs <- c("No oncogenic alt.", "Oncogenic alt.") 
      df_cln[[sel_event]] <- factor(df_cln[[sel_event]], levels=legend.labs)
    }

    if (tolower(args$tumor_types) %in% c("all", "dna_cohort")){
      # compute log hazard ratio, adjusted for tumor type
      covs <- c(sel_event, "Tumor_Type")
      df_dat <- df_cln[,covs]
      df_dat <- encode_discrete_exclusive(df_dat=df_dat, cov_name="Tumor_Type", name_sep="_", remove_first=T,
                                          remove_col=T, ref_level="LUAD")
      out <- train_survival(surv=surv, df=df_dat, name_model="coxph_standard", train_indices=NULL, train_foldid=NULL)
      df_out <- out$df_cov
      df_out <- df_out %>% mutate(Z=Coefficient/Standard_Error, pval=1-pchisq(Z^2, df=1))
      hr <- exp(df_out$Coefficient[1])
      se <- df_out$Standard_Error[1]
      z_alpha <- 1-qnorm(0.025)
      pval <- df_out$pval[1]
      hr_low <-  exp(df_out$Coefficient[1] - z_alpha*se)
      hr_high <- exp(df_out$Coefficient[1] + z_alpha*se)
      hr_annot <- paste("Adj. HR",  round(hr,2), paste0("[", round(hr_low, 2), ";", round(hr_high, 2), "]"))
    } else  {
      hr_annot <- NULL
    }

    if (tolower(args$tumor_types) %in% c("all", "dna_cohort")){
      title <- sel_event
    } else {
      title <- paste(sel_event, "-", paste0(args$tumor_types, collapse=" & "))
    }
    # tables.height <- 0.135 * length(legend.labs)
    tables.height <- 0.125 * length(legend.labs)
    fit <- eval(parse(text=(paste("survfit(Surv(Survival_Time, Survival_Status) ~", sel_event, ", data=df_cln)"))))
    output <- file.path(args$output, paste0(sel_event, ".pdf"))
    draw_km_curves(df_cln, fit, output=output, width=args$width, height=args$height, title=title,
                   xlim=88, break.time.by=12, tables.height=tables.height, conf.int=F, legend.labs=legend.labs,
                   pval.coord=c(12,0.95), ylab=NULL, annot=hr_annot, annot_x=12, annot_y=0.7)

    # pie plot
    mask_ok_1 <- tolower(args$tumor_types) %in% c("all", "dna_cohort")
    mask_ok_2 <- args$tumor_types=="LUAD" & sel_event=="STK11"
    mask_ok_3 <- args$tumor_types=="PRAD" & sel_event=="ERG"

    if (mask_ok_1 | mask_ok_2 | mask_ok_3){
      df_count <- df_cln %>% filter(.data[[sel_event]]!="No oncogenic alt.") %>% group_by(Tumor_Type) %>%
        dplyr::summarize(Count=n()) %>% rename(Label=Tumor_Type) %>% mutate(Percentage=Count/sum(Count))
      fig <- pie_count_plotly(fig, df_count, pal=NULL, values2colors=tumortypes2colors, textinfo='label', fontsize=8)
      output <- file.path(args$output, paste0(sel_event, "_pie.pdf"))
      save_image(fig, output, width=250, height=150)
    }
  }
}

# run ==================================================================================================================

if (getOption('run.main', default=TRUE)) {
  parser <- ArgumentParser(description='Draw KM curves to assess prognostic value of genetic alterations.')
  parser$add_argument("--cohort", type="character", help="Name of the cohort.", default="prism")
  parser$add_argument("--tumor_types", type="character", nargs="+", default="LUAD", 
                      help="List of tumor types to be considered.")
  parser$add_argument("--event_counts", type="character", help="Path to table of events counts.",
                      default="../../../results/combined_alterations/count/count_by_gene_DNA_RNA_annotated_prism.tsv")
  parser$add_argument("--sample_types", type="character", nargs="+", default=c("DNA_T", "RNA_T"),
                      help="Only patients with these samples will be considered.")
  parser$add_argument("--min_counts", type="integer", default=20,
        help="Threshold on the minimum number of events for the variable to be used for making 2 distinct populations.")
  parser$add_argument("--granularity", type="character", default="gene",
                      help="Choose 'gene' or 'alteration'.")
  parser$add_argument("--width", type="double", default=5, help="Width of output plot in inches.")
  parser$add_argument("--height", type="double", default=3, help="Height of output plot in inches.")
  parser$add_argument("--output", type="character", help="Path to output folder.",
                      default="../../../results/survival_analysis/km_curves/LUAD/by_gene")
  parser$add_argument("-l", "--log", type="character", help="Path where log file will be saved.")
  args <- parser$parse_args()

  # log file
  log <- file(args$log, open="wt")
  sink(log)
  sink(log, type="message")

  print(args)
  main(args)
}
