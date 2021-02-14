# eval_rnaseq_map.utils
# ~~~~~~~~~~~~~~~~~~~~~
#
# This module procided some utiilty functions for plotting

library(tidyverse)

source(here::here("scripts/utils/base.R"))


to_lab <- function(x) {
  LABS <- list(
    recall = "Recall",
    precision = "Precision",
    f1 = "F1 score",
    spearman = "Spearman's rho",
    nrmse = "NRMSE",
    auc = "AUC",
    spearman_ = "Spearman's rho (intersection)",
    nrmse_ = "NRMSE (intersection)",
    auc_ = "AUC (intersection)",
    all = "All biotypes",
    mrna = "mRNA",
    protein_coding = "mRNA",
    .protein_coding = "mRNA", # XXX
    lncrna = "lncRNA",
    .lncrna = "lncRNA", # XX
    other = "Others",
    others = "Others",
    length = "Read length (bases)",
    depth = "Library size (million reads)",
    reps = "Number of replicates",
    gencode = "GENCODE",
    gencode_basic = "GENCODE-Basic",
    gencode_refseq = "GENCODE_RefSeq",
    refseq = "RefSeq",
    noncode = "NONCODE"
  )

  tryCatch({
    return(
      x %>%
        map_chr(function(x)
          LABS[[x]])
    )
  },
  error = function(e) {
    message(e)
    return(x)
  })
}


scale_ <- function(group) {
  if (group == "abbr") {
    return(
      list(color = "tool",
                   shape = "upper_abbr",
                   linetype = "upper_abbr")
    )
  }

  list(color = group,
                   shape = group,
                   linetype = group)
}


to_color <- function(x) {
  COLORS <- c(
    "star" = "#00AFBB",
    "rsem" = "#00AFBB",
    "ebseq" = "#00AFBB",
    "hisat2" = "#E7B800",
    "stringtie" = "#E7B800",
    "ballgown" = "#E7B800",
    "tophat2" = "#FC4E07",
    "cuffdiff" = "#FC4E07",
    "kallisto" = "#2B046E",
    "sleuth" = "#2B046E",
    "sleuth_wt" = "#2B046E",
    "edger" = "#888888",
    "incomplete" = "#00AFBB",
    "complete" = "#E7B800",
    "complicated" = "#FC4E07",
    "mrna" = "#00AFBB",
    "lncrna" = "#E7B800",
    "other" = "#FC4E07",
    "others" = "#FC4E07",
    "all" = "#888888"
  )

  sapply(as.character(x), function(x)
    COLORS[[x]], USE.NAMES = FALSE)
}


add_condition_lab <-
  function(plot,
           lab_add = NULL,
           var_remove = NULL) {

    .condition <- list(length = "read length = 100 bases",
                       depth = "library size = 40 M reads",
                       reps = "n = 3")
    if (length(var_remove) > 0)
      .condition[[var_remove]] <- NULL
    if (length(lab_add) > 0)
      .condition[[length(.condition) + 1]] <- lab_add

    .condition_lab <-
      paste0("Condition: ", paste(unlist(.condition), collapse = ", "))

    plot %>% ggpubr::annotate_figure(bottom = ggpubr::text_grob(.condition_lab, size = 16))
  }


plot_line <-
  function(data_,
           group_,
           tilt = FALSE,
           facet.by = NULL,
           theme_ = NULL,
           ...) {

    .default <- list(
      x = "condition",
      y = "value",
      group = group_,
      color = group_,
      fill = group_
    )

    .kwargs <- c(...)

    .scale_ <- scale_(group_)

    .kwargs <- c(.kwargs, list(group = group_, color = group_, shape = group_))

    .default[names(.kwargs)] <- .kwargs

    .colors <- data_ %>%
      distinct(get(!!group_), get(!!.scale_$color)) %>%
      deframe() %>%
      tolower()

    .colors <- .colors %>%
      to_color() %>%
      set_names(names(.colors))

    .shapes <- data_ %>%
      distinct(get(!!group_), get(!!.scale_$shape)) %>%
      deframe() %>%
      factor()

    .theme <- theme(
      text = element_text(size = 18),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      strip.text.x = element_text(size = 16)
    )

    if (length(theme_) > 0) .theme <- theme_

    if (tilt) .theme <- .theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    g <- do.call(ggpubr::ggline, c(list(data_), .default))
    g <- g + scale_colour_manual(name = group_, values = .colors)
    g <- g + scale_shape_manual(name = group_, values = .shapes)
    g <- g + theme_bw() + .theme

    if (length(facet.by)) g <- ggpubr::facet(g + theme_bw() + .theme, facet.by = facet.by)

    g
  }


plot_box <-
  function(data_,
           group_,
           tilt = FALSE,
           facet.by = NULL,
           theme_ = NULL,
           ...) {

    .default <- list(
      x = "abbr",
      y = "value",
      group = group_,
      color = group_
    )

    .kwargs <- c(...)

    .scale_ <- scale_(group_)

    .kwargs <- c(.kwargs, list(group = group_, color = group_))

    .default[names(.kwargs)] <- .kwargs

    .colors <- data_ %>%
      distinct(get(!!group_), get(!!.scale_$color)) %>%
      deframe() %>%
      tolower()

    .colors <- .colors %>%
      to_color() %>%
      set_names(names(.colors))

    .shapes <- data_ %>%
      distinct(get(!!group_), get(!!.scale_$shape)) %>%
      deframe() %>%
      factor()

    .theme <- theme(
      text = element_text(size = 18),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      strip.text.x = element_text(size = 16)
    )

    if (length(theme_) > 0) .theme <- theme_

    if (tilt) .theme <- .theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

    g <- do.call(ggpubr::ggboxplot, c(list(data_), .default))
    g <- g + scale_colour_manual(name = group_, values = .colors)
    g <- g + theme_bw() + .theme

    if (length(facet.by)) g <- ggpubr::facet(g + theme_bw() + .theme, facet.by = facet.by)

    g
  }


plot_bar_ <- function(data_,
                      x = "bin",
                      y = "value",
                      facet.by,
                      group = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      ylim = NULL,
                      breaks = NULL,
                      labels = NULL,
                      tilt = FALSE) {
  theme_ <- theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

  if (tilt)
    theme_ <-
      theme_ + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  g <- data_ %>%
    ggbarplot(
      x,
      y,
      group = group,
      color = "darkgrey",
      fill = "darkgrey",
      label = FALSE,
      position = position_dodge(0.9),
      lab.size = 7,
      xlab = xlab,
      ylab = ylab
    )
  g <- g + scale_x_discrete(breaks = breaks, labels = labels)
  g <- g + theme_ + theme(legend.position = "none")
  g <- facet(g + theme_bw() + theme_, facet.by = facet.by)

  g
}


plot_bar <- function(data_,
                     group_ = NULL,
                     tilt = FALSE,
                     facet.by = NULL,
                     theme_ = NULL,
                     ...) {

  .default <- list(
    x = "condition",
    y = "value",
    group = group_,
    color = group_,
    fill = group_,
    position = position_dodge(0.9)
  )

  .kwargs <- c(...)

  .scale_ <- scale_(group_)

  .default[names(.kwargs)] <- .kwargs

  .colors <- data_ %>%
    distinct(get(!!group_), get(!!.scale_$color)) %>%
    deframe() %>%
    tolower()

  .colors <- .colors %>%
    to_color() %>%
    set_names(names(.colors))

  .theme <- theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

  if (length(theme_) > 0) .theme <- theme_

  if (tilt)
    .theme <-
    .theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  g <- do.call(ggpubr::ggbarplot, c(list(data_), .default))
  g <- g + scale_colour_manual(name = group_, values = .colors)
  g <- g + scale_fill_manual(name = group_, values = .colors)
  g <- g + .theme

  if (length(facet.by)) g <- ggpubr::facet(g + theme_bw() + .theme, facet.by = facet.by)

  g
}


plot_scatter_ <- function(data_,
                     group_ = NULL,
                     lims = NULL,
                     tilt = FALSE,
                     facet.by = NULL,
                     theme_ = NULL,
                     ...) {

  .default <- list(
    x = "x",
    y = "y",
    group = group_,
    color = group_,
    fill = group_,
    size = 3
  )

  .kwargs <- c(...)

  .scale_ <- scale_(group_)

  .default[names(.kwargs)] <- .kwargs

  .colors <- data_ %>%
    distinct(get(!!group_), get(!!.scale_$color)) %>%
    deframe() %>%
    tolower()

  .colors <- .colors %>%
    to_color() %>%
    set_names(names(.colors))

  .theme <- theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

  if (length(theme_) > 0) .theme <- theme_

  if (tilt)
    .theme <-
    .theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))

  g <- do.call(ggpubr::ggscatter, c(list(data_), .default))
  g <- g + scale_colour_manual(name = group_, values = .colors)
  g <- g + scale_fill_manual(name = group_, values = .colors)
  g <- g + theme_bw() + .theme

  if (length(facet.by)) g <- ggpubr::facet(g + theme_bw() + .theme, facet.by = facet.by)

  g
}


plot_scatter <-
  function(data_,
           labs,
           title,
           lims,
           group,
           plot.title.vjust = -6) {
    scale_ <- scale_(group)

    colors_ <- data %>%
      distinct(get(!!group), get(!!scale_$color)) %>%
      deframe() %>%
      tolower()

    colors_ <- colors_ %>%
      to_color() %>%
      set_names(names(colors_))

    shapes_ <- data %>%
      distinct(get(group), get(scale_$shape)) %>%
      deframe() %>%
      factor()

    theme_ <- theme(
      text = element_text(size = 18),
      plot.title = element_text(vjust = plot.title.vjust, hjust = 0.025),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 18)
    )

    lim_ <- lims[[tolower(title)]]

    g <- data_ %>% ggpubr::ggscatter(
      "x",
      "y",
      color = group,
      shape = group,
      size = 3,
      xlab = labs[1],
      ylab = labs[2],
      title = to_lab(title)
    )
    g <- g + scale_colour_manual(name = "group", values = colors_)
    g <- g + scale_shape_manual(name = "group", values = shapes_)
    g <-
      g + geom_abline(slope = 1,
                      intercept = 0,
                      linetype = "dotted")
    g <- g + xlim(lim_) + ylim(lim_) + coord_fixed()
    g <- g + theme_bw() + theme_
    g
  }


# NOTE: Plot ROC without color
plot_roc_ <- function(data_, lfc = 1, fdr = 0.05) {
  g <- data_ %>%
    ggplot2::ggplot(aes(d = D, m = M))
  g <- g + plotROC::geom_roc(n.cuts = FALSE,, size = 1)

  g
}


plot_roc <- function(data,
                     group,
                     facet.by,
                     lfc = 1,
                     fdr = 0.05) {
  scale_ <- scale_(group)

  colors_ <- data_ %>%
    distinct(get(!!group), get(!!scale_$color)) %>%
    deframe() %>%
    tolower()

  colors_ <- colors_ %>%
    to_color() %>%
    set_names(names(colors_))

  shapes_ <- data_ %>%
    distinct(get(group), get(scale_$shape)) %>%
    deframe() %>%
    factor()

  theme_ <- theme(
    text = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18)
  )

  g <- data_ %>%
    ggplot(aes(
      d = D,
      m = M,
      color = get(group),
      linetype = get(group)
    ))
  g <- g + plotROC::geom_roc(n.cuts = FALSE,, size = 1)
  g <-
    g + geom_abline(slope = 1,
                    intercept = 0,
                    linetype = "dotted")
  g <- g + scale_colour_manual(name = "group", values = colors_)
  g <- g + scale_shape_manual(name = "group", values = shapes_)
  g <- g + style_roc(theme = theme_bw) + theme_
  g <- g + xlim(0, 1) + ylim(0, 1) + coord_fixed()
  g <- g + xlab("False positive rate") + ylab("True positive rate")
  g <- facet(g + theme_bw() + theme_, facet.by = facet.by)

  g
}


plot_tile <-
  function(data_,
           x,
           y,
           group,
           title = NULL,
           xlab = NULL,
           ylab = NULL) {
    g <- data_ %>%
      ggplot(aes_string(x, y, group)) +
      geom_tile(aes(fill = count)) +
      geom_text(aes(fill = data_$count, label = round(data_$count, 1)))

    g <- g + labs(title = title, x = xlab, y = ylab)
    g <- g + scale_fill_gradient(low = "#ffffff", high = "#24417d")

    g
  }


plot_intersect <- function(data_) {
  theme_ <- theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 18),
    plot.margin = unit(c(0, 5, 0, 5), "lines"),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

  limits_ <- data_ %>% min_max_()

  g <- data_ %>%
    ggcorrplot::ggcorrplot(
      method = "circle",
      hc.order = TRUE,
      outline.color = "#ffffff",
      lab = FALSE,
      tl.cex = 16,
      show.legend = TRUE
    )
  g <-
    g + scale_fill_gradient(low = "#ffffff",
                            high = "#24417d",
                            limits = limits_)
  g <- g + theme_

  g
}


plot_intersect_ <- function(data_) {
  theme_ <- theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 18),
    plot.margin = unit(c(0, 5, 0, 5), "lines"),
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

  limits_ <- data_ %>% min_max_()

  g <- data_ %>%
    ggcorrplot::ggcorrplot(
      method = "square",
      hc.order = TRUE,
      outline.color = "#ffffff",
      lab = TRUE,
      lab_size = 12,
      tl.cex = 36,
      show.legend = FALSE
    )
  g <-
    g + scale_fill_gradient(low = "#ffffff",
                            high = "#24417d",
                            limits = limits_)
  g <- g + theme_

  g
}


save_plot <-
  function(plot,
           sub_name,
           name,
           width = 6.4 * 2.4,
           height = 3.2 * 2.8) {
    path_ <- here::here(output_dir, paste(PREFIX, sub_name, name, sep = "_"))
    message(path_)
    plot(plot)
    # NOTE: Specify cario_pdf to display unicode character
    ggsave(
      file = path_,
      plot = plot,
      dpi = 300,
      width = width,
      height = height,
      device = cairo_pdf
    )

    return(NULL)
  }


add_lab <- function(g, text) {
  if(is.null(g)) return(NULL)
  g %>% ggpubr::annotate_figure(fig.lab = text, fig.lab.size = 21, fig.lab.face = "bold")
}


rm_legend <- function(g) {
  if (is.null(g)) return(NULL)
  (g + guides(feature_type = FALSE)) + theme(legend.position = "none")
}
