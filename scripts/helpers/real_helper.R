# Utils
#
# This module procided some utiilty helpers for DE step

h.de <- new.env(); source(here::here("scripts/helpers/de_helper.R"), h.de)

library(tidyverse)


load_qpcr_result <- function(input_path) {
  # TODO: Eliminate NA value
  result <- input_path %>% data.table::fread(sep = "\t", stringsAsFactors = FALSE)

  columns <- c("id_ref", "value", "flag_detection")
  nsamples <- ncol(result) / length(columns)
  groups <- rep(c("A", "B", "C", "D"), each = 4)
  sample_ids <- paste0(groups, seq(1, 4))

  colnames_ <- c()
  for (s in sample_ids) {
    colnames_ <- c(colnames_, paste(s, columns, sep = "_"))
  }

  colnames(result) <- colnames_
  colsel <- grep("value", colnames(result))

  result <- result %>% select(colsel) %>% rowid_to_column(var = "ID")
  colnames(result) <- colnames(result) %>% strsplit("_") %>% sapply("[", 1)
  result <- result %>%
    pivot_longer(-1, names_to = "sample", values_to = "value") %>%
    mutate(group = substr(sample, 1, 1))

  return(result)
}


load_results <- function(paths, ballgown_coef = 1) {
  results_de <- list()
  results_de_sig <- list()

  for (i in seq_along(paths)) {
    names(paths)[i] %>% print
    tool_ <- names(paths)[i] %>% strsplit("-") %>% sapply(function(x) x[length(x)])

    load_func_ <- get(paste0("load_", tool_), envir = h.de)
    if (tool_ == "ballgown") load_func_ <- partial(load_func_, coef = ballgown_coef)

    results_de[[i]] <- load_func_(paths[i])
    results_de_sig[[i]] <- results_de[[i]] %>% h.de$extract_de(t_logfc = 1, t_pval = 0.05)
  }
  names(results_de) <- names(paths)
  names(results_de_sig) <- names(paths)

  return(
    list(
      de = results_de,
      sig = results_de_sig
    )
  )
}
