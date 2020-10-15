# utl
#
# This module procided some utiilty helpers for quant step

library(tidyverse)

utl <- new.env(); source(here::here("scripts/utils/base.R"), utl)


PATTERNS <- list(
  transcript = "eval_align.merged.sqlite$"
)

load_result <- function(path, table_name) {
  conn <- RSQLite::SQLite() %>% RSQLite::dbConnect(path, synchronous = "off")

  .data <- conn %>% RSQLite::dbGetQuery(paste0("select * from ", table_name, ";"))
  RSQLite::dbDisconnect(conn)

  .data
}


calc_metrics <- function(
  confusion_matrix,
  feature_ids = c()
) {

  .confusion_matrix <- confusion_matrix

  if (length(feature_ids) > 0) {
    .confusion_matrix <- .confusion_matrix %>% filter(qname %in% feature_ids)
  }

  .metrics <- .confusion_matrix %>%
    mutate(pct_tp = tp / total) %>%
    mutate(pct_fp = fp / total) %>%
    mutate(pct_tn = 0.0) %>%
    mutate(pct_fn = fn / total)

  .metrics <- .metrics %>%
    group_by(source_id) %>%
    summarize(tp = sum(pct_tp), fp = sum(pct_fp), tn = sum(pct_tn), fn = sum(pct_fn))

  .metrics <- .metrics %>%
    mutate(recall = tp / (tp + fn)) %>%
    mutate(precision = tp / (tp + fp)) %>%
    mutate(f1 = 2 * tp / (2 * tp + fp + fn)) %>%
    select(-tp, -fp, -tn, -fn)

  .metrics
}
