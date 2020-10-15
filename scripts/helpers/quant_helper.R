# utl
#
# This module procided some utiilty helpers for quant step

library(tidyverse)

utl <- new.env(); source(here::here("scripts/utils/base.R"), utl)


PATTERNS <- list(
  transcript = c(
    cuffdiff = "conv_cuffdiff_to_raw/count_matrix_transcript.tsv",
    any = "conv_any_to_raw/count_matrix_transcript.tsv"
  ),
  gene = c(
    cuffdiff = "conv_cuffdiff_to_raw/count_matrix_gene.tsv",
    any = "conv_any_to_raw/count_matrix_gene.tsv"
  )
)


load_groundtruth <- function(path) {
  message(path)

  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  .data[is.na(.data$lambda), ]$lambda <- 0
  .data[is.na(.data$fc), ]$fc <- 1
  .data <- .data %>%
    mutate(ctrl = lambda, case = lambda * fc) %>%
    select(feature_id, ctrl, case)

  .data
}


load_result <- function(path) {
  message(path)

  .data <- data.table::fread(path, sep = "\t", stringsAsFactors = FALSE)
  colnames(.data)[1] <- "feature_id"

  .data
}


to_cpm <- function(x, log = FALSE, pseudo_count = utl$PSEUDO_COUNT, f = edgeR::cpm) {
  if (log) f <- partial(f, log = TRUE, prior.count = pseudo_count)

  .cpm <- cbind(
    select(x, !where(is.numeric)),
    x %>%
    select(where(is.numeric)) %>% f
  )

  .cpm
}


to_tpm <- function(x, lengths, log = FALSE, pseudo_count = utl$PSEUDO_COUNT, f = utl$tpm) {
  if (log) f <- partial(f, log = TRUE, prior.count = pseudo_count)

  .tpm <- cbind(
    select(x, !where(is.numeric)),
    x %>%
    column_to_rownames(var = colnames(.)[1]) %>%
    select(where(is.numeric)) %>% f(lengths = lengths)
  )
  rownames(.tpm) <- NULL

  .tpm
}


fill_lack <- function(x, annotation) {
  .feature_ids <- union(
    pull(x, 1),
    annotation$feature_id
  )

  .x <- data.frame(
    feature_id = .feature_ids,
    stringsAsFactors = FALSE
  )

  .x <- .x %>%
    left_join(x, by = c("feature_id" = colnames(x)[1])) %>%
    mutate(across(-1, ~ replace_na(.x, 0.0)))

  .x
}


.join <- function(est, true, f = dplyr::left_join) {
  .true <- true %>%
    pivot_longer(-feature_id, names_to = "group", values_to = "value")

  .est <- est %>%
    pivot_longer(-feature_id, names_to = "group", values_to = "value")

  .joined <- .true %>%
    f(.est, by = c("feature_id", "group"))

  .joined
}


calc_metric <- function(
  est,
  true,
  method,
  feature_ids = c(),
  join_ = dplyr::left_join
) {
  # XXX: Exclude 0 value on estimated
  .joined <- compare(est, true, feature_ids, join_) %>%
    filter(value.x > 0) %>%
    filter(value.y > 0)

  tryCatch(
    {
      .metric <- method(.joined$value.x, .joined$value.y)
    },
    error = function(e) {
      message(e)
      "Return NaN." %>% message
      .metric <<- NaN
    }
  )

  .metric
}


# NOTE: Calculate metric for each group and feature
calc_metric2 <- function(
  est,
  true,
  method,
  feature_ids1 = c(),
  feature_ids2 = c(),
  join_ = dplyr::left_join
) {
  # XXX: Exclude 0 value on estimated
  .joined <- compare2(est, true, feature_ids1, feature_ids2, join_) %>%
    filter(value.x > 0) %>%
    filter(value.y > 0)

  tryCatch(
    {
      .metric <- method(.joined$value.x, .joined$value.y)
    },
    error = function(e) {
      message(e)
      "Return NaN." %>% message
      .metric <<- NaN
    }
  )

  .metric
}


compare <- function(
  est,
  true,
  feature_ids = c(),
  join_ = dplyr::left_join
) {
  if (length(feature_ids) < 1) feature_ids <- true$feature_id

  .joined <- .join(est, true, join_) %>%
    filter(feature_id %in% feature_ids)

  nrow(.joined) %>% message()

  .joined
}


compare2 <- function(
  est,
  true,
  feature_ids1 = c(),
  feature_ids2 = c(),
  join_ = dplyr::left_join
) {
  if (length(feature_ids1) < 1) feature_ids1 <- true$feature_id
  if (length(feature_ids2) < 1) feature_ids2 <- true$feature_id

  .joined <- .join(est, true, join_) %>%
    filter((group == "ctrl" & feature_id %in% feature_ids1) | (group == "case" & feature_id %in% feature_ids2))

  nrow(.joined) %>% message()

  .joined
}


get_comparisions <- function(x, ests, true) {
  f <- function(
    ids,
    ests,
    true,
    m = purrr::map,
    f
  ) {
    ests %>%
      mutate(
        value = m(
          data,
          ~ f(.x, true, ids)
        )
      ) %>% select(-data)
  }

  x %>%
    mutate(
      comparision = map(feature_ids, ~ f(.x, ests, true, m = purrr::map, f = compare))
    )
}


calc_metrics <- function(est, true, feature_ids = c()) {
  .metrics <- tibble(
    spearman = calc_spearman(est = est, true = true, feature_ids = feature_ids),
    nrmse = calc_nrmse(est = est, true = true, feature_ids = feature_ids)
  )

}


..calc_metrics.. <- function(x, ests, true) {
  f <- function(
    ids,
    ests,
    true,
    m = purrr::map,
    f
  ) {
    ests %>%
      mutate(
        value = m(
          data,
          ~ f(.x, true, ids)
        )
      ) %>% select(-data)
  }

  x %>%
    mutate(
      spearman = map(feature_ids, ~ f(.x, ests, true, m = purrr::map_dbl, f = calc_spearman))
    ) %>%
    mutate(
      nrmse = map(feature_ids, ~ f(.x, ests, true, m = purrr::map_dbl, f = calc_nrmse))
    )
}


calc_metrics2 <- function(est, true, feature_ids1 = c(), feature_ids2 = c()) {
  .metrics <- tibble(
    spearman = calc_spearman2(est = est, true = true, feature_ids1 = feature_ids1, feature_ids2 = feature_ids2),
    nrmse = calc_nrmse2(est = est, true = true, feature_ids1 = feature_ids1, feature_ids2 = feature_ids2)
  )

}


..calc_metrics2.. <- function(x, ests, true) {
  f <- function(
    ids1,
    ids2,
    ests,
    true,
    f,
    m = purrr::map
  ) {
    ests %>%
      mutate(
        value = m(
          data,
          ~ f(.x, true, ids1, ids2)
        )
      ) %>% select(-data)
  }

  x %>%
    mutate(
      spearman = map2(feature_ids1, feature_ids2, ~ f(.x, .y, ests, true, f = calc_spearman2, m = purrr::map_dbl))
    ) %>%
    mutate(
      nrmse = map2(feature_ids1, feature_ids2, ~ f(.x, .y, ests, true, f = calc_nrmse2, m = purrr::map_dbl))
    )
}


calc_spearman <- partial(calc_metric, method = partial(cor, method = "spearman", use = "complete.obs"))


calc_spearman2 <- partial(calc_metric2, method = partial(cor, method = "spearman", use = "complete.obs"))


calc_nrmse <- partial(calc_metric, method = partial(hydroGOF::nrmse, na.rm = TRUE))


calc_nrmse2 <- partial(calc_metric2, method = partial(hydroGOF::nrmse, na.rm = TRUE))
