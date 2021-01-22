# eval_rnaseq_map.helpers
# ~~~~~~~~~~~~~~~~~~~~~~~
#
# This module procided some utiilty helpers for DE step

library(tidyverse)


PATTERNS <- list(
  transcript = c(
    cuffdiff = "de_cuffdiff/isoform_exp.diff$",
    ebseq = "de_ebseq/transcript/result.tsv$",
    ballgown = "de_ballgown/result_transcript.tsv$",
    sleuth = "de_sleuth/result_transcript_wt.tsv",
    edger = "de_edger/transcript/result_.*.tsv$"
  ),
  gene = c(
    cuffdiff = "de_cuffdiff/gene_exp.diff$",
    ebseq = "de_ebseq/gene/result.tsv$",
    ballgown = "de_ballgown/result_gene.tsv$",
    sleuth = "de_sleuth/result_gene_wt.tsv",
    edger = "de_edger/gene/result_.*.tsv$"
  )
)


load_groundtruth <- function(path) {
  .data <- path %>% data.table::fread(header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  .data <- .data %>% mutate(logfc = log2(fc), pval_adj = NA)
  .data <- .data %>% select(feature_id, logfc, pval_adj)

  .data
}


load_edger <- function(path) {
  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  .data <- .data %>% select(GeneID, logFC, FDR)
  colnames(.data) <- c("feature_id", "logfc", "pval_adj")

  .data
}


load_ebseq <- function(path) {
  .data <- data.table::fread(path, header = FALSE, skip = 1, sep = "\t", stringsAsFactors = FALSE)

  .data <- .data %>% mutate(logfc = -1 * log2(V4))
  .data <- .data %>% select(V1, logfc, V2)
  colnames(.data) <- c("feature_id", "logfc", "pval_adj")

  .data
}


load_ballgown <- function(path, coef = 1) {
  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  .data <- .data %>% mutate(logfc = coef * log2(as.numeric(fc)))

  tryCatch(
    {
      .data <- .data %>% select(transcriptIDs, logfc, qval)
    },
    error = function(e) {
      .data <<- .data %>% select(`id`, logfc, qval)
    }
  )

  colnames(.data) <- c("feature_id", "logfc", "pval_adj")

  .data
}


load_cuffdiff <- function(path) {
  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  .data <- .data %>% filter(status == "OK")

  .data <- .data %>% select(test_id, `log2(fold_change)`, q_value)
  colnames(.data) <- c("feature_id", "logfc", "pval_adj")

  .max <- (.data %>% filter(is.finite(logfc)))$logfc %>% max()
  .min <- (.data %>% filter(is.finite(logfc)))$logfc %>% min()

  tryCatch(
    {
      .data[.data$logfc == Inf, ]$logfc <- .max
    },
    error = function(e) {
      .data <<- .data
    }
  )
  tryCatch(
    {
      .data[.data$logfc == -Inf, ]$logfc <- .min
    },
    error = function(e) {
      .data <<- .data
    }
  )

  .data %>% filter(is.finite(logfc))
}



load_sleuth <- function(path) {
  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  .data <- .data %>% select(target_id, qval)
  colnames(.data) <- c("feature_id", "pval_adj")

  .data
}


load_sleuth_wt <- function(path, coef = -1) {
  .data <- data.table::fread(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  tryCatch({
    .data <- .data %>% select(target_id, b, qval)
    colnames(.data) <- c("feature_id", "logfc", "pval_adj")

    .data <- .data %>% mutate(logfc = coef * logfc)
  },
  error = function(e) {
    message(paste0(e, "Return df filled logfc column with NA."))
    .data <<- .data %>%
      mutate(b = NA) %>%
      select(target_id, b, qval)
    colnames(.data) <<- c("feature_id", "logfc", "pval_adj")
  })

  .data
}


extract_de <- function(x, t_logfc, t_pval_adj) {
  if (any(colnames(x) %in% "logfc")) {
    return(x %>% filter(abs(logfc) >= t_logfc & pval_adj < t_pval_adj))
  }

  x %>% filter(pval_adj < t_pval_adj)
}


fill_lack <- function(x, annotation) {
  .feature_ids <- union(
    x$feature_id,
    annotation$feature_id
  )
  .x <- data.frame(
    feature_id = .feature_ids,
    stringsAsFactors = FALSE
  )

  .x <- .x %>%
    left_join(x, by = "feature_id") %>%
    mutate(is_lack = FALSE)

  .is_lack <- is.na(.x$logfc)
  if (sum(.is_lack) > 0) {
    .x[.is_lack, ]$is_lack <- TRUE
    .x[.is_lack, ]$logfc <- 0
    .x[.is_lack, ]$pval_adj <- 1
  }

  .x
}



.join_fill <- function(est, true, feature_ids, join_ = dplyr::left_join) {
  if (length(feature_ids) < 1) feature_ids <- true$feature_id

  .joined <- true %>%
    join_(est, by = "feature_id")

  if (length(feature_ids) > 0) .joined <- .joined %>% filter(feature_id %in% feature_ids)
  nrow(.joined) %>% message()

  .is_lack <- is.na(.joined$logfc.y)
  if (sum(.is_lack) > 0) {
    .joined[.is_lack, ]$logfc.y <- 0

    tryCatch(
      {
        .joined[.is_lack, ]$pval_adj.y <- 1
      },
      error = function(e) {
        paste0("Error occured: ", e) %>% message()
      }
    )
  }

  .joined
}


.calc_metric <- function(
  est,
  true,
  method,
  join_ = dplyr::left_join,
  feature_ids = c()
) {

  .joined <- .join_fill(est, true, feature_ids, join_)

  tryCatch(
    {
      .metric <- method(.joined$logfc.x, .joined$logfc.y)
    },
    error = function(e) {
      paste0("Error occured: ", e) %>% message()
      .metric <<- NA
    }
  )

  .metric
}


counts_n_tested <- function(
  est,
  feature_ids = c()
) {

  n <- 0

  tryCatch(
    {
      n <- est %>%
        filter(feature_id %in% feature_ids & !is_lack) %>%
        nrow
    },
    error = function(e) {
      n <<- est %>%
        filter(feature_id %in% feature_ids) %>%
        nrow

    }
  )

  n
}


calc_spearman <- partial(.calc_metric, method = partial(cor, method = "spearman", use = "complete.obs"))


calc_pearson <- partial(.calc_metric, method = partial(cor, method = "pearson", use = "complete.obs"))


calc_nrmse <- partial(.calc_metric, method = partial(hydroGOF::nrmse, na.rm = TRUE))


calc_roc_df <- function(est, true, feature_ids = c(), cutoff_logfc = 1, join_ = dplyr::left_join) {
  .true <- true
  .true$D <- (abs(.true$logfc) >= cutoff_logfc) %>% as.integer()

  .joined <- .join_fill(est, .true, feature_ids, join_)
  .joined <- .joined %>% mutate(m = pval_adj.y)

  # NOTE: Set 1 for records that are not DE
  tryCatch(
    {
      .joined$m[!(abs(.joined$logfc.y) >= cutoff_logfc)] <- 1
    },
    error = function(e) {
      "Lack of `logfc` column. skipped." %>% message()
    }
  )
  .joined$m[is.na(.joined$m)] <- 1

  .joined %>%
    mutate(M = 1 - m) %>%
    select(D, M)
}


calc_auc <- function(x) {
  x %>%
    utl$plot_roc_() %>%
    plotROC::calc_auc() %>%
    rename(value = AUC) %>%
    select(value) %>%
    tibble %>%
    .$value
}


calc_metrics <- function(est, true, feature_ids = c()) {
  .est <- est
  if (length(feature_ids) > 0) {
    .est <- .est %>% filter(feature_id %in% feature_ids)
  }

  .metrics <- tibble(
    spearman = calc_spearman(est = est, true = true, feature_ids = feature_ids),
    nrmse = calc_nrmse(est = est, true = true, feature_ids = feature_ids)
  )

  .metrics <- .metrics %>%
    mutate(roc = list(calc_roc_df(est = est, true = true, feature_ids = feature_ids)))

  .metrics <- .metrics %>%
    mutate(auc = map_dbl(roc, ~ calc_auc(.x)) %>% unlist) %>%
    select(-roc)
}


calc_intersects <- function(x, rel = TRUE) {
  .sets <- x$feature_ids

  .n_intersect <- crossing(target = .sets, query = .sets) %>%
    pmap(function(target, query) intersect(target, query)) %>%
    map(length) %>%
    unlist %>%
    matrix(ncol = nrow(x)) %>% data.frame

  .names <- x %>% pull(abbr)

  .abs2rel_ <- function(x) {
    # NOTE: The query itself always equal max
    .denominator <- max(x)

    .rel <- x / .denominator
    .rel[is.nan(.rel)] <- 0

    .rel
  }

  if (rel) {
    .n_intersect <- .n_intersect %>% apply(1, .abs2rel) %>% data.frame
  }

  .rownames(.n_intersect) <- colnames(.n_intersect) <- .names
  .n_intersect
}
