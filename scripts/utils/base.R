# eval_rnaseq_map.utils
# ~~~~~~~~~~~~~~~~~~~~~
#
# This module procided some utiilty functions

library(tidyverse)


PSEUDO_COUNT <- 0.01
NUM_BINS <- 4
NUM_BINS2 <- 81


load_gtf <- function(path, tags = c("gene_id", "gene_name", "gene_type", "transcript_id", "transcript_name", "transcript_type"), types = c("transcript")) {
  rtracklayer::readGFF(path, version = 2L, tags = tags, filter = list(type = types)) %>%
    rename(feature = type)
}

# NOTE: Load with manually parsing
load_gtf_ <- function(path, types = c("gene", "transcript", "exon")) {
  .attributes <- path %>%
    rtracklayer::readGFF(nrows = 100,
            tags = character(0)) %>% .$attributes

  attribute_keys <- .attributes %>%
    gsub("; ", ";", .) %>%
    strsplit(";") %>%
    map(trimws) %>%
    map(~strsplit(., " ")) %>%
    unlist %>%
    .[c(TRUE, FALSE)] %>% unique

  tags <- c("gene_id", "gene_name", "transcript_id", "transcript_name")

  # NOTE: This Order is important
  extra_keys <- c("gene_type", "class", "geneSuperClass", "tcat", "transcript_type")

  appended_key <- NA

  for (k in extra_keys) {
    if (any(grepl(k, attribute_keys))) {
      appended_key <- k
      tags <- c(tags, k)
    }
  }

  paste0("tags: ", paste(tags, collapse = ", ")) %>% message

  gtf <- rtracklayer::readGFF(path, version = 2L, tags = tags, filter = list(type = types)) %>%
    rename(feature = type)

  # NOTE: NONCODE has no transcript type because all records are lncRNA
  if (is.na(appended_key)) {
    "No key appended. Fill transcript type with `lncRNA`." %>% message
    gtf <- gtf %>%
      mutate(transcript_type = "lncRNA")

    appended_key <- "transcript_type"
  }

  if (grepl("lncRNA_", path) || grepl("long_noncoding_RNAs", path)) {
    gtf$gene_type <- "lncRNA"
    gtf$transcript_type <- "lncRNA"
  }

  if (grepl("pc_", path) || grepl(".pc.", path)) {
    gtf$gene_type <- "protein_coding"
    gtf$transcript_type <- "protein_coding"
  }

  gtf <- gtf %>%
    rename(transcript_type = !!appended_key)

  # NOTE: If lack of gene_type,
  # fill gene type with concatenation of transcript_type(s) (protein coding > lncRNA > others)

  .filled_fantom <- function(gtf) {
    gtf$gene_type <- "others"

    if (any(gtf$geneSuperClass == "all_lncRNA")) gtf[gtf$geneSuperClass == "all_lncRNA",]$gene_type <- "lncRNA"
    if (any(gtf$geneSuperClass == "all_mRNA")) gtf[gtf$geneSuperClass == "all_mRNA",]$gene_type <- "protein_coding"

    gtf
  }

  if (any(grepl("geneSuperClass", colnames(gtf)))) {
    gtf <- .filled_fantom(gtf)
  }

  .filled_generic <- function(gtf) {
    .genes <- gtf %>%
      group_by(gene_id) %>%
      summarize(transcript_types = toString(sort(unique(transcript_type))))

    .gene_ids_lncrna <- .genes %>%
      filter(grepl("lncRNA", transcript_types)) %>% .$gene_id

    .gene_ids_pc <- .genes %>%
      filter(grepl("protein_coding", transcript_types) | grepl("mRNA", transcript_types)) %>% .$gene_id

    gtf$gene_type <- "others"
    if (length(.gene_ids_lncrna) > 0) gtf[gtf$gene_id %in% .gene_ids_lncrna,]$gene_type <- "lncRNA"
    if (length(.gene_ids_pc) > 0) gtf[gtf$gene_id %in% .gene_ids_pc,]$gene_type <- "protein_coding"

    gtf
  }

  if (!any(grepl("gene_type", colnames(gtf)))) {
    gtf <- .filled_generic(gtf)
  }

  gtf
}


load_annotations <- function(types = c("transcript"), feature_id_ = "transcript_id") {
  gencode <- here::here("shared/assets/references/grch38/annotations/gencode/gencode.v31.annotation.gtf") %>% load_gtf(types = types)
  gencode_basic <- here::here("shared/assets/references/grch38/annotations/gencode/gencode.v31.basic.annotation.gtf") %>% load_gtf(types = types)
  gencode_refseq <- here::here("shared/assets/references/grch38/annotations/gencode_refseq/gencode.v31_refseq.v109.20190607.gtf") %>% load_gtf(types = c("exon"))

  # NOTE: Append general biotype to comprehensive annotation
  gencode_lncrna <- here::here("shared/assets/references/grch38/annotations/gencode/gencode.v31.long_noncoding_RNAs.gtf") %>% load_gtf(types = types)

  gencode$biotype <- "others"
  transcript_ids_lncrna <- gencode_lncrna$transcript_id[!is.na(gencode_lncrna$transcript_id)] %>% unique()
  gencode$biotype[gencode$transcript_id %in% transcript_ids_lncrna] <- "lncrna"
  gencode$biotype[gencode$transcript_type == "protein_coding"] <- "mrna"

  gencode_basic$biotype <- "others"
  gencode_basic$biotype[gencode_basic$transcript_id %in% transcript_ids_lncrna] <- "lncrna"
  gencode_basic$biotype[gencode_basic$transcript_type == "protein_coding"] <- "mrna"

  # HACK: To be simple
  gencode_refseq$biotype <- "others"
  gencode_refseq$biotype[gencode_refseq$transcript_type %in% c("protein_coding", "mRNA")] <- "mrna"
  gencode_refseq$biotype[gencode_refseq$transcript_id %in% transcript_ids_lncrna] <- "lncrna"
  transcript_ids_refseq <- setdiff(gencode_refseq$transcript_id, gencode$transcript_id)
  gencode_refseq$biotype[(gencode_refseq$transcript_id %in% transcript_ids_refseq) & (gencode_refseq$transcript_type == "lncRNA")] <- "lncrna"

  gtfs <- list(
    gencode = gencode,
    gencode_basic = gencode_basic,
    gencode_refseq = gencode_refseq
  ) %>%
  map(~.x %>% mutate(feature_id = get(feature_id_)))
}


feature_length <- function(gtf, feature_ = "exon", key_feature_ = "transcript_id") {
  gtf %>%
    filter(feature == feature_) %>%
    group_by(.dots = key_feature_) %>%
    summarize(sum(end - start)) %>%
    data.frame
}


to_factor <- function(x) {
  if (any(x %in% c("mrna", "lncrna", "others", "all"))) {
    return(
      factor(x,
        labels = c("mRNA", "lncRNA", "Others", "All"),
        levels = c("mrna", "lncrna", "others", "all")
      )
    )
  }
  if (any(x %in% c("incomplete", "complete", "complicated"))) {
    return(
      factor(x,
        labels = c("Incomplete", "Complete", "Complicated"),
        levels = c("incomplete", "complete", "complicated")
      )
    )
  }
  if (any(x %in% c("gencode_basic", "gencode", "gencode_refseq"))) {
    return(
      factor(x,
        labels = c("GENCODE-Basic", "GENCODE", "GENCODE+Refseq"),
        levels = c("gencode_basic", "gencode", "gencode_refseq")
      )
    )
  }
}


to_combination <- function(x) {
  # FIXME: Workaround
  suffix <- ifelse(endsWith(x, "_wt.tsv"), "_wt", "")

  x <- dirname(x)

  if (!is.na(strsplit(x, "gencode")[[1]][2])) {
    x <- paste0("gencode", strsplit(x, "gencode")[[1]][2])
  }

  replacements <-
    c(
      "/conv_sam_to_bam",
      "/conv_rsem_to_ebseq_matrix",
      "/prep_prepde",
      "/align_",
      "/quant_",
      "/de_",
      "/transcript",
      "/gene",
      "/summary.txt",
      "/conv_stringtie_to_raw",
      "/conv_any_to_raw",
      "/conv_any2_raw",
      "/conv_any2_raw_tximport",
      "/conv_rsem_to_matrix",
      "/conv_cuffdiff_to_raw",
      "_nofilter"
    )

  for (r in replacements) {
    x <- gsub(r, "/", x)
  }

  x <- gsub("//", "/", x)
  x <- gsub("/", "-", x)
  x <- gsub("-$", "", x)
  x <- paste0(x, suffix)

  x
}


to_xlab <- function(x) {
  XLABS <- list(
    length = "Read length (bases)",
    depth = "Library size (million reads)",
    reps = "Number of replicates"
  )

  XLABS[[x]]
}


find_paths <- function(input_dir, patterns) {
  .all <- input_dir %>% list.files(recursive = TRUE)

  .paths <- patterns %>%
    map(function(p) .all[grepl(p, .all)]) %>%
    map(function(x) file.path(input_dir, x))

  .paths
}


annotation_used <- function(x) {
  PATTERNS <- c("gencode_refseq", "gencode_basic", "gencode")
  for (p in PATTERNS) {
    if (grepl(p, x)) {
      return(p)
    }
  }

  NULL
}


replaces <- function(x, replacements) {
  for (r in names(replacements)) {
    x <- gsub(r, replacements[[r]], x)
  }

  x
}


to_scenario <- function(x, factorize = TRUE, capitalize = FALSE) {
  SCENARIOS <- c(
    "gencode_refseq" = "complicated",
    "gencode_basic" = "incomplete",
    "gencode" = "complete"
  )

  scenario <- replaces(x, SCENARIOS)
  if (capitalize) scenario <- to_capital(scenario)
  if (factorize) scenario <- to_factor(scenario)

  scenario
}


to_abbr <- function(x, to_scenario = TRUE) {
  ABBRS <- c(
    "star" = "Sr",
    "rsem" = "Rs",
    "ebseq" = "Eb",
    "hisat2" = "Hs",
    "stringtie" = "St",
    "ballgown" = "Ba",
    "tophat2" = "Th",
    "cuffdiff" = "Cu",
    "kallisto" = "Ka",
    "sleuth_lrt" = "Sl(LRT)",
    "sleuth_wt" = "Sl",
    "sleuth" = "Sl",
    "edger" = "Ed"
  )

  if (to_scenario) {
    ABBRS <- c(
      ABBRS,
      c(
        "gencode_refseq" = "Complicated",
        "gencode_basic" = "Incomplete",
        "gencode" = "Complete"
      )
    )
  } else {
    ABBRS <- c(
      ABBRS,
      c(
        "gencode_refseq" = "GENCODE+RefSeq",
        "gencode_basic" = "GENCODE.Basic",
        "gencode" = "GENCODE"
      )
    )
  }


  replaces(x, ABBRS)
}


to_tool <- function(x) {
  TOOLS <- c(
    "Sr" = "star",
    "Rs" = "rsem",
    "Eb" = "ebseq",
    "Hs" = "hisat2",
    "St" = "stringtie",
    "Ba" = "ballgown",
    "Th" = "tophat2",
    "Cu" = "cuffdiff",
    "Ka" = "kallisto",
    "Sl(LRT)" = "sleuth_lrt",
    "Sl(WT)" = "sleuth_wt",
    "Sl" = "sleuth",
    "Ed" = "edger"
  )

  replaces(x, TOOLS)
}


transcript_id_to_biotype <- function(x, annotation) {
  if (is.na(!any(ls(envir = utils) %>% charmatch("biotypes")))) {
    "Call" %>% print()
    biotypes <- annotation %>%
      select(transcript_id, biotype) %>%
      deframe()
  }

  biotypes[[x]]
}


mean_by_group <- function(x, prefixes = c("ctrl", "case"), ordered = FALSE) {
  if (ordered) {
    n_rep <- ncol(x) / 2
    groups <- c(rep("ctrl", n_rep), rep("case", n_rep))
    colnames(x) <- paste0(groups, "_", seq(n_rep) - 1)
    prefixes <- c("ctrl", "case")
  }

  axis_ <- x %>% select(!where(is.numeric))
  mean_ <- data.frame(
    ctrl = x[, grepl(prefixes[1], colnames(x))] %>% apply(1, mean),
    case = x[, grepl(prefixes[2], colnames(x))] %>% apply(1, mean)
  )

  cbind(axis_, mean_)
}


tpm <- function(x, lengths, log = FALSE, prior.count = PSEUDO_COUNT) {
  if (log) x <- x + prior.count
  rownames(lengths) <- pull(lengths, 1)
  .length <- lengths[rownames(x),]
  .tpm <- t(t((x / pull(.length, 2))) * 1e6 / colSums((x / pull(.length, 2))))
  if (log) .tpm <- log2(.tpm)

  .tpm
}


.condition <- function(x) {
  x %>%
    strsplit("\\.") %>%
    sapply("[", 1) %>%
    as.numeric()
}


.annotation <- function(x) {
  x %>%
    strsplit("\\.") %>%
    sapply("[", 2) %>%
    strsplit("-") %>%
    sapply("[", 1)
}


.combination <- function(x) {
  x %>%
    strsplit("\\.") %>%
    sapply("[", 2) %>%
    strsplit("-") %>%
    lapply("[", -1) %>%
    sapply(paste, collapse = "-")
}


.last <- function(x, sep = "-") {
  if (sep == ".") sep <- paste0("\\", sep)

  x %>%
    strsplit(sep) %>%
    sapply(function(x) x[length(x)])
}


.except_last <- function(x, sep = "-") {
  collapse_ <- sep
  if (sep == ".") sep <- paste0("\\", sep)

  x %>%
    strsplit(sep) %>%
    lapply(function(x) x[-length(x)]) %>%
    sapply(paste, collapse = collapse_)
}


tool <- .last


upper <- .except_last


decorate_metrics <- function(metrics) {
  metrics %>%
    mutate(condition = .condition(name)) %>%
    mutate(annotation = .annotation(name)) %>%
    mutate(scenario = to_scenario(annotation)) %>%
    mutate(combination = .combination(name)) %>%
    mutate(tool = tool(name)) %>%
    mutate(abbr = to_abbr(combination)) %>%
    mutate(upper_abbr = upper(abbr)) %>%
    mutate(upper_abbr = ifelse(upper_abbr == "", abbr, upper_abbr))
}


split_tibble <- function(tibble, col) {
  tibble %>% split(., .[, col])
}


to_capital <- function(x) {
  s <- strsplit(x, " ")[[1]]

  paste(
    toupper(substring(s, 1, 1)),
    substring(s, 2),
    sep = "", collapse = " "
  )
}


biotype_to_feature_ids <- function(annotation, biotype_ = "all") {
  if (biotype_ == "all") biotype_ <- TRUE

  annotation %>%
    filter(biotype %in% biotype_) %>%
    .$feature_id
}


multi_join <- function(list_df, by, method = "inner_join", cols = c()) {
  rename_ <- function(x, y) {
    i <- grep(paste(by, collapse = "|"), colnames(x))
    colnames(x)[-i] <- paste(colnames(x)[-i], y, sep = ".")

    x
  }

  # NOTE: Add postfix as value.a, value.b, value.c...
  list_df_renamed_ <- map2(
    list_df,
    letters[1:length(list_df)],
    rename_
  )

  method_ <- get(method)
  joined_ <- Reduce(function(x, y) method_(x, y, by = by), list_df_renamed_)
  if (length(cols) > 1) joined_ <- joined_ %>% select(matches(paste(cols, collapse = "|")))

  joined_
}


compare_metrics <- function(x, key, by) {
  x %>%
    split_tibble(key) %>%
    multi_join(by = by) %>%
    select(by, matches(paste(key, "value", "name", sep = "|")))
}


filter_by_condition <- function(x, condition = "main") {
  CONDITIONS <- list(
    main = c(100, 20, 3, 0)
  )

  .cond <- CONDITIONS[[condition]]
  x %>% filter(condition %in% .cond)
}


filter_by_combination <- function(x, combination = "main", col = "combination") {
  COMBINATIONS <- list(
    main = c(
      "star",
      "star-rsem",
      "star-rsem-ebseq",
      "hisat2",
      "hisat2-stringtie",
      "hisat2-stringtie-ballgown",
      "tophat2",
      "tophat2-cuffdiff",
      "kallisto",
      "kallisto-sleuth",
      "kallisto-sleuth_wt"
    )
  )

  .comb <- COMBINATIONS[[combination]]

  .names <- x %>% pull(col)

  x[.names %in% .comb,]
}

mode <- function(x) {
  .y <- table(x) %>%
    enframe %>% filter(value == max(value)) %>%
    select(2, 1) %>%
    deframe

  .names <- names(.y)

  .y <- as.numeric(.y)
  names(.y) <- .names

  .y
}


min_max <- function(x) {
  c(min(x), max(x))
}


set_biotype <- function(df, annotation) {
  biotypes_ <- annotation %>% distinct(feature_id, biotype)
  df %>%
    left_join(
      biotypes_,
      by = "feature_id"
    )
}


set_bin <- function(x, col, num_bins = NUM_BINS, exclude = FALSE) {
  .breaks <- function(x, num_bins) {
    x %>%
      pull(col) %>%
      quantile(seq(0, 1, 1 / num_bins), na.rm = TRUE)
  }

  .labels <- function(num_bins) {
    .ends <- seq(0, 1, 1 / num_bins) * 100
    .starts <- .ends
    .starts <- .starts[-length(.starts)]
    .starts[1] <- 0

    paste0("\u2265 ", .starts, "%")
  }

  .x <- x
  if (exclude) {
    "mode was excluded." %>% message
    paste0(mode(x[col]), ":", names(mode(x[col])))  %>% message

    .x <- x %>% filter(!(UQ(rlang::sym(col)) %in% mode(UQ(rlang::sym(col))) | is.na(UQ(rlang::sym(col)))))
  }

  breaks <- .breaks(.x, num_bins)
  paste0(col, " bins:") %>% message
  breaks %>% paste0(sep = " ") %>% message

  # NOTE: Overflow values replace to max bin
  breaks[length(breaks)] <- +Inf
  .max <- as.character(round(max(pull(.x, col)), 3))
  .x$bin <- cut(pull(.x, col), breaks, include.lowest = TRUE, right = FALSE) %>% fct_relabel(~gsub("Inf", .max, .x))

  if (exclude) {
    .x <- .x %>% bind_rows(
      x %>%
        filter(UQ(rlang::sym(col)) %in% mode(UQ(rlang::sym(col))) | is.na(UQ(rlang::sym(col)))) %>%
        mutate(bin = NA)
    )
  }

  .x
}


remove_spikein <- function(mat, prefix = 'ERCC-') {
  keep <- !grepl(prefix, rownames(mat))
  return(mat[keep,])
}


calc_slope <- function(x) {
  .x <- x %>% mutate(bin = as.numeric(factor(bin)))

  return(lm(.x$value ~ .x$bin)$coefficients[[2]])
}

rm_private_vars <- function(){
  rm(list = setdiff(ls(envir = .GlobalEnv, all.names = TRUE), ls(envir = .GlobalEnv)), envir = .GlobalEnv)
}
