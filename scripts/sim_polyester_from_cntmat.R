#
# Generate simulation data from kallisto"s estimated abundance using polyester
#
# Usage:
#   Rscript this.R <transcripts.fasta> <count_matrix.tsv> [readlen] [output_dir]
#


# requirements
options(stringAsFactors = FALSE)

library(tidyverse)
library(readr)
library(data.table)
library(tools)

library(Biostrings)
library(ballgown)
library(polyester)


# set wd to ../../ from this script
root_path <- function() {
  tryCatch({
        return_path <- rstudioapi::getSourceEditorContext()$path %>% dirname %>% dirname
      }, error = function(e) {
        return_path <<- getwd()
      }, finally = {
        return(return_path)
      }
  )
}

root_path() %>% setwd

argv <- commandArgs(trailingOnly = TRUE)

fasta_path <- normalizePath(argv[1])
count_path <- normalizePath(argv[2])
readlen    <- as.numeric(argv[3])
output_dir <- argv[4]

message("Argments:")
print(argv)

if (is.na(readlen)) {
  readlen <- 100
}

if (is.na(output_dir)) {
  output_dir <- "."
}

#
# Load transcript sequences
#
# CAUTION: to avoid polyester"s bug, don"t use gtf
# https://github.com/alyssafrazee/polyester/issues/24
fasta <- fasta_path %>% readDNAStringSet

#
# Load counts and shape to matrix
#
counts <- fread(count_path, header = TRUE, sep = "\t")
counts <- data.frame(feature_id = names(fasta)) %>% left_join(counts, by = "feature_id")
counts[is.na(counts)] <- 0

readmat <- counts[, -1] %>% as.matrix
readmat <- readmat %>% round
rownames(readmat) <- counts[, 1]

#
# Generate reads
#
PAIRED          <- TRUE
READLEN         <- readlen # default
# FRAGLEN       <- 250 # default
ERROR_RATE      <- 0.005 # default
BIAS            <- "none" # CHANGE: Due to bias options include bug, change from "rnaf" to "none"
FRAG_GC_BIAS    <- "none" # default
STRAND_SPECIFIC <- TRUE # CAUTION: this parameter will generate FR directed reads (NOT RF)
SEED            <- 12345
OUTDIR          <- output_dir

if (!file.exists(OUTDIR)) dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
counts %>% write_tsv(paste(OUTDIR, "readmat.tsv", sep = "/"))

simulate_experiment_countmat(fasta = fasta_path,
                            readmat = readmat, paired = PAIRED,
                            readlen = READLEN, error_rate = ERROR_RATE,
                            bias = BIAS, strand_specific = STRAND_SPECIFIC,
                            seed = SEED, outdir = OUTDIR)
