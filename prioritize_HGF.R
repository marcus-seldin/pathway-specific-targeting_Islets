#!/usr/bin/env Rscript

# Usage: Rscript prioritize_HGF.R DATASET_DIR NEW_OUTPUT_DIR [discovery|go]
#
# Best tested approach: the 25 pancreatic genes most strongly correlated
# with visceral adipose HGF in discovery donors. This panel placed HGF
# #1,791 of 220,810 rankable gene–tissue pairs in held-out donors.
# Correlation signs agreed for 24 of 25 transcripts. Leading discovery
# genes were TIGD5, ALDH1L1, and TREX1. No GO term passed the enrichment
# threshold, so there is no supported core GO pathway for this panel.
#
# Optional GO mode tests two literature-defined pancreatic processes:
#   GO:0032024 — positive regulation of insulin secretion.
#     This 45-gene panel placed subcutaneous adipose HGF #3,934 of
#     486,382 pairs across all donors. Islet HGF increased insulin
#     production and improved islet function in mouse studies.
#   GO:0044342, GO:0061469, GO:1904692 — pancreatic beta-cell
#     proliferation and its regulation. This 18-gene panel placed
#     subcutaneous adipose HGF #5,513 overall across all donors.
# See PMID 10625667 and 11723058. Use HGF's all-candidate rank because
# the historical UniProt secretion list omits HGF.

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) == 3L) args[[3L]] else "discovery"
if (!length(args) %in% c(2L, 3L) || !mode %in% c("discovery", "go")) {
  stop("Usage: Rscript prioritize_HGF.R DATASET_DIR NEW_OUTPUT_DIR [discovery|go]")
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
engine <- if (mode == "go") "pathway_prior_screen.R" else "reverse_pancreas.R"
workflow <- new.env()
workflow$screen_anchor <- "HGF"
workflow$screen_args <- args[1:2]
sys.source(file.path(dirname(normalizePath(script)), engine), envir = workflow)
