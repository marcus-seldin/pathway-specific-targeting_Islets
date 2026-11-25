#!/usr/bin/env Rscript

# Usage: Rscript prioritize_SERPINB1.R DATASET_DIR NEW_OUTPUT_DIR
#
# Main pancreatic target: beta-cell proliferation and its regulation.
#   GO:0044342, GO:0061469, GO:1904692 — type B pancreatic cell
#     proliferation, its regulation, and its positive regulation.
#     The combined 18-gene panel placed liver SERPINB1 #353 of 35,737
#     secreted pairs overall and #38 among secreted liver pairs. This is
#     consistent with experiments in which liver-derived SerpinB1
#     promotes beta-cell proliferation (PMID 26701651).
#   GO:0010951 — negative regulation of endopeptidase activity.
#     This 53-gene panel reflects SerpinB1's protease-inhibitor activity
#     and placed liver SERPINB1 #560 of 35,737 secreted pairs.
#
# These ranks use all donors. Only 14 matched liver–pancreas donors
# were available in the held-out split, below the 25-pair threshold.

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
workflow <- new.env()
workflow$screen_anchor <- "SERPINB1"
sys.source(file.path(dirname(normalizePath(script)), "pathway_prior_screen.R"), envir = workflow)
