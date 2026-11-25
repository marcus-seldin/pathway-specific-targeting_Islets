#!/usr/bin/env Rscript

# Usage: Rscript prioritize_FGF21.R DATASET_DIR NEW_OUTPUT_DIR
#
# Test two pancreatic pathways motivated by experimental evidence that
# FGF21 protects beta cells during lipid stress (PMID 16936195, 31654570).
#   GO:2001243 — negative regulation of intrinsic apoptotic signaling.
#     This 25-gene panel gave the best alternative rank in the tested
#     screen: visceral adipose FGF21 was #3,724 of 35,737 secreted pairs.
#   GO:0070542 and GO:0071398 — response and cellular response to fatty
#     acid. The combined 41-gene panel placed liver FGF21 #3,815 of 35,737.
#
# These ranks use all donors. The original insulin-secretion panel placed
# liver FGF21 #3,639. Neither stress-response panel improved that rank.

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
workflow <- new.env()
workflow$screen_anchor <- "FGF21"
sys.source(file.path(dirname(normalizePath(script)), "pathway_prior_screen.R"), envir = workflow)
