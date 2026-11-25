#!/usr/bin/env Rscript

# Rank GTEx gene–tissue pairs by their association with pancreatic
# insulin-secretion gene expression, then plot the rankings and compare
# one candidate's cross-tissue and within-tissue scores.
#
# From the repository root:
#   Rscript "script for crosstissue correlation with pancreas pathways.R" /path/to/datasets /path/to/results
#
# Defaults: datasets ./datasets, output ./results, WGCNA threads 1.
# Set SLURM_CPUS_PER_TASK (provided by Slurm) to use allocated cores.
# Install WGCNA, ggplot2, and forcats before running.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript <script> [dataset_dir] [output_dir]")
dataset_dir <- if (length(args) >= 1L) args[[1L]] else "datasets"
output_dir <- if (length(args) >= 2L) args[[2L]] else "results"

input_files <- file.path(dataset_dir, c(
  expression = "GTEx NA included env.RData",
  phenotypes = "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
  secreted = "uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab",
  ontology = "uniprot-human-genes and goterms mapping.tab"
))
names(input_files) <- c("expression", "phenotypes", "secreted", "ontology")
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files)) {
  stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
}

library(ggplot2)
threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(threads) || threads < 1L) stop("SLURM_CPUS_PER_TASK must be a positive integer.")
if (threads > 1L) WGCNA::allowWGCNAThreads(nThreads = threads)

# 1. Load the manuscript's prefiltered GTEx v8 expression data.
#    The saved table has one row per gene–tissue pair. Transpose it so
#    that rows are donors and columns are gene–tissue pairs.
workspace <- new.env()
load(input_files[["expression"]], envir = workspace)
expression <- workspace$GTEx_subfiltered
if (is.null(expression$gene_tissue)) {
  stop("The GTEx workspace must contain GTEx_subfiltered with a gene_tissue column.")
}
rownames(expression) <- expression$gene_tissue
expression$gene_tissue <- NULL
expression <- as.data.frame(t(expression))
rm(workspace)

# The original workflow checks donor sex counts for context. Donor sex
# is neither a filter nor a covariate in the correlation analysis.
phenotypes <- read.delim(input_files[["phenotypes"]])
donor_ids <- sub("^GTEX-", "", phenotypes$SUBJID)
donor_sex <- ifelse(phenotypes$SEX == 1, "M", "F")
message("Sex counts for donors matched to the expression table:")
print(table(donor_sex[donor_ids %in% rownames(expression)]))

# 2. Define the pancreas target genes and annotate secreted candidates.
#    Use the historical UniProt files to reproduce the manuscript ranking.
#    Text matching includes GO descriptions such as regulation of insulin secretion.
ontology <- read.delim(input_files[["ontology"]])
pathway_genes <- ontology$Gene.names...primary..[
  grepl("insulin secretion", ontology$Gene.ontology..biological.process.)
]
secreted <- read.delim(input_files[["secreted"]], check.names = FALSE)
secreted_genes <- secreted[["Gene names  (primary )"]]

source_tissues <- c(
  "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)",
  "Small Intestine - Terminal Ileum", "Stomach", "Thyroid",
  "Muscle - Skeletal", "Pituitary", "Liver", "Kidney - Cortex",
  "Heart - Left Ventricle", "Colon - Sigmoid", "Adrenal Gland", "Artery - Aorta"
)
gene_symbol <- function(x) sub("_.*", "", x)
tissue_name <- function(x) sub(".*_", "", x)
target_tissue <- "Pancreas"

# 3. Calculate Ssec for each source gene–tissue pair.
#    WGCNA uses biweight midcorrelation with pairwise complete observations.
#    Ssec is the mean of -log10(p) over the target pathway genes. The
#    original workflow replaces missing and zero p-values with 0.5.
#    The score summarizes association strength across both directions.
score_pathway <- function(source, target) {
  if (!ncol(source) || !ncol(target)) {
    stop("No expression columns matched the source tissue or target pathway.")
  }
  p_values <- WGCNA::bicorAndPvalue(
    source, target, use = "pairwise.complete.obs"
  )$p
  p_values[is.na(p_values) | p_values == 0] <- 0.5
  scores <- data.frame(
    Ssec_score = rowMeans(-log10(p_values)),
    gene_tissue = colnames(source)
  )
  scores[order(scores$Ssec_score, decreasing = TRUE), , drop = FALSE]
}

expression_tissues <- tissue_name(colnames(expression))
pathway_columns <- gene_symbol(colnames(expression)) %in% pathway_genes
pancreas_pathway <- expression[, expression_tissues == target_tissue & pathway_columns,
                               drop = FALSE]
message(nrow(expression), " GTEx donor rows and ", ncol(pancreas_pathway),
        " pancreatic pathway expression columns.")
scores <- score_pathway(
  expression[, expression_tissues %in% source_tissues, drop = FALSE],
  pancreas_pathway
)

# 4. Label and export the rankings. The original significance flag uses
#    a descriptive cutoff of three standard deviations above the mean Ssec.
scores$gene_symbol <- gene_symbol(scores$gene_tissue)
scores$secreted <- ifelse(scores$gene_symbol %in% secreted_genes, "Secreted", "Non-secreted")
cutoff <- mean(scores$Ssec_score + 3 * sd(scores$Ssec_score))
scores$Sig_P1e3 <- ifelse(scores$Ssec_score > cutoff, "Significant", "NS")
scores$rank_overall <- seq_len(nrow(scores))
secreted_scores <- scores[scores$secreted == "Secreted", , drop = FALSE]
secreted_scores$rank_secreted <- seq_len(nrow(secreted_scores))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.table(scores, file.path(output_dir, "all_scores.tsv"), sep = "\t", row.names = FALSE)
write.table(secreted_scores, file.path(output_dir, "secreted_scores.tsv"),
            sep = "\t", row.names = FALSE)

# 5. Recreate the ranking plots. Use one color per secretion annotation.
plot_scores <- function(data, title) {
  ggplot(data, aes(x = forcats::fct_inorder(gene_tissue), y = Ssec_score,
                   fill = secreted)) +
    geom_col() +
    scale_fill_manual(values = c("Secreted" = "darkorange1", "Non-secreted" = "darkorchid")) +
    labs(title = title, x = "Gene and source tissue", y = "GTEx cross-tissue Ssec score",
         fill = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5, hjust = 1))
}

pdf(file.path(output_dir, "plots.pdf"), width = 12, height = 7)
significant <- scores[scores$Sig_P1e3 == "Significant", , drop = FALSE]
if (nrow(significant)) {
  pie(table(tissue_name(significant$gene_tissue)), main = "Tissue distribution of significant scores")
  pie(table(tissue_name(head(significant$gene_tissue, 100))), main = "Tissue distribution of top 100")
  print(plot_scores(head(significant, 50), "Top associations with pancreatic insulin secretion"))
}
if (nrow(secreted_scores)) {
  print(plot_scores(head(secreted_scores, 40), "Top secreted candidates"))
}

# 6. Compare within- and cross-tissue Ssec for one candidate.
#    The original script chose the fourth-ranked secreted gene–tissue pair.
#    For each comparison, divide Ssec by its mean across all source genes
#    before retrieving the selected gene's score.
comparison_rank <- 4L
if (nrow(secreted_scores) >= comparison_rank) {
  candidate <- secreted_scores[comparison_rank, ]
  origin <- tissue_name(candidate$gene_tissue)
  origin_columns <- expression_tissues == origin
  pancreas_columns <- expression_tissues == target_tissue
  comparisons <- list(
    scores[tissue_name(scores$gene_tissue) == origin, , drop = FALSE],
    score_pathway(expression[, pancreas_columns, drop = FALSE], pancreas_pathway),
    score_pathway(expression[, origin_columns, drop = FALSE],
                  expression[, origin_columns & pathway_columns, drop = FALSE])
  )
  names(comparisons) <- c(paste(origin, target_tissue, sep = " to "),
                          paste(target_tissue, target_tissue, sep = " to "),
                          paste(origin, origin, sep = " to "))
  comparison <- do.call(rbind, lapply(names(comparisons), function(label) {
    data <- comparisons[[label]]
    normalized <- data$Ssec_score / mean(data$Ssec_score)
    keep <- gene_symbol(data$gene_tissue) == candidate$gene_symbol
    data.frame(comparison = rep(label, sum(keep)), gene_tissue = data$gene_tissue[keep],
               normalized_Ssec = normalized[keep])
  }))
  print(ggplot(comparison, aes(comparison, normalized_Ssec, fill = comparison)) +
          geom_col(show.legend = FALSE) +
          labs(title = paste("Within- and cross-tissue scores for", candidate$gene_symbol),
               x = NULL, y = "Normalized Ssec score") +
          theme_classic())
}
invisible(dev.off())
print(sessionInfo())
