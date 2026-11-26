# Pathway-specific targeting of pancreatic function

These R scripts use filtered GTEx v8 expression data to prioritize candidate regulators of pancreatic processes. The original analysis ranks genes from 13 source tissues by how strongly their expression correlates across donors with insulin-secretion transcripts in bulk pancreas. The companion scripts apply the same scoring approach to pancreatic transcript panels for HGF, FGF21, and SERPINB1.

## Datasets

Download the analysis inputs from the [Google Drive dataset folder](https://drive.google.com/drive/folders/1ESiTRAQA--LecraBZzAVdWK7iKkliMGF?usp=sharing). Put the four files in `datasets/` and keep their filenames, or pass the path to your dataset directory when running a script.

| File | Role |
| --- | --- |
| `GTEx NA included env.RData` | Filtered GTEx v8 workspace containing `GTEx_subfiltered`. The analysis uses 310 donor rows after transposition. Missing observations are retained, so matched donor counts vary by tissue pair. |
| `GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt` | GTEx subject annotations. The original script reports donor sex counts for context and does not use sex as a filter or covariate. |
| `uniprot-human-genes and goterms mapping.tab` | Historical UniProt mapping of human primary gene symbols to Gene Ontology (GO) terms. Used for target selection and enrichment. |
| `uniprot-secreted-filtered-organism__Homo+sapiens+(Human)+[9606]_.tab` | Historical UniProt list used to annotate candidates and generate the secreted-only ranking. |

## Scripts and pancreatic target panels

| Script | Target panel and result |
| --- | --- |
| `script for crosstissue correlation with pancreas pathways.R` | Selects pancreatic genes annotated to GO biological processes containing insulin secretion. Pituitary HPR ranked **#3 of 35,737** secreted geneâ€“tissue pairs overall and **#1** among secreted pituitary candidates. |
| `prioritize_HGF.R` | By focusing on positive regulation of insulin secretion and pancreatic B-cell proliferation, HGF emerges (top 1%) as a top-ranked candidate |
| `prioritize_FGF21.R` | Tests fatty-acid response and negative regulation of intrinsic apoptotic signaling, two processes relevant to B-cell survival during lipid stress. The apoptosis-control panel placed visceral adipose FGF21 is prioritized in top 1%. The original insulin-secretion screen ranked liver FGF21 slightly higher. |
| `prioritize_SERPINB1.R` | Tests pancreatic B-cell proliferation and negative regulation of endopeptidase activity. The 18-gene proliferation panel placed liver SERPINB1 in the top 1% using all donors. |

Keep `reverse_pancreas.R` and `pathway_prior_screen.R` alongside the three `prioritize_*.R` scripts. The first provides the donor-split workflow and GO overrepresentation analysis. The second scores fixed, literature-motivated GO panels using all available donors. Each protein script documents the exact GO terms and the supporting experimental studies.

## Run the analyses

Install the R packages `WGCNA`, `ggplot2`, and `forcats`. The verified GTEx runs used R 4.1.2 and WGCNA 1.71.

The commands below assume the annotated original script is at the repository root and the companion scripts are in `protein_scripts/`. Adjust the directory paths to match your layout. Use a fresh output directory for each companion run.

```bash
Rscript "script for crosstissue correlation with pancreas pathways.R" datasets results/hpr
Rscript protein_scripts/prioritize_HGF.R datasets results/hgf
Rscript protein_scripts/prioritize_HGF.R datasets results/hgf_go go
Rscript protein_scripts/prioritize_FGF21.R datasets results/fgf21
Rscript protein_scripts/prioritize_SERPINB1.R datasets results/serpinb1
```

The original script writes `all_scores.tsv`, `secreted_scores.tsv`, and `plots.pdf`. The companion scripts write the pancreas panel membership in `panels.tsv`, complete per-panel rankings, and an `anchor_ranks.tsv` summary. The donor-split workflow also exports discovery and held-out correlations and GO enrichment tables. The fixed GO workflow records the exact pathway definitions in `panel_definitions.tsv`.

## Scoring and interpretation

For each source geneâ€“tissue pair, the scripts calculate biweight midcorrelations against the target transcripts in pancreas. Ssec is the mean of âˆ’logâ‚â‚€(p) across the target panel. The companion workflows require at least 25 paired donor observations and detectable expression for both transcripts in at least 10 paired donors. The HGF donor-split result uses donors reserved for validation. The fixed GO results are descriptive rankings from the full cohort.

The rankings summarize coexpression with bulk pancreas. They can prioritize candidates for follow-up, but do not establish direct effects of circulating proteins or identify the responding cell type. The secreted rankings use the supplied historical UniProt snapshot. HGF is absent from that list, so its reported rank uses all candidates. To explore other organs or processes, update the target tissue and target-gene selection for the relevant workflow.

## References

- Velez et al. [Genetic variation of putative myokine signaling is dominated by biological sex and sex hormones](https://elifesciences.org/articles/76887).
- Zhou, Tamburini, Van et al. [Leveraging inter-individual transcriptional correlation structure to infer discrete signaling mechanisms across metabolic tissues](https://elifesciences.org/articles/88863).

For questions about the data or analysis, contact [mseldin@uci.edu](mailto:mseldin@uci.edu).


If any datasets or insight would help, please contact mseldin@uci.edu

All datasets used for the analysis proposed are provided here:
https://drive.google.com/drive/folders/1ESiTRAQA--LecraBZzAVdWK7iKkliMGF?usp=sharing
