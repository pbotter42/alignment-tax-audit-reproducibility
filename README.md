# Alignment Tax Reproducibility Package

This package contains the minimal code, data, manuscript source, and generated outputs needed to reproduce the archived `stat.ML` submission.

## Contents

- `code/run_alignment_tax_audit.R`
  Main analysis script for the item-level audit.
- `data/hugging_item_level_data_final_april8_2024.csv`
  Item-level Open LLM Leaderboard snapshot used in the paper.
- `data/test_level_data_july23_2024_Levenshtein_20_N_641_factors.csv`
  Model filter and metadata table used to restrict the audit sample.
- `data/jackknife_results_checkpoint.csv`
  Precomputed leave-one-out structural impact values used in the manuscript workflow.
- `outputs/`
  Generated figures and summary tables used by the manuscript.
- `manuscript/alignment_tax_archive.tex`
  Submission source.
- `manuscript/references.bib`
  Bibliography source.
- `manuscript/llm_article_arxiv_v1.pdf`
  Compiled submission PDF.

## Quick Run

From this folder:

```bash
Rscript code/run_alignment_tax_audit.R
```

That writes the regenerated figures and summary tables into `outputs/`.

## Notes

- The default workflow uses the checkpointed jackknife file because recomputing all leave-one-out PCA fits is much slower.
- To force a full jackknife recomputation, run:

```bash
Rscript code/run_alignment_tax_audit.R --recompute_jackknife=true
```

- The manuscript currently references the figure files in `outputs/`.
