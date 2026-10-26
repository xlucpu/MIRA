# MIRA

An R package for the frozen two-step classification of renal cell carcinoma expression profiles into C1a, C1b and C2. The package includes the fixed model, prediction functions, an artificial example and usage documentation. Prediction uses base R and its standard stats/utils packages only; `glmnet` is not required.

## Install and run

```r
install.packages("remotes") # only if needed
remotes::install_github("xlucpu/MIRA")

library(MIRA)
example <- mira_example()
prediction <- predict_mira(example$expression, batch_id = example$batch_id)
prediction[, c("sample_id", "predicted_label", "p_C1", "p_C1a_given_C1")]
```

The included example is entirely synthetic. Its expected labels and probabilities are available as `example$expected`.

For your own data:

```r
expression <- read_mira_expression("expression.tsv")
prediction <- predict_mira(expression, batch_id = "my_study")
write.table(prediction, "predictions.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(attr(prediction, "unlabeled_batch_calibration"), "batch_calibration.rds")
```

`predict_MIRA()`, `predict_MIRA_v2()` and `predict_MIRA_final()` are aliases with the same default rule. Help is available with `?predict_mira`.

## Input

- A numeric **gene-by-sample** expression matrix, with unique gene symbols and sample IDs.
- Every gene in the fixed **18,686-gene background** must be present. Extra genes are ignored. Obtain the background with `mira_model()$gene_universe`.
- Missing background genes, duplicate identifiers and nonfinite required expression values are rejected. No silent gene mapping or imputation is performed.
- Supply processed, comparable gene-level expression within each assay batch. The package does not convert raw counts to TPM.
- A marker-only measurement of the 16 feature genes is insufficient: ranks require the full background. Do not supply MCP-counter abundance scores.

## Frozen model and rule

The first component uses individual within-sample percentile ranks of eight MHC-II genes. The second uses the arithmetic mean of the ranks of eight B-lineage genes. `mira_model()` exposes the exact genes, coefficients and transformations.

Within each sample, ranks are `(rank(expression, ties.method="average") - 0.5) / 18686`. Features are centered and scaled by the mean and sample SD of the complete unlabeled target batch, followed by frozen component-specific training transformations. A feature SD below `1e-10` is replaced by one.

1. `P(C1) < 0.5` assigns **C2**.
2. Otherwise, `P(C1a | C1) >= 0.5` assigns **C1a**; the remaining cases are **C1b**.

Both original elastic-net logistic components used `alpha=0.5` and `lambda=0.003`. The final coefficients were fitted using 539 development samples for the first step and 248 original C1 samples for the second. The package does not train or tune models.

The output includes both component probabilities and derived joint scores. The final label is determined by the sequential thresholds, **not** the largest joint score. An explicit `decision_rule="legacy_joint_argmax"` option is retained for historical comparisons.

## Batch requirements

Provide the complete eligible study/assay population, defined without selecting samples by original/predicted subtype or outcome. Both components are calibrated using this full batch, including samples subsequently assigned C2. Do not recalculate calibration in the predicted C1 subset.

For multiple assay batches, provide one batch ID per expression column. A named vector is matched to sample IDs:

```r
batch <- setNames(sample_annotation$assay_batch, sample_annotation$sample_id)
prediction <- predict_mira(expression, batch_id = batch)
```

Changing batch composition can change predictions. The minimum of two samples per batch is a computational requirement, not a validated effective minimum. Isolated single-sample application has not been validated. Model scores are not validated clinical confidence or treatment-benefit estimates.

## Verification

The packaged implementation reproduces the original fixed model without refitting. Separate retrospective leave-one-cohort-out models reproduce the 539 predictions underlying Figure S11; the distributed model is the final pooled fit, not a substitute for those held-out predictions. Aggregate results are included in `inst/extdata/S11_reproduction.tsv` and `pooled_model_equivalence.tsv`.

| Held-out cohort | Three-class balanced accuracy | Second-step AUC within original C1 |
|---|---:|---:|
| BIONIKK | 0.758 | 0.625 |
| JAVELIN | 0.744 | 0.837 |
| CheckMate | 0.808 | 0.929 |

Performance estimates reflect retrospective agreement with the original cohort-specific labels. No patient-level measurements, original training objects or labels are distributed.

## Development and license

This package is generated from the maintainer's unified `revision.R` analysis source; the distributed source contains prediction code only. The package tests use synthetic data. License: MIT.

