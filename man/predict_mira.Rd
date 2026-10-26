\name{predict_mira}
\alias{predict_mira}
\alias{predict_MIRA}
\alias{predict_MIRA_v2}
\alias{predict_MIRA_final}
\title{Apply the Frozen Two-Step MIRA Classifier}
\usage{
predict_mira(expression, batch_id = NULL, model = NULL,
  decision_rule = c("two_step", "legacy_joint_argmax"))
predict_MIRA(expression, batch_id = NULL, model = NULL,
  decision_rule = c("two_step", "legacy_joint_argmax"))
predict_MIRA_v2(expression, batch_id = NULL, model = NULL,
  decision_rule = c("two_step", "legacy_joint_argmax"))
predict_MIRA_final(expression, batch_id = NULL, model = NULL,
  decision_rule = c("two_step", "legacy_joint_argmax"))
}
\arguments{
\item{expression}{Numeric genes-by-samples matrix with unique gene symbols and sample IDs. All 18,686 fixed background genes are required. Additional genes are ignored.}
\item{batch_id}{Study or assay batch identifier, either one value for all columns or one value per sample. A named vector is matched to the sample IDs. Each batch must contain at least two samples.}
\item{model}{Optional compact model list with the same frozen specification; normally leave \code{NULL} to use the packaged model. No training is performed.}
\item{decision_rule}{The default \code{"two_step"} applies sequential 0.5 thresholds. \code{"legacy_joint_argmax"} explicitly reproduces the historical joint-score maximum rule and may give different labels.}
}
\description{Predicts C1a, C1b or C2 using frozen logistic coefficients and target-batch calibration. All exported aliases use the same final default rule.}
\details{
Expression is ranked within each sample against the fixed background. Eight separate MHC-II gene ranks enter the first component; the mean rank of eight B-lineage markers enters the second component. Within each complete batch, each feature is centered and divided by its sample standard deviation. Standard deviations below \code{1e-10} are replaced by one. Frozen component-specific training transformations and coefficients are then applied.

If \eqn{P(C1)<0.5}, the label is C2. Otherwise, \eqn{P(C1a|C1)\ge0.5} gives C1a and the remaining cases give C1b. Do not replace this rule with the maximum joint score.

Provide the entire eligible study or assay population defined without subtype or outcome selection. Both components use calibration from this entire batch, including samples subsequently assigned C2. Changing batch composition can change predictions. Two samples is a computational floor, not a validated effective minimum; isolated single-sample use is not validated. Missing background genes, duplicate IDs and nonfinite required expression values cause errors. The function does not convert raw counts to TPM or map gene names.
}
\value{A data frame containing sample IDs, predicted labels, the component probabilities \code{p_C1} and \code{p_C1a_given_C1}, derived joint scores, the rule and batch metadata. The \code{unlabeled_batch_calibration} attribute records means, standard deviations and batch sizes. Probabilities are model scores, not validated clinical confidence or treatment-benefit estimates.}
\examples{
example <- mira_example()
prediction <- predict_mira(example$expression, example$batch_id)
prediction[, c("sample_id", "predicted_label", "p_C1", "p_C1a_given_C1")]
}
\seealso{\code{\link{mira_model}}, \code{\link{mira_example}}, \code{\link{read_mira_expression}}}

