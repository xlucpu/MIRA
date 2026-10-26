\name{mira_example}
\alias{mira_example}
\title{Load a Synthetic Expression Example}
\usage{mira_example()}
\description{Provides an artificial example for testing installation and learning the prediction interface. It is not patient data or clinical validation data.}
\value{A list with a numeric \code{expression} matrix, \code{batch_id} and a data frame of \code{expected} predictions. The matrix has 18,686 gene rows and 12 synthetic sample columns.}
\examples{
example <- mira_example()
dim(example$expression)
example$expected[, c("sample_id", "predicted_label")]
}
\seealso{\code{\link{predict_mira}}}

