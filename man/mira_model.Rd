\name{mira_model}
\alias{mira_model}
\title{Load the Packaged Frozen Numeric Model}
\usage{mira_model()}
\description{Loads and validates the packaged fixed coefficients, training transformations and gene definitions.}
\value{A compact model list. It contains numeric parameters and aggregate training information, with no patient-level expression, training labels or original fitted objects.}
\examples{length(mira_model()$gene_universe)}
\seealso{\code{\link{predict_mira}}}

