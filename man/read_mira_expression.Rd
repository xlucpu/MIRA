\name{read_mira_expression}
\alias{read_mira_expression}
\title{Read a Gene Expression Matrix}
\usage{read_mira_expression(path)}
\arguments{\item{path}{Path to an RDS matrix or data frame, or a tab-delimited text file. Compressed TSV files are supported.}}
\description{Reads expression for subsequent validation and prediction. In TSV files, the first column contains gene symbols and remaining named columns contain numeric sample expression.}
\value{A genes-by-samples numeric matrix. Gene coverage and identifiers are checked when passed to \code{predict_mira}.}
\seealso{\code{\link{predict_mira}}}

