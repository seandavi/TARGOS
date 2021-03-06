% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/TargetOsteo.R
\name{getmiRNA}
\alias{getmiRNA}
\title{Get TARGET miRNA data (only 89 discovery samples) as an ExpressionSet}
\usage{
getmiRNA(DCC = options("DCC"))
}
\arguments{
\item{DCC}{The string location of the DCC directory}
}
\description{
This function loads the miRNA data and does sample matching with
the clinical data (loaded using \code{\link{getClinical()}}). The
actual file name used is "miRNA dCt values.txt". No further normalization
is performed.
}
\examples{
mirna = getmiRNA()
mirna
}

