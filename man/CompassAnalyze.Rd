% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{CompassAnalyze}
\alias{CompassAnalyze}
\title{Normalization and ratio assessment}
\usage{
CompassAnalyze(
  object,
  adjust = TRUE,
  prop.top.enrich = 1,
  decreasing = TRUE,
  pseudo.count = 1,
  enone.ruv.factor = NULL,
  ratio.shrinkage = TRUE
)
}
\arguments{
\item{object}{Compass object}

\item{adjust}{Whether to perform linear regression-based adjustment, default: TRUE}

\item{prop.top.enrich}{Proportion of top-enriched genes to use for
adjustment of non-specific enrichment, by default all genes are used.}

\item{decreasing}{Whether using decreasing enriched genes for
the adjustment of non-specific enrichment, default: TRUE (decreasing)}

\item{pseudo.count}{A numeric scalar of pseudo-counts to be added to each gene, default: 1.}

\item{enone.ruv.factor}{List of RUV factor from \code{enONE}, only normalization
factor generated by RUV method can be passed.}

\item{ratio.shrinkage}{Whether to perform shrinkage of ratio, default TRUE.}
}
\value{
Compass object
}
\description{
Normalization and ratio assessment
}
