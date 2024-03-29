% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{calcNCIN}
\alias{calcNCIN}
\title{Compute NCIN ratio}
\usage{
calcNCIN(
  data,
  spike.in.prefix = NULL,
  enrich.group,
  input.id = "Input",
  enrich.id = "Enrich",
  scale.factor,
  adjust.factor,
  pseudo.count = 1
)
}
\arguments{
\item{data}{A un-normalized count data matrix of shape n x p, where n is the
number of samples and p is the number of features.}

\item{spike.in.prefix}{A character specify the prefix of spike-in id, e.g.,
"^FB" stands for fly spike-in id, default: NULL.}

\item{enrich.group}{Vector of enrichment group, e.g., c("Input","Enrich","Input","Enrich").}

\item{input.id}{Input library id, must be consistent with \code{enrich.group}, e.g., "Input".}

\item{enrich.id}{Enrich library id, must be consistent with \code{enrich.group}, e.g., "Enrich".}

\item{scale.factor}{Vector of scale factors generated from \code{calcScaleFactor}.}

\item{adjust.factor}{Vector of adjust factors generated from \code{calcAdjustFactor}.}

\item{pseudo.count}{A numeric scalar of pseudo-counts to be added to each gene, default: 1.}
}
\value{
List containing the unfiltered and filtered ratio of each gene from each sample.
}
\description{
Compute NCIN ratio
}
