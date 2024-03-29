% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\name{getRatio}
\alias{getRatio}
\alias{getRatio,Compass,character-method}
\title{Accessor of ratio}
\usage{
getRatio(object, slot = c("sample", "spike_in"), filter = FALSE)

\S4method{getRatio}{Compass,character}(object, slot = c("sample", "spike_in"), filter = FALSE)
}
\arguments{
\item{object}{Compass.}

\item{slot}{Which slot to get, one of \code{sample} or \code{spike_in}.}

\item{filter}{Whether to get the filtered ratio, default FALSE.}
}
\value{
list of ratio table
}
\description{
Get all or filtered ratio
}
