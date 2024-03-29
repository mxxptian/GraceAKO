% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhq_threshold.R
\name{bhq_threshold}
\alias{bhq_threshold}
\title{Function for applying the Standard Benjamini-Hochberg for controlling False discovery rate}
\usage{
bhq_threshold(pvals, fdr = 0.1)
}
\arguments{
\item{pvals}{a vector of the feature statistics, AGCs}

\item{fdr}{the pre-specified false discovery rate level}
}
\value{
The threshold value corresponding to the Standard Benjamini-Hochberg procedure
}
\description{
Function for applying the Standard Benjamini-Hochberg for controlling False discovery rate
}
