% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{ves_characteristics_bin}
\alias{ves_characteristics_bin}
\title{Distance bins of vessel characteristics}
\usage{
ves_characteristics_bin(vsf, binsize = 200, diam = "Dcircle")
}
\arguments{
\item{vsf}{\code{sf} object of vessel polygons}

\item{binsize}{Size of distance bins to away from cambium.}
}
\value{
A data.frame with characteristics of vessels included in each distance bin.
Includes vessel grouping indices, mean and median vessel area, kurtosis, skewness,
hydraulically-weighted diameter (Dh), total theoretical conductivity (Kh_total),
and mean theoretical conductivity (Kh_mean).
}
\description{
Get characteristics of vessels by binned distance from cambium.
}
