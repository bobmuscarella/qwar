% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{cam_dist}
\alias{cam_dist}
\title{Measure distance along rays to cambium}
\usage{
cam_dist(csf, rsf, vsf, parallel = F, cores = detectCores(logical = FALSE))
}
\arguments{
\item{csf}{\code{sf} object of cambium line}

\item{rsf}{\code{sf} object of ray line(s)}

\item{vsf}{\code{sf} object of vessel polygons}

\item{parallel}{Logical if processing should be done in parallel (not yet implemented).}

\item{cores}{Number of cores to use for parallel processing.}
}
\value{
A vector of distances (in microns) from each feature to the
cambium line, along the closest ray.
}
\description{
Measure distance along rays to cambium line of features.
}
