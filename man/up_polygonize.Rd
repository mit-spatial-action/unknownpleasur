% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unknownpleasur.R
\name{up_polygonize}
\alias{up_polygonize}
\title{Section Lines to Polygons}
\usage{
up_polygonize(x, baselines, id_col, mode = "planar")
}
\arguments{
\item{x}{\code{sf} object containing sectional lines.}

\item{baselines}{\code{sf} object containing baselines.}

\item{id_col}{Name of id column.}

\item{mode}{If \code{planar}, results will be \code{"POLYGON"}s. If \code{xyz}, results will be `"LINESTRINGS"``.}
}
\value{
An \code{sf} object.
}
\description{
Section Lines to Polygons
}
