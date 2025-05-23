% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unknownpleasur.R
\name{up_elevate}
\alias{up_elevate}
\title{Elevate Lines Based on Raster Values}
\usage{
up_elevate(
  x,
  interval,
  raster,
  angle = 90,
  scale = "interval",
  factor = 1,
  mode = "xyz"
)
}
\arguments{
\item{x}{\code{sf} object containing \code{LINESTRING}s.}

\item{interval}{Line spacing. Can be a \code{units} object.}

\item{raster}{\code{terra} \code{Spatraster} from which to extract elevations.}

\item{angle}{Angle of lines in degrees.}

\item{scale}{One of "actual" or "interval." If "actual", use actual raster values. If "interval," base scale on interval between lines.}

\item{factor}{Z exaggeration.}

\item{mode}{One of \code{"xyz"}, or \code{"planar"}. If \code{"xyz"} returns linestrings
with Z component. If \code{"planar"} returns lines offset by value.}
}
\value{
\code{sf} object containing \code{LINESTRING}s.
}
\description{
Elevate Lines Based on Raster Values
}
