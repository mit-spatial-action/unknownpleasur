% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unknownpleasur.R
\name{up_unknown_pleasures}
\alias{up_unknown_pleasures}
\title{Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.}
\usage{
up_unknown_pleasures(
  x,
  raster,
  interval,
  line_angle,
  scale,
  factor,
  elev_angle = line_angle,
  mode = "planar",
  mask = TRUE,
  polygon = TRUE
)
}
\arguments{
\item{x}{\code{sf} object that will be treated as an extent..}

\item{raster}{\code{raster} object.}

\item{interval}{Line and sample spacing. Can be \code{units} object. Otherwise will be in unit of CRS.}

\item{line_angle}{Angle (in degrees) at which to create lines.}

\item{scale}{One of "actual" or "interval." If "actual", use actual raster values. If "interval," base scale on interval between lines.}

\item{factor}{Z exaggeration, essentially.}

\item{elev_angle}{Angle at which to offset points.}

\item{mode}{If \code{planar}, results will be planar offset lines. If \code{xyz}, lines will be offset on \code{LINESTRING} z axis.}

\item{mask}{If \code{true}, output lines will be clipped to edge of input \code{x}.}

\item{polygon}{If \code{TRUE}, outputs polygons (or closed linestrings if paired with \code{mode = "xyz"}. If \code{FALSE}, outputs lines with no baseline.}
}
\value{
\code{sf} object  containing "Unknown Pleasures" features.
}
\description{
Transforms regularly-spaced lines into an "Unknown Pleasures"-esque set of regular section cuts based on raster value.
}
