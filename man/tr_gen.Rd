% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{tr_gen}
\alias{tr_gen}
\title{Creates a list containing the trajectories information.}
\usage{
tr_gen(data, projection = CRS("+proj=longlat"), origin = "1970-01-01 UTC")
}
\arguments{
\item{projection}{Projection string of class CRS-class.}

\item{origin}{Optional. Origin of the date-time. Only needed in the internal process to create an object of type POSIXct.}
}
\value{
Every component of the list has the following subcomponents:
#' \itemize{
  \item An object of class Track.
  \item An object of class SpatialPoints.
  \item Subject identifier
  \item Trip identifier.
}
}
\description{
Creates a list containing the trajectories information.
}
\keyword{internal}
