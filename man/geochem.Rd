\name{geochem}
\alias{geochem}
\docType{data}
\title{Geochemical Measurements}
\description{
  Measurements of 11 geochemical constituents on 1000 samples, used to
  illustrate \code{\link{chipairs}}. The values are the ranks of the
  measurements within each variable, with tied values given their average rank.
}
\usage{geochem}
\format{
  A data frame with 1000 rows and 11 numeric columns: \code{TiO2},
  \code{Cr2O3}, \code{FeO}, \code{MnO}, \code{MgO}, \code{CaO}, \code{Ni},
  \code{Zn}, \code{Ga}, \code{Y} and \code{Zr}.
}
\source{
  ValueMetrics Australia. Included with the permission of the copyright holder.
}
\examples{
summary(geochem[, 1:4])
chipairs(geochem[, c("Cr2O3", "MgO", "Ni")])
}
\keyword{datasets}
