\name{chidepth}
\alias{chidepth}
\title{Data Depth for Multivariate Chi-Plots}
\description{
  Calculates the data depth of each row of a data matrix with respect to the
  whole sample. These depths are the basis of the lambda-values in
  \code{\link{chimult}}.
}
\usage{
chidepth(X, depth = c("RP", "TD", "PD", "MhD"), Nsim = 100)
}
\arguments{
  \item{X}{a numeric matrix without missing values.}
  \item{depth}{the depth function: one of the names below, or a function.}
  \item{Nsim}{the number of random directions for \code{"RP"} and \code{"PD"}.}
}
\details{
  All the depth functions are implemented in base \R:
  \describe{
    \item{\code{"RP"}}{random projection depth: the average, over \code{Nsim}
      random directions, of the univariate halfspace depth of the projected
      data. The directions are generated as in \code{fda.usc::mdepth.RP}, so the
      results are identical for the same random number seed. This is the
      default.}
    \item{\code{"TD"}}{the average of the univariate halfspace (Tukey) depths of
      the coordinates, as \code{fda.usc::mdepth.TD}. Note that this is not the
      exact multivariate Tukey depth.}
    \item{\code{"PD"}}{projection depth (Zuo and Serfling, 2000),
      \eqn{1 / (1 + O(x))} with outlyingness
      \eqn{O(x) = \max_u |u^T x - \mathrm{med}(u^T X)| / \mathrm{MAD}(u^T X)}{O(x) = max_u |u'x - med(u'X)| / MAD(u'X)},
      the maximum being taken over \code{Nsim} random directions.}
    \item{\code{"MhD"}}{Mahalanobis depth, \eqn{1 / (1 + d^2)} with \eqn{d} the
      Mahalanobis distance from the mean, as \code{fda.usc::mdepth.MhD}. This is
      very fast but parametric.}
  }
  The names \code{"mdepth.RP"}, \code{"mdepth.TD"}, \code{"mdepth.MhD"} and
  \code{"Projection"} are also accepted.

  A user-specified \code{depth} must be a function of the data matrix that
  returns one depth per row, either as a numeric vector or as a list with
  component \code{dep}. For example, \code{function(X) fda.usc::mdepth.RP(X,
  proj = 500)$dep}, or a function calling \code{depth()} in package
  \pkg{DepthProc}.
}
\value{A numeric vector of depths, one for each row of \code{X}.}
\references{
  Zuo, Y. and Serfling, R. (2000). General notions of statistical depth function.
  \emph{The Annals of Statistics}, \bold{28}, 461--482.
}
\seealso{\code{\link{chimult}}}
\examples{
set.seed(1)
X <- matrix(rnorm(200), 100, 2)
d <- cbind(RP = chidepth(X, "RP"), TD = chidepth(X, "TD"),
           PD = chidepth(X, "PD"), MhD = chidepth(X, "MhD"))
round(cor(d, method = "spearman"), 2)

## a user-supplied depth: inverse distance from the coordinatewise median
chidepth(X[1:5, ], function(X) 1 / (1 + sqrt(rowSums(t(t(X) - apply(X, 2, median))^2))))
}
\keyword{multivariate}
