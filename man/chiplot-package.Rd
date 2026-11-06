\name{chiplot-package}
\alias{chiplot-package}
\docType{package}
\title{Chi-Plots for Assessing Bivariate and Multivariate Dependence}
\description{
  The chi-plot (Fisher and Switzer, 1985, 2001) is a graphical aid to
  interpreting a scatterplot of two variables: under independence the points
  form a horizontal band, and different forms of association produce
  characteristically different patterns. Fisher and Lee (2026) extend the idea
  to multivariate data, to explore whether a multivariate distribution factorizes
  into independent blocks of variables.
}
\details{
  The main functions are:
  \describe{
    \item{\code{\link{chiplot}}}{a chi-plot for two variables;}
    \item{\code{\link{chipairs}}}{a scatterplot matrix with the corresponding
      chi-plots above the diagonal;}
    \item{\code{\link{chimult}}}{chi- and lambda-values, permutation tests and the
      multivariate chi-plot for a given partition of the variables into blocks;}
    \item{\code{\link{chipart}}}{a search over all partitions of the variables,
      keeping those that show a significant departure from independence;}
    \item{\code{\link{chidepth}}}{the data depth functions used for the
      lambda-values of multivariate chi-plots.}
  }
  The calculations of chi-values are carried out in compiled C code. All other
  code uses base \R only.
}
\author{
  Nicholas I. Fisher and Alan J. Lee.
  Maintainer: Nicholas I. Fisher <nif@valuemetrics.com.au>
}
\references{
  Fisher, N. I. and Switzer, P. (1985). Chi-plots for assessing dependence.
  \emph{Biometrika}, \bold{72}, 253--265. \doi{10.1093/biomet/72.2.253}

  Fisher, N. I. and Switzer, P. (2001). Graphical assessment of dependence: is a
  picture worth a hundred tests? \emph{The American Statistician}, \bold{55},
  233--239. \doi{10.1198/000313001317098248}

  Fisher, N. I. and Lee, A. J. (2026). Chi-plots for assessing multivariate
  dependence. Submitted to \emph{Journal of Computational and Graphical
  Statistics}.
}
\keyword{package}
\keyword{hplot}
