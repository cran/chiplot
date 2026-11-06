\name{chiplot}
\alias{chiplot}
\alias{chivals}
\alias{plot.chiplot}
\alias{print.chiplot}
\title{Chi-Plot for Assessing Dependence of Two Variables}
\description{
  \code{chivals} calculates the chi- and lambda-values of Fisher and Switzer
  (1985) for a bivariate sample; \code{chiplot} calculates them and draws the
  chi-plot.
}
\usage{
chiplot(x, y = NULL, ctrl.lines = 2, jitter = NULL, jitter.amount = 0.25,
        cex.points = NULL, pch.type = 19, dot.col = "black", main = "",
        cex.axis = 1, cex.lam = 1.2, plot = TRUE, \dots)

chivals(x, y = NULL, jitter = NULL, jitter.amount = 0.25)

\method{plot}{chiplot}(x, ctrl.lines = 2, cex.points = NULL, pch.type = 19,
     dot.col = "black", main = "", cex.axis = 1, cex.lam = 1.2,
     xlab = expression(lambda), ylab = expression(chi),
     col.ctrl = "red", lty.ctrl = 1, \dots)

\method{print}{chiplot}(x, \dots)
}
\arguments{
  \item{x}{a numeric vector or ordered factor; or, if \code{y} is \code{NULL}, a
    matrix or data frame whose first two columns are used. For the methods, an
    object of class \code{"chiplot"}.}
  \item{y}{a numeric vector or ordered factor of the same length as \code{x}.}
  \item{ctrl.lines}{control lines enclosing, on average, the given percentage of
    points under independence: \code{0} (none), \code{1} (90\%), \code{2} (95\%)
    or \code{3} (99\%). The values \code{0.90}, \code{0.95} and \code{0.99} are
    also accepted.}
  \item{jitter}{which variables to jitter, to break ties (for example, for
    ordered categorical variables): \code{TRUE} (both), \code{FALSE} (neither),
    a logical vector, or column numbers or names. The default \code{NULL}
    jitters ordered factors only.}
  \item{jitter.amount}{size of the jitter, as a fraction (less than 0.5) of the
    smallest gap between distinct values, so that the order of distinct values
    is preserved.}
  \item{cex.points}{size of the plotted points. The default \code{NULL} chooses
    a size, as large as possible without crowding, from the number of points.}
  \item{pch.type, dot.col}{plotting symbol and colour of the points.}
  \item{main}{title of the plot.}
  \item{cex.axis}{size of the axis annotation.}
  \item{cex.lam}{size of the axis labels \eqn{\lambda} and \eqn{\chi}.}
  \item{plot}{logical: if \code{FALSE} the values are calculated but not plotted.}
  \item{xlab, ylab}{axis labels.}
  \item{col.ctrl, lty.ctrl}{colour and line type of the control lines.}
  \item{\dots}{further graphical parameters passed to \code{\link{plot}}.}
}
\details{
  Observations with a missing value in either variable are removed. For each
  point \eqn{(x_i, y_i)}, let \eqn{F_i}, \eqn{G_i} and \eqn{H_i} be the
  proportions of the other \eqn{n-1} points with \eqn{x_j \le x_i},
  \eqn{y_j \le y_i}, and both, respectively. Then
  \deqn{\chi_i = (H_i - F_i G_i) / \{F_i (1 - F_i) G_i (1 - G_i)\}^{1/2},}{chi_i = (H_i - F_i G_i) / sqrt(F_i (1 - F_i) G_i (1 - G_i)),}
  \deqn{\lambda_i = 4\,\mathrm{sign}(\tilde F_i \tilde G_i)\,
        \max(\tilde F_i^2, \tilde G_i^2),}{lambda_i = 4 sign(F~_i G~_i) max(F~_i^2, G~_i^2),}
  where \eqn{\tilde F_i = F_i - 1/2}{F~_i = F_i - 1/2} and
  \eqn{\tilde G_i = G_i - 1/2}{G~_i = G_i - 1/2}. Points with
  \eqn{|\lambda_i| \ge 4\{1/(n-1) - 1/2\}^2}, at the edges of the data, are not
  plotted. Under independence the points lie in a horizontal band about
  \eqn{\chi = 0}, approximately the stated percentage of them between the
  control lines at \eqn{\pm c_p / \sqrt{n}}{+/- c_p / sqrt(n)}, with
  \eqn{c_p} = 1.54, 1.78 and 2.18 for 90\%, 95\% and 99\%.

  As in the original program, the sine transformation of the chi-values
  mentioned by Fisher and Switzer (1985) is not used. The chi- and
  lambda-values depend only on the ranks of the data.
}
\value{
  An object of class \code{"chiplot"} (returned invisibly by \code{chiplot}): a
  list with components
  \item{lambda, chi}{the lambda- and chi-values of the plotted points;}
  \item{index}{the row numbers, in the original data, of the plotted points;}
  \item{n}{the number of complete pairs;}
  \item{npts}{the number of points plotted;}
  \item{n.missing}{the number of pairs removed because of missing values;}
  \item{names}{the names of the two variables.}
}
\references{
  Fisher, N. I. and Switzer, P. (1985). Chi-plots for assessing dependence.
  \emph{Biometrika}, \bold{72}, 253--265. \doi{10.1093/biomet/72.2.253}

  Fisher, N. I. and Switzer, P. (2001). Graphical assessment of dependence: is a
  picture worth a hundred tests? \emph{The American Statistician}, \bold{55},
  233--239. \doi{10.1198/000313001317098248}
}
\seealso{\code{\link{chipairs}} for all pairs of variables;
  \code{\link{chimult}} for blocks of variables.}
\examples{
set.seed(1)
op <- par(mfrow = c(2, 2))
## Independent normal variables
x <- rnorm(100); y <- rnorm(100)
chiplot(x, y, main = "Independence")
## Bivariate normal, correlation 0.7
y <- 0.7 * x + sqrt(1 - 0.7^2) * rnorm(100)
chiplot(x, y, main = "Correlation 0.7")
## Non-monotone association: points near a half circle
x <- rnorm(2000); y <- rnorm(2000)
keep <- y >= 0 & abs(x^2 + y^2 - 1) <= 0.5
chiplot(x[keep][1:100], y[keep][1:100], main = "Half ring")
## Ordered categorical variable: jittered automatically
u <- rnorm(100)
grade <- cut(u + rnorm(100), c(-Inf, -1, 0, 1, Inf), ordered_result = TRUE)
res <- chiplot(u, grade, main = "Ordered categories")
par(op)
res
}
\keyword{hplot}
\keyword{nonparametric}
