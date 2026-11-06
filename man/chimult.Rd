\name{chimult}
\alias{chimult}
\alias{plot.chimult}
\alias{print.chimult}
\title{Multivariate Chi-Plot for a Partition of the Variables}
\description{
  For a multivariate sample and a partition of the variables into blocks,
  \code{chimult} calculates the chi- and lambda-values of Fisher and Lee (2026),
  permutation tests of independence of the blocks, and the material for the
  multivariate chi-plot, which is drawn by \code{plot}.
}
\usage{
chimult(Data, D, nperm = 100, Nsim = 100, depth = "RP", jitter = NULL,
        jitter.amount = 0.25, plot = FALSE, verbose = TRUE, \dots)

\method{plot}{chimult}(x, pc.level = 0.95, main = "", col1 = "red",
     col2 = "blue", pt = 12, box = TRUE, cex = NULL, mar = c(4, 4, 2, 2),
     guide.option = 1, del = 0.2, \dots)

\method{print}{chimult}(x, \dots)
}
\arguments{
  \item{Data}{a numeric matrix or data frame (\eqn{n \times p}{n x p}), each row
    an observation of a \eqn{p}-variate random variable. Ordered factors are
    allowed.}
  \item{D}{the block sizes, in column order: for example, if \eqn{p = 5},
    \code{D = c(2, 3)} compares the first two columns with the last three. At
    least two blocks are needed. Reorder the columns of \code{Data} first for
    other partitions.}
  \item{nperm}{the number of random permutations used for the p-values and
    null distributions.}
  \item{Nsim}{the number of random directions for the \code{"RP"} and
    \code{"PD"} data depths.}
  \item{depth}{the data depth used for the lambda-values: \code{"RP"} (random
    projection depth, the default), \code{"TD"}, \code{"PD"}, \code{"MhD"}, or a
    function. See \code{\link{chidepth}}.}
  \item{jitter, jitter.amount}{jittering to break ties, for example for ordered
    categorical variables: \code{TRUE}, \code{FALSE}, a logical vector, or column
    numbers or names; see \code{\link{chiplot}}. The default \code{NULL} jitters
    ordered factors only.}
  \item{plot}{logical: draw the chi-plot?}
  \item{verbose}{logical: report removed observations and jittered variables?}
  \item{\dots}{for \code{chimult}, arguments passed to \code{plot.chimult}; for
    \code{plot.chimult}, graphical parameters passed to \code{\link{plot}}.}
  \item{x}{an object of class \code{"chimult"}.}
  \item{pc.level}{probability level for the guidelines.}
  \item{main}{title of the chi-plot.}
  \item{col1, col2}{colours of the horizontal guidelines and of the vertical
    reference lines.}
  \item{pt}{point size of the text listing the groups.}
  \item{box}{logical: draw a box around the list of groups?}
  \item{cex}{size of the points; by default chosen from the sample size.}
  \item{mar}{plot margins, as in \code{\link{par}}.}
  \item{guide.option}{the horizontal guidelines: \code{1}, at the normal
    percentiles for \code{pc.level}; \code{2}, at the percentiles of the
    permutation distribution of the chi-values; \code{3}, both; \code{4}, option
    1 if the two sets are within \code{del} of each other, otherwise option 2.}
  \item{del}{tolerance used by \code{guide.option = 4}.}
}
\details{
  Observations with missing values are removed. Each column is centred at its
  median and divided by its interquartile range.

  For each observation \eqn{u_i}, the generalized orthant is the set of
  \eqn{v} with \eqn{v_j \le u_{ij}} if \eqn{u_{ij} > 0} and \eqn{v_j > u_{ij}}
  otherwise. With \eqn{H_i} the proportion of the sample in this orthant, and
  \eqn{H_{1i}, \dots, H_{ki}} the proportions in the corresponding orthants for
  the variables of each block, the chi-value is
  \deqn{\chi_i = \sqrt{n}\,(H_i - P_i) / \{P_i (1 - P_i S_i)\}^{1/2},}{chi_i = sqrt(n) (H_i - P_i) / sqrt(P_i (1 - P_i S_i)),}
  where \eqn{P_i = \prod_l H_{li}} and
  \eqn{S_i = 1 + \sum_l (1 - H_{li}) / H_{li}}. Under independence of the blocks
  the chi-values are approximately standard normal. The lambda-value of each
  observation is its data depth, transformed by the empirical distribution
  function of the depths.

  The null distributions are estimated by permuting the rows of blocks
  \eqn{1, \dots, k-1}, keeping the rows within each block together.
  \code{sigP.chi} is the proportion of permutations whose mean squared
  chi-value is at least that of the data; \code{sigP.lam} is the p-value of a
  two-sample Kolmogorov--Smirnov test comparing the data depths with those of
  the permuted data.

  The plot shows the variables in each group, the chi-plot, and, if
  \code{sigP.lam <= 1 - pc.level}, a Q-Q plot of the data depths against their
  permutation distribution.

  With the same random number seed, \code{chimult} gives the same results as
  the authors' original \code{chimult.lam} function with
  \code{fda.usc::mdepth.RP}.
}
\value{
  An object of class \code{"chimult"}, a list with components
  \item{chi}{the chi-values of the data.}
  \item{chisq.actual}{the mean of the squared chi-values.}
  \item{chi.values.perm}{the chi-values of the permuted data sets.}
  \item{chisq.values.perm}{the mean squared chi-values of the permuted data sets.}
  \item{sigP.chi}{the permutation p-value of the mean squared chi-value.}
  \item{lambda.actual}{the data depths of the observations.}
  \item{lambda}{the lambda-values: the data depths transformed by their empirical
    distribution function.}
  \item{lambda.values.perm}{the data depths of the permuted data sets.}
  \item{sigP.lam}{the p-value of the Kolmogorov--Smirnov test for the depths.}
  \item{X}{the centred and scaled data.}
  \item{D, k}{the block sizes and the number of blocks.}
  \item{depth, nperm, Nsim}{the settings used.}
  \item{n.removed, rows}{the number of observations removed because of missing
    values, and the rows of \code{Data} used.}
  \item{jittered}{the names of jittered variables.}
  \item{call}{the call.}
}
\references{
  Fisher, N. I. and Lee, A. J. (2026). Chi-plots for assessing multivariate
  dependence. Submitted to \emph{Journal of Computational and Graphical
  Statistics}.
}
\seealso{\code{\link{chipart}} to search over partitions;
  \code{\link{chidepth}}; \code{\link{chiplot}}.}
\examples{
## Two independent blocks of two correlated normal variables
set.seed(161225)
n <- 200
z1 <- matrix(rnorm(2 * n), n) \%*\% chol(matrix(c(1, 0.5, 0.5, 1), 2))
z2 <- matrix(rnorm(2 * n), n) \%*\% chol(matrix(c(1, 0.7, 0.7, 1), 2))
Data <- data.frame(z1, z2)
colnames(Data) <- paste0("Var.", 1:4)
res <- chimult(Data, D = c(2, 2), nperm = 100)
res
plot(res)

## Dependence between the blocks
Data$Var.3 <- Data$Var.3 + 0.8 * Data$Var.1
res2 <- chimult(Data, D = c(2, 2), nperm = 100)
plot(res2, guide.option = 3)

\donttest{
## Five variables, mutually uncorrelated but with a hole in the data
set.seed(280926)
z <- matrix(rnorm(20000 * 5), ncol = 5)
z <- z[sqrt(rowSums((z - 1)^2)) >= 2.5, ][1:200, ]
res3 <- chimult(z, D = rep(1, 5), nperm = 500)
plot(res3, guide.option = 4, pt = 10)
}
}
\keyword{hplot}
\keyword{multivariate}
\keyword{nonparametric}
