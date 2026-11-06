\name{chipart}
\alias{chipart}
\alias{chi.partitions}
\alias{plot.chipart}
\alias{print.chipart}
\title{Search All Partitions of the Variables with Multivariate Chi-Plots}
\description{
  \code{chipart} enumerates the partitions of the \eqn{p} variables into two or
  more blocks, runs \code{\link{chimult}} for each, and keeps the results for
  partitions where either test statistic has p-value less than \code{SigProb}.
  \code{chi.partitions} lists the partitions examined.
}
\usage{
chipart(Data, kfix = 0, SigProb = 0.05, nperm = 100, Nsim = 100,
        depth = "RP", nblocks = NULL, nperm.screen = NULL,
        screen.level = 0.2, jitter = NULL, jitter.amount = 0.25,
        verbose = TRUE)

chi.partitions(p, kfix = 0, nblocks = NULL)

\method{plot}{chipart}(x, which = NULL, ask = NULL, \dots)

\method{print}{chipart}(x, \dots)
}
\arguments{
  \item{Data}{a numeric matrix or data frame (\eqn{n \times p}{n x p}).}
  \item{kfix}{if positive, the first \code{kfix} variables are always kept
    together as one block, and only the remaining variables are partitioned.
    This reduces the number of partitions from \eqn{B_p - 1} to
    \eqn{B_{p - kfix}}{B_(p - kfix)}, where \eqn{B_m} is the \eqn{m}th Bell
    number.}
  \item{SigProb}{the significance level used to select interesting partitions.}
  \item{nperm, Nsim, depth}{as in \code{\link{chimult}}.}
  \item{nblocks}{if given, only partitions with these numbers of blocks are
    examined; for example, \code{nblocks = 2} compares each subset of the
    variables with its complement.}
  \item{nperm.screen}{if given, a smaller number of permutations used to screen
    all partitions; partitions with either p-value at most \code{screen.level}
    are then re-run with \code{nperm} permutations.}
  \item{screen.level}{the screening threshold used with \code{nperm.screen}.}
  \item{jitter, jitter.amount}{jittering to break ties, applied once to the
    data before the search, as in \code{\link{chimult}}.}
  \item{verbose}{logical: report progress?}
  \item{p}{the number of variables.}
  \item{x}{an object of class \code{"chipart"}.}
  \item{which}{partition numbers to plot; by default all interesting partitions.}
  \item{ask}{logical: pause before each plot? By default, on interactive devices
    when there is more than one plot.}
  \item{\dots}{arguments passed to \code{\link{plot.chimult}}.}
}
\details{
  The number of partitions grows quickly with \eqn{p}: the Bell numbers are 1,
  2, 5, 15, 52, 203, 877, 4140, 21147, 115975, \dots. The options
  \code{kfix}, \code{nblocks}, \code{nperm.screen} and a fast depth such as
  \code{depth = "MhD"} reduce the computing time.

  For each partition the variables are reordered so that the blocks are
  contiguous, keeping the original order within each block, and
  \code{\link{chimult}} is called. A partition is \dQuote{interesting} if
  \code{sigP.chi} or \code{sigP.lam} is less than \code{SigProb}. No allowance
  is made for multiple testing.
}
\value{
  \code{chipart} returns an object of class \code{"chipart"}, a list with
  components
  \item{summary}{a data frame with one row per partition: the partition number,
    the groups, the number of blocks, \code{sigP.chi}, \code{sigP.lam}, the
    number of permutations used, and whether the partition is interesting.}
  \item{results}{a list of \code{"chimult"} objects for the interesting
    partitions, named by partition number.}
  \item{partitions}{the partitions: a matrix with one column per partition
    giving the block number of each variable.}
  \item{variables, SigProb, kfix, nperm, nperm.screen, Nsim, depth, n,
    n.removed, call}{details of the search.}

  \code{chi.partitions} returns the matrix of partitions.
}
\references{
  Fisher, N. I. and Lee, A. J. (2026). Chi-plots for assessing multivariate
  dependence. Submitted to \emph{Journal of Computational and Graphical
  Statistics}.
}
\seealso{\code{\link{chimult}}}
\examples{
## the partitions examined for 4 variables, and with the first kept fixed
chi.partitions(4)
chi.partitions(4, kfix = 1)

## four variables, of which the second and fourth are dependent
set.seed(5)
W <- matrix(rnorm(100 * 4), 100, 4)
W[, 4] <- W[, 4] + W[, 2]
colnames(W) <- paste0("V", 1:4)
res <- chipart(W, nperm = 20, verbose = FALSE)
res
plot(res, which = 3)

\donttest{
## screening with fewer permutations
res2 <- chipart(W, nperm = 200, nperm.screen = 20)
res2$summary
}
}
\keyword{multivariate}
\keyword{hplot}
