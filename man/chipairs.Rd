\name{chipairs}
\alias{chipairs}
\title{Scatterplot Matrix with Chi-Plots}
\description{
  A scatterplot matrix, like \code{\link{pairs}}, with the scatterplots below
  the diagonal and the corresponding chi-plots above the diagonal.
}
\usage{
chipairs(Data, labels = NULL, ctrl.lines = 2, jitter = NULL,
         jitter.amount = 0.25, cex.points = NULL, pch.type = 19,
         dot.col = "black", cex.label = NULL, cex.axis = 0.7,
         plgap = 0.2, axes = TRUE, main = NULL, \dots)
}
\arguments{
  \item{Data}{a numeric matrix or data frame (\eqn{n \times p}{n x p}), each
    column a variable. Ordered factors are allowed.}
  \item{labels}{variable names for the diagonal; by default the column names.}
  \item{ctrl.lines}{control lines on the chi-plots, as in \code{\link{chiplot}}.}
  \item{jitter, jitter.amount}{jittering of variables to break ties, as in
    \code{\link{chiplot}}; by default ordered factors are jittered.}
  \item{cex.points}{size of the points. By default it is chosen from the number
    of points and the number of panels.}
  \item{pch.type, dot.col}{plotting symbol and colour of the points.}
  \item{cex.label}{size of the variable names on the diagonal.}
  \item{cex.axis}{size of the axis annotation.}
  \item{plgap}{space between panels, in lines of text.}
  \item{axes}{logical: draw axes on the outer panels? The scatterplots are
    annotated on the bottom and left, the chi-plots (on the scale -1 to 1) on the
    top and right.}
  \item{main}{an overall title.}
  \item{\dots}{further graphical parameters passed to \code{\link{plot}}.}
}
\details{
  As in \code{\link{pairs}}, the panel in row \eqn{i} and column \eqn{j} has
  variable \eqn{j} on the horizontal axis and variable \eqn{i} on the vertical
  axis. Each pair of variables uses all the observations for which both are
  present. The chi-plot is symmetric in the two variables.
}
\value{\code{NULL}, invisibly.}
\seealso{\code{\link{chiplot}}, \code{\link{pairs}}}
\examples{
chipairs(geochem[, 1:5], main = "geochem")

## a smaller sample, with the points emphasised
chipairs(geochem[1:300, c("Cr2O3", "MgO", "CaO", "Ni")], cex.points = 0.5)
}
\keyword{hplot}
