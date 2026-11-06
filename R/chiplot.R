## Bivariate chi-plots: Fisher & Switzer (1985, 2001).

## Control-line constants for ctrl.lines = 1, 2, 3 (90%, 95%, 99%);
## the lines are drawn at +/- clims / sqrt(n).
.chi_clims <- c(1.54, 1.78, 2.18)

.ctrl_index <- function(ctrl.lines) {
    if (is.null(ctrl.lines) || length(ctrl.lines) != 1L || is.na(ctrl.lines))
        stop("'ctrl.lines' must be 0, 1, 2 or 3 (or 0.90, 0.95, 0.99)", call. = FALSE)
    if (ctrl.lines %in% 0:3) return(as.integer(ctrl.lines))
    m <- match(round(ctrl.lines, 2), c(0.90, 0.95, 0.99))
    if (is.na(m))
        stop("'ctrl.lines' must be 0, 1, 2 or 3 (or 0.90, 0.95, 0.99)", call. = FALSE)
    m
}

## Two numeric vectors from (x, y), or from the first two columns of x.
.xy_data <- function(x, y, xname, yname) {
    if (is.null(y)) {
        if (is.null(dim(x)) || ncol(x) < 2L)
            stop("supply 'x' and 'y', or a matrix or data frame with two columns as 'x'",
                 call. = FALSE)
        X <- .as_numeric_matrix(x[, 1:2, drop = FALSE])
    } else {
        X <- .as_numeric_matrix(data.frame(x = x, y = y))
        colnames(X) <- c(xname, yname)
    }
    X
}

chivals <- function(x, y = NULL, jitter = NULL, jitter.amount = 0.25) {
    X <- .xy_data(x, y, deparse1(substitute(x)), deparse1(substitute(y)))
    ord <- attr(X, "ordered")
    ok <- stats::complete.cases(X)
    X <- X[ok, , drop = FALSE]
    n <- nrow(X)
    if (n < 5L)
        stop("insufficient data for a chi-plot: fewer than 5 complete pairs",
             call. = FALSE)
    jcols <- .jitter_cols(jitter, X)
    if (!length(jcols) && is.null(jitter) && length(ord)) jcols <- ord
    X <- .jitter_columns(X, jcols, jitter.amount)
    r <- .chivals2_native(X[, 1L], X[, 2L])
    structure(list(lambda = r$lambda,
                   chi = r$chi,
                   index = which(ok)[r$index],
                   n = n,
                   npts = length(r$chi),
                   n.missing = sum(!ok),
                   names = colnames(X)),
              class = "chiplot")
}

chiplot <- function(x, y = NULL, ctrl.lines = 2, jitter = NULL, jitter.amount = 0.25,
                    cex.points = NULL, pch.type = 19, dot.col = "black", main = "",
                    cex.axis = 1, cex.lam = 1.2, plot = TRUE, ...) {
    if (!is.null(y)) {
        X <- data.frame(x, y)
        names(X) <- c(deparse1(substitute(x)), deparse1(substitute(y)))
        res <- chivals(X, jitter = jitter, jitter.amount = jitter.amount)
    } else {
        res <- chivals(x, jitter = jitter, jitter.amount = jitter.amount)
    }
    if (plot)
        graphics::plot(res, ctrl.lines = ctrl.lines, cex.points = cex.points,
                       pch.type = pch.type, dot.col = dot.col, main = main,
                       cex.axis = cex.axis, cex.lam = cex.lam, ...)
    invisible(res)
}

## Draw one chi-plot in the current plot region (used by plot.chiplot and chipairs).
.draw_chiplot <- function(res, ctrl.lines = 2, cex.points = NULL, pch.type = 19,
                          dot.col = "black", col.ctrl = "red", lty.ctrl = 1,
                          axes = TRUE, cex.axis = 1, ...) {
    if (is.null(cex.points)) cex.points <- .auto_cex(res$npts)
    graphics::plot(res$lambda, res$chi, xlim = c(-1, 1), ylim = c(-1, 1), type = "n",
                   xlab = "", ylab = "", axes = FALSE, ...)
    graphics::box()
    if (axes) {
        graphics::axis(1, cex.axis = cex.axis)
        graphics::axis(2, cex.axis = cex.axis, las = 1)
    }
    graphics::abline(h = 0, v = 0)
    ic <- .ctrl_index(ctrl.lines)
    if (ic > 0L) {
        cl <- .chi_clims[ic] / sqrt(res$n)
        graphics::abline(h = c(-cl, cl), col = col.ctrl, lty = lty.ctrl)
    }
    graphics::points(res$lambda, res$chi, cex = cex.points, pch = pch.type, col = dot.col)
    invisible(NULL)
}

plot.chiplot <- function(x, ctrl.lines = 2, cex.points = NULL, pch.type = 19,
                         dot.col = "black", main = "", cex.axis = 1, cex.lam = 1.2,
                         xlab = expression(lambda), ylab = expression(chi),
                         col.ctrl = "red", lty.ctrl = 1, ...) {
    op <- graphics::par(pty = "s")
    on.exit(graphics::par(op))
    .draw_chiplot(x, ctrl.lines = ctrl.lines, cex.points = cex.points,
                  pch.type = pch.type, dot.col = dot.col, col.ctrl = col.ctrl,
                  lty.ctrl = lty.ctrl, cex.axis = cex.axis, ...)
    graphics::title(xlab = xlab, ylab = ylab, cex.lab = cex.lam, line = 2.5)
    if (nzchar(main)) graphics::title(main = main)
    invisible(x)
}

print.chiplot <- function(x, ...) {
    cat("Chi-plot of", x$names[2L], "against", x$names[1L], "\n")
    cat("  ", x$n, " complete pairs", sep = "")
    if (x$n.missing) cat(" (", x$n.missing, " with missing values removed)", sep = "")
    cat("; ", x$npts, " points plotted\n", sep = "")
    cat("  Proportion of |chi| beyond the 95% control lines: ",
        formatC(mean(abs(x$chi) > .chi_clims[2L] / sqrt(x$n)), format = "f", digits = 3),
        "\n", sep = "")
    invisible(x)
}

## Scatterplot matrix with chi-plots above the diagonal.
chipairs <- function(Data, labels = NULL, ctrl.lines = 2, jitter = NULL,
                     jitter.amount = 0.25, cex.points = NULL, pch.type = 19,
                     dot.col = "black", cex.label = NULL, cex.axis = 0.7,
                     plgap = 0.2, axes = TRUE, main = NULL, ...) {
    X <- .as_numeric_matrix(Data)
    p <- ncol(X)
    if (p < 2L) stop("'Data' must have at least two columns", call. = FALSE)
    if (is.null(labels)) labels <- colnames(X)
    if (length(labels) != p)
        stop("'labels' must have one element per column of 'Data'", call. = FALSE)
    jcols <- .jitter_cols(jitter, X)
    if (!length(jcols) && is.null(jitter)) jcols <- attr(X, "ordered")
    X <- .jitter_columns(X, jcols, jitter.amount)

    pscale <- max(0.35, min(1, 2.5 / p))
    if (is.null(cex.label)) cex.label <- max(0.6, min(1.5, 6 / p))
    outer.mar <- if (axes) 2.5 else 1
    op <- graphics::par(mfrow = c(p, p), mar = rep(plgap, 4),
                        oma = c(outer.mar, outer.mar,
                                outer.mar + if (is.null(main)) 0 else 1.5, outer.mar),
                        mgp = c(2, 0.4, 0), tcl = -0.2, pty = "m")
    on.exit(graphics::par(op))

    for (i in seq_len(p)) {
        for (j in seq_len(p)) {
            if (i == j) {
                graphics::plot.new()
                graphics::plot.window(c(0, 1), c(0, 1))
                graphics::box()
                graphics::text(0.5, 0.5, labels[i], cex = cex.label)
                next
            }
            ok <- stats::complete.cases(X[, i], X[, j])
            ## as in pairs(): column j on the horizontal axis, row i on the vertical
            x <- X[ok, j]
            y <- X[ok, i]
            if (sum(ok) < 5L) {
                graphics::plot.new()
                graphics::box()
                next
            }
            if (i < j) {
                res <- structure(c(.chivals2_native(x, y)[c("lambda", "chi")],
                                   list(n = sum(ok))), class = "chiplot")
                res$npts <- length(res$chi)
                cexp <- if (is.null(cex.points)) .auto_cex(res$npts, pscale) else cex.points
                .draw_chiplot(res, ctrl.lines = ctrl.lines, cex.points = cexp,
                              pch.type = pch.type, dot.col = dot.col, axes = FALSE,
                              lty.ctrl = 2, col.ctrl = "red", ...)
                if (axes && i == 1L)
                    graphics::axis(3, at = c(-1, 0, 1), cex.axis = cex.axis)
                if (axes && j == p)
                    graphics::axis(4, at = c(-1, 0, 1), cex.axis = cex.axis, las = 1)
            } else {
                cexp <- if (is.null(cex.points)) .auto_cex(sum(ok), pscale) else cex.points
                graphics::plot(x, y, xlab = "", ylab = "", axes = FALSE,
                               cex = cexp, pch = pch.type, col = dot.col, ...)
                graphics::box()
                if (axes && i == p) graphics::axis(1, cex.axis = cex.axis)
                if (axes && j == 1L) graphics::axis(2, cex.axis = cex.axis, las = 1)
            }
        }
    }
    if (!is.null(main))
        graphics::mtext(main, side = 3, outer = TRUE, line = outer.mar + 0.3,
                        cex = 1.2, font = 2)
    invisible(NULL)
}
