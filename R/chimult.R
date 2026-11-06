## Multivariate chi-plots: Fisher & Lee (2026).

## Prepare the data matrix for chimult(): numeric, complete rows, optional jitter.
.chimult_data <- function(Data, jitter, jitter.amount, verbose) {
    X <- .as_numeric_matrix(Data)
    ord <- attr(X, "ordered")
    ok <- stats::complete.cases(X)
    n.removed <- sum(!ok)
    X <- X[ok, , drop = FALSE]
    if (n.removed && verbose)
        message(n.removed, " observation(s) with missing values removed")
    jcols <- .jitter_cols(jitter, X)
    if (!length(jcols) && is.null(jitter) && length(ord)) {
        jcols <- ord
        if (verbose)
            message("jittering ordered factor(s): ",
                    paste(colnames(X)[ord], collapse = ", "))
    }
    X <- .jitter_columns(X, jcols, jitter.amount)
    attr(X, "ordered") <- NULL
    list(X = X, n.removed = n.removed, rows = which(ok), jittered = colnames(X)[jcols])
}

## Centre each column at its median and rescale by its interquartile range
## (by the standard deviation, or not at all, if the IQR is zero).
.centre_scale <- function(X) {
    X <- t(t(X) - apply(X, 2, stats::median))
    s <- apply(X, 2, stats::IQR)
    zero <- !(s > 0)
    if (any(zero)) {
        s[zero] <- apply(X[, zero, drop = FALSE], 2, stats::sd)
        s[!(s > 0)] <- 1
    }
    t(t(X) / s)
}

chimult <- function(Data, D, nperm = 100, Nsim = 100, depth = "RP",
                    jitter = NULL, jitter.amount = 0.25, plot = FALSE,
                    verbose = TRUE, ...) {
    cl <- match.call()
    prep <- .chimult_data(Data, jitter, jitter.amount, verbose)
    X <- prep$X
    n <- nrow(X); p <- ncol(X)
    D <- .check_D(D, p)
    k <- length(D)
    cD <- c(0L, cumsum(D))
    if (n < 5L)
        stop("too few complete observations (", n, ") for a chi-plot", call. = FALSE)
    nperm <- as.integer(nperm)
    if (is.na(nperm) || nperm < 1L)
        stop("'nperm' must be a positive whole number", call. = FALSE)
    depthfun <- .depth_function(depth, Nsim)

    ## Centre each column at its median, rescale using its IQR
    X <- .centre_scale(X)

    ## Basic chi and lambda values
    chi <- .chimult_native(X, D)
    lambda.actual <- depthfun(X)
    lambda <- stats::ecdf(lambda.actual)(lambda.actual)    # EDF of the lambdas
    chisq.actual <- mean(chi^2)

    ## Permutation distributions: permute the rows of blocks 1, ..., k-1
    chi.values.perm <- matrix(0, n, nperm)
    lambda.temp.perm <- matrix(0, n, nperm)
    chisq.values.perm <- numeric(nperm)
    data.perm <- X
    for (iperm in seq_len(nperm)) {
        for (k0 in seq_len(k - 1L)) {
            krange <- (cD[k0] + 1L):cD[k0 + 1L]
            data.perm[, krange] <- X[sample(n), krange]
        }
        chi.values.perm[, iperm] <- .chimult_native(data.perm, D)
        chisq.values.perm[iperm] <- mean(chi.values.perm[, iperm]^2)
        lambda.temp.perm[, iperm] <- depthfun(data.perm)
    }

    sigP.chi <- mean(chisq.actual <= chisq.values.perm)
    sigP.lam <- suppressWarnings(
        stats::ks.test(lambda.actual, as.vector(lambda.temp.perm))$p.value)

    res <- list(
        chi = chi,
        chisq.actual = chisq.actual,
        chi.values.perm = as.vector(chi.values.perm),
        chisq.values.perm = chisq.values.perm,
        sigP.chi = sigP.chi,
        lambda.actual = lambda.actual,
        lambda = lambda,
        lambda.values.perm = as.vector(lambda.temp.perm),
        sigP.lam = sigP.lam,
        X = X,
        D = D,
        k = k,
        depth = if (is.function(depth)) "user-supplied" else .depth_name(depth[1L]),
        nperm = nperm,
        Nsim = as.integer(Nsim),
        n.removed = prep$n.removed,
        rows = prep$rows,
        jittered = prep$jittered,
        call = cl)
    class(res) <- "chimult"
    if (plot) graphics::plot(res, ...)
    res
}

print.chimult <- function(x, ...) {
    cat("Multivariate chi-plot analysis\n")
    gt <- groupText(x)
    for (i in seq_along(gt)) cat("  Group ", i, ": ", gt[i], "\n", sep = "")
    cat("  n = ", nrow(x$X), sep = "")
    if (x$n.removed) cat(" (", x$n.removed, " with missing values removed)", sep = "")
    cat("; depth = ", x$depth, sep = "")
    if (x$depth %in% c("RP", "PD")) cat(" (Nsim = ", x$Nsim, ")", sep = "")
    cat("; nperm = ", x$nperm, "\n", sep = "")
    if (length(x$jittered))
        cat("  Jittered: ", paste(x$jittered, collapse = ", "), "\n", sep = "")
    cat("  Mean squared chi-value: ", formatC(x$chisq.actual, digits = 4, format = "g"),
        ";  permutation p-value = ", .format_p(x$sigP.chi), "\n", sep = "")
    cat("  Data depth (lambda) KS test p-value = ", .format_p(x$sigP.lam), "\n", sep = "")
    invisible(x)
}

## Basic (lambda, chi) chi-plot, plus a Q-Q plot of the data depths if the KS
## test detects a difference between the permuted and unpermuted lambda values.
##   guide.option = 1  horizontal guidelines at normal percentiles for pc.level
##   guide.option = 2  guidelines at the permutation-distribution percentiles
##   guide.option = 3  both sets of guidelines
##   guide.option = 4  option 1 if the two sets are within 'del' of each other,
##                     otherwise option 2
plot.chimult <- function(x, pc.level = 0.95, main = "", col1 = "red", col2 = "blue",
                         pt = 12, box = TRUE, cex = NULL, mar = c(4, 4, 2, 2),
                         guide.option = 1, del = 0.2, ...) {
    res <- x
    if (!guide.option %in% 1:4)
        stop("'guide.option' must be 1, 2, 3 or 4", call. = FALSE)
    lambda <- res$lambda
    chi <- res$chi
    if (is.null(cex)) cex <- 0.75 * .auto_cex(length(chi), 1.3)

    ## guideline parameters
    pc.twosided <- pc.level + 0.5 * (1 - pc.level)
    pc.guideline <- stats::qnorm(pc.twosided)
    y.upper <- max(chi, 1.2 * pc.guideline, na.rm = TRUE)
    y.lower <- min(chi, -1.2 * pc.guideline, na.rm = TRUE)
    close <- TRUE
    if (guide.option > 1) {
        cvp.quants <- stats::quantile(res$chi.values.perm,
                                      c(0.5 * (1 - pc.level), 1 - 0.5 * (1 - pc.level)),
                                      names = FALSE)
        pd.lower <- cvp.quants[1L]; pd.upper <- cvp.quants[2L]
        crit.diff <- max(abs(pd.lower - (-pc.guideline)), abs(pd.upper - pc.guideline))
        close <- (crit.diff < del)
        y.upper <- max(chi, 1.2 * pc.guideline, pd.upper, na.rm = TRUE)
        y.lower <- min(chi, -1.2 * pc.guideline, pd.lower, na.rm = TRUE)
    }
    ylims <- c(y.lower, y.upper)

    xlab <- expression(paste("Transformed ", lambda))
    ylab <- expression(chi)
    mtxt.chi <- paste0("Estimated p-value = ", .format_p(res$sigP.chi))
    show.lam <- isTRUE(res$sigP.lam <= 1 - pc.level)
    op <- graphics::par(mfrow = if (show.lam) c(2, 2) else c(1, 2), pty = "s", mar = mar)
    on.exit(graphics::par(op))

    ## panel listing the groups, then the basic chi-plot
    writeGroupNames(res, pt, box)
    graphics::plot(lambda, chi, pch = 19, xlim = c(0, 1), ylim = ylims,
                   xlab = xlab, ylab = ylab, cex = cex, ...)
    if (nzchar(main)) graphics::title(main = main)
    graphics::mtext(mtxt.chi, cex = 0.75)
    graphics::abline(h = 0)
    normal.lines <- guide.option %in% c(1, 3) || (guide.option == 4 && close)
    perm.lines <- guide.option %in% c(2, 3) || (guide.option == 4 && !close)
    if (normal.lines)
        graphics::abline(h = c(-pc.guideline, pc.guideline), lty = 1, col = col1)
    if (perm.lines)
        graphics::abline(h = c(pd.lower, pd.upper), lty = 4, col = col1)
    graphics::abline(v = (1:9) / 10, lty = 3, col = col2)

    if (show.lam) {
        ## The data depth test: compare raw lambda values with their permutation distribution
        lambda.raw <- res$lambda.actual
        lambda.perm <- res$lambda.values.perm
        lims <- c(0, max(0.5, lambda.raw, lambda.perm))
        xlab.lam <- expression(paste("Null distribution of ", lambda))
        ylab.lam <- expression(paste("Actual distribution of ", lambda))
        mtxt.lam <- paste0("Estimated p-value = ", .format_p(res$sigP.lam))
        stats::qqplot(lambda.perm, lambda.raw, plot.it = TRUE, xlab = xlab.lam,
                      ylab = ylab.lam, xlim = lims, ylim = lims, pch = 19, cex = cex)
        graphics::mtext(mtxt.lam, cex = 0.75)
        graphics::abline(0, 1)
    }
    invisible(res)
}
