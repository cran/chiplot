## Data depth functions used for the lambda-values of multivariate chi-plots.
## All are written in base R. "RP", "TD" and "MhD" reproduce
## fda.usc::mdepth.RP, fda.usc::mdepth.TD and fda.usc::mdepth.MhD; for "RP"
## the random directions are generated in the same way, so results agree
## exactly for the same random number seed.

chidepth <- function(X, depth = c("RP", "TD", "PD", "MhD"), Nsim = 100) {
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    if (anyNA(X)) stop("'X' must not contain missing values", call. = FALSE)
    fun <- .depth_function(depth, Nsim)
    fun(X)
}

## Resolve the 'depth' argument to a function of the data matrix that returns
## one depth value per row.
.depth_function <- function(depth, Nsim) {
    if (is.function(depth)) {
        f <- depth
        return(function(X) {
            d <- f(X)
            if (is.list(d) && !is.null(d$dep)) d <- d$dep
            d <- as.numeric(d)
            if (length(d) != nrow(X) || anyNA(d))
                stop("the 'depth' function must return one non-missing value per row",
                     call. = FALSE)
            d
        })
    }
    if (!is.character(depth))
        stop("'depth' must be one of \"RP\", \"TD\", \"PD\", \"MhD\" or a function",
             call. = FALSE)
    name <- .depth_name(depth[1L])
    Nsim <- as.integer(Nsim)
    if (name %in% c("RP", "PD") && (is.na(Nsim) || Nsim < 1L))
        stop("'Nsim' must be a positive whole number", call. = FALSE)
    switch(name,
           RP  = function(X) .depth_RP(X, Nsim),
           TD  = function(X) .depth_TD(X),
           PD  = function(X) .depth_PD(X, Nsim),
           MhD = function(X) .depth_MhD(X))
}

.depth_name <- function(depth) {
    key <- tolower(gsub("[^A-Za-z]", "", depth))
    map <- c(rp = "RP", mdepthrp = "RP", randomprojection = "RP",
             td = "TD", mdepthtd = "TD", tukey = "TD",
             pd = "PD", projection = "PD",
             mhd = "MhD", mdepthmhd = "MhD", mahalanobis = "MhD")
    if (!key %in% names(map))
        stop("unknown 'depth': ", depth,
             " (use \"RP\", \"TD\", \"PD\", \"MhD\" or a function)", call. = FALSE)
    unname(map[key])
}

## Univariate halfspace (Tukey) depth of each z[i] in the sample z:
## min(Fn(z), 1 - Fn(z - 1e-15)), with Fn the empirical d.f. (as fda.usc::mdepth.TD1)
.depth_TD1 <- function(z) {
    s <- sort(z)
    n <- length(z)
    pmin(findInterval(z, s) / n, 1 - findInterval(z - 1e-15, s) / n)
}

## Random projection depth (as fda.usc::mdepth.RP): average of the univariate
## halfspace depths of the data projected onto Nsim random directions.
.depth_RP <- function(X, Nsim) {
    d <- ncol(X)
    u <- matrix(stats::runif(d * Nsim, -1, 1), Nsim, d)
    proj <- u / sqrt(rowSums(u * u))
    z <- proj %*% t(X)
    pdep <- matrix(NA_real_, nrow = nrow(X), ncol = Nsim)
    for (i in seq_len(Nsim)) pdep[, i] <- .depth_TD1(z[i, ])
    apply(pdep, 1, mean)
}

## Coordinate-wise halfspace depth (as fda.usc::mdepth.TD): the average of the
## univariate halfspace depths of the coordinates. This is not the exact
## multivariate Tukey depth, which is very costly to compute in high dimensions.
.depth_TD <- function(X) {
    pdep <- apply(X, 2, .depth_TD1)
    if (!is.matrix(pdep)) pdep <- matrix(pdep, nrow = nrow(X))
    apply(pdep, 1, mean)
}

## Projection depth (Zuo & Serfling 2000), approximated using Nsim random
## directions: 1 / (1 + max_u |u'x - med(u'X)| / mad(u'X)).
.depth_PD <- function(X, Nsim) {
    d <- ncol(X)
    u <- matrix(stats::rnorm(d * Nsim), Nsim, d)
    proj <- u / sqrt(rowSums(u * u))
    z <- proj %*% t(X)
    out <- numeric(nrow(X))
    for (i in seq_len(Nsim)) {
        zi <- z[i, ]
        s <- stats::mad(zi)
        if (s > 0) out <- pmax(out, abs(zi - stats::median(zi)) / s)
    }
    1 / (1 + out)
}

## Mahalanobis depth (as fda.usc::mdepth.MhD): 1 / (1 + squared Mahalanobis distance).
.depth_MhD <- function(X) {
    mu <- colMeans(X)
    S <- stats::cov(X)
    Sinv <- tryCatch(solve(S), error = function(e) NULL)
    if (is.null(Sinv)) {
        sv <- svd(S)
        pos <- sv$d > max(dim(S)) * max(sv$d) * .Machine$double.eps
        Sinv <- sv$v[, pos, drop = FALSE] %*% (t(sv$u[, pos, drop = FALSE]) / sv$d[pos])
    }
    Xc <- sweep(X, 2, mu)
    1 / (1 + rowSums((Xc %*% Sinv) * Xc))
}
