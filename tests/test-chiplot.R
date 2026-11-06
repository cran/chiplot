## Tests for package chiplot (base R only; run by R CMD check)
library(chiplot)

## ---- bivariate chi-plot -------------------------------------------------------
set.seed(1)
x <- rnorm(100); y <- 0.5 * x + rnorm(100)
a <- chivals(x, y)
b <- chivals(y, x)
stopifnot(inherits(a, "chiplot"),
          all.equal(a$chi, b$chi), all.equal(a$lambda, b$lambda),   # symmetric
          all(abs(a$chi) <= 1), all(abs(a$lambda) <= 1),
          a$n == 100, a$npts == length(a$chi))

## reference values from the original Fortran algorithm (chisub.f), n = 6
xr <- c(1, 2, 3, 4, 5, 6); yr <- c(2, 1, 4, 3, 6, 5)
r <- chivals(xr, yr)
F <- G <- (rank(xr) - 1) / 5
G <- (rank(yr) - 1) / 5
H <- sapply(1:6, function(i) sum(xr[-i] <= xr[i] & yr[-i] <= yr[i]) / 5)
lam <- 4 * ifelse((F - 0.5) * (G - 0.5) < 0, -1, 1) * pmax((F - 0.5)^2, (G - 0.5)^2)
keep <- abs(lam) < 4 * (1 / 5 - 0.5)^2 & F * G * (1 - F) * (1 - G) != 0
chi <- (H - F * G) / sqrt(F * G * (1 - F) * (1 - G))
stopifnot(all.equal(r$lambda, lam[keep]), all.equal(r$chi, chi[keep]),
          identical(r$index, which(keep)))

## missing values
x2 <- x; x2[c(3, 50)] <- NA
d <- chivals(x2, y)
stopifnot(d$n == 98, d$n.missing == 2, !any(d$index %in% c(3, 50)))
stopifnot(inherits(try(chivals(1:4, 1:4), silent = TRUE), "try-error"))

## plotting runs
pdf(NULL)
chiplot(x, y, main = "test")
chipairs(cbind(x, y, z = rnorm(100)))
invisible(dev.off())

## ---- multivariate chi-values: reference implementation of Hval.cpp ---------
ref_chimult <- function(X, D) {
    n <- nrow(X); p <- ncol(X); cD <- c(0, cumsum(D))
    clamp <- function(h) { h[h < 1e-8] <- 1 / n; h[h > 1 - 1e-8] <- 1 - 1 / n; h }
    sapply(seq_len(n), function(i) {
        xi <- X[i, ]
        IN <- t(t(X) <= xi)
        IN[, xi <= 0] <- !IN[, xi <= 0]
        h <- clamp(mean(rowSums(IN) == p))
        hb <- sapply(seq_along(D), function(l) {
            cols <- (cD[l] + 1):cD[l + 1]
            clamp(mean(rowSums(IN[, cols, drop = FALSE]) == length(cols)))
        })
        P <- prod(hb)
        sqrt(n) * (h - P) / sqrt(P * (1 - P * (1 + sum((1 - hb) / hb))))
    })
}
set.seed(2)
for (n in c(30, 64, 65, 129)) {
    X <- matrix(rnorm(n * 4), n, 4)
    X[, 2] <- round(X[, 2])                      # ties
    X <- t(t(X) - apply(X, 2, median))
    for (D in list(c(2, 2), c(1, 3), c(1, 1, 2))) {
        stopifnot(isTRUE(all.equal(chiplot:::.chimult_native(X, D), ref_chimult(X, D),
                                   tolerance = 1e-12)))
    }
}

## ---- chimult ------------------------------------------------------------------
set.seed(4)
Z <- matrix(rnorm(400), 100, 4); Z[, 3] <- Z[, 3] + Z[, 1]
res <- chimult(Z, c(2, 2), nperm = 50, verbose = FALSE)
stopifnot(inherits(res, "chimult"), res$sigP.chi < 0.05,
          length(res$chi) == 100, length(res$chi.values.perm) == 100 * 50,
          all(res$lambda > 0 & res$lambda <= 1))
set.seed(4); res.b <- chimult(Z, c(2, 2), nperm = 50, verbose = FALSE)
stopifnot(identical(res$sigP.chi, res.b$sigP.chi))                   # reproducible
for (dp in c("TD", "PD", "MhD"))
    stopifnot(length(chimult(Z, c(2, 2), nperm = 5, depth = dp, verbose = FALSE)$lambda) == 100)
Zm <- Z; Zm[5, 2] <- NA
stopifnot(nrow(chimult(Zm, c(2, 2), nperm = 5, verbose = FALSE)$X) == 99)
stopifnot(inherits(try(chimult(Z, c(2, 3), verbose = FALSE), silent = TRUE), "try-error"))
pdf(NULL); plot(res, guide.option = 3); plot(res, guide.option = 4); invisible(dev.off())

## random projection depth agrees with fda.usc::mdepth.RP for the same seed
if (requireNamespace("fda.usc", quietly = TRUE)) {
    set.seed(9); d1 <- fda.usc::mdepth.RP(Z, proj = 50)$dep
    set.seed(9); d2 <- chidepth(Z, "RP", Nsim = 50)
    stopifnot(isTRUE(all.equal(unname(d1), d2)))
    stopifnot(isTRUE(all.equal(unname(fda.usc::mdepth.TD(Z)$dep), chidepth(Z, "TD"))))
}

## ---- partitions and chipart -----------------------------------------------------
stopifnot(ncol(chi.partitions(4)) == 14, ncol(chi.partitions(6)) == 202,
          ncol(chi.partitions(5, kfix = 2)) == 5, ncol(chi.partitions(5, nblocks = 2)) == 15,
          !anyDuplicated(apply(chi.partitions(6), 2, paste, collapse = "")))
set.seed(5)
W <- matrix(rnorm(80 * 4), 80, 4); W[, 4] <- W[, 4] + W[, 2]
cp <- chipart(W, nperm = 30, verbose = FALSE)
stopifnot(inherits(cp, "chipart"), nrow(cp$summary) == 14,
          length(cp$results) == sum(cp$summary$Interesting))
## partitions separating variables 2 and 4 should be the interesting ones
sep24 <- apply(cp$partitions, 2, function(g) g[2] != g[4])
stopifnot(all(cp$summary$Interesting[sep24]))
pdf(NULL); plot(cp, ask = FALSE); invisible(dev.off())
cat("All tests passed\n")
