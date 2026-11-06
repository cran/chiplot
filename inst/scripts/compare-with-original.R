## Compare package chiplot with the original code of Fisher & Lee, using the same
## random number seed. Needs the original files 'Chiplot Functions.21Oct26.R' and
## 'Hval.cpp' in 'origDir', and packages Rcpp, fda.usc, magrittr and stringr.

origDir <- "."            # ***** folder containing the original files
library(chiplot)
library(Rcpp); library(fda.usc); library(magrittr); library(stringr)
old <- new.env()
sourceCpp(file.path(origDir, "Hval.cpp"), env = old)
sys.source(file.path(origDir, "Chiplot Functions.21Oct26.R"), envir = old)

rmvn <- function(n, sigma) matrix(rnorm(n * ncol(sigma)), n) %*% chol(sigma)
compare <- function(label, Data, D, seed, nperm = 200) {
    set.seed(seed); o <- old$chimult.lam(Data, D, nperm = nperm, Nsim = 100)
    set.seed(seed); n <- chimult(Data, D, nperm = nperm, Nsim = 100, verbose = FALSE)
    els <- c("chi", "chisq.actual", "chi.values.perm", "chisq.values.perm", "sigP.chi",
             "lambda.actual", "lambda", "lambda.values.perm", "sigP.lam", "X")
    d <- sapply(els, function(e) max(abs(unname(o[[e]]) - unname(n[[e]]))))
    t.old <- system.time({ set.seed(seed); old$chimult.lam(Data, D, nperm = nperm, Nsim = 100) })[[3]]
    t.new <- system.time({ set.seed(seed); chimult(Data, D, nperm = nperm, Nsim = 100, verbose = FALSE) })[[3]]
    cat(sprintf("\n%s\n  p-values old: %.3f, %.4g   new: %.3f, %.4g\n  time old %.1f s, new %.1f s\n",
                label, o$sigP.chi, o$sigP.lam, n$sigP.chi, n$sigP.lam, t.old, t.new))
    cat("  largest differences:\n"); print(signif(d, 3))
}

S <- rbind(cbind(matrix(c(1, .5, .5, 1), 2), matrix(0, 2, 2)),
           cbind(matrix(0, 2, 2), matrix(c(1, .7, .7, 1), 2)))
set.seed(161225); Ex1 <- data.frame(rmvn(200, S)); names(Ex1) <- paste0("Var.", 1:4)
compare("Example 3.1", Ex1, c(2, 2), 1)
S[1:2, 3:4] <- S[3:4, 1:2] <- 0.4
set.seed(181225); Ex2 <- data.frame(rmvn(200, S)); names(Ex2) <- paste0("Var.", 1:4)
compare("Example 3.2", Ex2, c(2, 2), 2)
set.seed(3); Ex3 <- data.frame(matrix(rnorm(1000 * 6), 1000, 6))
compare("n = 1000, p = 6, three blocks", Ex3, c(2, 3, 1), 3)
