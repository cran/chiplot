## Search all partitions of the variables for departures from independence.

## Bell numbers B(0), ..., B(m) via the Bell triangle.
.bell <- function(m) {
    out <- numeric(m + 1L)
    out[1L] <- 1
    row <- 1
    for (i in seq_len(m)) {
        new <- numeric(i + 1L)
        new[1L] <- row[i]
        for (j in seq_len(i)) new[j + 1L] <- new[j] + row[j]
        row <- new
        out[i + 1L] <- row[1L]
    }
    out
}

## All set partitions of p items, as a p x B(p) integer matrix of group labels
## (restricted growth strings, in lexicographic order; the first column puts
## all items in one group).
.setparts <- function(p) {
    p <- as.integer(p)
    if (p < 1L) stop("'p' must be at least 1", call. = FALSE)
    nb <- .bell(p)[p + 1L]
    out <- matrix(0L, p, nb)
    a <- rep(1L, p)
    m <- rep(1L, p)          # m[i] = max(a[1:i])
    col <- 0L
    repeat {
        col <- col + 1L
        out[, col] <- a
        i <- p
        while (i > 1L && a[i] > m[i - 1L]) i <- i - 1L
        if (i == 1L) break
        a[i] <- a[i] + 1L
        m[i] <- max(m[i - 1L], a[i])
        if (i < p) {
            a[(i + 1L):p] <- 1L
            m[(i + 1L):p] <- m[i]
        }
    }
    out
}

## Partitions to be examined: all those with at least two groups, or, if
## kfix > 0, those in which the first kfix variables form one group and the
## remaining variables are partitioned in every possible way.
chi.partitions <- function(p, kfix = 0, nblocks = NULL) {
    p <- as.integer(p); kfix <- as.integer(kfix)
    if (is.na(kfix) || kfix < 0L || kfix >= p)
        stop("'kfix' must be a whole number between 0 and p - 1", call. = FALSE)
    if (kfix == 0L) {
        parts <- .setparts(p)[, -1L, drop = FALSE]
    } else {
        rest <- .setparts(p - kfix) + 1L
        parts <- rbind(matrix(1L, kfix, ncol(rest)), rest)
    }
    if (!is.null(nblocks)) {
        k <- apply(parts, 2, max)
        parts <- parts[, k %in% nblocks, drop = FALSE]
    }
    parts
}

chipart <- function(Data, kfix = 0, SigProb = 0.05, nperm = 100, Nsim = 100,
                    depth = "RP", nblocks = NULL, nperm.screen = NULL,
                    screen.level = 0.2, jitter = NULL, jitter.amount = 0.25,
                    verbose = TRUE) {
    cl <- match.call()
    prep <- .chimult_data(Data, jitter, jitter.amount, verbose)
    X <- prep$X
    p <- ncol(X)
    if (p < 2L) stop("'Data' must have at least two columns", call. = FALSE)
    parts <- chi.partitions(p, kfix, nblocks)
    npart <- ncol(parts)
    if (npart == 0L) stop("no partitions satisfy the restrictions", call. = FALSE)
    if (!is.null(nperm.screen)) {
        nperm.screen <- as.integer(nperm.screen)
        if (is.na(nperm.screen) || nperm.screen < 1L)
            stop("'nperm.screen' must be a positive whole number", call. = FALSE)
    }
    if (verbose) {
        message("Number of partitions = ", npart)
        if (npart > 2000L)
            message("This is a large search: consider kfix, nblocks or nperm.screen.")
    }

    run_one <- function(jpart, np) {
        g <- parts[, jpart]
        ord <- order(g)
        chimult(X[, ord, drop = FALSE], tabulate(g), nperm = np, Nsim = Nsim,
                depth = depth, jitter = FALSE, verbose = FALSE)
    }

    sigP.chi <- sigP.lam <- rep(NA_real_, npart)
    groups <- character(npart)
    results <- vector("list", npart)
    nperm.used <- rep(as.integer(nperm), npart)
    first.np <- if (is.null(nperm.screen)) nperm else nperm.screen
    for (jpart in seq_len(npart)) {
        if (verbose && (jpart %% max(1L, npart %/% 20L) == 0L || jpart == npart))
            message("  partition ", jpart, " of ", npart)
        res <- run_one(jpart, first.np)
        nperm.used[jpart] <- res$nperm
        if (!is.null(nperm.screen) &&
            min(res$sigP.chi, res$sigP.lam) <= screen.level && nperm > nperm.screen) {
            res <- run_one(jpart, nperm)
            nperm.used[jpart] <- res$nperm
        }
        sigP.chi[jpart] <- res$sigP.chi
        sigP.lam[jpart] <- res$sigP.lam
        groups[jpart] <- groupLabels(res)
        if (min(res$sigP.chi, res$sigP.lam) < SigProb) results[[jpart]] <- res
    }

    interesting <- pmin(sigP.chi, sigP.lam) < SigProb
    summary <- data.frame(Partition = seq_len(npart), Groups = groups,
                          Blocks = apply(parts, 2, max),
                          sigP.chi = sigP.chi, sigP.lam = sigP.lam,
                          nperm = nperm.used, Interesting = interesting,
                          stringsAsFactors = FALSE)
    names(results) <- as.character(seq_len(npart))
    results <- results[interesting]

    out <- list(summary = summary, results = results, partitions = parts,
                variables = colnames(X), SigProb = SigProb, kfix = as.integer(kfix),
                nperm = as.integer(nperm), nperm.screen = nperm.screen, Nsim = Nsim,
                depth = if (is.function(depth)) "user-supplied" else .depth_name(depth[1L]),
                n = nrow(X), n.removed = prep$n.removed, call = cl)
    class(out) <- "chipart"
    out
}

print.chipart <- function(x, ...) {
    s <- x$summary
    cat("Search of partitions for departures from independence\n")
    cat("  ", length(x$variables), " variables, n = ", x$n, "; ", nrow(s),
        " partitions examined", sep = "")
    if (x$kfix > 0L)
        cat(" (first ", x$kfix, " variable(s) kept together)", sep = "")
    cat("\n  depth = ", x$depth, "; nperm = ", x$nperm, sep = "")
    if (!is.null(x$nperm.screen)) cat(" (screening with ", x$nperm.screen, ")", sep = "")
    cat("\n\n")
    if (!any(s$Interesting)) {
        cat("No interesting partitions (both p-values >= ", x$SigProb, ")\n", sep = "")
    } else {
        cat("Interesting partitions (chi or lambda p-value < ", x$SigProb, "):\n\n", sep = "")
        tab <- s[s$Interesting, c("Partition", "Groups", "sigP.chi", "sigP.lam")]
        tab <- tab[order(pmin(tab$sigP.chi, tab$sigP.lam), tab$sigP.chi), ]
        names(tab) <- c("Partition", "Groups", "P(chi)", "P(lambda)")
        print(tab, row.names = FALSE, digits = 3)
    }
    nint <- sum(!s$Interesting)
    if (nint) cat("\n", nint, " partition(s) not significant\n", sep = "")
    invisible(x)
}

plot.chipart <- function(x, which = NULL, ask = NULL, ...) {
    if (!length(x$results)) {
        message("No interesting partitions to plot")
        return(invisible(x))
    }
    ids <- names(x$results)
    if (!is.null(which)) {
        ids <- intersect(ids, as.character(which))
        if (!length(ids)) stop("none of 'which' is an interesting partition", call. = FALSE)
    }
    if (is.null(ask)) ask <- grDevices::dev.interactive() && length(ids) > 1L
    if (ask) {
        oask <- grDevices::devAskNewPage(TRUE)
        on.exit(grDevices::devAskNewPage(oask))
    }
    for (id in ids) graphics::plot(x$results[[id]], ...)
    invisible(x)
}
