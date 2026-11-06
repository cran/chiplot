## Internal helpers shared by chiplot(), chipairs(), chimult() and chipart().

## Interfaces to the compiled routines in src/ (see src/chicore.c)
.chivals2_native <- function(x, y) .Call(C_chivals2, as.double(x), as.double(y))
.chimult_native <- function(X, D) .Call(C_chimult, X, as.integer(D))

## Convert Data (matrix or data frame) to a numeric matrix with column names.
## Ordered factors become their integer codes; their column numbers are
## returned in attribute "ordered" so that they can be jittered.
.as_numeric_matrix <- function(Data, prefix = "Var.") {
    ordered.cols <- integer(0)
    if (is.data.frame(Data)) {
        ok <- vapply(Data, function(v) is.numeric(v) || is.ordered(v) || is.logical(v),
                     logical(1))
        if (!all(ok))
            stop("all columns of 'Data' must be numeric or ordered factors; not: ",
                 paste(names(Data)[!ok], collapse = ", "), call. = FALSE)
        ordered.cols <- which(vapply(Data, is.ordered, logical(1)))
        Data[] <- lapply(Data, function(v) if (is.ordered(v)) as.numeric(as.integer(v))
                                           else as.numeric(v))
        X <- as.matrix(Data)
    } else {
        X <- as.matrix(Data)
        if (!is.numeric(X) && !is.logical(X))
            stop("'Data' must be a numeric matrix or a data frame", call. = FALSE)
    }
    storage.mode(X) <- "double"
    if (is.null(colnames(X)))
        colnames(X) <- paste0(prefix, seq_len(ncol(X)))
    attr(X, "ordered") <- unname(ordered.cols)
    X
}

## Columns selected by a 'jitter' specification: TRUE/FALSE, a logical vector,
## column numbers or column names.
.jitter_cols <- function(jitter, X) {
    p <- ncol(X)
    if (is.null(jitter) || identical(jitter, FALSE)) return(integer(0))
    if (isTRUE(jitter)) return(seq_len(p))
    if (is.logical(jitter)) {
        if (length(jitter) != p)
            stop("a logical 'jitter' must have one element per variable", call. = FALSE)
        return(which(jitter))
    }
    if (is.character(jitter)) {
        m <- match(jitter, colnames(X))
        if (anyNA(m))
            stop("unknown variable(s) in 'jitter': ",
                 paste(jitter[is.na(m)], collapse = ", "), call. = FALSE)
        return(m)
    }
    j <- as.integer(jitter)
    if (any(is.na(j) | j < 1L | j > p))
        stop("'jitter' column numbers must lie between 1 and ", p, call. = FALSE)
    unique(j)
}

## Break ties in the selected columns by adding uniform noise of at most
## 'amount' times the smallest gap between distinct values, so that the
## ordering of distinct values (e.g. ordered categories) is preserved.
.jitter_columns <- function(X, cols, amount = 0.25) {
    if (!length(cols)) return(X)
    if (!is.numeric(amount) || length(amount) != 1L || amount <= 0 || amount >= 0.5)
        stop("'jitter.amount' must be a number in (0, 0.5)", call. = FALSE)
    n <- nrow(X)
    for (j in cols) {
        u <- sort(unique(X[, j]))
        gap <- if (length(u) > 1L) min(diff(u)) else 1
        X[, j] <- X[, j] + stats::runif(n, -amount * gap, amount * gap)
    }
    X
}

## A default point size: as large as possible without crowding.
.auto_cex <- function(n, scale = 1) {
    scale * min(1, max(0.25, sqrt(150 / max(n, 1))))
}

## P-values printed to three decimal places, as in the original code.
.format_p <- function(p) {
    if (is.na(p)) return("NA")
    formatC(p, format = "f", digits = 3)
}

## Validate a vector of block sizes D for p variables.
.check_D <- function(D, p) {
    if (!is.numeric(D) || anyNA(D) || any(D < 1) || any(D != round(D)))
        stop("'D' must be a vector of positive whole numbers", call. = FALSE)
    D <- as.integer(D)
    if (length(D) < 2L)
        stop("'D' must define at least two blocks of variables", call. = FALSE)
    if (sum(D) != p)
        stop("the block sizes in 'D' sum to ", sum(D), " but 'Data' has ", p,
             " columns", call. = FALSE)
    D
}

## Variable names in each block: one string per block, e.g. "Var.1, Var.2".
groupText <- function(res, sep = ", ") {
    cD <- c(0L, cumsum(res$D))
    vapply(seq_along(res$D), function(k0) {
        paste(colnames(res$X)[(cD[k0] + 1L):cD[k0 + 1L]], collapse = sep)
    }, character(1))
}

## All blocks in one string, e.g. "Var.1,Var.2 Versus Var.3,Var.4".
groupLabels <- function(res) {
    paste(groupText(res, sep = ","), collapse = " Versus ")
}

## Split a string into lines of at most 'width' characters, breaking after
## commas where possible (as in the original writeGroupNames()).
.wrap_at_commas <- function(s, width) {
    out <- character(0)
    while (nchar(s) > 0L) {
        if (nchar(s) < width) {
            out <- c(out, s)
            break
        }
        commas <- gregexpr(",", s, fixed = TRUE)[[1L]]
        commas <- commas[commas > 0L]
        fits <- commas[commas <= width]
        end <- if (length(fits)) max(fits)
               else if (length(commas)) min(commas)
               else nchar(s)
        out <- c(out, substr(s, 1L, end))
        s <- substring(s, end + 1L)
    }
    out
}

## Draw the panel listing the variables in each group (top left of a chi-plot).
writeGroupNames <- function(res, pt = 12, box = TRUE) {
    k <- length(res$D)
    group.text <- paste("Group", seq_len(k), ": ", groupText(res), sep = "")
    ptSize <- match(pt, c(12, 11, 10, 9))
    if (is.na(ptSize)) ptSize <- 3L
    lineLength <- c(20, 22, 25, 29)[ptSize]
    lines <- lapply(group.text, .wrap_at_commas, width = lineLength)
    cex <- pt / graphics::par("ps")
    collapseChrs <- if (k < 8) "\n\n" else "\n"
    myText <- paste(vapply(lines, paste, character(1), collapse = "\n"),
                    collapse = collapseChrs)
    graphics::plot(c(0, 1), c(0, 1), axes = FALSE, type = "n", xlab = "", ylab = "")
    if (box) graphics::box()
    graphics::text(0.05, 0.95, myText, adj = c(0, 1), family = "mono", cex = cex)
    invisible(myText)
}
