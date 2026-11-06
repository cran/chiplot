## Examples from N. I. Fisher & A. J. Lee, "Chi-plots for assessing multivariate
## dependence", using package chiplot.
## The simulated examples need no data files. For the real-data examples, set
## 'filesDir' to the folder containing LASERI.csv, Fisher.dji.d.csv,
## Fisher.nkx.d.csv, Fisher.spxd.d.csv, large.pva.raw.csv and abalone.csv.

library(chiplot)
filesDir <- "."           # ***** REDEFINE THIS DIRECTORY APPROPRIATELY

## multivariate normal sample (base R replacement for mvtnorm::rmvnorm)
rmvn <- function(n, sigma) matrix(rnorm(n * ncol(sigma)), n) %*% chol(sigma)

######################
## EXAMPLE 3.1  Two independent bivariate blocks
set.seed(161225)
rho1 <- 0.5; rho2 <- 0; rho3 <- 0.7
S1 <- matrix(c(1, rho1, rho1, 1), 2, 2); S2 <- matrix(rho2, 2, 2)
S3 <- matrix(c(1, rho3, rho3, 1), 2, 2)
sigma <- rbind(cbind(S1, S2), cbind(S2, S3))
Data <- data.frame(rmvn(200, sigma)); colnames(Data) <- paste0("Var.", 1:4)
resEx1 <- chimult(Data, D = c(2, 2), nperm = 1000, Nsim = 100)
plot(resEx1, pt = 12, guide.option = 1, mar = c(4, 4, 2, 2))

## EXAMPLE 3.2  Two correlated bivariate blocks
set.seed(181225)
S2 <- matrix(0.4, 2, 2)
sigma <- rbind(cbind(S1, S2), cbind(S2, S3))
Data <- data.frame(rmvn(200, sigma)); colnames(Data) <- paste0("Var.", 1:4)
resEx2 <- chimult(Data, D = c(2, 2), nperm = 1000, Nsim = 100)
plot(resEx2, pt = 12, guide.option = 4, mar = c(4, 4, 2, 2))

## EXAMPLE 3.3A  Large off-centre hole in the data
set.seed(280926)
NN <- 500000; P <- 5
dum <- matrix(rnorm(NN * P), NN, P)
keep <- sqrt(rowSums((dum - 1)^2)) >= 5          # outside radius 5 of (1, ..., 1)
data <- dum[keep, ]
Data <- data.frame(data[sample(nrow(data), 200), ]); colnames(Data) <- paste0("Var.", 1:5)
resEx3A <- chimult(Data, D = rep(1, P), nperm = 1000, Nsim = 100)
plot(resEx3A, pt = 10, guide.option = 4, mar = c(4, 4, 2, 2))

## EXAMPLE 3.3B  Spherical central hole, p = 10
sample_outside_sphere <- function(n, p, D) {
    log_Sd <- pchisq(D^2, df = p, lower.tail = FALSE, log.p = TRUE)
    R <- sqrt(qchisq(log_Sd + log(runif(n)), df = p, lower.tail = FALSE, log.p = TRUE))
    Z <- matrix(rnorm(n * p), n, p)
    Z / sqrt(rowSums(Z^2)) * R
}
set.seed(201026)
Data <- data.frame(sample_outside_sphere(500, 10, 2.5)); colnames(Data) <- paste0("Var.", 1:10)
resEx3B <- chimult(Data, D = rep(1, 10), nperm = 1000, Nsim = 100)
plot(resEx3B, pt = 12, guide.option = 3, mar = c(4, 4, 2, 2))

## EXAMPLE 3.4  Three groups, one correlation between groups 1 and 3
set.seed(201225)
rho1 <- 0.5; rho2 <- 0; rho3 <- 0.7; rho4 <- 0.4; rho5 <- 0.5
S11 <- matrix(rho1, 3, 3); diag(S11) <- 1
S22 <- matrix(c(1, rho3, rho3, 1), 2, 2); S33 <- matrix(c(1, rho4, rho4, 1), 2, 2)
sigma <- rbind(cbind(S11, matrix(rho2, 3, 2), matrix(rho2, 3, 2)),
               cbind(matrix(rho2, 2, 3), S22, matrix(rho2, 2, 2)),
               cbind(matrix(rho2, 2, 3), matrix(rho2, 2, 2), S33))
sigma[2, 6] <- sigma[6, 2] <- rho5
Data <- data.frame(rmvn(200, sigma)); colnames(Data) <- paste0("Var.", 1:7)
resEx4 <- chimult(Data, D = c(3, 2, 2), nperm = 1000, Nsim = 100)
plot(resEx4, pt = 12, guide.option = 4, mar = c(4, 4, 2, 2))

######################
## Real data examples

## LASERI data
LASERI.full <- read.csv(file.path(filesDir, "LASERI.csv"))
LASERI <- LASERI.full[, c(2, 3, 4, 5, 6, 9, 10, 11, 14, 15, 18, 29, 22, 23, 8)]
res5.1 <- chimult(LASERI, D = c(7, 8), nperm = 1000, Nsim = 100)
plot(res5.1, pt = 9, guide.option = 3, mar = c(4, 4, 2, 2))

## Stock market indices
select <- c(1:5, 41:45, 81:85, 121:125, 161:165, 201:204)
weekdays5 <- function(file, prefix) {
    d <- read.csv(file.path(filesDir, file))[select, ]
    out <- data.frame(matrix(unlist(d), ncol = 5, byrow = TRUE))
    colnames(out) <- paste0(prefix, c("Mon", "Tue", "Wed", "Thu", "Fri"))
    out
}
All.df <- data.frame(weekdays5("Fisher.dji.d.csv", "DJIA"),
                     weekdays5("Fisher.spxd.d.csv", "SP500"),
                     weekdays5("Fisher.nkx.d.csv", "NIKAI"))
res2A <- chimult(All.df, D = c(5, 5, 5), nperm = 1000, Nsim = 100)
plot(res2A, pt = 10, guide.option = 3, mar = c(4, 4, 2, 2))
res2B <- chimult(All.df, D = c(10, 5), nperm = 1000, Nsim = 100)
plot(res2B, pt = 10, guide.option = 3, mar = c(4, 4, 2, 2))

## Survey data: ordered categorical responses, jittered to break ties
large.pva.raw <- read.csv(file.path(filesDir, "large.pva.raw.csv"), header = FALSE)
names(large.pva.raw) <- paste0("A", seq_len(ncol(large.pva.raw)))
Data <- large.pva.raw[, c(50:52, 40:45)]
res3 <- chimult(Data, D = c(3, 5, 1), nperm = 1000, Nsim = 100, jitter = TRUE)
plot(res3, pt = 12, guide.option = 4, mar = c(4, 4, 2, 2))
## (The original script added 0.01 * rnorm() noise after median-centring;
##  jitter = TRUE adds uniform noise within +/- 0.25 of the gap between categories.)

## Abalone data: females, one of the more interesting partitions
abalone.raw <- read.csv(file.path(filesDir, "abalone.csv"))
females.raw <- abalone.raw[abalone.raw[, 1] == "F", -1]
Data <- females.raw[, c(1, 8, 2, 3, 4)]
res4 <- chimult(Data, D = c(2, 1, 1, 1), nperm = 1000, Nsim = 100)
plot(res4, pt = 12, guide.option = 4, mar = c(4, 4, 2, 2))

## All partitions of the five abalone variables (51 partitions)
search <- chipart(Data, nperm = 200, nperm.screen = 50)
search
plot(search)
