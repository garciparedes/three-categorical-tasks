## Author: Sergio García Prado

rm(list = ls())

library(magrittr)

k <- 3
n <- 1000
p <- c(0.3, 0.4)
lambda <- c(1, 5, 10)
x1 <- rpois(n * p[1], lambda[1])
x2 <- rpois(n * p[2], lambda[2])
x3 <- rpois(n * (1 - sum(p)), lambda[3])
x <- c(x1, x2, x3)

LogLikeliHood <- function(p, lambda, x) {
  sum(log(p * dpois(x, lambda[1:(k - 1)]) + (1 - sum(p)) * dpois(x, lambda[k])))
}

LogLikeliHoodTheta <- function(theta, x) {
  p <- theta[1:(k - 1)]
  lambda <- theta[k:(2 * k - 1)]
  logL <- LogLikeliHood(p, lambda, x)
  return(logL)
}

NegativeLogLikelihoodTheta <- function(...) {
  - LogLikeliHoodTheta(...)
}

# theta.initial <- c(0.3, 0.4, 1, 50, 100)
# opt <- optim(theta.initial, NegativeLogLikelihoodTheta, x = x, method= "L-BFGS-B",
#              lower = rep(10e-4, (k - 1) + k), upper = c(rep(1 - 10e-4, k -1), rep(Inf, k)))
#
# (theta.hat <- opt$par)


LogLikeliHoodPoisson <- function(lambda, x) {
  # n <- length(x)
  # logL <- sum(- n * log(lambda) - sum(x) / lambda)
  logL <- sum(dpois(x, lambda, log = TRUE))
  return(logL)
}

NegativeLogLikeliHoodPoisson <- function(...) {
  - LogLikeliHoodPoisson(...)
}


OptimizeLambda <- function(x) {
  # opt <- optim(runif(1), NegativeLogLikeliHoodPoisson, x = x,
  #              lower = 10e-4, upper = Inf, method = "L-BFGS-B")
  # lambda.hat <- opt$par
  
  lambda.hat <- mean(x)
  return(lambda.hat)
}

CalculateB <- function(x, lambda) {
  sapply(x, function(obs) {
    dpois(obs, lambda, log = TRUE) %>%
    which.max()
  })
}


OptimizePoissonMixtureEM <- function(x, k) {
  lambda <- sort(runif(k, min = min(x), max = max(x)))
  b <- CalculateB(x, lambda)
  b.old <- rep(0, n)
  i <- 0

  # Fijamos criterio de convergencia en la estabilización de etiquetas de la mixtura.
  # Ya que cuando esta se estabiliza, los EMV de lambda no variarán.
  while (any(b != b.old)) {
    i <- i + 1
    lambda <- sapply(1:k, function(j) {
      OptimizeLambda(x[b == j])
    })
    b.old <- b
    b <- CalculateB(x, lambda)
  }
  result <- list(lambda = lambda,
                 group = b,
                 steps = i,
                 ratios = prop.table(table(b)))
  return(result)
}

## EM Algorithm
iters <- 100

cat("Running EM...\n")
lambda <- matrix(rep(0, iters * k), iters, k)
ratios <- matrix(rep(0, iters * k), iters, k)
for (r in 1:iters) {
  opt <- OptimizePoissonMixtureEM(x, k)
  lambda[r, ] <- opt$lambda
  ratios[r, ] <- opt$ratios
  cat('Iteration: ',r, ', Steps: ', opt$steps, ', Lambda: ',  lambda[r, ], ', p: ', ratios[r, ],'\n')
}

cat("Finished EM. \n")
cat("EMV: \n")
cat("Lambda: ",colMeans(lambda, na.rm = TRUE), '\n')
cat("p: ", colMeans(ratios, na.rm = TRUE), '\n')
