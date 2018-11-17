## Author: Sergio García Prado

rm(list = ls())

library(magrittr)

k <- 5
n <- 1000
p <- c(0.2, 0.2, 0.2, 0.2, 0.2)  # por simplicidad sum(p) = 1
lambda <- c(1, 5, 10, 15, 20)
x <- (1:k) %>%
  sapply(function(i) {
    rpois(n * p[i], lambda[i])
  }) %>%
  unlist()


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

  p <- (1:k) %>%
    sapply(function(i) {
      mean(b == i)
    })

  result <- list(lambda = lambda,
                 group = b,
                 steps = i,
                 ratios = p)
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
