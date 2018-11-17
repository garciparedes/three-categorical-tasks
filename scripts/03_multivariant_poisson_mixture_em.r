## Author: Sergio García Prado

rm(list = ls())

library(magrittr)

# Dimensions
d <- 2

# Mixtures
k <- 3

# Cases
n <- 1000

p <- matrix(c(0.1, 0.5,
              0.6, 0.2,
              0.3, 0.3),
            3, 2, byrow = TRUE)

lambda <- matrix(c( 1, 10,
                    5, 50,
                   10, 70),
                 3, 2, byrow = TRUE)

x <- sapply(1:d, function(i) {
    sapply(1:k, function(j){
      rpois(n * p[j, i], lambda[j, i])
    }) %>%
    unlist()
  }) %>%
  unlist() %>%
  matrix(n, d)


LogLikeliHoodPoisson <- function(lambda, x) {
 logL <- sapply(1:length(lambda), function(i) {
      sum(dpois(x[, i], lambda[i], log = TRUE))
    }) %>%
    sum()
 return(logL)
}

NegativeLogLikeliHoodPoisson <- function(...) {
 - LogLikeliHoodPoisson(...)
}

OptimizeLambda <- function(x) {
  # opt <- optim(runif(2), NegativeLogLikeliHoodPoisson, x = x,
  #              lower = 10e-4, upper = Inf, method = "L-BFGS-B")
  # lambda.hat <- opt$par
  lambda.hat <- colMeans(x)
  return(lambda.hat)
}

CalculateB <- function(x, lambda) {
  apply(x, 1, function(obs) {
    sapply(1:nrow(lambda), function(i) {
         sum(dpois(obs, lambda[i, ], log = TRUE))
      }) %>% which.max()
  })
}


OptimizeMultiVariantPoissonMixtureEM <- function(x, k, d) {
  lambda <- (1:d) %>%
    sapply(function(i) {
      sort(runif(k, min = min(x[, i]), max = max(x[, i])))
    })

  b <- CalculateB(x, lambda)
  b.old <- rep(0, n)

  i <- 0
  while (any(b != b.old)) {
    # Fijamos criterio de convergencia en la estabilización de etiquetas de la mixtura.
    # Ya que cuando esta se estabiliza, los EMV de lambda no variarán.

    i <- i + 1
    lambda <- sapply(1:k, function(j) {
      OptimizeLambda(x[b == j, ])
    }) %>% t()

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
                 ratios = p
               )
  return(result)
}

a <- OptimizeMultiVariantPoissonMixtureEM(x, k, d)


plot(x, col = a$group + 1)
