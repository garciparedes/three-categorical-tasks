## Author: Sergio García Prado

rm(list = ls())

library(magrittr)
library(mvtnorm)
# Dimensions
d <- 2

# Mixtures
k <- 3

# Cases
n <- 1000

p <- c(0.1,
       0.6,
       0.3)

theta <- list(list(mean = c(3,  4),
                sigma = matrix(c(2, 1,
                                 1, 2),  2, 2, byrow = TRUE)),
           list(mean = c(8,  4),
                sigma = matrix(c( 2, -1,
                               -1, 2), 2, 2, byrow = TRUE)),
           list(mean = c(9, 10),
                sigma = matrix(c(1, 0,
                               0, 1), 2, 2, byrow = TRUE)))

x <- sapply(1:k, function(j){
      rmvnorm(n * p[j], theta[[j]]$mean, theta[[j]]$sigma)
  }) %>%
  { do.call(rbind, .) }

plot(x)


OptimizeTheta <- function(x) {
  theta.hat <- list()
  theta.hat$mean <- colMeans(x)
  theta.hat$sigma <- cov(x)
  return(theta.hat)
}


CalculateB <- function(x, theta) {
  apply(x, 1, function(obs) {
    sapply(1:length(theta), function(i) {
      dmvnorm(obs, theta[[i]]$mean, theta[[i]]$sigma, log = TRUE)
    }) %>% which.max()
  })
}


RandomTheta<- function(x, k, d) {
  means <- (1:d) %>%
    sapply(function(i) {
      sort(runif(k, min = min(x[, i]), max = max(x[, i])))
    })

  theta <- apply(means,1, function(m) {
    list(mean = m, sigma = diag(d))
  })

  return(theta)
}

OptimizeMultiVariantNormalMixtureEM <- function(x, k, d) {
  theta <- RandomTheta(x, k, d)
  b <- CalculateB(x, theta)
  b.old <- rep(0, n)

  i <- 0
  while (any(b != b.old)) {
    # Fijamos criterio de convergencia en la estabilización de etiquetas de la mixtura.
    # Ya que cuando esta se estabiliza, los EMV de theta no variarán.

    i <- i + 1
    theta <- lapply(1:k, function(j) {
      OptimizeTheta(x[b == j, ])
    })

    b.old <- b
    b <- CalculateB(x, theta)
  }

  p <- (1:k) %>%
   sapply(function(i) {
     mean(b == i)
   })

  result <- list(theta = theta,
                 group = b,
                 steps = i,
                 ratios = p)
  return(result)
}

(opt.em <- OptimizeMultiVariantNormalMixtureEM(x, k, d))

plot(x, col = opt.em$group + 1)
