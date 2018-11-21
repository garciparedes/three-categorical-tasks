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

theta <- matrix(c(3,  4,
                  8,  4,
                  9, 10),
                k, d, byrow = TRUE)

x <- sapply(1:k, function(i){
      rmvnorm(n * p[i], theta[i, ])
  }) %>%
  { do.call(rbind, .) }

plot(x)


OptimizeTheta <- function(x) {
  colMeans(x)
}

CalculateB <- function(x, theta) {
  apply(x, 1, function(obs) {
    sapply(1:nrow(theta), function(i) {
      - sqrt(sum((obs - theta[i,]) ^ 2))
    }) %>%
    which.max()
  })
}

RandomTheta<- function(x, k, d) {
  means <- (1:d) %>%
    sapply(function(i) {
      sort(runif(k, min = min(x[, i]), max = max(x[, i])))
    })

  return(means)
}

OptimizeKMeansEM <- function(x, k, d) {
  theta <- RandomTheta(x, k, d)
  b <- CalculateB(x, theta)
  b.old <- rep(0, n)

  i <- 0
  while (any(b != b.old)) {
    # Fijamos criterio de convergencia en la estabilización de etiquetas de la mixtura.
    # Ya que cuando esta se estabiliza, los EMV de theta no variarán.

    i <- i + 1
    theta <- (1:k) %>%
      sapply(function(j) {
        OptimizeTheta(x[b == j, ])
      }) %>%
      t()

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

(opt.em <- OptimizeKMeansEM(x, k, d))

plot(x, col = opt.em$group + 1)
