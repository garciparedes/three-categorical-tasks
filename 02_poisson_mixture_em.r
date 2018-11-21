##
## Title: Mixtura de k distribuciones de Poisson.
##
## Author: Sergio García Prado
##
## Date: Noviembre de 2018
##

rm(list = ls())

library(magrittr)

# Número de Mixturas.
k <- 5

# Tamaño muestral.
n <- 1000

# Pesos de la Mixtura.
p <- c(0.2, 0.2, 0.2, 0.2, 0.2)  # por simplicidad sum(p) = 1

# Parámetros Distribucionales.
lambda <- c(1, 5, 10, 15, 20)


# Muestra.
x <- (1:k) %>%
  sapply(function(i) {
    rpois(n * p[i], lambda[i])
  }) %>%
  unlist()


LogLikeliHoodPoisson <- function(lambda, x) {
  # Método manual.
  # n <- length(x)
  # logL <- sum(- n * log(lambda) - sum(x) / lambda)

  # Método basado en la función de densidad implementada en R.
  logL <- sum(dpois(x, lambda, log = TRUE))
  return(logL)
}


NegativeLogLikeliHoodPoisson <- function(...) {
  # Opuesto de la LogVerosimilitud (para poder optimizarse con optim()).

  - LogLikeliHoodPoisson(...)
}


OptimizeLambda <- function(x) {
  # Optimización del EMV para los parámetros de la distribución.

  # Búsqueda iterativa del EMV.
  # opt <- optim(runif(1), NegativeLogLikeliHoodPoisson, x = x,
  #              lower = 10e-4, upper = Inf, method = "L-BFGS-B")
  # lambda.hat <- opt$par

  # Búsqueda analítica.
  lambda.hat <- mean(x)
  return(lambda.hat)
}


CalculateB <- function(x, lambda) {
  # Indica para qué mixtura se maximiza la verosimilitud perfil en cada
  # observación.

  sapply(x, function(obs) {
    dpois(obs, lambda, log = TRUE) %>%
    which.max()
  })
}


OptimizePoissonMixtureEM <- function(x, k) {
  # Calcula los estimadores máximos verosímiles mediante el Algorimo EM.

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

## Algoritmo EM.

# Iteraciones.
iters <- 10

cat("Running EM...\n")
lambda <- matrix(rep(0, iters * k), iters, k)
ratios <- matrix(rep(0, iters * k), iters, k)
for (r in 1:iters) {
  opt.em <- OptimizePoissonMixtureEM(x, k)
  lambda[r, ] <- opt.em$lambda
  ratios[r, ] <- opt.em$ratios
  cat('Iteration: ',r, ', Steps: ', opt.em$steps, ', Lambda: ',  lambda[r, ], ', p: ', ratios[r, ],'\n')
}

## Resultados.

cat("Finished EM. \n")
cat("EMV: \n")
cat("Lambda: ",colMeans(lambda, na.rm = TRUE), '\n')
cat("p: ", colMeans(ratios, na.rm = TRUE), '\n')
