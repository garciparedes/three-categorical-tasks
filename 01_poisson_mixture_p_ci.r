##
## Title: Intervalo de Confianza basado en el estadístico Razón de
##        Verosimilitud para una Mixtura de 2 distribuciones de Poisson.
##
## Author: Sergio García Prado
##
## Date: Noviembre de 2018
##

rm(list = ls())


# Tamaño muestral.
n <- 1000

# Proporción de la mixtura.
p <- 0.7

# Parámetros de Poisson.
lambda1 <- 2
lambda2 <- 5

# Muestra.
x1 <- rpois(n * p, lambda1)
x2 <- rpois(n * (1 - p), lambda2)
x <- c(x1, x2)


LogLikeliHood <- function(p, lambda1, lambda2, x) {
  sum(log(p * dpois(x, lambda1) + (1 - p) * dpois(x, lambda2)))
}


LogLikeliHoodTheta <- function(theta, x) {
  p <- theta[1]
  lambda1 <- theta[2]
  lambda2 <- theta[3]

  logL <- LogLikeliHood(p, lambda1, lambda2, x)
  return(logL)
}


NegativeLogLikelihoodTheta <- function(...) {
  # Opuesto de la LogVerosimilitud (para poder optimizarse con optim()).

  - LogLikeliHoodTheta(...)
}


theta.initial <- c(runif(1), runif(2) * 10)
if (theta.initial[3] < theta.initial[2]) {
  # Para resolver la identificabilidad forzamos a que el valor lambda de la
  # primera mixtura tenga un valor menor que el de la segunda.
  temp <- theta.initial[2]
  theta.initial[2] <- theta.initial[3]
  theta.initial[3] <- temp
}


# Búsqueda de los EMVs de manera iterativa.
opt <- optim(theta.initial, NegativeLogLikelihoodTheta, x = x, method= "L-BFGS-B",
             lower = c(10e-4, 10e-4, 10e-4), upper = c(1 - 10e-4, Inf, Inf))


# EMVs.
(theta.hat <- opt$par)


f <- function(p, theta, x, alpha) {
  # Función auxiliar para el cálculo del punto donde se iguala el valor
  # chi-cuadrado a nivel 1 - alpha con el valor observado.
  2 * (LogLikeliHoodTheta(theta, x) - LogLikeliHood(p, theta[2], theta[3], x)) - qchisq(1 - alpha, df = 1)
}

# Nivel Alpha.
alpha <- 0.05

# Confidence Intervals
c(uniroot(f, c(0, theta.hat[1]), theta = theta.hat, x = x, alpha = alpha)$root,
  uniroot(f, c(theta.hat[1], 1), theta = theta.hat, x = x, alpha = alpha)$root)
