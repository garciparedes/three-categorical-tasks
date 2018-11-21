##
## Title: Mixtura de k distribuciones de multivariantes con d dimensiones,
##        siguiendo distribuciones de Poisson en las marginales.
##
## Author: Sergio García Prado
##
## Date: Noviembre de 2018
##
## Notas: Se ha elegido fijar la dimensionalidad en 2 para poder hacer
##        representaciones gráficas de manera sencilla. Sin embargo,
##        es posible hacer probar el funcionamiento en más dimensiones
##        variando dichos parámetros.
##

rm(list = ls())

library(magrittr)


# Número de Dimensiones.
d <- 2

# Número de Mixturas.
k <- 3

# Número de Observaciones en la muestra.
n <- 1000

# Pesos de la Mixtura.
p <- c(0.1, 0.6, 0.3)  # por simplicidad sum(p) = 1

# Parámetros de las distribuciones de Poisson.
lambda <- matrix(c( 1, 10,
                    5, 50,
                   10, 70),
                 3, 2, byrow = TRUE)

# Muestra.
x <-
  sapply(1:d, function(j) {
    sapply(1:k, function(i){
      rpois(n * p[i], lambda[i, j])
    }) %>%
    unlist()
  }) %>%
  unlist() %>%
  matrix(n, d)

# Diagrama de la Muestra
plot(x)


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
  # Optimización del EMV para los parámetros de la distribución.

  # Búsqueda Iterativa.
  # opt <- optim(runif(2), NegativeLogLikeliHoodPoisson, x = x,
  #              lower = 10e-4, upper = Inf, method = "L-BFGS-B")
  # lambda.hat <- opt$par

  # Búsqueda analítica.
  lambda.hat <- colMeans(x)
  return(lambda.hat)
}

CalculateB <- function(x, lambda) {
  # Indica para qué mixtura se maximiza la verosimilitud perfil en cada
  # observación.

  apply(x, 1, function(obs) {
    sapply(1:nrow(lambda), function(i) {
         sum(dpois(obs, lambda[i, ], log = TRUE))
      }) %>% which.max()
  })
}


OptimizeMultiVariantPoissonMixtureEM <- function(x, k, d) {
  # Calcula los estimadores máximos verosímiles mediante el Algorimo EM.

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
                 ratios = p)
  return(result)
}


# Buscamos los EMVs.
# Nota: Este método puede llamarse de manera iterativa, para conseguir
#       distintas soluciones y después promediar los resultados. Sin embargo,
#       en los casos de prueba ha bastado con una única iteración para
#       encontrar el óptimo por lo que se ha decidido no implementar.
(opt.em <- OptimizeMultiVariantPoissonMixtureEM(x, k, d))

# Diagrama de la muestra con las categorías obtenidas por el EMV asignadas.
plot(x, col = opt.em$group + 1)
