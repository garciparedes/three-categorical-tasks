##
## Title: Mixtura de k distribuciones de Normales multivariantes con d
##        dimensiones.
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
library(mvtnorm)


# Número de Dimensiones.
d <- 2

# Número de Mixturas.
k <- 3

# Tamaño Muestral.
n <- 1000

# Pesos de la Mixtura.
p <- c(0.1, 0.6, 0.3)  # por simplicidad sum(p) = 1

# Parámetro de la muestra.
theta <- list(list(mean = c(3,  4),
                   sigma = matrix(c(2, 1,
                                    1, 2),  2, 2, byrow = TRUE)),
              list(mean = c(8,  4),
                  sigma = matrix(c( 2, -1,
                                   -1, 2), 2, 2, byrow = TRUE)),
              list(mean = c(9, 10),
                  sigma = matrix(c(1, 0,
                                   0, 1), 2, 2, byrow = TRUE)))

# Muestra.
x <- sapply(1:k, function(j){
      rmvnorm(n * p[j], theta[[j]]$mean, theta[[j]]$sigma)
  }) %>%
  { do.call(rbind, .) }

# Diagrama de Dispersión de la muestra.
plot(x)


OptimizeTheta <- function(x) {
  # Optimización del EMV para los parámetros de la distribución.

  # Búsqueda Analítica
  theta.hat <- list(mean = colMeans(x),
                    sigma = cov(x))
  return(theta.hat)
}


CalculateB <- function(x, theta) {
  # Indica para qué mixtura se maximiza la verosimilitud perfil en cada
  # observación.
  apply(x, 1, function(obs) {
    sapply(1:length(theta), function(i) {
      dmvnorm(obs, theta[[i]]$mean, theta[[i]]$sigma, log = TRUE)
    }) %>% which.max()
  })
}


RandomTheta<- function(x, k, d) {
  # Generación de parámetros iniciales de las k distribuciones.

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
  # Calcula los estimadores máximos verosímiles mediante el Algorimo EM.

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


# Buscamos los EMVs.
# Nota: Este método puede llamarse de manera iterativa, para conseguir
#       distintas soluciones y después promediar los resultados. Sin embargo,
#       en los casos de prueba ha bastado con una única iteración para
#       encontrar el óptimo por lo que se ha decidido no implementar.
(opt.em <- OptimizeMultiVariantNormalMixtureEM(x, k, d))

# Diagrama de la muestra con las categorías obtenidas por el EMV asignadas.
plot(x, col = opt.em$group + 1)
