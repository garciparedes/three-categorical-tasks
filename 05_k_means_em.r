##
## Title: Implementación del Algorimo EM para la resolución del problema de
##        las K-Medias.
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

# Número de Observaciones.
n <- 1000

# Número de Mixturas.
k <- 3

# Proporción de las mixturas.
p <- c(0.1, 0.6, 0.3)  # por simplicidad sum(p) = 1

# Parámetros para la generación de la muestra.
theta <- matrix(c(3,  4,
                  8,  4,
                  9, 10),
                k, d, byrow = TRUE)

# Muestra.
x <- sapply(1:k, function(i){
      rmvnorm(n * p[i], theta[i, ])
  }) %>%
  { do.call(rbind, .) }

# Diagrama de Dispersión de la muestra.
plot(x)


OptimizeTheta <- function(x) {
  # Búsqueda del punto cerntral de los datos de entrada.
  colMeans(x)
}


CalculateB <- function(x, theta) {
  # Indica para qué mixtura se maximiza la verosimilitud perfil en cada
  # observación.
  apply(x, 1, function(obs) {
    sapply(1:nrow(theta), function(i) {
      - sqrt(sum((obs - theta[i,]) ^ 2))
    }) %>%
    which.max()
  })
}


RandomTheta<- function(x, k, d) {
  # Generación de k puntos aleatoriamente en el rango de los datos.
  # Se utiliza como puntos de partida para resolver el problema de
  # las K-Medias.

  means <- (1:d) %>%
    sapply(function(i) {
      sort(runif(k, min = min(x[, i]), max = max(x[, i])))
    })

  return(means)
}


OptimizeKMeansEM <- function(x, k, d) {
  # Calcula los estimadores máximos verosímiles mediante el Algorimo EM.

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

# Buscamos los EMVs.
# Nota: Este método puede llamarse de manera iterativa, para conseguir
#       distintas soluciones y después promediar los resultados. Sin embargo,
#       en los casos de prueba ha bastado con una única iteración para
#       encontrar el óptimo por lo que se ha decidido no implementar.
(opt.em <- OptimizeKMeansEM(x, k, d))

# Diagrama de la muestra con las categorías obtenidas por el EMV asignadas.
plot(x, col = opt.em$group + 1)
