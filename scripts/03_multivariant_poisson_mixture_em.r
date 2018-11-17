## Author: Sergio García Prado

rm(list = ls())

library(magrittr)

# Dimensions
d <- 2

# Mixtures
k <- 3

# Cases
n <- 1000

#
p.x <- c(0.3, 0.4)
lambda.x <- c(1, 5, 10)
x1 <- rpois(n * p.x[1], lambda.x[1])
x2 <- rpois(n * p.x[2], lambda.x[2])
x3 <- rpois(n * (1 - sum(p.x)), lambda.x[3])

p.y <- c(0.1, 0.6)
lambda.y <- c(10, 50, 70)

y1 <- rpois(n * p.y[1], lambda.y[1])
y2 <- rpois(n * p.y[2], lambda.y[2])
y3 <- rpois(n * (1 - sum(p.y)), lambda.y[3])



data <- matrix(c(x1, x2, x3, y1, y2, y3),
               n, d)

data
