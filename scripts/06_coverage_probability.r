## Author: Sergio García Prado

rm(list = ls())

library(magrittr)

WaldConfidenceInterval <- function(y, n, alpha) {
  p.hat <- y / n
  z <- qnorm(1 - alpha / 2)
  confidence.interval <- p.hat + c(-1, 1) * z * sqrt(p.hat * (1 - p.hat) /n )

  return(confidence.interval)
}


ScoreConfidenceInterval <- function(y, n, alpha) {
  p.hat <- y / n
  z <- qnorm(1 - alpha / 2)

  pivot <- (p.hat + (z ^ 2) / (2 * n)) / (1 + (z ^ 2) / n)
  bound <- z / (1 + (z ^ 2) / n) * sqrt((p.hat * (1 - p.hat) / n) + ((z ^ 2) / (4 * n ^ 2)))

  confidence.interval <- pivot + c(-1, 1) * bound
  return(confidence.interval)
}


LogLikelihood <- function(p, y, n) {
  y * log(p) + (n - y) * log(1 - p)
}


LikelihoodRatioConfidenceInterval <- function(y, n, alpha) {
  confidence.interval <- rep(0, 2)

  p.hat <- y / n

  f <- function(p) {
    2 * (LogLikelihood(p.hat, y, n) - LogLikelihood(p, y, n)) - qchisq(1 - alpha, df = 1)
  }

  # obligatory tryCatch() block to fix y = 0 or y = n cases.
  confidence.interval <- tryCatch({
    c(uniroot(f, c(0,   p.hat))$root,
      uniroot(f, c(p.hat, 1))$root)
  }, error=function(e) {rep(NA, 2)})

  return(confidence.interval)
}


Amplitude <- function(ci) {
  mean(ci[, 2] - ci[, 1], na.rm = TRUE)
}


Coverage <- function(ci, p) {
  mean(ci[, 1] <= p & p <= ci[, 2], na.rm = TRUE)
}


alpha <- 0.05
p <- 0.1
n <- 100

iters <- 10000

wald.ci <- matrix(rep(0, iters * 2), iters, 2)
score.ci <- matrix(rep(0, iters * 2), iters, 2)
likelihoodratio.ci <- matrix(rep(iters * 2), iters, 2)

for (i in 1:iters) {
  y <- rbinom(1, size = n, prob = p)

  wald.ci[i, ] <- WaldConfidenceInterval(y, n, alpha)
  score.ci[i, ] <- ScoreConfidenceInterval(y, n, alpha)
  likelihoodratio.ci[i, ] <- LikelihoodRatioConfidenceInterval(y, n, alpha)
}
confidence.intervals <- list(wald=wald.ci, score=score.ci, likelihood=likelihoodratio.ci)
summarised.ci <- lapply(confidence.intervals, function(ci) {
    data.frame(left=mean(ci[, 1], na.rm = TRUE), right=mean(ci[, 2], na.rm =TRUE),
               amplitude=Amplitude(ci), coverage=Coverage(ci, p))
  }) %>%
  { do.call(rbind, .) }

summarised.ci
