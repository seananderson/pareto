library(dplyr)
library(rstan) # >= 2.7.0-1
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)
dpareto <- function(x, xm, alpha)
  ifelse(x > xm, alpha * xm^alpha / (x^(alpha+1)), 0)
qpareto <- function(p, xm, alpha)
  ifelse(p < 0 | p > 1, NaN, xm * (1 - p)^(-1 / alpha))
rpareto <- function(n, xm, alpha)
  qpareto(runif(n), xm, alpha)

## hierarchical distribution on alphas
set.seed(1234)
cv_alpha <- 0.05
mu_alpha <- 1.2
n_groups <- 30
xmins <- runif(n_groups, 2, 20)
N <- 50
alphas <- rgamma(n_groups, shape = 1/cv_alpha^2, rate = 1/(mu_alpha * cv_alpha^2))
d <- data.frame(xm = rep(xmins, each = N), alpha = rep(alphas, each = N), n = 1)
d <- plyr::mdply(d, rpareto) %>% rename(y = V1)
d$group <- rep(seq_len(n_groups), each = N)
ggplot(d, aes(y)) + geom_histogram() + facet_wrap(~group)

m <- stan("pareto-hier.stan", iter = 2000,
  data = list(
    N = length(d$y),
    y = d$y,
    y_mins = xmins,
    J = n_groups,
    group_index = d$group))

m2 <- stan("pareto-hier2.stan", iter = 2000,
  data = list(
    N = length(d$y),
    y = d$y,
    y_mins = xmins,
    J = n_groups,
    group_index = d$group))

a <- extract(m)$alphas
a2 <- extract(m2)$alphas
med <- apply(a, 2, median)
med2 <- apply(a2, 2, median)
jit <- jitter(rep(0, n_groups), 0.3)

par(mfrow = c(2, 2))

plot(c(jit-0.5, jit+0.5), c(med2, med), xlim = c(-0.5, 0.6), log = "y")
segments(jit-0.5, med2, jit+0.5, med)
points(jit + 0.55, alphas, col = "red")

ci <- apply(a, 2, quantile, probs = c(0.1, 0.9))

plot(med, log = "y", ylab = "alpha", ylim = range(ci), pch = 4)
segments(1:length(med), ci[1, ], 1:length(med), ci[2, ], col = "#00000030")
points(1:length(med), alphas, col = "red", pch = 20)
points(1:length(med2), med2, col = "blue", pch = 4)

plot(density(extract(m)$cv_alpha))
abline(v = cv_alpha)
plot(density(extract(m)$mu_alpha))
abline(v = mu_alpha)
