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
set.seed(1236)
cv_alpha <- 0.15
cv_alpha2 <- 0.07
mu_alpha <- 1.0
n_groups <- 10
n_groups2 <- 10
xmins <- runif(n_groups * n_groups2, 2, 15)
N <- 50
alphas <- rgamma(n_groups, shape = 1/cv_alpha^2, rate = 1/(mu_alpha * cv_alpha^2))
alphas2 <- lapply(alphas, function(x)
  rgamma(n_groups2, shape = 1/cv_alpha2^2, rate = 1/(x * cv_alpha2^2))) %>%
  unlist()
plot(alphas2);points(rep(alphas, each = n_groups2), pch = 3)
ei <- data.frame(alphas2 = alphas2, site_index = 1:length(alphas2), eco_index = rep(1:n_groups, each = n_groups2))
d <- data.frame(xm = rep(xmins, each = N), alpha = rep(alphas2, each = N), n = 1)
d <- plyr::mdply(d, rpareto) %>% rename(y = V1)
d <- data.frame(d, site = rep(ei$site, each = N))
# d$group <- rep(seq_len(n_groups), each = N)
ggplot(d, aes(y)) + geom_histogram() + facet_wrap(~site)

m <- stan("pareto-hier-2-level.stan", iter = 1000,
  data = list(
    N = length(d$y),
    J = n_groups2 * n_groups,
    K = n_groups,
    site_index = d$site,
    eco_index = ei$eco_index,
    y = d$y,
    y_mins = xmins
    ))

a <- extract(m)$alphas_site
a2 <- extract(m)$alphas_eco
cv_eco <- extract(m)$cv_alpha_eco
cv_site <- extract(m)$cv_alpha_site

med <- apply(a, 2, median)
med2 <- apply(a2, 2, median)

# jit <- jitter(rep(0, n_groups), 0.3)
pdf("heir-2-level-sim.pdf")
plot(alphas2, med)
plot(alphas, med2)

plot(alphas2);points(rep(alphas, each = n_groups2), pch = 3)
points(med, col = "red")
hist(cv_eco);abline(v = cv_alpha, col = "red")
hist(cv_site);abline(v = cv_alpha2, col = "red")

dev.off()

# # pdf("heir-alphas-sim.pdf", width = 6, height = 6)
# par(mfrow = c(2, 2))
#
# plot(c(jit-0.5, jit+0.5), c(med2, med), xlim = c(-0.5, 0.6), log = "y")
# segments(jit-0.5, med2, jit+0.5, med)
# points(jit + 0.55, alphas, col = "red")
#
# ci <- apply(a, 2, quantile, probs = c(0.1, 0.9))
#
# plot(med, log = "y", ylab = "alpha", ylim = range(ci), pch = 4)
# segments(1:length(med), ci[1, ], 1:length(med), ci[2, ], col = "#00000030")
# points(1:length(med), alphas, col = "red", pch = 20)
# points(1:length(med2), med2, col = "blue", pch = 4)
#
# plot(density(extract(m)$cv_alpha))
# abline(v = cv_alpha)
# plot(density(extract(m)$mu_alpha))
# abline(v = mu_alpha)
# # dev.off()
