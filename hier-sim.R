library(dplyr)
library(rstan)
library(ggplot2)

# set.seed(1)
alpha <- 0.2
sd_beta <- 0.2
mu_beta <- 1
n_groups <- 40
sigma <- 4
betas <- rnorm(n_groups, mu_beta, sd_beta)
N <- 100
x <- seq(-2, 2, length.out = N)
d <- data.frame(alpha = rep(alpha, each = N), beta = rep(betas, each = N), n = 1,
  sigma = sigma, x = rep(x, n_groups))
d <- mutate(d,
  mean = alpha + beta * x,
  y = rnorm(N * n_groups, mean, sigma))
d$group <- rep(seq_len(n_groups), each = N)
med2 <- plyr::daply(d, "group", function(x) {
  m2 <- lm(y ~ x, data = x)
  coef(m2)[[2]]
})

# ggplot(d, aes(x, y)) + geom_point() + facet_wrap(~group)

library(lme4)
m <- lmer(y ~ x + (x | group), data = d)

med <- coef(m)$group[,2]
# med2 <- as.numeric(coef(m2)[-c(1)])

par(mfrow = c(1, 1))
jit <- jitter(rep(0, n_groups), 0.3)
plot(c(jit-0.5, jit+0.5), c(med2, med), xlim = c(-0.5, 0.6))

segments(jit-0.5, med2, jit+0.5, med)
points(jit + 0.55, betas, col = "red")
