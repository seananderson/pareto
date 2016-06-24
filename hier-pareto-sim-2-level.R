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
simulation_test_pareto <- function(
  cv_alpha = 0.45,
  cv_alpha2 = 0.27,
  mu_alpha = 1.0,
  n_groups = 10,
  n_groups2 = 7,
  N = 40
) {
  xmins <- runif(n_groups * n_groups2, 2, 15)
  # xmins <- runif(n_groups * n_groups2, 8, 8)
  alphas_eco <- rgamma(n_groups, shape = 1/cv_alpha^2, 
    rate = 1/(mu_alpha * cv_alpha^2))
  alphas_site <- lapply(alphas_eco, function(x)
    rgamma(n_groups2, shape = 1/cv_alpha2^2, rate = 1/(x * cv_alpha2^2))) %>%
      unlist()
    # plot(alphas_site);points(rep(alphas_eco, each = n_groups2), pch = 3)
    ei <- data.frame(alphas_site = alphas_site, site_index = 1:length(alphas_site),
      eco_index = rep(1:n_groups, each = n_groups2))
    d <- data.frame(xm = rep(xmins, each = N), alpha = rep(alphas_site, each = N),
      n = 1)
    d <- plyr::mdply(d, rpareto) %>% rename(y = V1)
    d <- data.frame(d, site = rep(ei$site, each = N))
    # ggplot(d, aes(y)) + geom_histogram() + facet_wrap(~site)


#    browser()
    m <- stan("pareto-hier-2-level.stan", iter = 2000, chains = 4,
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
    cv_eco <- extract(m)$cv_alpha_eco %>% median()
    cv_site <- extract(m)$cv_alpha_site %>% median()
    # xmin <- extract(m)$xmin %>% median()

    med <- apply(a, 2, median)
    med2 <- apply(a2, 2, median)

    max_rhat <- max(summary(m)$summary[, "Rhat"])
    min_neff <- min(summary(m)$summary[, "n_eff"])

    list(a_site = alphas_site, a_eco = alphas_eco, a_site_hat = med, a_eco_hat = med2,
      scalars = data.frame(cv_site = cv_alpha2, cv_eco = cv_alpha, 
      cv_site_hat = cv_site, cv_eco_hat = cv_eco, #xmin = min(xmins), xmin_hat = xmin,
      max_rhat = max_rhat, min_neff = min_neff))

}

set.seed(123)
output <- list()
for (i in 1:1)   {
  output[[i]] <- simulation_test_pareto()
  output[[i]]$rep <- i
}

a_site_sim <- plyr::ldply(output, function(x) {
  data.frame(rep = x$rep, a_site = x$a_site, a_site_hat = x$a_site_hat)})
a_eco_sim <- plyr::ldply(output, function(x) {
  data.frame(rep = x$rep, a_eco = x$a_eco, a_eco_hat = x$a_eco_hat)})
scalars_sim <- plyr::ldply(output, function(x) {
  data.frame(rep = x$rep, x$scalars)})

ggplot(a_site_sim, aes(a_site, a_site_hat)) + geom_point() + facet_wrap(~rep) +
 geom_abline(intercept = 0, slope = 1)

ggplot(a_eco_sim, aes(a_eco, a_eco_hat)) + geom_point() + facet_wrap(~rep) +
 geom_abline(intercept = 0, slope = 1)

scalars_sim
x <- seq(0, 25, length.out = 100)
plot(x, dnorm(x, 0, 2), type = "l", ylim = c(0, 0.26))
lines(x, dcauchy(x, 0, 5), lty = 2)

plot(scalars_sim$cv_site_hat);abline(h = scalars_sim$cv_site, col = "red")
plot(scalars_sim$cv_eco_hat);abline(h = scalars_sim$cv_eco, col = "red")
plot(scalars_sim$cv_eco_hat);abline(h = scalars_sim$cv_eco, col = "red")

source("load_dat.R")
## fish.l
# head(fish.l)
# nrow(fish.l)
set.seed(1) # downsample to experiment:
ids <- as.character(unique(fish.l$Ecoregion))[1:12]
d <- fish.l[fish.l$Ecoregion %in% ids, ]
ids <- as.character(unique(d$SiteCode))[1:100]
d <- d[d$SiteCode %in% ids, ]
nrow(fish.l)
nrow(d)
ec <- select(d, Ecoregion)  %>% unique()  %>% mutate(eco_id = 1:n())
si <- select(d, SiteCode)  %>% unique()  %>% mutate(site_id = 1:n())
d <- inner_join(d, ec) %>% inner_join(si)
xm <- d %>% group_by(site_id)  %>% 
  summarize(xmin = min(pcbm))
ec2 <- select(d, site_id, eco_id) %>% unique()

m <- stan("pareto-hier-2-level.stan", iter = 500, chains = 4,
  data = list(
    N = length(d$pcbm),
    J = max(d$site_id),
    K = max(d$eco_id),
    site_index = d$site_id,
    eco_index = ec2$eco_id,
    y = d$pcbm,
    y_mins = xm$xmin))
b <- broom::tidy(m, estimate.method = "median", rhat = TRUE, ess = TRUE)
# plot(scalars_sim$cv_site_hat);abline(h = scalars_sim$cv_site, col = "red")
# plot(scalars_sim$xmin_hat);abline(h = scalars_sim$xmin, col = "red")
# scalars_sim <- scalars_sim %>% mutate(proportional_error = gcc0)
# Probably best to fix the minimum value 
# Models are very sensitive to this value being wrong 
# The data is estimated almost exactly as being the minimum observed value always 
# Can be very bad if this minimum parameter isn't allowed to get big enough 
# And very bad computationally if you don't limit the parameter to the minimum observed 






# jit <- jitter(rep(0, n_groups), 0.3)
pdf("heir-2-level-sim.pdf")
plot(alphas2, med)
plot(alphas, med2)

plot(alphas2);points(rep(alphas, each = n_groups2), pch = 3)
points(med, col = "red")
hist(cv_eco);abline(v = cv_alpha, col = "red")
hist(cv_site);abline(v = cv_alpha2, col = "red")

plot(m, pars = c("cv_alpha_eco", "cv_alpha_site"), show_density = TRUE)
traceplot(m, pars = c("cv_alpha_eco", "cv_alpha_site"))

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
