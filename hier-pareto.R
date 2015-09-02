library(dplyr)
library(rstan) # >= 2.7.0-1
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)

source("load_dat.R") # creates fish.l
set.seed(1) # downsample to experiment:
d <- fish.l[sample(seq_len(nrow(fish.l)), 50000), ] # downsample for speed for now
xmins <- group_by(d, SiteCode) %>% summarise(xmin = min(pcbm), n = n()) %>%
  filter(n >= 5)
xmins$site_index <- 1:nrow(xmins)
n_groups <- max(xmins$site_index)
d <- inner_join(d, xmins)
d <- d %>% arrange(site_index, pcbm)

# by site:
m <- stan("pareto-hier.stan", iter = 500,
  data = list(
    N = length(d$pcbm),
    y = d$pcbm,
    y_mins = xmins$xmin,
    J = n_groups,
    group_index = d$site_index),
  pars = c("alphas", "cv_alpha", "mu_alpha"))

m_summ <- summary(m)$summary
max(m_summ[, "Rhat"])
min(m_summ[, "n_eff"])
saveRDS(m, file = "m.rds")

# by ecoregion:
xmins <- group_by(d, Ecoregion) %>% summarise(xmin = min(pcbm), n = n())
d$xmin <- NULL
d$n <- NULL
d$ecoregion_index <- NULL
xmins$ecoregion_index <- 1:nrow(xmins)
n_groups <- max(xmins$ecoregion_index)
d <- inner_join(d, xmins)
d <- d %>% arrange(ecoregion_index, pcbm)

m2 <- stan("pareto-hier.stan", iter = 500,
  data = list(
    N = length(d$pcbm),
    y = d$pcbm,
    y_mins = xmins$xmin,
    J = n_groups,
    group_index = d$ecoregion_index),
  pars = c("alphas", "cv_alpha", "mu_alpha"))
saveRDS(m2, file = "m2.rds")

m2_summ <- summary(m2)$summary
max(m2_summ[, "Rhat"])
min(m2_summ[, "n_eff"])

# look at posteriors:
a <- extract(m)$alphas
a2 <- extract(m2)$alphas
med <- apply(a, 2, median)
med2 <- apply(a2, 2, median)
jit <- jitter(rep(0, n_groups), 0.3)

pdf("hier-alphas.pdf", width = 7, height = 10)
par(mfrow = c(4, 2))

ord <- order(med)
ci <- apply(a, 2, quantile, probs = c(0.1, 0.9))
plot(med[ord], log = "y", ylab = "alpha", ylim = range(ci), pch = 20, cex = 0.3)
segments(1:length(med), ci[1, ord], 1:length(med), ci[2, ord], col = "#00000010")

ord <- order(med2)
ci <- apply(a2, 2, quantile, probs = c(0.1, 0.9))
plot(med2[ord], log = "y", ylab = "alpha", ylim = range(ci), pch = 20, cex = 0.5)
segments(1:length(med2), ci[1, ord], 1:length(med2), ci[2, ord], col = "#00000070")

hist(extract(m)$cv_alpha)
hist(extract(m2)$cv_alpha)
hist(extract(m)$mu_alpha)
hist(extract(m2)$mu_alpha)

q <- data.frame(med = med, site_index = 1:length(med))
q2 <- data.frame(med2 = med2, ecoregion_index = 1:length(med2))
d <- inner_join(d, q)
d <- inner_join(d, q2)
dd <- d[,c("med", "med2")]
dd <- unique(dd)

plot(dd$med, dd$med2, col = "#00000030")
dev.off()
