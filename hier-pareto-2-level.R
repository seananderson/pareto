library(dplyr)
library(rstan) # >= 2.7.0-1
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
library(ggplot2)

source("load_dat.R") # creates fish.l
set.seed(1) # downsample to experiment:
d <- fish.l[sample(seq_len(nrow(fish.l)), 10000), ] # downsample for speed for now

xmins <- group_by(d, SiteCode, Ecoregion) %>% summarise(xmin = min(pcbm), n = n()) %>%
  filter(n >= 5)
d <- inner_join(d, xmins)
d$eco_index <- as.numeric(as.factor(as.character(d$Ecoregion)))
d <- d %>% arrange(eco_index)
ei <- select(d, xmin, SiteCode, eco_index) %>% unique() %>%
  mutate(site_index = 1:length(SiteCode))
d <- inner_join(d, ei)

n_site <- max(d$site_index)
n_eco <- max(d$eco_index)

d <- d %>% arrange(eco_index, site_index)

m <- stan("pareto-hier-2-level.stan", iter = 100,
  data = list(
    N = length(d$pcbm),
    J = n_site,
    K = n_eco,
    site_index = d$site_index,
    eco_index = ei$eco_index,
    y = d$pcbm,
    y_mins = ei$xmin
    ),
  pars = c("alphas_site", "cv_alpha_site", "mu_alpha", "cv_alpha_eco", "alphas_eco"))

m_summ <- summary(m)$summary
max(m_summ[, "Rhat"])
min(m_summ[, "n_eff"])
saveRDS(m, file = "m.rds")

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
