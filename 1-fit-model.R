library(dplyr)
library(ggplot2)
library(rstan) 
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source("load_dat.R")
## fish.l
fish.l <- fish.l %>% group_by(SiteCode) %>% mutate(n = n()) %>% filter(n>200) %>% as_data_frame() %>% 
  mutate(SiteCode = as.character(SiteCode), Ecoregion = as.character(Ecoregion))
n_eco <- length(unique(fish.l$Ecoregion))
n_sample <- 20

set.seed(1) # downsample to experiment:
ids <- as.character(unique(fish.l$Ecoregion))[sample(1:n_eco, n_sample)]
d <- fish.l[fish.l$Ecoregion %in% ids, ]
# ids <- as.character(unique(d$SiteCode))[1:100]
# d <- d[d$SiteCode %in% ids, ]
nrow(fish.l)
nrow(d)
ec <- select(d, Ecoregion)  %>% unique()  %>% mutate(eco_id = 1:n())
si <- select(d, SiteCode)  %>% unique()  %>% mutate(site_id = 1:n())
d <- inner_join(d, ec) %>% inner_join(si)
xm <- d %>% group_by(site_id)  %>% 
  summarize(xmin = min(pcbm))
ec2 <- group_by(d, site_id) %>% 
  summarize(eco_id = eco_id[1]) # CHECK! that these are unique and only 1
# filter(xm,xmin>500 ) 

# m <- stan("pareto-hier-2-level.stan", iter = 400, chains = 4,
#   data = list(
#     N = length(d$pcbm),
#     J = max(d$site_id),
#     K = max(d$eco_id),
#     site_index = d$site_id,
#     eco_index = ec2$eco_id,
#     y = d$pcbm,
#     y_mins = xm$xmin))
# m
# rstan::traceplot(m, pars = "lp__")
# 
# b <- broom::tidy(m, estimate.method = "median", rhat = TRUE, ess = TRUE,
#   conf.int = TRUE)
   
# ggplot(b, aes(estimate, term)) +
#   geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
#   geom_point()

hist(xm$xmin)
unique(xm$xmin)
# x <- seq(0, 15, length.out = 100)
# plot(x, dnorm(x, 0, 10), type = "l", ylim = c(0, 0.26))
# lines(x, dcauchy(x, 0, 10), lty = 2)

m2 <- stan("pareto-hier-2-level-xmin.stan", iter = 400, chains = 4,
  data = list(
    N = length(d$pcbm),
    J = max(d$site_id),
    K = max(d$eco_id),
    site_index = d$site_id,
    eco_index = ec2$eco_id,
    y = d$pcbm,
    ymin = min(c(xm$xmin, 32))))
m2                               
saveRDS(m2, file = "pareto-hier-xmin.rds")
rstan::traceplot(m2, pars = "lp__")

b2 <- broom::tidy(m2, estimate.method = "median", rhat = TRUE, ess = TRUE,
  conf.int = TRUE)
   
b2 <- b2 %>% mutate(term_group = gsub("\\[[0-9]*\\]","",term)) 
b2 <- b2 %>% mutate(id = gsub("[a-z_]+\\[([0-9]*)\\]","\\1",term)) 

ggplot(filter(b2, term != "xmin"), aes(estimate, term, color = term_group)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_point()
ggsave("pareto-hier-xmin.pdf", width = 10, height = 40)

ids <- select(d, site_id, eco_id, Ecoregion) %>% unique()
as <- filter(b2, term_group == "alphas_site") %>% 
  mutate(site_id = as.numeric(id)) %>% 
    inner_join(ids)
as <- group_by(as, eco_id) %>% mutate(site_id2 = 1:n()) %>% 
  as_data_frame()
ae <- filter(b2, term_group == "alphas_eco") %>% 
  mutate(eco_id = as.numeric(id))  %>% 
    select(eco_id, estimate) %>% 
      rename(group_estimate = estimate)
as <- inner_join(as, ae)

ggplot(as, aes(estimate, site_id2)) +       
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0) +
  geom_point() +
  facet_wrap(~Ecoregion) +
  scale_x_log10() +
  geom_vline(aes(xintercept = group_estimate), lty = 2)
ggsave("pareto-hier-xmin-eco.pdf", width = 9, height = 9)

# hist(extract(m2)$xmin, breaks=25)
cp <- function(x, xm, alpha) (xm / x) ^ alpha
xmin <- filter(b2, term_group == "xmin")$estimate
output <- plyr::adply(filter(b2, term_group == "alphas_site"), 1, function(x) {
  id <- as.numeric(x[["id"]])
  sizes <- dplyr::filter(d, site_id == id)$pcbm %>% sort()
  s <- cp(sizes, xm = xmin, alpha = x[["estimate"]])
  P <- ecdf(sizes)
  o <- data.frame(site_id = id, sizes = sizes, cdf = 1 - P(sizes), fit = s)
  o
})
output <- inner_join(output, ids)
output <- group_by(output, eco_id) %>% mutate(site_id2 = site_id - min(site_id) + 1) %>% 
  as_data_frame()
p <- ggplot(output, aes(sizes, cdf, colour = as.factor(site_id2))) +
  geom_point(alpha = 0.01) +
  facet_wrap(~Ecoregion) + 
  scale_x_log10() +
  geom_line(aes(y = fit)) +
    theme(legend.position = "none")
  print(p)
ggsave("pareto-hier-xmin-cdf.pdf", width = 10, height = 9)
