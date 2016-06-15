data {
  int<lower=0> N; // rows of data
  int<lower=0> J; // number of sites
  int<lower=0> K; // number of ecoregions
  int site_index[N];
  int eco_index[J];
  vector[N] y; // vector to hold observations
  real y_mins[J]; // minimum observations
}
parameters {
  real<lower=0.01> mu_alpha;
  real<lower=0.01> cv_alpha_eco;
  real<lower=0.01> alphas_eco[K];
  real<lower=0.01> cv_alpha_site;
  real<lower=0.01> alphas_site[J];
}
model {
  // containers to hold calculations:
  vector[N] y_min;
  vector[N] alpha;
  vector[J] den;

  // priors:
  cv_alpha_eco ~ cauchy(0, 2);
  cv_alpha_site ~ cauchy(0, 2);
  mu_alpha ~ cauchy(0, 5);
  alphas_site ~ cauchy(0, 10);
  alphas_eco ~ cauchy(0, 10);

  // ecoregion level alphas:
  alphas_eco ~ gamma(1/cv_alpha_eco^2, 1/(mu_alpha * cv_alpha_eco^2));

  for (j in 1:J) {
    den[j] <- 1/(alphas_eco[eco_index[j]] * cv_alpha_site^2);
  }
  // site level alphas:
  alphas_site ~ gamma(1/cv_alpha_site^2, den);

  for (i in 1:N) {
    y_min[i] <- y_mins[site_index[i]];
    alpha[i] <- alphas_site[site_index[i]];
  }

  y ~ pareto(y_min, alpha); // vectorized
}
