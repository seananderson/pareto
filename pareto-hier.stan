data {
  int<lower=0> N; // rows of data
  int<lower=0> J; // group #
  int group_index[N];
  vector[N] y; // vector to hold observations
  real y_mins[J];
}
parameters {
  real<lower=0.01> mu_alpha;
  real<lower=0.01> cv_alpha;
  real<lower=0.01> alphas[J];
}
model {
  vector[N] y_min;
  vector[N] alpha;
  cv_alpha ~ cauchy(0, 2.5);
  mu_alpha ~ cauchy(0, 5);
  alphas ~ cauchy(0, 10);
  alphas ~ gamma(1/cv_alpha^2, 1/(mu_alpha * cv_alpha^2)); // vectorized
  for (i in 1:N) {
    y_min[i] <- y_mins[group_index[i]];
    alpha[i] <- alphas[group_index[i]];
  }
  y ~ pareto(y_min, alpha); // vectorized
}
