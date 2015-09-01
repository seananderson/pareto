data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
}
transformed data {
  real<lower=0> y_min;
  y_min <- min(y);
}
parameters {
  real<lower=0.01> alpha;
}
model {
  alpha ~ cauchy(0, 3);
  y ~ pareto(y_min, alpha);
}
