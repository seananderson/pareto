data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> ymin_upper;
}
parameters {
  real<lower=0.01> alpha;
  real<lower=0.01, upper=ymin_upper> y_min;
}
model {
  alpha ~ cauchy(0, 3);
  y ~ pareto(y_min, alpha);
}
