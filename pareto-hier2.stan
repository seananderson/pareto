data {
  int<lower=0> N; // rows of data
  int<lower=0> J; // group #
  int group_index[N];
  vector[N] y; // vector to hold observations
  real y_mins[J];
}
parameters {
  real<lower=0.01, upper=100> alphas[J];
}
transformed parameters {
  vector[N] y_min;
  vector[N] alpha;
  for (i in 1:N) {
    y_min[i] <- y_mins[group_index[i]];
    alpha[i] <- alphas[group_index[i]];
  }
}
model {
  y ~ pareto(y_min, alpha);
}
