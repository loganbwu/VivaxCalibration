// These models are for testing a variable-length chain of delays in an ODE model

functions {
  vector my_ode(real time, vector y, vector theta) {
    real alpha = theta[1];
    real dummy = theta[2];
    
    vector[num_elements(y)] dydt;
    for (i in 1:num_elements(y)) {
      if (i == 1) {
        dydt[i] = -alpha * y[i];
      } else if (i < num_elements(y)) {
        dydt[i] = -alpha * y[i] + alpha * y[i-1];
      } else {
        dydt[i] = alpha * y[i-1];
      }
    }
    return dydt;
  }
}

data {
  int<lower=0>            T;
  real                    t0;
  array[T] real<lower=t0> ts;
  int<lower=1>            n_delays;
  vector[n_delays+1]      y0;
  array[T] real           cases;
  real dummy;
}

parameters {
  real<lower=0> sigma; // noise
  real alpha;
}

transformed parameters {
  vector[2] theta;
  theta[1] = alpha;
  theta[2] = dummy;
  // real y[T, n_delays+1];
  array[T] vector[n_delays+1] y = ode_bdf(my_ode, y0, t0, ts, theta);
}

model {
  sigma ~ cauchy(0, 5);
  theta[1] ~ uniform(0, 1);
  for (t in 1:T) {
    cases[t] ~ normal(y[t][n_delays+1], sigma);
  }
}
