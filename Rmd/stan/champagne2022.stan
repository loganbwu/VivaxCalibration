// structure based on https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

functions {
  real[] champagne(real t, real[] y, real[] theta, 
  real[] x_r, int[] x_i) {
    
    real Il = y[1];
    real I0 = y[2];
    real Sl = y[3];
    real S0 = y[4];
    real CumulativeInfections = y[5];
    
    // int N = x_i[1]; # unused
    real r = x_r[1];
    real gammal = x_r[2];
    real f = x_r[3];
    real alpha = x_r[4];
    real beta = x_r[5];
    real rho = x_r[6];
    real delta = x_r[7];
    
    real lambda = theta[1];
    // real delta = theta[2];
    
    real dIl_dt = (1-alpha)*(lambda*(Il+I0)+delta)*(S0+Sl) + (lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gammal*Il - r*Il;
    real dI0_dt = -(lambda*(Il+I0)+delta)*I0 + gammal*Il - r*I0;
    real dSl_dt = -(1-alpha*(1-beta))*(lambda*(Il+I0)+delta+f)*Sl + alpha*(1-beta)*(lambda*(Il+I0)+delta)*S0 - gammal*Sl + r*Il;
    real dS0_dt = -(1-alpha*beta)*(lambda*(Il+I0)+delta)*S0 + (lambda*(I0+Il)+delta)*alpha*beta*Sl + alpha*beta*f*Sl + gammal*Sl + r*I0;
    real dCumulativeInfections = (lambda*(Il+I0)+delta)*(S0+Sl);
    
    return {dIl_dt, dI0_dt, dSl_dt, dS0_dt, dCumulativeInfections};
  }
}

data {
  int<lower=1> n_times;
  real<lower=0> y0[5];
  real t0;
  real ts[n_times];
  int N;
  int cases[n_times];
  real<lower=0> r;
  real<lower=0> gammal;
  real<lower=0> f;
  real<lower=0, upper=1> alpha;
  real<lower=0, upper=1> beta;
  real<lower=0, upper=1> rho;
  real<lower=0> delta;
}

transformed data {
  real x_r[7] = {
    r,
    gammal,
    f,
    alpha,
    beta,
    rho,
    delta
  };
  int x_i[1] = { N };
}

parameters {
  real<lower=0, upper=9> lambda;
  // real<lower=0> delta;
  
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_times, 5];
  real incidence[n_times - 1];
  real phi = 1. / phi_inv;
  {
    real theta[1];
    theta[1] = lambda;
    // theta[2] = delta;
    
    y = integrate_ode_bdf(champagne, y0, t0, ts, theta, x_r, x_i);
  }
  
  for (i in 1:n_times-1)
    incidence[i] = fmax(1e-12, (y[i+1, 5] - y[i, 5]) * N); # incorrect, just for testing
}

model {
  //priors
  lambda ~ normal(0, 1e4);
  delta ~ normal(0, 1e4);
  
  phi_inv ~ exponential(5);
  
  //sampling distribution
  for (i in 1:(n_times - 1))
    cases[i] ~ neg_binomial_2(incidence[i], phi);
}

generated quantities {
  real sim_cases[n_times-1];
  
  for (i in 1:(n_times - 1))
    sim_cases[i] = neg_binomial_2_rng(incidence[i], phi);
}
