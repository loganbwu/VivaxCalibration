functions {
  real[] sir(real t, real[] y, real[] theta, 
  real[] x_r, int[] x_i) {
    
    real Il = y[1];
    real I0 = y[2];
    real Sl = y[3];
    real S0 = y[4];
    // int N = x_i[1]; # unused
    
    real lambda = theta[1];
    real r = theta[2];
    real gammal = theta[3];
    real f = theta[4];
    real delta = theta[5];
    real alpha = theta[6];
    real beta = theta[7];
    real rho = theta[8];
    
    real dIl_dt = (1-alpha)*(lambda*(Il+I0)+delta)*(S0+Sl) + (lambda*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gammal*Il - r*Il;
    real dI0_dt = -(lambda*(Il+I0)+delta)*I0 + gammal*Il - r*I0;
    real dSl_dt = -(1-alpha*(1-beta))*(lambda*(Il+I0)+delta+f)*Sl;
    real dS0_dt = -(1-alpha*beta)*(lambda*(Il+I0)+delta)*S0 + (lambda*(I0+Il)+delta)*alpha*beta*Sl + alpha*beta*f*Sl + gammal*Sl + r*I0;
    
    return {dIl_dt, dI0_dt, dSl_dt, dS0_dt};
  }
}

data {
  int<lower=1> n_times;
  real<lower=0, upper=1> y0[4];
  real t0;
  real ts[n_times];
  int N;
  int cases[n_times];
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  real<lower=0> lambda;
  real<lower=0> r;
  real<lower=0> gammal;
  real<lower=0> f;
  real<lower=0> delta;
  real<lower=0, upper=1> alpha;
  real<lower=0, upper=1> beta;
  real<lower=0, upper=1> rho;
  
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_times, 4];
  real incidence[n_times - 1];
  real phi = 1. / phi_inv;
  {
    real theta[8];
    theta[1] = lambda;
    theta[2] = 1./60; #r;
    theta[3] = 1./223; #gammal;
    theta[4] = 1./72; #f;
    theta[5] = delta;
    theta[6] = 0.21; #alpha;
    theta[7] = 0.66; #beta;
    theta[8] = 0.21; #rho;
    
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }
  
  for (i in 1:n_times-1)
  // for (i in 2:n_times)
    incidence[i] = y[i, 1] - y[i+1, 1]; # incorrect, just for testing
}

model {
  //priors
  lambda ~ normal(0, 1e6);
  delta ~ normal(0, 1e6);
  
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
  // cases ~ neg_binomial_2(col(to_matrix(y), 1) + col(to_matrix(y), 2), phi);
  cases[1:(n_times-1)] ~ neg_binomial_2(incidence, phi);
}

generated quantities {
  real pred_cases[n_times-1];
  pred_cases = neg_binomial_2_rng(incidence, phi);
}
