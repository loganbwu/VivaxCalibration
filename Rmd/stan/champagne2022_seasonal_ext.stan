// structure based on https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

// In this version, we allow the seasonality parameters to be estimated

functions {
  real suitability(real t, real eps, real kappa, real phase) {
    return(eps + (1-eps)*pi()/beta(0.5, kappa+0.5)*((1+sin(2*pi()/365.25*(t-phase))))^kappa); # no normalising beta function
  }
  
  real[] champagne(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    
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
    real delta = x_r[6];
    real phase = x_r[7];
    
    real omega = suitability(t, theta[2], theta[3], phase);
    real foi = theta[1] * omega;
    
    real dIl_dt = (1-alpha)*(foi*(Il+I0)+delta)*(S0+Sl) + (foi*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gammal*Il - r*Il;
    real dI0_dt = -(foi*(Il+I0)+delta)*I0 + gammal*Il - r*I0;
    real dSl_dt = -(1-alpha*(1-beta))*(foi*(Il+I0)+delta+f)*Sl + alpha*(1-beta)*(foi*(Il+I0)+delta)*S0 - gammal*Sl + r*Il;
    real dS0_dt = -(1-alpha*beta)*(foi*(Il+I0)+delta)*S0 + (foi*(I0+Il)+delta)*alpha*beta*Sl + alpha*beta*f*Sl + gammal*Sl + r*I0;
    real dCumulativeInfections = (foi*(Il+I0)+delta)*(S0+Sl);
    
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
  real<lower=0> delta;
  real<lower=0> phase;
}

transformed data {
  real x_r[7] = {
    r,
    gammal,
    f,
    alpha,
    beta,
    delta,
    phase
  };
  int x_i[1] = { N };
}

parameters {
  real<lower=0, upper=0.05> lambda;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> eps;
  real<lower=0, upper=30> kappa;
  //real<lower=0, upper=365.25> phase; // note: may be an issue with this being a cyclical variable
}

transformed parameters{
  real y[n_times, 5];
  real<lower=0> incidence[n_times - 1];
  real phi = 1. / phi_inv;
  {
    real theta[3];
    theta[1] = lambda;
    theta[2] = eps;
    theta[3] = kappa;
    //theta[4] = phase;
    
    y = integrate_ode_bdf(champagne, y0, t0, ts, theta, x_r, x_i);
  }
  
  for (i in 1:n_times-1)
    incidence[i] = fmax(1e-12, (y[i+1, 5] - y[i, 5]) * N); # incorrect, just for testing
}

model {
  //priors
  lambda ~ normal(0, 1e4);
  eps ~ uniform(0, 1);
  kappa ~ exponential(0.1);
  
  phi_inv ~ exponential(0.1);
  
  //sampling distribution
  for (i in 1:(n_times - 1))
    cases[i] ~ neg_binomial_2(incidence[i], phi);
}

generated quantities {
  real sim_cases[n_times-1];
  real susceptible[n_times-1];
  real R0[n_times-1];
  real Rc[n_times-1];
  
  for (i in 1:(n_times - 1)) {
    sim_cases[i] = neg_binomial_2_rng(incidence[i], phi);
    susceptible[i] = y[i, 4];
    real foi = lambda * suitability(ts[i], eps, kappa, phase);
    R0[i] = foi/r + foi * f / (gammal * (f + gammal + r));
    Rc[i] = foi * (1-alpha) * (gammal+r) * (f + gammal) / (r * (gammal * (f + gammal + r) + alpha*f * (beta*(r + gammal) - gammal)));
  }
}
