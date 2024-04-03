// structure based on https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html

// Changes from v1: Phase is now a model parameter

functions {
  real suitability(real t, real eps, real kappa, real phase) {
    // To speed up the ODE solver, we use omega=1 to get to equilibrium and only begin seasonality soon before the measurement period
    real omega;
    if (t < (-20*365.25)) {
      omega = 1;
    } else {
      omega = eps + (1-eps)*pi()/beta(0.5, kappa+0.5)*((1+sin(2*pi()*(t - phase)/365.25))/2)^kappa;
    }
    return(omega);
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
    
    real eps = theta[2];
    real kappa = theta[3];
    real phase = theta[4];
    real omega = suitability(t, eps, kappa, phase);
    real foi = theta[1] * omega;
    
    real dIl_dt = (1-alpha)*(foi*(Il+I0)+delta)*(S0+Sl) + (foi*(Il+I0)+delta)*I0 + (1-alpha)*f*Sl - gammal*Il - r*Il;
    real dI0_dt = -(foi*(Il+I0)+delta)*I0 + gammal*Il - r*I0;
    real dSl_dt = -(1-alpha*(1-beta))*(foi*(Il+I0)+delta+f)*Sl + alpha*(1-beta)*(foi*(Il+I0)+delta)*S0 - gammal*Sl + r*Il;
    real dS0_dt = -(1-alpha*beta)*(foi*(Il+I0)+delta)*S0 + (foi*(I0+Il)+delta)*alpha*beta*Sl + alpha*beta*f*Sl + gammal*Sl + r*I0;
    real dCumulativeInfections = (foi*(Il+I0)+delta)*(S0+Sl) + f*Sl;
    
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
}

transformed data {
  real x_r[6] = {
    r,
    gammal,
    f,
    alpha,
    beta,
    delta
  };
  int x_i[1] = { N };
  // We need to add an extra timepoint before to make difference calculations valid
  // Note that i+1 along ts_extended and y is the ith component of ts
  real ts_extended[n_times+1];
  real dt = ts[2] - ts[1];
  ts_extended[1] = ts[1] - dt;
  for (i in 1:n_times) {
    ts_extended[i+1] = ts[i];
  }
}

parameters {
  real<lower=0, upper=0.9> lambda;
  real<lower=0, upper=1> eps;
  real<lower=0> kappa;
  real<lower=0, upper=365.25> phase;
  real<lower=0> phi_inv;
}

transformed parameters{
  real y[n_times+1, 5];
  real incidence[n_times];
  real phi = 1. / phi_inv;
  {
    real theta[4];
    theta[1] = lambda;
    theta[2] = eps;
    theta[3] = kappa;
    theta[4] = phase;
    
    y = integrate_ode_bdf(champagne, y0, t0, ts_extended, theta, x_r, x_i);
  }
  
  for (i in 1:n_times) {
    incidence[i] = fmax(1e-12, (y[i+1, 5] - y[i, 5]) * N * alpha);
  }
}

model {
  //priors
  lambda ~ exponential(1);
  // eps ~ uniform(0, 1);
  kappa ~ exponential(0.1);
  // phase ~ uniform(0, 365.25);
  phi_inv ~ exponential(5);
  
  //sampling distribution
  for (i in 1:n_times) {
    cases[i] ~ neg_binomial_2(incidence[i], phi);
  }
}

generated quantities {
  real sim_cases[n_times];
  real susceptible[n_times];
  real infectious[n_times];
  real latent[n_times];
  real R0[n_times];
  real Rc[n_times];
  real foi[n_times];
  
  for (i in 1:n_times) {
    susceptible[i] = y[i+1, 4];
    infectious[i] = y[i+1, 1] + y[i+1, 2];
    latent[i] = y[i+1, 3];
    sim_cases[i] = neg_binomial_2_rng(incidence[i], phi);
    foi[i] = lambda * suitability((ts_extended[i+1]+ts_extended[i])/2, eps, kappa, phase);
    R0[i] = foi[i]/r + foi[i] * f / (gammal * (f + gammal + r));
    Rc[i] = foi[i] * (1-alpha) * (gammal+r) * (f + gammal) / (r * (gammal * (f + gammal + r) + alpha*f * (beta*(r + gammal) - gammal)));
  }
}
