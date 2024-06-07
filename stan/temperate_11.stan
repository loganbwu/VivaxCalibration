// These models are for testing a variable-length chain of delays in an ODE model

// Changes from v2: omega is now 'suitability' function and equations are verified against the clinical pathways
// Changes from v3: changed observation model for cases/relapses
// Changes from v5: observation model is negative binomial not poisson
// Changes from v6: previous model fit parameters were just lambda and phi_inv. Now we add relapse_clinical_immunity
// Changes from v7: now includes all primaries/relapses, not just clinical primaries/relapses
// Changes from v8 (TODO): seasonal parameters are now estimated.
// Changes from v9: relapse clinical immunity is now a fixed data variable
// Changes from v9_noRCI: seasonality is now fixed but alpha and beta are fit with strongly informative priors
// Changes from v10: Setup for overwinter as per draft chapter. Alpha and beta have strongly informative priors. Seasonality is still estimated.

functions {
  real suitability(real t, real eps, real kappa, real phase) {
    real omega = eps + (1-eps)*pi()/beta(0.5, kappa+0.5)*((1+sin(2*pi()*(t - phase)/365.25))/2)^kappa;
    return(omega);
  }
  
  vector my_ode(real time, vector y, vector theta, real[] x_r, int[] x_i) {
    
    // Unpack theta
    real lambda = theta[1];
    real alpha = theta[2];
    real beta = theta[3];
    real relapse_clinical_immunity = theta[4];
    real p_long = theta[5];
    real p_silent = theta[6];
    // real eps = theta[7];
    real kappa = theta[7];
    real phase = theta[8];
    
    
    // Unpack x_r
    real N = x_r[1];
    real gamma_d = x_r[2];
    real gamma_l = x_r[3];
    real delta = x_r[4];
    real f = x_r[5];
    real r = x_r[6];
    real eps = x_r[7];
    
    // Unpack x_i
    int n_dormant = x_i[1];
    
    // Convenience quantities
    int n_stages = n_dormant + 1;
    int n_compartments = 2 + 3*n_stages;
    int active = n_dormant + 1;
    real recover = r;
    real silent = p_silent;
    real primary = 1 - p_silent;
    real treatedprimary = alpha;
    real untreatedprimary = 1-treatedprimary;
    real treatedsilent = 0;
    real untreatedsilent = 1-treatedsilent;
    real treatedrelapse = treatedprimary * (1-relapse_clinical_immunity); // This should change and/or be a fit variable
    real untreatedrelapse = 1-treatedrelapse;
    real complete = beta;
    real incomplete = 1-beta;
    real clear_d = gamma_d;
    real clear_l = gamma_l;
    real advance = delta*n_dormant;
    real relapse = f;
    real long_hyp = p_long;
    real short_hyp = 1-long_hyp;
    
    // Give reasonable names to variables
    real S0 = y[1];
    real I0 = y[2];
    array[n_stages] real Sl;
    array[n_stages] real Scl;
    array[n_stages] real Icl;
    for (i in 1:n_stages) {
      Sl[i] = y[2+i];
      Scl[i] = y[2+n_stages+i];
      Icl[i] = y[2+2*n_stages+i];
    }
    
    // Force of infection
    // real infect = lambda * suitability(time, eps, kappa, phase) * (I0 + sum(Icl)) + phi_2;
    real infect = lambda * suitability(time, eps, kappa, phase) * (I0 + sum(Icl));
    
    // Compute derivatives
    // S0
    real refactor_0 = infect*(short_hyp*treatedprimary*complete + long_hyp*(primary*treatedprimary*complete + silent*treatedsilent*complete));
    real dS0dt = -S0*infect*(short_hyp*(treatedprimary*incomplete+untreatedprimary) + long_hyp*(primary*(treatedprimary*incomplete + untreatedprimary) + silent*(treatedsilent*incomplete + untreatedsilent))) +
    sum(Sl)*refactor_0 +
    Sl[active]*relapse*treatedprimary*complete +
    sum(Sl[1:n_dormant])*clear_d +
    Sl[active]*clear_l +
    sum(Scl)*refactor_0 +
    Scl[active]*relapse*treatedrelapse*complete +
    sum(Scl[1:n_dormant])*clear_d +
    Scl[active]*clear_l +
    I0*recover;
    
    // Sl
    real dSldt[n_stages];
    real refactor_1 = infect*(short_hyp + long_hyp*(primary + silent*treatedsilent*complete)) + advance + clear_d;
    dSldt[1] = -Sl[1]*refactor_1 +
    S0*infect*long_hyp*silent*(treatedsilent*incomplete + untreatedsilent);
    for (i in 2:n_stages-1) {
      dSldt[i] = -Sl[i]*refactor_1 +
      Sl[i-1]*advance;
    }
    dSldt[n_stages] = -Sl[n_stages]*(infect*(short_hyp + long_hyp*(primary + silent*treatedsilent*complete)) + relapse + clear_l) +
    Sl[n_stages-1]*advance;
    
    // Scl
    real dScldt[n_stages];
    real refactor_2 = infect*(short_hyp + long_hyp*(primary*(treatedprimary*complete + untreatedprimary) + silent*treatedsilent*complete)) + advance + clear_d;
    real refactor_21 = infect*long_hyp*primary*treatedprimary*incomplete;
    dScldt[1] = -Scl[1]*refactor_2 +
    S0*refactor_21 +
    Sl[1]*refactor_21 +
    Icl[1]*recover;
    for (i in 2:n_stages-1) {
      dScldt[i] = -Scl[i]*refactor_2 +
      Scl[i-1]*advance +
      Sl[i]*refactor_21 +
      Icl[i]*recover;
    }
    dScldt[n_stages] = -Scl[n_stages]*(infect*(short_hyp*(treatedprimary*complete + untreatedprimary) + long_hyp*(primary*(treatedprimary*complete + untreatedprimary) + silent*treatedsilent*complete)) + relapse*(treatedrelapse*complete + untreatedrelapse) + clear_l) +
    S0*infect*short_hyp*treatedprimary*incomplete +
    sum(Sl)*infect*short_hyp*treatedprimary*incomplete +
    Sl[n_stages]*refactor_21 +
    Sl[n_stages]*relapse*treatedprimary*incomplete +
    sum(Scl[1:n_dormant])*infect*short_hyp*treatedprimary*incomplete +
    Scl[n_stages-1]*advance +
    Icl[n_stages]*recover;
    
    // I0
    real dI0dt = -I0*(infect + recover) +
    sum(Icl[1:n_dormant])*clear_d +
    Icl[active]*clear_l;
    
    // Icl
    real dIcldt[n_stages];
    real refactor_3 = infect*short_hyp + advance + clear_d + recover;
    real refactor_31 = infect*long_hyp*primary*untreatedprimary;
    dIcldt[1] = -Icl[1]*refactor_3 +
    S0*refactor_31 +
    Sl[1]*refactor_31 +
    Scl[1]*refactor_31 +
    I0*infect*long_hyp;
    for (i in 2:n_stages-1) {
      dIcldt[i] = -Icl[i]*refactor_3 +
      Sl[i]*refactor_31 +
      Scl[i]*refactor_31 +
      Icl[i-1]*advance;
    }
    dIcldt[n_stages] = -Icl[n_stages]*(clear_l + recover) +
    S0*infect*short_hyp*untreatedprimary +
    sum(Sl)*infect*short_hyp*untreatedprimary +
    Sl[n_stages]*refactor_31 +
    Sl[n_stages]*relapse*untreatedprimary +
    sum(Scl)*infect*short_hyp*untreatedprimary +
    Scl[n_stages]*refactor_31 +
    Scl[n_stages]*relapse*untreatedrelapse +
    I0*infect*short_hyp +
    sum(Icl[1:n_dormant])*infect*short_hyp +
    Icl[n_stages-1]*advance;
    
    // A clinical case is defined by the number of treatments. It is not affected by the outcome of treatment.
    // We also count the first 'relapse' with no clinical immunity (i.e., cryptic primary infection) as a primary infection because it presents as a primary.
    real dAllPrimary = (S0*infect*(short_hyp + long_hyp*primary) +
    sum(Sl)*infect*(short_hyp + long_hyp*primary) +
    sum(Scl)*infect*(short_hyp + long_hyp*primary) +
    Sl[active]*relapse) * N;
    
    real dAllRelapse = Scl[active]*relapse * N;
    
    real dClinicalPrimary = dAllPrimary * treatedprimary;
    
    real dClinicalRelapse = dAllRelapse * treatedrelapse;
    
    // Assign derivatives
    vector[num_elements(y)] dydt;
    dydt[1] = dS0dt;
    dydt[2] = dI0dt;
    for (i in 1:n_stages) {
      dydt[2+i] = dSldt[i];
      dydt[2+n_stages+i] = dScldt[i];
      dydt[2+2*n_stages+i] = dIcldt[i];
    }
    dydt[n_compartments+1] = dAllPrimary;
    dydt[n_compartments+2] = dAllRelapse;
    dydt[n_compartments+3] = dClinicalPrimary;
    dydt[n_compartments+4] = dClinicalRelapse;
    
    return dydt;
  }
}

data {
  int<lower=1> n_times;
  real t0;
  real ts[n_times];
  int cases[n_times];
  
  real<lower=0> N;
  real<lower=0> gamma_d;
  real<lower=0> gamma_l;
  real<lower=0> delta;
  real<lower=0> f;
  real<lower=0> r;
  real<lower=0, upper=1> eps;
  real<lower=0> alpha_shape1;
  real<lower=0> alpha_shape2;
  real<lower=0> relapse_clinical_immunity_shape1;
  real<lower=0> relapse_clinical_immunity_shape2;
  real<lower=0> p_silent_shape1;
  real<lower=0> p_silent_shape2;
  
  int<lower=1> n_dormant;
  
  vector[2 + 3*(n_dormant+1) + 4] y0; # last elements are incidence related
}

transformed data {
  real x_r[7] = {
    N,
    gamma_d,
    gamma_l,
    delta,
    f,
    r,
    eps
  };
  
  int x_i[1] = {
    n_dormant
  };
  
  int n_stages = n_dormant + 1;
  int n_compartments = 2 + 3*n_stages;
  int len_y = n_compartments + 4;
  
  real ts_extended[n_times+1];
  real dt = ts[2] - ts[1];
  ts_extended[1] = ts[1] - dt;
  for (i in 1:n_times) {
    ts_extended[i+1] = ts[i];
  }
}

parameters {
  real<lower=0, upper=0.5> lambda;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> alpha;
  real<lower=0, upper=1> beta;
  real<lower=0, upper=1> relapse_clinical_immunity;
  real<lower=0, upper=1> p_long;
  real<lower=0, upper=1> p_silent;
  // real<lower=0, upper=1> eps;
  real<lower=1, upper=10> kappa;
  real<lower=-91.3125, upper=273.9375> phase; // Approx boundary adjustment as peak timing is 365.25/4 + phase
}

transformed parameters {
  vector[8] theta;
  theta[1] = lambda;
  theta[2] = alpha;
  theta[3] = beta;
  theta[4] = relapse_clinical_immunity;
  theta[5] = p_long;
  theta[6] = p_silent;
  theta[7] = kappa;
  theta[8] = phase;
  real phi = 1. / phi_inv;
  
  real incidence[n_times];
  real incidence_primary[n_times];
  real incidence_relapse[n_times];
  
  // real suitability_eval[12];
  // for (i in 1:12) {
    //   suitability_eval[i] = suitability(i * 365.25 / 12, eps, kappa, phase);
    // }
    // array[n_times+1] vector[len_y] y_extended = ode_bdf(my_ode, y0, t0, ts_extended, theta, suitability_eval, x_r, x_i);
    array[n_times+1] vector[len_y] y_extended = ode_bdf(my_ode, y0, t0, ts_extended, theta, x_r, x_i);
    for (i in 1:n_times) {
      incidence[i] = y_extended[i+1][n_compartments+3] - y_extended[i][n_compartments+3] +
      y_extended[i+1][n_compartments+4] - y_extended[i][n_compartments+4];
      incidence_primary[i] = y_extended[i+1][n_compartments+3] - y_extended[i][n_compartments+3];
      incidence_relapse[i] = y_extended[i+1][n_compartments+4] - y_extended[i][n_compartments+4];
    }
    array[n_times] vector[len_y] y;
    for (i in 1:n_times) {
      y[i] = y_extended[i+1];
    }
}

model {
  lambda ~ exponential(5);
  phi_inv ~ exponential(5);
  alpha ~ beta(alpha_shape1, alpha_shape2);
  beta ~ beta(90, 10);
  kappa ~ exponential(0.1);
  relapse_clinical_immunity ~ beta(relapse_clinical_immunity_shape1, relapse_clinical_immunity_shape2); // very little clinical immunity
  p_long ~ beta(100, 1); // almost complete long relapse
  p_silent ~ beta(p_silent_shape1, p_silent_shape2);
  
  cases ~ neg_binomial_2(incidence, phi);
  
}

generated quantities {
  int sim_cases[n_times] = neg_binomial_2_rng(incidence, phi);
}
