// These models are for testing a variable-length chain of delays in an ODE model

functions {
  vector my_ode(real time, vector y, vector theta, real[] x_r, int[] x_i) {
    
    // Unpack theta
    real lambda = theta[1];
    
    // Unpack x_r
    real alpha = x_r[1];
    real beta = x_r[2];
    real relapse_clinical_immunity = x_r[3];
    real gamma_d = x_r[4];
    real gamma_l = x_r[5];
    real delta = x_r[6];
    real phi = x_r[7];
    real f = x_r[8];
    real r = x_r[9];
    real p_long = x_r[10];
    real p_silent = x_r[11];
    real population_size = x_r[12];
    
    // Unpack x_i
    int n_dormant = x_i[1];
    
    // Convenience quantities
    int n_stages = n_dormant + 1;
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
    // real clearance = c(rep(gamma_d, n_dormant), gamma_l);
    real clear_d = gamma_d;
    real clear_l = gamma_l;
    real advance = delta*n_dormant;
    real relapse = f;
    real long_hyp = p_long;
    real short_hyp = 1-long_hyp;
    // array[n_dormant] int inactive = 1:n_dormant;
    int inactive[n_dormant];
    for (i in 1:n_dormant) inactive[i] = i;
    int active = n_dormant + 1;
    
    // Give reasonable names to variables
    real S0 = y[1];
    real I0 = y[2];
    array[n_stages] real Sl;
    for (i in 1:n_stages) {
      Sl[i] = y[2+i];
    }
    array[n_stages] real Scl;
    for (i in 1:n_stages) {
      Scl[i] = y[2+n_stages+i];
    }
    array[n_stages] real Icl;
    for (i in 1:n_stages) {
      Icl[i] = y[2+2*n_stages+i];
    }
    real Incidence = y[2+3*n_stages+1];
    
    real FOI = lambda * (I0 + sum(Icl));
    real infect = FOI;
    
    vector[num_elements(y)] dydt;
    
    // Compute derivatives
    dydt[1] = -S0*infect*(short_hyp*(treatedprimary*incomplete*incomplete+untreatedprimary) + long_hyp*primary*(treatedprimary*incomplete+untreatedprimary)) +
    sum(Sl)*(infect*(short_hyp*treatedprimary*complete+long_hyp*(primary*treatedprimary*complete+silent*treatedsilent*complete))) + Sl[n_stages]*relapse*treatedprimary*complete + sum(Sl[1:n_dormant])*clear_d + Sl[active]*clear_l +
    sum(Scl)*infect*(short_hyp*treatedprimary + long_hyp*(primary*treatedprimary*complete + silent*treatedsilent*complete)) + Scl[n_stages]*relapse*treatedrelapse*complete + sum(Scl[1:n_dormant])*clear_d + Scl[active]*clear_l;
    
    for (i in 1:n_stages) {
      
    }
    
    return dydt;
  }
}

data {
  int<lower=1> n_times;
  real t0;
  real ts[n_times];
  int cases[n_times];
  
  real alpha;
  real beta;
  real relapse_clinical_immunity;
  real gamma_d;
  real gamma_l;
  real delta;
  real phi;
  real f;
  real r;
  real p_long;
  real p_silent;
  real population_size;
  
  int n_dormant;
}

transformed data {
  real x_r[13] = {
    alpha,
    beta,
    relapse_clinical_immunity,
    gamma_d,
    gamma_l,
    delta,
    phi,
    f,
    r,
    p_long,
    p_silent,
    population_size
  };
  
  int x_i[1] = {
    n_dormant
  };
  
  int n_stages = n_dormant + 1;
  int n_compartments = 2 + 3*n_stages;
  // real<lower=0, upper=1> y0[n_compartments];
  vector[n_compartments + 1] y0; # last element is cumulative incidence
}

parameters {
  real lambda;
}

transformed parameters {
  vector[1] theta;
  theta[1] = lambda;
  
  // int n_compartments = 2 + 3*n_stages;
  array[n_times] vector[n_compartments] y = ode_bdf(my_ode, y0, t0, ts, theta, x_r, x_i);
}

model {
  // sigma ~ cauchy(0, 5);
  theta[1] ~ uniform(0, 1);
  for (i in 1:n_times-1) {
    real incidence = y[i+1][num_elements(y[i+1])] - y[i][num_elements(y[i])];
    cases[i] ~ poisson(incidence);
  }
}
