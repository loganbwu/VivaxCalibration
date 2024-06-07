# Relies on VivaxODE/stan/temperate_v11.stan
stancode = "
functions {
  // Unused except for resimulation
  real suitability(real t, real eps, real kappa, real phase) {
    real omega = eps + (1-eps)*pi()/exp(lbeta(0.5, kappa+0.5))*((1+sin(2*pi()*(t - phase)/365.25))/2)^kappa;
    return(omega);
  }
  
  // v11
  vector my_ode(real time, vector y, vector theta, real[] x_r, int[] x_i) {
    
    // Unpack theta
    real alpha = theta[1];
    real beta = theta[2];
    real lambda = theta[3];
    // real phi = theta[4]; // Unused in mean ODE
    real kappa = theta[5];
    real phase = theta[6];
    real p_long = theta[7];
    real p_silent = theta[8];
    real p_RCI = theta[9];
    
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
    real treatedrelapse = treatedprimary * (1-p_RCI); // This should change and/or be a fit variable
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
    real infect = lambda * (eps + (1-eps)*pi()/exp(lbeta(0.5, kappa+0.5))*((1+sin(2*pi()*(time - phase)/365.25))/2)^kappa) * (I0 + sum(Icl));
    
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
"
rstan::expose_stan_functions(stanc(model_code = stancode))

make_log_prior = function(.data) {
  log_prior_func = function(x) {
    prior_alpha = dbeta(x[["alpha"]], .data$alpha_shape1, .data$alpha_shape2, log=T)
    prior_beta = dbeta(x[["beta"]], 90, 10, log=T)
    prior_lambda = dexp(x[["lambda"]], 5, log=T)
    prior_phi = dexp(1/x[["phi"]], 5, log=T)
    prior_kappa = dexp(x[["kappa"]], 0.1, log=T)
    # prior_phase = 0
    prior_p_long = dbeta(x[["p_long"]], .data$p_long_shape1, .data$p_long_shape2, log=T)
    prior_p_silent = dbeta(x[["p_silent"]], .data$p_silent_shape1, .data$p_silent_shape2, log=T)
    prior_p_RCI = dbeta(x[["p_RCI"]], .data$p_RCI_shape1, .data$p_RCI_shape2, log=T)
    return(prior_alpha + prior_beta + prior_lambda + prior_phi + prior_kappa + prior_p_long + prior_p_silent + prior_p_RCI)
  }
  return(log_prior_func)
}

make_log_likelihood = function(.data) {
  x_r = with(.data, c(N=N, gamma_d=gamma_d, gamma_l=gamma_l, delta=delta, f=f, r=r, eps=eps))
  x_i = with(.data, c(n_dormant=n_dormant))
  
  #' Stan version of ODE
  run_my_ode = function(t, state, parms) {
    return(list(my_ode(t, state, parms, x_r, x_i)))
  }
  
  # Define parameter bounds
  bounds = list(
    lower = c(alpha=0, beta=0, lambda=0, phi=0, kappa=0, phase=0, p_long=0, p_silent=0, p_RCI=0), 
    upper = c(alpha=1, beta=1, lambda=Inf, phi=Inf, kappa=Inf, phase=365.25, p_long=1, p_silent=1, p_RCI=1))
  
  log_likelihood_func = function(x) {
    if (any(x < bounds$lower) | any(x > bounds$upper)) {
      return(-Inf)
    }
    # Solve the ODE
    mean_ode = deSolve::ode(.data$y0, c(.data$t0, 0, .data$ts), run_my_ode, parms=x, method="bdf")
    clinical_incidence = mean_ode[,"ClinicalPrimary"] + mean_ode[,"ClinicalRelapse"]
    mean_cases = pmax(0, (clinical_incidence - lag(clinical_incidence))[-c(1, 2)])
    return(sum(dnbinom(.data$cases, mu = mean_cases, size = x[["phi"]], log = T)))
  }
  return(log_likelihood_func)
}

make_model = function(.data) {
  model = list(
    log_prior = make_log_prior(.data),
    log_likelihood = make_log_likelihood(.data)
  )
  return(model)
}