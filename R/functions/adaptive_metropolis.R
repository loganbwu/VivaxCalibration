source(here::here("R/functions/safe_mclapply.R"))

#' Run the adaptive Metropolis algorithm
#'
#' @param model list containing a prior and likelihood function
#' @param init vector of initial parameter values OR data frame
#' @param init_sd vector of initial parameter standard deviations before adaptation
#' @param init_covar initial covariance matrix of initial parameter proposals before adaptation
#' @param data for the model
#' @param n_iter number of sampling iterations to run
#' @param n_burnin
#' @param n_adapt
#' @param n_chains
#' @param thin proportion of sampling iterations to return
#' @param time_limit if provided in hours, estimates the number of iterations that can be performed, and the sampling stages will be scaled to fit
#' @param pars does nothing, for compatibility
#' @param include does nothing, for compatibility
metropolis_sampling = function(model, init, init_covar, data, n_iter, n_burnin=NULL, n_adapt=NULL, n_chains=1, thin=1, time_limit=NULL, previous_samples=NULL, pars=NULL, include=F, threading=TRUE) {
  
  n_samples = round(n_iter * thin)
  message("Generating samples across ", n_chains, " chains")
  
  if (n_chains > 1 & threading) {
    lapply_func = safe_mclapply
  } else {
    lapply_func = pblapply
    if (!is.null(time_limit)) {
      time_limit = time_limit / n_chains
    }
  }
  
  # Estimate required time
  n_time_estimate = 10
  start_eval = Sys.time()
  if (is.data.frame(init)) {
    # Just use the first one if we've given a dataframe of potential initial values
    time_estimate_init = head(init, 1)
  } else {
    # Otherwise use the vector
    time_estimate_init = init
  }
  for (i in seq_len(n_time_estimate)) {
    .UNUSED = model$log_prior(time_estimate_init) + model$log_likelihood(time_estimate_init)
  }
  duration_eval = as.numeric(Sys.time() - start_eval) / n_time_estimate
  if (!is.null(time_limit)) {
    total_specified_iterations = n_burnin + n_adapt + n_iter
    available_iterations = time_limit * 3600 / duration_eval
    time_to_spare = available_iterations / total_specified_iterations
    message("Time limit specified; multiplying specified iterations by ", round(time_to_spare, 1))
    n_burnin = ceiling(n_burnin * time_to_spare)
    n_adapt = ceiling(n_adapt * time_to_spare)
    n_iter = ceiling(n_iter * time_to_spare)
  }
  if (n_adapt > 0 & n_adapt < 100) {
    message("Potentially not enough adaptation iterations for an accurate covariance estimate")
  }
  
  # Benchmark to calculate estimated duration
  message("One evaluation took ", format_time(duration_eval), ". ", n_burnin + n_adapt + n_iter, " iterations per chain will take up to ", format_time(duration_eval*(n_burnin+n_adapt+n_iter)))
  
  # Perform burnin
  # proposal_covar = diag(length(init_sd)) * init_sd^2
  proposal_covar = init_covar
  if (is.data.frame(init)) {
    # If we've provided a dataframe of potential initial values, resample it for the number of chains then format as a list
    init_list = lapply_func(seq_len(n_chains), function(i) {
      init %>%
        slice_sample(n=1) %>%
        as.numeric() %>%
        setNames(names(init))
    })
  } else {
    init_list = rep(list(init), n_chains)
  }
  
  if (n_burnin > 0) {
    # Run chains
    message("Running ", n_burnin, " burnin iterations per chain...")
    burnin_output = lapply_func(init_list, function(init) {
      run_one_chain(model, init, proposal_covar, n_iter=n_burnin, thin=0.01)
    })
    
    # Get burnin samples
    burnin_samples = lapply_func(burnin_output, function(x) {
      x$sim %>%
        mutate(ess_weight = mean(x$ESS)) # just take the average ESS across all variables
    }) %>%
      bind_rows()
    
    # print(head(burnin_samples))
    message("Total ESS during burnin was ", sum(burnin_samples$ess_weight))
    
    # If burnin was unsuccessful (extremely low ESS) then terminate
    if (sum(burnin_samples$ess_weight) < 1) {
      warning("Burnin failed and resulted in an equivalent sample size close to zero; returning NULL.")
      return(NULL)
    }
    # Otherwise, update initial values using the best chains
    init_list = lapply_func(seq_len(n_chains), function(i) {
      burnin_samples %>%
        slice_sample(n=1, weight_by=ess_weight) %>%
        select(-ess_weight) %>%
        as.numeric() %>%
        setNames(names(init))
    })
  }
  
  # Run adaptation
  if (n_adapt > 0) {
    # Run chains
    message("Running ", n_adapt, " adaptation iterations per chain...")
    adapt_output = lapply_func(init_list, function(init) {
      run_one_chain(model, init, proposal_covar, n_iter=n_adapt, thin=min(1, 1000 / n_adapt))
    })
    
    # Process adaptation period
    adapt_samples = lapply(adapt_output, function(x) {
      x$sim
    }) %>%
      bind_rows()
    
    # Update proposal covariance
    # proposal_covar = 2.38^2 * cov(adapt_samples) / length(init)
    samples_covar = lapply(adapt_output, function(x) {
      x$sim
    }) %>%
      bind_rows() %>%
      cov()
    proposal_covar = 2.38^2 * samples_covar / length(init)
    
    # Update initial values for each chain
    init_list = lapply(adapt_output, function(x) {
      x$current_x
    })
  }
  
  # Run samples with new initial values and proposal covarian
  message("Running ", n_iter, " sampling iterations per chain...")
  sample_outputs = lapply_func(init_list, function(init) {
    run_one_chain(model, init, proposal_covar, n_iter=n_iter, thin=thin)
  })
  
  # message("Stopping to debug")
  
  if (any(sapply(sample_outputs, is.null))) {
    warning("Warning")
  }
  
  # Zip chain outputs together for nicer outputs
  zipped = list(data = data)
  for (element_name in names(sample_outputs[[1]])) {
    zipped[[element_name]] = lapply(sample_outputs, function(o) {
      o[[element_name]]
    })
  }
  
  
  return(structure(zipped, class = "metropolisfit"))
}

run_one_chain = function(model, init, covar, n_iter, thin=1) {
  
  n_samples = ceiling(n_iter * thin)
  
  current_x = init
  current_sample_row = 1
  samples = matrix(nrow=n_samples, ncol=length(init), dimnames=list(NULL, names(init)))
  samples_diagnostics = matrix(nrow=n_samples, ncol=4, dimnames=list(NULL, c("iteration", "accept", "ll", "lpp")))
  current_log_likelihood = model$log_likelihood(current_x)
  current_log_posterior = model$log_prior(current_x) + current_log_likelihood
  
  # Determine indices of iterations to sample
  samples_ix = round(seq(0, n_iter, length.out=n_samples+1))[-1]
  
  pb = progress::progress_bar$new(total=n_iter, clear=F, format="[:bar] :percent Elapsed: :elapsed ETA: :eta")
  for (i in seq_len(n_iter)) {
    # Make proposal
    proposal_x = MASS::mvrnorm(n=1, mu=current_x, Sigma=covar)
    
    # TESTING - Remove later
    # proposal_x[["phi"]] = 999 # why tf is this here again
    proposal_log_likelihood = model$log_likelihood(proposal_x)
    proposal_log_posterior = model$log_prior(proposal_x) + proposal_log_likelihood
    
    # Evaluate transition probability
    accept = exp(proposal_log_posterior - current_log_posterior) > runif(1)
    
    if (accept) {
      current_x = proposal_x
      current_log_likelihood = proposal_log_likelihood
      current_log_posterior = proposal_log_posterior
    }
    
    # Store result
    if (i == samples_ix[current_sample_row]) {
      samples[current_sample_row,] = current_x
      samples_diagnostics[current_sample_row,] = c(i, 1*accept, current_log_likelihood, current_log_posterior)
      current_sample_row = current_sample_row + 1
    }
    pb$tick()
  }
  samples = as_tibble(samples)
  
  # Restrict length of each chain
  if (nrow(samples) > 10000) {
    samples = samples %>%
      sample_n(10000)
  }
  samples_diagnostics = as_tibble(samples_diagnostics)
  if (nrow(samples_diagnostics) > 10000) {
    samples_diagnostics = samples_diagnostics %>%
      sample_n(10000)
  }
  if (n_samples > 1) {
    ESS = sapply(samples,
                 function(x) {
                   unname(coda::effectiveSize(x))
                 })
    # print(ESS)
  } else {
    ESS = 0 # not valid if there are no samples
  }
  return(list(sim = samples, # called `sim` to align with stan
              sim_diagnostics = samples_diagnostics,
              accept = mean(samples_diagnostics$accept),
              ESS = ESS,
              current_x = current_x))
}


extract.metropolisfit = function(fit, type="incidence", n_samples=100, alpha=0.05, threading=TRUE) {
  if (type == "incidence") {
    all_sims = bind_rows(fit$sim)
    # Address transformed variables
    if ("s" %in% names(all_sims)) {
      all_sims$p_long = 1 - all_sims$s
      all_sims$s = NULL
    }
    if ("p" %in% names(all_sims)) {
      all_sims$p_silent = 1 - all_sims$p
      all_sims$p = NULL
    }
    sample_ix = sample(1:nrow(all_sims), min(n_samples, nrow(all_sims)))
    
    x_r = with(fit$data, c(N=N, gamma_d=gamma_d, gamma_l=gamma_l, delta=delta, f=f, r=r, eps=eps))
    x_i = with(fit$data, c(n_dormant=n_dormant))
    
    run_my_ode = function(t, state, parms) {
      return(list(my_ode(t, state, parms, x_r, x_i)))
    }
    
    if (threading) {
      lapply_func = pbmclapply
    } else {
      lapply_func = pblapply
    }
    
    incidence = lapply_func(sample_ix, function(ix) {
      x = unlist(all_sims[ix,])
      mean_ode = deSolve::ode(fit$data$y0, c(fit$data$t0, 0, fit$data$ts), run_my_ode, parms=x, maxsteps=1e6)
      clinical_relapse = (mean_ode[,"ClinicalRelapse"] - lag(mean_ode[,"ClinicalRelapse"]))[-c(1,2)]
      if ("ClinicalPrimary" %in% colnames(mean_ode)) {
        # Run for v11
        clinical_primary = (mean_ode[,"ClinicalPrimary"] - lag(mean_ode[,"ClinicalPrimary"]))[-c(1,2)]
        clinical_incidence = clinical_primary + clinical_relapse
        
        data = tibble(time = fit$data$ts,
               clinical_primary = clinical_primary,
               clinical_relapse = clinical_relapse,
               clinical_incidence = clinical_incidence,
               # cases_lower = qpois(alpha/2, lambda=clinical_incidence),
               # cases_upper = qpois(1-(alpha/2), lambda=clinical_incidence))
               cases_lower = qnbinom(alpha/2, mu=clinical_incidence, size=x[["phi"]]),
               cases_upper = qnbinom(1-(alpha/2), mu=clinical_incidence, size=x[["phi"]]))
      } else {
        # Run for v12
        clinical_immediateprimary = (mean_ode[,"ClinicalImmediatePrimary"] - lag(mean_ode[,"ClinicalImmediatePrimary"]))[-c(1,2)]
        clinical_delayedprimary = (mean_ode[,"ClinicalDelayedPrimary"] - lag(mean_ode[,"ClinicalDelayedPrimary"]))[-c(1,2)]
        clinical_incidence = clinical_immediateprimary + clinical_delayedprimary + clinical_relapse
        
        data = tibble(time = fit$data$ts,
               clinical_immediateprimary = clinical_immediateprimary,
               clinical_delayedprimary = clinical_delayedprimary,
               clinical_relapse = clinical_relapse,
               clinical_incidence = clinical_incidence,
               # cases_lower = qpois(alpha/2, lambda=clinical_incidence),
               # cases_upper = qpois(1-(alpha/2), lambda=clinical_incidence))
               cases_lower = qnbinom(alpha/2, mu=clinical_incidence, size=x[["phi"]]),
               cases_upper = qnbinom(1-(alpha/2), mu=clinical_incidence, size=x[["phi"]]))
      }
      
      return(data)
    }) %>%
      bind_rows(.id = "ix")
    return(incidence)
  }
}

extract = function(...) UseMethod("extract")

rhat.metropolisfit = function(fit) {
  sim_matrices = lapply(names(fit$sim[[1]]) %>% setNames({.}),
                        function(param) {
                          x = lapply(fit$sim, function(s) {
                            s[[param]]
                          }) %>%
                            setNames(seq_along(fit$sim)) %>%
                            bind_cols() %>%
                            as.matrix() %>%
                            rstan::Rhat()
                        })
}

rhat = function(...) UseMethod("rhat")

#' @param dt time difference in seconds
format_time = function(dt) {
  if (dt < 1) {
    formatted = paste(round(dt*1000), "ms")
  } else if (dt < 100) {
    formatted = paste(round(dt, 1), "sec")
  } else if (dt < (3600)) {
    formatted = paste(round(dt/60, 1), "min")
  } else if (dt < (3600*24)) {
    formatted = paste(round(dt/3600, 1), "hr")
  } else {
    formatted = paste(round(dt/(3600*24), 1), "d")
  }
  return(formatted)
}