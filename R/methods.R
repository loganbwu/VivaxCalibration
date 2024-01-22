model_champagne2022 = "stan/champagne2022.stan"
stan_model_champagne2022 = stan_model(model_champagne2022)
model_champagne2022_poisson = "stan/champagne2022_poisson.stan"
stan_model_champagne2022_poisson = stan_model(model_champagne2022_poisson)

model_champagne2022_seasonal = "stan/champagne2022_seasonal.stan"
stan_model_champagne2022_seasonal = stan_model(model_champagne2022_seasonal)
model_champagne2022_seasonal_poisson = "stan/champagne2022_seasonal_poisson.stan"
stan_model_champagne2022_seasonal_poisson = stan_model(model_champagne2022_seasonal_poisson)

model_champagne2022_seasonal_ext = "stan/champagne2022_seasonal_ext2.stan"
stan_model_champagne2022_seasonal_ext = stan_model(model_champagne2022_seasonal_ext)

#' For all Stan solution functions, we allow it to initialise at the true parameters and avoid it getting stuck.
poisson_nonseasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, true_lambda=0.01) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$omega = .omega
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # .data_agg = aggregate_data(.data)
  
  optim_ext = optimizing(stan_model_champagne2022_poisson, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022_poisson,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda %>% as.numeric())
  return(.fit)
}

nonseasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, true_lambda=0.01, true_phi_inv=0.1) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$omega = .omega
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # .data_agg = aggregate_data(.data)
  
  optim_ext = optimizing(stan_model_champagne2022, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

poisson_seasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, .eps, true_lambda=0.01, true_phi_inv=0.1) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$omega = .omega
  .data$eps = .eps
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  optim_ext = optimizing(stan_model_champagne2022_seasonal_poisson, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022_seasonal_poisson,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

seasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, .eps, true_lambda=0.01, true_phi_inv=0.1) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  # .data$omega = .omega
  .data$eps = .eps
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  
  optim_ext = optimizing(stan_model_champagne2022_seasonal, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022_seasonal,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

extended_seasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  # .data$omega = .omega
  # .data$eps = true_eps
  # .data$phase = true_phase
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # Try to get a good initial value
  optim_ext = optimizing(stan_model_champagne2022_seasonal_ext, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv", "eps", "kappa")])
  
  .fit = sampling(stan_model_champagne2022_seasonal_ext,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

methods = tibble(
  method = c(
    "lambda_nonseasonal_poisson",
    "lambda_nonseasonal_negbin",
    "lambda_seasonal_poisson",
    "lambda_seasonal_negbin",
    "lambda_seasonal_negbin_ext")
)

# direct_sol = function(.cases, .population_size, .alpha, .beta, .omega) {
#   x = tibble(h = mean(.cases) / 30.4 / .population_size, # daily incidence per person
#              alpha = .alpha,
#              beta = .beta,
#              rho = .alpha,
#              omega = .omega,
#              prop_import = 0) %>%
#     calibrate_vivax_equilibrium(f=data_consts$f, gamma=data_consts$gammal, r=data_consts$r, return.all = TRUE)
#   
#   if (x$lambda == -2) {
#     return(NA_real_)
#   } else {
#     return(x$lambda)
#   }
# }

#' @param i index
run_scenario_method = function(i) {
  out_path = file.path("../run_scenario_method", paste0("row_", i, ".rds"))
  dir.create("../run_scenario_method")
  if (file.exists(out_path)) {
    result = read_rds(out_path)
    if (!is.null(result$fit)) {
      return(result)
    }
  }
  start = Sys.time()
  .method = data_scenarios_long$method[i]
  row = data_scenarios_long[i,]
  fit = withTimeout({
    if (.method == "lambda_nonseasonal_poisson") {
      poisson_nonseasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$transmission_rates)
    }
    else if (.method == "lambda_nonseasonal_negbin") {
      nonseasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, transmission_rates)
    }
    else if (.method == "lambda_seasonal_poisson") {
      poisson_seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    }
    else if (.method == "lambda_seasonal_negbin") {
      seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    }
    else if (.method == "lambda_seasonal_negbin_ext") {
      extended_seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    } else {
      stop("Method invalid")
    }
  }, timeout=timelimit_per_run, onTimeout="warning")
  end = Sys.time()
  
  result = list(fit = fit, time = end - start)
  write_rds(result, out_path, compress="gz")
  
  return(result)
}

extract_lambda = function(.result) {
  .fit = .result$fit
  param = rstan::extract(.fit, c("lambda"))[[1]]
  return(param)
}

extract_other = function(.result, var) {
  .fit = .result$fit
  param = rstan::extract(.fit, var)[[1]]
  return(param)
}

extract_phi_inv = function(.result) {
  .fit = .result$fit
  if (.fit@mode == 0) {
    if ("phi_inv" %in% names(.fit)) {
      param = rstan::extract(.fit, c("phi_inv"))[[1]]
    } else {
      param=NULL
    }
  } else {
    param = NULL
  }
  return(param)
}

extract_incidence = function(.result) {
  .fit = .result$fit
  if (.fit@mode != 0) {
    return(NULL)
  }
  incidence = rstan::extract(.fit, "incidence")[[1]]
  ix = sample(seq_len(dim(incidence)[1]), n_traces, replace=T)
  ts_sample = as_tibble(t(incidence[ix,])) %>%
    mutate(j = row_number()) %>%
    pivot_longer(-j, names_to = "trace", values_to = "incidence")
  return(ts_sample)
}

extract_runtime = function(.result) {
  return(.result$time)
}