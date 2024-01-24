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
poisson_nonseasonal_sol = function(data_consts, .cases, .population_size, .alpha, .beta, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0, refresh=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # .data_agg = aggregate_data(.data)
  
  # optim_ext = optimizing(stan_model_champagne2022_poisson, data=.data)
  # .theta_init = as.list(optim_ext$par[c("lambda")])
  .theta_init = list(lambda = true_lambda)
  
  .fit = sampling(stan_model_champagne2022_poisson,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = refresh)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda %>% as.numeric())
  return(.fit)
}

nonseasonal_sol = function(data_consts, .cases, .population_size, .alpha, .beta, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0, refresh=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # .data_agg = aggregate_data(.data)
  
  # optim_ext = optimizing(stan_model_champagne2022, data=.data)
  # .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  .theta_init = list(lambda = true_lambda,
                     phi_inv = true_phi_inv)
  
  .fit = sampling(stan_model_champagne2022,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = refresh)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

poisson_seasonal_sol = function(data_consts, .cases, .population_size, .alpha, .beta, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0, refresh=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$eps = true_eps
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # optim_ext = optimizing(stan_model_champagne2022_seasonal_poisson, data=.data)
  # .theta_init = as.list(optim_ext$par[c("lambda")])
  .theta_init = list(lambda = true_lambda)
  
  .fit = sampling(stan_model_champagne2022_seasonal_poisson,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = refresh)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

seasonal_sol = function(data_consts, .cases, .population_size, .alpha, .beta, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0, refresh=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$eps = true_eps
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  
  # optim_ext = optimizing(stan_model_champagne2022_seasonal, data=.data)
  # .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  .theta_init = list(lambda = true_lambda,
                     phi_inv = true_phi_inv)
  
  .fit = sampling(stan_model_champagne2022_seasonal,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = refresh)
  
  # return(rstan::extract(.fit, c("lambda"))$lambda)
  return(.fit)
}

extended_seasonal_sol = function(data_consts, .cases, .population_size, .alpha, .beta, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0, refresh=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  # .data$eps = true_eps
  # .data$phase = true_phase
  .data$n_times = length(.data$cases)
  .data$ts = .data$ts[seq_len(.data$n_times)]
  
  # Try to get a good initial value
  # optim_ext = optimizing(stan_model_champagne2022_seasonal_ext, data=.data)
  # .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv", "eps", "kappa")])
  .theta_init = list(lambda = true_lambda,
                     phi_inv = true_phi_inv,
                     eps = true_eps,
                     kappa = true_kappa)
  
  .fit = sampling(stan_model_champagne2022_seasonal_ext,
                  data = .data,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = refresh)
  
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
run_scenario_method = function(i, force=F, refresh=0) {
  # folder = "../run_scenario_method"
  # out_path = file.path(folder, paste0("row_", i, ".rds"))
  # if (!dir.exists(folder)) {
  #   dir.create(folder)
  # }
  # if (file.exists(out_path) & !force) {
  #   result = read_rds(out_path)
  #   if (!is.null(result$fit)) {
  #     return(result)
  #   }
  # }
  start = Sys.time()
  .method = data_scenarios_long$method[i]
  row = data_scenarios_long[i,]
  fit = withTimeout({
    if (.method == "lambda_nonseasonal_poisson") {
      f = poisson_nonseasonal_sol
    }
    else if (.method == "lambda_nonseasonal_negbin") {
      f = nonseasonal_sol
    }
    else if (.method == "lambda_seasonal_poisson") {
      f = poisson_seasonal_sol
    }
    else if (.method == "lambda_seasonal_negbin") {
      f = seasonal_sol
    }
    else if (.method == "lambda_seasonal_negbin_ext") {
      f = extended_seasonal_sol
    } else {
      stop("Method invalid")
    }
    result = f(data_consts = data_consts,
               .cases = row$cases[[1]],
               .population_size = row$population_size,
               .alpha = row$ascertainment_rates,
               .beta = row$radical_cure_rates,
               true_eps = row$seasonality_ratio,
               true_lambda = row$transmission_rates,
               true_phi_inv = 0.1,
               true_kappa = 1,
               true_phase = 0,
               refresh = 0)
  }, timeout=timelimit_per_run, onTimeout="warning")
  end = Sys.time()
  
  result = list(fit = fit, time = end - start)
  # write_rds(result, out_path, compress="gz")
  
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