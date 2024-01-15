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
  
  .data_agg = aggregate_data(.data)
  
  optim_ext = optimizing(stan_model_champagne2022_poisson, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022_poisson,
                  data = .data_agg,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  return(rstan::extract(.fit, c("lambda"))$lambda %>% as.numeric())
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
  
  .data_agg = aggregate_data(.data)
  
  optim_ext = optimizing(stan_model_champagne2022, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022,
                  data = .data_agg,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
                  control = list(max_treedepth = 4),
                  refresh = 0)
  
  return(rstan::extract(.fit, c("lambda"))$lambda)
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
  
  return(rstan::extract(.fit, c("lambda"))$lambda)
}

seasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, .eps, true_lambda=0.01, true_phi_inv=0.1) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$omega = .omega
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
  
  return(rstan::extract(.fit, c("lambda"))$lambda)
}

extended_seasonal_sol = function(.cases, .population_size, .alpha, .beta, .omega, true_eps, true_lambda=0.01, true_phi_inv=0.1, true_kappa=1, true_phase=0) {
  .data = data_consts
  .data$cases = .cases
  .data$population_size = .population_size
  .data$alpha = .alpha
  .data$beta = .beta
  .data$omega = .omega
  .data$eps = true_eps
  .data$phase = true_phase
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
  
  return(rstan::extract(.fit, c("lambda"))$lambda)
}

methods = tibble(
  method = c(
    "lambda_nonseasonal_poisson",
    "lambda_nonseasonal",
    "lambda_seasonal_poisson",
    "lambda_seasonal",
    "lambda_seasonal_ext")
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