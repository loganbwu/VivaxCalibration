library(tidyverse)
library(rstan)
library(rstansim) # devtools::install_github("ewan-keith/rstansim")
library(parallel)
library(patchwork)
library(pbmcapply)
library(pbapply)
source("R/constants.R")
n_cores = parallelly::availableCores()
options(mc.cores = n_cores)
message("Running on ", n_cores, " cores")
rstan_options(auto_write = TRUE)

n_years = 5
n_iter = 250 # should be at least 500
n_chains = 2
n_repetitions = 10 # how many times to duplicate each scenario
cores_per_sampler = 1 # set to n_chains if not running lots of scenarios
limit_runs = Inf # set to a finite number for testing, or Inf to run all

model_champagne2022 = "Rmd/stan/champagne2022.stan"
stan_model_champagne2022 = stan_model(model_champagne2022)
model_champagne2022_poisson = "Rmd/stan/champagne2022_poisson.stan"
stan_model_champagne2022_poisson = stan_model(model_champagne2022_poisson)

model_champagne2022_seasonal = "Rmd/stan/champagne2022_seasonal.stan"
stan_model_champagne2022_seasonal = stan_model(model_champagne2022_seasonal)
model_champagne2022_seasonal_poisson = "Rmd/stan/champagne2022_seasonal_poisson.stan"
stan_model_champagne2022_seasonal_poisson = stan_model(model_champagne2022_seasonal_poisson)

model_champagne2022_seasonal_ext = "Rmd/stan/champagne2022_seasonal_ext.stan"
stan_model_champagne2022_seasonal_ext = stan_model(model_champagne2022_seasonal_ext)

# perform simulation study
dt = years/annual_subdivisions
t0 = -50*years
t = seq_len(n_years*annual_subdivisions) * dt
n_times = length(t)
N = 1000 # population size

#initial conditions
I_init = 0.01
y0 = c(Il=0, I0=I_init, Sl=0, S0=1-I_init, CumulativeInfections=0)

# constants for Stan
data_consts = list(n_times = n_times+1,
                   y0 = y0,
                   t0 = t0,
                   ts = seq_len(n_times+1) * dt,
                   N = N,
                   cases = rep(99999, n_times+1),
                   r = 1./60, # r
                   gammal = 1./223, # gammal
                   f = 1./72, # f
                   alpha = 0.21, # alpha
                   beta = 0.66, # beta
                   rho = 0.21, # rho
                   delta = 0,
                   eps = 0,
                   kappa = 1,
                   phase = 0
)

# ascertainment_rates = c(0.25, 0.5, 0.75, 1)
# radical_cure_rates = seq(0, 1, by=0.2)
radical_cure_rates = 0.66
ascertainment_rates = data_consts$alpha
seasonality_ratio = seq(0, 1, length.out=3)
# radical_cure_rates = data_consts$beta
transmission_rates = seq(0.01, 0.02, by=0.005)
importation_rate = 0 # because constant importation makes less sense in seasonal transmission
population_size = N

data_scenarios = expand_grid(
  ascertainment_rates,
  radical_cure_rates,
  transmission_rates,
  importation_rate,
  seasonality_ratio,
  population_size
) %>%
  mutate(ID = LETTERS[row_number()], .before=0)

.simulate_cases = function(alpha, beta, lambda, delta, eps, N) {
  data = data_consts
  data$alpha = alpha
  data$beta = beta
  data$delta = delta
  data$eps = eps # 0=full seasonality, 1=no seasonality
  data$N = N
  
  real_params = list(lambda=lambda, phi_inv=0.1)
  synth_data = simulate_data(
    file = model_champagne2022_seasonal,
    data_name = "dummy_data",
    input_data = data,
    param_values = real_params,
    vars = c("ts", "sim_cases", "susceptible")
  )
  synth_data_rds = readRDS(synth_data$datasets[1])
  indx <- sapply(synth_data_rds, length)
  synth_df = lapply(synth_data_rds, function(x) {length(x) = max(indx); x}) %>%
    as.data.frame() %>%
    drop_na()
}
# Add cases onto a dataframe of scenarios based on its parameter columns
simulate_cases = function(.scenarios) {
  cases_scenarios = mclapply(seq_len(nrow(.scenarios)), function(i) {
    dat = .scenarios[i,]
    x = .simulate_cases(dat$ascertainment_rates, dat$radical_cure_rates, dat$transmission_rates, dat$importation_rate, dat$seasonality_ratio, dat$population_size)
  })
  
  .scenarios$cases = lapply(cases_scenarios, function(x) {x$cases})
  .scenarios$ts = lapply(cases_scenarios, function(x) {x$ts})
  .scenarios
}
message("Running ", nrow(data_scenarios), " scenarios for ", n_repetitions, " repetitions")

data_scenarios = data_scenarios %>%
  slice(rep(1:n(), each = n_repetitions)) %>%
  simulate_cases()

message("Simulated cases for ", nrow(data_scenarios), " instances")


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
  
  if (true_lambda == 0.01) {
    return(NULL)
  }
  
  optim_ext = optimizing(stan_model_champagne2022_poisson, data=.data)
  .theta_init = as.list(optim_ext$par[c("lambda", "phi_inv")])
  
  .fit = sampling(stan_model_champagne2022_poisson,
                  data = .data_agg,
                  iter = n_iter,
                  chains = n_chains,
                  init = rep(list(.theta_init), n_chains), # Start from MLE solution
                  cores = cores_per_sampler,
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

data_scenarios_long = data_scenarios %>%
  tidyr::crossing(methods) %>%
  head(limit_runs)

message("Executing a total of ", nrow(data_scenarios_long), " fits")

run_scenario_method = function(i) {
  start = Sys.time()
  .method = data_scenarios_long$method[i]
  est = with(data_scenarios_long, {
    if (.method == "lambda_nonseasonal_poisson") {
      poisson_nonseasonal_sol(cases[[i]], population_size[i], ascertainment_rates[i], radical_cure_rates[i], 1, transmission_rates[i])
    }
    else if (.method == "lambda_nonseasonal") {
      nonseasonal_sol(cases[[i]], population_size[i], ascertainment_rates[i], radical_cure_rates[i], 1, transmission_rates[i])
    }
    else if (.method == "lambda_seasonal_poisson") {
      poisson_seasonal_sol(cases[[i]], population_size[i], ascertainment_rates[i], radical_cure_rates[i], 1, seasonality_ratio[i], transmission_rates[i])
    }
    else if (.method == "lambda_seasonal") {
      seasonal_sol(cases[[i]], population_size[i], ascertainment_rates[i], radical_cure_rates[i], 1, seasonality_ratio[i], transmission_rates[i])
    }
    else if (.method == "lambda_seasonal_ext") {
      extended_seasonal_sol(cases[[i]], population_size[i], ascertainment_rates[i], radical_cure_rates[i], 1, seasonality_ratio[i], transmission_rates[i])
    } else {
      stop("Method invalid")
    }
  })
  end = Sys.time()
  
  return(list(estimate = est, time = end - start))
}

tictoc::tic()
estimates_all = pbmclapply(seq_len(nrow(data_scenarios_long)), run_scenario_method)
data_scenarios_long$estimate = lapply(estimates_all, function(x) {x$estimate})
data_scenarios_long$time = lapply(estimates_all, function(x) {x$time})
tictoc::toc()

# data_scenarios = run_all(data_scenarios)
workspace_filename = format(Sys.time(), "%Y%M%d %H%M%S.Rdata")
save.image(workspace_filename)

message("Done.")
