library(R.utils)
library(tidyverse)
library(rstan)
library(rstansim) # devtools::install_github("ewan-keith/rstansim")
library(parallel)
library(patchwork)
library(pbmcapply)
library(pbapply)
library(readr)
source("R/constants.R")
source("R/methods.R")
n_cores = parallelly::availableCores()
options(mc.cores = n_cores)
message("Running on ", n_cores, " cores")
rstan_options(auto_write = TRUE)

n_years = 5
n_iter = 500 # should be at least 500
n_chains = 2
n_repetitions = 5 # how many times to duplicate each scenario
cores_per_sampler = 1 # set to n_chains if not running lots of scenarios
limit_runs = Inf # set to a finite number for testing, or Inf to run all
timelimit_per_run = 60*0.1#60 * 30

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

.simulate_cases = function(alpha=0.5, beta=0.5, lambda=1, delta=0, eps=0, N=100, index=0) {
  data = data_consts
  data$alpha = alpha
  data$beta = beta
  data$delta = delta
  data$eps = eps # 0=full seasonality, 1=no seasonality
  data$N = N
  
  real_params = list(lambda=lambda, phi_inv=0.1)
  
  if (index == 0) {
    index = sample.int(999999999, 1)
  }
  synth_data = suppressMessages(simulate_data(
    file = model_champagne2022_seasonal,
    path = "sim_data",
    data_name = paste0("data_", index),
    input_data = data,
    param_values = real_params,
    vars = c("ts", "sim_cases", "susceptible")
  ))
  # print(synth_data$datasets[1])
  synth_data_rds = readRDS(synth_data$datasets[1])
  file.remove(synth_data$datasets[1])
  indx <- sapply(synth_data_rds, length)
  synth_df = lapply(synth_data_rds, function(x) {length(x) = max(indx); x}) %>%
    as.data.frame() %>%
    drop_na()
}
# Add cases onto a dataframe of scenarios based on its parameter columns
simulate_cases = function(.scenarios) {
  cases_scenarios = mclapply(seq_len(nrow(.scenarios)), function(i) {
    dat = .scenarios[i,]
    x = .simulate_cases(dat$ascertainment_rates, dat$radical_cure_rates, dat$transmission_rates, dat$importation_rate, dat$seasonality_ratio, dat$population_size, index=i)
  })
  
  .scenarios$cases = lapply(cases_scenarios, function(x) {x$cases})
  .scenarios$ts = lapply(cases_scenarios, function(x) {x$ts})
  .scenarios
}
message("Simulating ", nrow(data_scenarios), " scenarios for ", n_repetitions, " repetitions")

data_scenarios = data_scenarios %>%
  slice(rep(1:n(), each = n_repetitions)) %>%
  simulate_cases()
message("Simulated cases for ", nrow(data_scenarios), " instances")

data_scenarios_long = data_scenarios %>%
  tidyr::crossing(methods) %>%
  head(limit_runs)

message("Executing a total of ", nrow(data_scenarios_long), " fits")

#' @param i index
run_scenario_method = function(i) {
  start = Sys.time()
  .method = data_scenarios_long$method[i]
  row = data_scenarios_long[i,]
  est = withTimeout({
    if (.method == "lambda_nonseasonal_poisson") {
      poisson_nonseasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$transmission_rates)
    }
    else if (.method == "lambda_nonseasonal") {
      nonseasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, transmission_rates)
    }
    else if (.method == "lambda_seasonal_poisson") {
      poisson_seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    }
    else if (.method == "lambda_seasonal") {
      seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    }
    else if (.method == "lambda_seasonal_ext") {
      extended_seasonal_sol(row$cases[[1]], row$population_size, row$ascertainment_rates, row$radical_cure_rates, 1, row$seasonality_ratio, row$transmission_rates)
    } else {
      stop("Method invalid")
    }
  }, timeout=timelimit_per_run, onTimeout="warning")
  end = Sys.time()
  
  result = list(estimate = est, time = end - start)
  if (!dir.exists("run_scenario_method")) dir.create("run_scenario_method")
  out_path = file.path("run_scenario_method", paste0("row_", i, ".rds"))
  write_rds(result, out_path, compress="gz")
  
  return(result)
}

tictoc::tic()
# for (i in seq_len(nrow(data_scenarios_long))) {
#   print(i)
#   x = run_scenario_method(i)
# }
estimates_all = pbmclapply(seq_len(nrow(data_scenarios_long)), run_scenario_method)
print(lapply(estimates_all, class))
data_scenarios_long$estimate = lapply(estimates_all, function(x) {x$estimate})
data_scenarios_long$time = lapply(estimates_all, function(x) {x$time})
tictoc::toc()

# data_scenarios = run_all(data_scenarios)
workspace_filename = format(Sys.time(), "%Y%M%d %H%M%S.Rdata")
save.image(workspace_filename)

message("Done.")
