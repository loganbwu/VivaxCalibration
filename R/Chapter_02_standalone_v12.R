start_time = Sys.time()

library(optparse)

option_list <- list(
  make_option(c("-h", "--hours"), default = 8,
              help = "Target total runtime in hours"),
  make_option(c("-i", "--index"), default = "2",
              help = "Scenario index/indices to run (e.g., '3' or '1,3,5')"),
  make_option(c("-c", "--chains"), default = 1,
              help = "Number of chains/cores")
)

# This gracefully handles both command line and interactive use
opt <- parse_args(OptionParser(option_list = option_list))
max_hours = as.numeric(opt$hours)
# Parse scenario indices
if (!is.null(opt$index)) {
  scenario_indices <- as.integer(strsplit(opt$index, ",")[[1]])
} else {
  scenario_indices <- NA
}
chains = as.integer(opt$chains)


start_time = Sys.time()

library(R.utils)
library(tidyverse)
library(rstan)
library(parallel)
library(pbmcapply)
library(pbapply)
library(memoise)
library(here)
library(RColorBrewer)
source(here("R", "constants.R"))
source(here("R", "load_functions.R"))

rstan_options(auto_write = TRUE, threads_per_chain = 1)

# n_cores = parallelly::availableCores() - 1
n_cores = chains
options(mc.cores = n_cores)
message("Running on ", n_cores, " cores")

# Load model and metropolis algorithm (copied from VivaxODE project folder)
source(here("R", "models/temperate_v12.R"))
source(here("R", "functions/adaptive_metropolis.R"))
# Load priors for each scenario
source(here("R", "priors.R"))
my_state_init = state_init_4
start_time = Sys.time()

china_selections = tribble(
  ~Region, ~min, ~max,
  "Dengzhou", "2004-01-01", "2009-01-01",
  "Guantang", "1977-01-01", "1982-01-01",
  "Huangchuan", NA, NA,
  "Xiayi", NA, NA
) %>%
  mutate(min = as.Date(min),
         max = as.Date(max))

china_data = read_rds(here("Rmd", "china_data.rds"))

data_baseline = list(
  # t0 = -30*years,
  t0 = -5*years,
  gamma_d = 1/434.,
  gamma_l = 1/223,
  f = 1/72,
  r = 1/60,
  eps = 0
)
# Add each scenario's data on
data_scenarios = lapply(seq_len(nrow(scenarios)), function(i) {
  data_scenario = data_baseline
  scenario_specific = scenarios[i,]
  for (name in names(scenario_specific)) {
    if (!is.character(scenario_specific[[name]])) {
      data_scenario[[name]] = scenario_specific[[name]]
    }
  }
  data_scenario$y0 = my_state_init(data_scenario, I0=0.01)
  
  # Add case data
  region_name = scenario_specific$region
  ts_region = china_data %>%
    filter(Region == region_name)
  first_year = min(year(ts_region$Date))
  ts_region = ts_region %>%
    mutate(ts = as.numeric(Date - as.Date(paste0(first_year, "-01-01"))))
  data_scenario$ts = ts_region$ts
  data_scenario$cases = ts_region$Cases
  return(data_scenario)
}) %>%
  setNames(scenarios$name)


# CHANGE THIS
# data_scenarios = head(data_scenarios, length(data_scenarios/2))

# Define inits - start at the mean of all scenarios and the mean sd of all chains
init = c(
  alpha = 0.173,
  beta = 0.899,
  lambda = 0.0304,
  phi = 1.42,
  kappa = 2.45,
  phase = 214,
  p_long = 0.775,
  p_silent = 0.245,
  p_RCI = 0.635
)

init_sd = c(alpha = 0.0293,
            beta = 0.0281,
            lambda = 0.0162,
            phi = 0.277,
            kappa = 0.510,
            phase = 6.21,
            p_long = 0.0841,
            p_silent = 0.0599,
            p_RCI = 0.0411)
init_covar = diag(length(init_sd)) * init_sd^2

init_covar = structure(
  c(0.001104, -3.805e-05, -8.047e-05, 0.01212, 0.002304, 
    0.06224, -5.32e-05, -0.0002501, -0.001002, -3.805e-05, 0.0006094, 
    7.558e-06, -0.00134, -0.001085, 0.1463, -1.584e-05, 0.0002763, 
    -0.0002124, -8.047e-05, 7.558e-06, 1.838e-05, -0.001252, -0.000429, 
    0.02418, 6.135e-05, 0.0001365, 6.5e-05, 0.01212, -0.00134, -0.001252, 
    0.8752, 0.4019, -40.36, -0.01619, -0.03659, 0.05608, 0.002304, 
    -0.001085, -0.000429, 0.4019, 0.9305, -50.06, -0.01586, -0.02713, 
    0.03045, 0.06224, 0.1463, 0.02418, -40.36, -50.06, 4815, 1.711, 
    3.533, -5.482, -5.32e-05, -1.584e-05, 6.135e-05, -0.01619, -0.01586, 
    1.711, 0.005007, 0.001365, -0.001304, -0.0002501, 0.0002763, 
    0.0001365, -0.03659, -0.02713, 3.533, 0.001365, 0.005868, -0.00446, 
    -0.001002, -0.0002124, 6.5e-05, 0.05608, 0.03045, -5.482, -0.001304, 
    -0.00446, 0.01553),
  dim = c(9L, 9L),
  dimnames = list(c("alpha", "beta", "lambda", "phi", "kappa", "phase", "p_long", "p_silent", "p_RCI"),
                  c("alpha", "beta", "lambda", "phi", "kappa", "phase", "p_long", "p_silent", "p_RCI")))

# Now, can we define a function that provides the inits based on the scenario?
# the stored model stats file will store a small sample from the previous run,
# to be used as initial values (if it exists). Therefore we can skip the burn-in
# period.
stored_model_stats_file = here("Rmd", "model_stats_2026-02-12.rds")
scenario_init = function(name, init, init_sd) {
  if (file.exists(stored_model_stats_file)) {
    # Get initial values and covariance from the saved file
    stored_model_stats = read_rds(stored_model_stats_file)
    if (name %in% names(stored_model_stats)) {
      if (is.data.frame(stored_model_stats[[name]]$init)) {
        init_df = stored_model_stats[[name]]$init
        init_row = head(init_df, 1)
        init_vec = as.numeric(init_row) %>% setNames(names(init_row))
        stored_model_stats[[name]]$init = init_vec
      }
      return(stored_model_stats[[name]])
    }
  }
  return(list(init = init, init_covar = init_covar, use_defaults=TRUE))
}

samp_results = rep(list(NULL), length(data_scenarios)) %>%
  setNames(names(data_scenarios))

models = lapply(data_scenarios, make_model)

use_threading = TRUE
num_scenarios = ifelse(any(is.na(scenario_indices)), length(data_scenarios), length(scenario_indices))
max_hours_per_chain = max_hours / num_scenarios / ifelse(use_threading, 1, n_cores)
chains_per_scenario = n_cores

do_scenario = function(i, force_initialisation = FALSE) {
  model_name = names(models)[i]
  model_stats = scenario_init(model_name, init, init_sd)
  if ("use_defaults" %in% names(model_stats) | force_initialisation) {
    # We're using the default initial values and covariance so need to perform adaptation
    n_burnin = 400
    n_adapt = 100
  } else if (is.data.frame(model_stats$init)) {
    # We've provided a range of initial values so burnin is not required
    n_burnin = 0
    n_adapt = 0
  } else {
    # We've provided a mean initial value so we need to burnin to get some variation
    n_burnin = 100
    n_adapt = 0
  }
  n_burnin = 0
  n_adapt = 0
  samp = metropolis_sampling(models[[i]],
                             init = model_stats$init,
                             init_covar = model_stats$init_covar,
                             data = data_scenarios[[i]],
                             n_iter = 1000,
                             n_burnin = n_burnin,
                             n_adapt = n_adapt,
                             n_chains = chains_per_scenario,
                             time_limit = max_hours_per_chain,
                             threading = use_threading)
}


count_successful_chains = function(x) {
  if (is.null(x)) {
    return(0)
  } else {
    null_chains = sapply(x$sim, is.null)
    return(sum(!null_chains))
  }
}

calculate_lpp = function(x) {
  if (is.null(x)) {
    return(NA)
  } else {
    lpp = bind_rows(x$sim_diagnostics) %>%
      pull(lpp) %>%
      mean()
    return(lpp)
  }
}

calculate_accept = function(x) {
  if (is.null(x)) {
    return(NA)
  } else {
    accept_ratio = mean(unlist(x$accept))
    return(accept_ratio)
  }
}

calculate_ess = function(x) {
  if (is.null(x)) {
    return(NA)
  } else {
    return(mean(unlist(x$ESS)))
  }
}

calculate_acceptance = function(x) {
  if (is.null(x)) {
    return(NA)
  } else {
    return(mean(unlist(x$accept)))
  }
}

# Run iterations
for (i in seq_len(length(data_scenarios))) {
  if (!any(is.na(scenario_indices)) & !i %in% scenario_indices) {
    # If specified, only run specific scenarios - useful for a distributed workload
    next
  }
  message("- Processing scenario ", i)
  samp_result = do_scenario(i, force_initialisation=TRUE)
  successful_chains = count_successful_chains(samp_result)
  ess = calculate_ess(samp_result)
  
  if (successful_chains >= chains_per_scenario & ess > 1) {
    # Assess goodness of convergence
    lpp = calculate_lpp(samp_result)
    
    
    # Save results
    samp_results[[i]] = samp_result
    rds_filename = here::here("Rmd", paste0("workspaces/samp_results[[", i, "]].rds"))
    message("Successful chains: ", successful_chains, ". ESS: ", ess)
    message("Saving results to ", rds_filename)
    write_rds(samp_results[[i]], rds_filename)
  }
  else {
    message("! Encountered an issue with sampler results so did not save. Successful chains: ", successful_chains, ". ESS: ", ess)
  }
}

# basic trace plot
scenario_to_plot = 1
samp_results[[scenario_to_plot]]$sim %>%
  bind_rows(.id = "chain") %>%
  mutate(i = row_number(),
         chain = as.integer(chain)) %>%
  sample_n(min(n(), 1000)) %>%
  filter(i > nrow({.})/10) %>%
  pivot_longer(-c(i, chain)) %>%
  ggplot(aes(x=i, y=value, color=as.character(chain))) +
  geom_step() +
  facet_wrap(vars(name), scales="free_y") +
  theme(legend.position = "none")

# basic density plot
samp_results[[scenario_to_plot]]$sim %>%
  bind_rows(.id = "chain") %>%
  sample_n(min(n(), 1000)) %>%
  mutate(i = row_number(),
         chain = as.integer(chain)) %>%
  filter(i > nrow({.})/10) %>%
  pivot_longer(-c(i, chain)) %>%
  ggplot(aes(x=value, color=as.character(chain))) +
  geom_density() +
  facet_wrap(vars(name), scales="free") +
  theme(legend.position = "none")

successful_chains = sapply(samp_results, function(x) {
  if (is.null(x)) {
    return(0)
  } else {
    null_chains = sapply(x$sim, is.null)
    return(sum(!null_chains))
  }
})

# Assess goodness of convergence
lpps = sapply(samp_results, calculate_lpp)
ess = sapply(samp_results, calculate_ess)
acceptance_ratios = sapply(samp_results, calculate_acceptance)

# Now use samp_results to calculate the optimal inits
get_inits = function(results) {
  random_sample = results$sim %>%
    bind_rows() %>%
    sample_n(min(100, nrow(.))) # Just get 100 to make the data size smaller.
  new_covar = 2.38^2 * cov(bind_rows(results$sim)) / length(results$current_x[[1]])
  return(list(init=random_sample, init_covar=new_covar))
}
new_inits = lapply(samp_results, get_inits)

# If a scenario's log posterior probability seems extremely bad, it indicates a failure to converge and we don't want to write those inits. Instead, we will generate inits from the other scenarios.

# Similarly, if a scenario's covariance is zero on the diagonals, there is also a failure to converge.

# Calculate which scenarios are extremely bad
lpps_mean = mean(lpps, na.rm=T)
lpps_sd = sd(lpps, na.rm=T)
lpps_pnorm = pnorm(lpps, mean=lpps_mean, sd=lpps_sd)
okay_scenarios = names(samp_results)[which(lpps_pnorm <= 0.95 & ess > 1)]
extremely_bad_scenarios = names(samp_results)[!names(samp_results) %in% okay_scenarios]

# Create some inits to use
okay_init = lapply(new_inits[okay_scenarios], function(x) {
  x$init
}) %>%
  bind_rows() %>%
  sample_n(min(100, nrow(.)))

# Calculate the root-mean-square variances across the other scenarios
okay_init_covar = (Reduce(
  '+',
  lapply(new_inits[okay_scenarios], function(x) {
    sqrt(x$init_covar * diag(nrow=ncol(x$init)))
  })
) / ncol(okay_init)) ^ 2

# Replace the bad inits with something generic so hopefully the next run won't get stuck in a bad space
if (length(extremely_bad_scenarios) > 0) {
  warning("There were ", length(extremely_bad_scenarios), " scenarios where the log posterior probability was unusually poor compared to the others. Potential failure to fit to data.")
}
for (scenario in extremely_bad_scenarios) {
  new_inits[[scenario]] = list(
    init = okay_init,
    init_covar = okay_init_covar
  )
}

# write_rds(new_inits, stored_model_stats_file)
results_files = list.files(here("Rmd/workspaces"), "samp_results.*\\.rds$")
for (x in results_files) {
  i = x %>% str_extract("[0-9]+") %>% as.integer()
  samp_results[[i]] = read_rds(here("Rmd/workspaces", paste0("samp_results[[", i, "]].rds")))
  names(samp_results)[i] = names(data_scenarios)[i]
}

# Resimulate trajectories for use in plotting (so we don't need to compile the model later)
resim = pbmclapply(samp_results, function(samp) {
  if (is.null(samp)) {
    return(NULL)
  }
  resim = extract(samp, "incidence", n_samples=100, threading=F)
}) %>%
  bind_rows(.id = "Scenario") %>%
  mutate(Scenario = fct_inorder(Scenario)) %>%
  left_join(scenarios, by=c("Scenario" = "name")) %>%
  mutate(name_short = name_short %>% str_replace_all(", ", ",\n")) %>%
  pivot_longer(matches("^clinical"), names_to="metric") %>%
  mutate(metric = metric %>% case_match(
    "clinical_incidence" ~ "Total",
    "clinical_relapse" ~ "Relapse",
    "clinical_primary" ~ "Primary"
  ))

# Do some diagnostic plots
posterior_seasonal_1 = pblapply(samp_results, function(samp) {
  if (is.null(samp)) {
    return(NULL)
  }
  bind_cols(
    bind_rows(samp$sim),
    bind_rows(samp$sim_diagnostics, .id = "chain")
  ) %>%
    slice_sample(n = 1000) %>%
    mutate(phase = phase + years/4) %>%
    pivot_longer(-c(iteration, accept, lpp, ll, chain), names_to = "parameter", values_to = "value")
}) %>%
  bind_rows(.id = "Scenario") %>%
  mutate(Scenario = fct_inorder(Scenario))#
posterior_seasonal_2 = posterior_seasonal_1 %>%
  group_by(Scenario, parameter)
x1 = scenarios %>%
  select(name, name_short, region)
posterior_seasonal = posterior_seasonal_2 %>%
  left_join(scenarios, by=c("Scenario" = "name")) %>%
  mutate(name_short = name_short %>% str_replace_all(", ", ",\n"))

plot_original_data = lapply(data_scenarios, function(x) {
  tibble(time = x$ts,
         cases = x$cases)
}) %>%
  bind_rows(.id = "Scenario") %>%
  left_join(scenarios, by=c("Scenario" = "name")) %>%
  mutate(Scenario = fct_inorder(Scenario))

resim_seasonality = pblapply(samp_results[1:2], function(samp) {
  if (is.null(samp)) {
    return(NULL)
  }
  t = seq(0, years, length.out=500)
  samp_sim = bind_rows(samp$sim)
  samp_rand = sample.int(nrow(samp_sim), min(nrow(samp_sim), 500))
  suitability_traces = lapply(samp_rand, function(ix) {
    samp_suitability = tibble(
      time = t,
      suitability = sapply(t, function(tt) {
        omega = suitability(tt, samp$data$eps, samp_sim[[ix, "kappa"]], samp_sim[[ix, "phase"]])
      })
    )
  }) %>%
    bind_rows(.id = "trace")
}) %>%
  bind_rows(.id = "Scenario") %>%
  left_join(scenarios, by=c("Scenario" = "name")) %>%
  mutate(Scenario = fct_inorder(Scenario))

# If baseline is not present, get the first name_shortest and plot these. Otherwise the plots break.
# if ("Baseline" %in% posterior_seasonal$name_shortest) {
#   first_name_shortest = "Baseline"
# } else {
#   first_name_shortest = first(posterior_seasonal$name_shortest)
# }
# 
# trace_plot = posterior_seasonal %>%
#   filter(name_shortest == first_name_shortest) %>%
#   ggplot(aes(x = iteration, y = value, color=Scenario, group=interaction(Scenario, chain))) +
#   geom_step(alpha=0.25) +
#   coord_cartesian(xlim = c(0, NA)) +
#   scale_y_log10() +
#   facet_wrap(vars(parameter), scales="free", labeller=plot_labeller_novar) +
#   labs(subtitle = "Parameter traces")
# trace_plot
# filename = paste0("../plots/china_trace.png")
# ggsave(filename, width=8, height=8)

end_time = Sys.time()
print(end_time - start_time)
print(paste("Time elapsed:", end_time - start_time))

# workspace_filename = paste0("workspaces/Chapter_02_china_metropolis_", Sys.Date(), ".RData")
workspace_filename = here("Rmd", "workspaces/Chapter_02_china_metropolis_latest.RData")
save.image(workspace_filename)

