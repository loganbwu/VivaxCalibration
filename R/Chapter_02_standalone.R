args <- commandArgs(trailingOnly = TRUE)
# Check that there is exactly one argument. If not, provide a script usage statement.
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)>1) {
  stop("More than one argument specified.n", call.=FALSE)
}
max_hours = as.numeric(args[1])

library(R.utils)
library(tidyverse)
library(rstan)
library(parallel)
library(patchwork)
library(pbmcapply)
library(pbapply)
library(memoise)
library(RColorBrewer)
library(ggtext)
setwd("Rmd")
source("../R/constants.R")
source("../R/load_functions.R")

rstan_options(auto_write = TRUE, threads_per_chain = 1)

n_cores = parallelly::availableCores()
options(mc.cores = n_cores)
message("Running on ", n_cores, " cores")

# Load model and metropolis algorithm (copied from VivaxODE project folder)
source(file = "../R/models/temperate_v11.R")
source(file = "../R/functions/adaptive_metropolis.R")
# Load priors for each scenario
source(file = "../R/priors.R")
my_state_init = state_init_3
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

china_data = read_rds("china_data.rds")

data_baseline = list(
  t0 = -30*years,
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
  phase = 119,
  p_long = 0.775,
  p_silent = 0.245,
  p_RCI = 0.123
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

# Now, can we define a function that provides the inits based on the scenario?
stored_model_stats_file = "model_stats.rds"
scenario_init = function(name, init, init_sd) {
  if (file.exists(stored_model_stats_file)) {
    stored_model_stats = read_rds(stored_model_stats_file)
    if (name %in% names(stored_model_stats)) {
      return(stored_model_stats[[name]])
    }
  }
  return(list(init = init, init_covar = init_covar, use_defaults=TRUE))
}

samp_results = rep(list(NULL), length(data_scenarios)) %>%
  setNames(names(data_scenarios))

models = lapply(data_scenarios, make_model)

chains_per_scenario = max(1, floor(n_cores / length(data_scenarios)))
# max_hours_per_scenario = max_hours / length(data_scenarios)
max_hours_per_scenario = max_hours # run all scenarios concurrently
do_scenario = function(i) {
  model_name = names(models)[i]
  model_stats = scenario_init(model_name, init, init_sd)
  if ("use_defaults" %in% names(model_stats)) {
    n_burnin = 400
    n_adapt = 100
  } else {
    # Assume we've already burned in and adapted by using the initial values
    n_burnin = 100
    n_adapt = 0
  }
  samp = metropolis_sampling(models[[i]],
                             init = model_stats$init,
                             init_covar = model_stats$init_covar,
                             data = data_scenarios[[i]],
                             n_iter = 1000,
                             n_burnin = n_burnin,
                             n_adapt = n_adapt,
                             n_chains = chains_per_scenario,
                             time_limit = max_hours_per_scenario)
}
samp_results_lapply = pbmclapply(seq_len(length(data_scenarios)), do_scenario)

for (i in seq_len(length(data_scenarios))) {
  samp_results[[i]] = samp_results_lapply[[i]]
}
rm(samp_results_lapply)

# Now use samp_results to calculate the optimal inits
get_inits = function(results) {
  new_init = bind_rows(results$current_x) %>%
    colMeans()
  new_covar = 2.38^2 * cov(bind_rows(results$sim)) / length(results$current_x[[1]])
  return(list(init=init, init_covar=new_covar))
}
new_inits = lapply(samp_results, get_inits)
write_rds(new_inits, stored_model_stats_file)


# Resimulate trajectories
resim = pbmclapply(samp_results, function(samp) {
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

end_time = Sys.time()
print(end_time)
print(end_time - start_time)

workspace_filename = paste0("Chapter_02_china_metropolis_", Sys.Date(), ".RData")
save.image(workspace_filename)

# Do some diagnostic plots
posterior_seasonal_1 = pblapply(samp_results, function(samp) {
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

# If baseline is not present, get the first name_shortest and plot these. Otherwise the plots break.
if ("Baseline" %in% posterior_seasonal$name_shortest) {
  first_name_shortest = "Baseline"
} else {
  first_name_shortest = first(posterior_seasonal$name_shortest)
}

trace_plot = posterior_seasonal %>%
  filter(name_shortest == first_name_shortest) %>%
  ggplot(aes(x = iteration, y = value, color=Scenario, group=interaction(Scenario, chain))) +
  geom_step(alpha=0.25) +
  coord_cartesian(xlim = c(0, NA)) +
  scale_y_log10() +
  facet_wrap(vars(parameter), scales="free", labeller=plot_labeller_novar) +
  labs(subtitle = "Parameter traces")
trace_plot
filename = paste0("../plots/china_trace.png")
ggsave(filename, width=8, height=8)

