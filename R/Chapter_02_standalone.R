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

if (!require(MalariaData)) {
  china_data = read_rds("china_data.rds")
} else {
  library(MalariaData)
  regions = c("Xiayi", "Guantang", "Dengzhou", "Huangchuan") %>% setNames({.})
  china_data_all = lapply(regions, function(x) {load_region(x, species = "Vivax", source="Local")}) %>%
    bind_rows(.id = "Region") %>%
    select(Region, Date, Cases) %>%
    left_join(china_selections, by="Region") %>%
    mutate(Date = as.Date(paste(year(Date), month(Date), "01", sep="-")) + months(1) - days(1)) %>% # Align to last day of the month
    mutate(include = ifelse(Date >= min & Date < max, "grey20", "grey"),
           include = replace_na(include, "grey"))
  
  china_data = china_data_all %>%
    filter(include == "grey20") %>%
    select(-include)
  
  write_rds(china_data, "china_data.rds")
}

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
data_scenarios = head(data_scenarios, length(data_scenarios/2))

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

samp_results = rep(list(NULL), length(data_scenarios)) %>%
  setNames(names(data_scenarios))

max_hours = 24
models = lapply(data_scenarios, make_model)

chains_per_scenario = max(1, floor(n_cores / length(data_scenarios)))
max_hours_per_scenario = max_hours# / length(data_scenarios)
do_scenario = function(i) {
  samp = metropolis_sampling(models[[i]],
                             init = init,
                             init_sd = init_sd,
                             data = data_scenarios[[i]],
                             n_iter = 1000,
                             n_burnin = 400,
                             n_adapt = 100,
                             n_chains = chains_per_scenario,
                             time_limit = max_hours_per_scenario)
}
samp_results_lapply = pbmclapply(seq_len(length(data_scenarios)), do_scenario)

for (i in seq_len(length(data_scenarios))) {
  samp_results[[i]] = samp_results_lapply[[i]]
}

end_time = Sys.time()
print(end_time)
print(end_time - start_time)
# workspace_filename = paste0("workspaces/Chapter_02_china_metropolis_", Sys.Date(), ".RData")
rds_filename = paste0("samp_results/Chapter_02_china_metropolis_", Sys.Date(), ".rds")
# save.image(workspace_filename)
write_rds(samp_results, file=rds_filename,compress="xz")

