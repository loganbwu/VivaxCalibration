
#' Aggregate time in twelve month blocks
aggregate_time = function(ts) {
  every_twelfth = seq(annual_subdivisions, length(ts), annual_subdivisions) # c(1, 13, 25, ...)
  return(ts[every_twelfth])
}

aggregate_cases = function(cases) {
  every_twelfth = seq(annual_subdivisions, length(cases), annual_subdivisions)
  cases = sapply(every_twelfth, function(mo) {
    ixs = seq_along(cases)
    sum(cases[ixs >= (mo - annual_subdivisions) & ixs < mo])
  })
  return(cases)
}

aggregate_eps = function(eps=NULL) {
  1
}

#' Aggregate in twelve month blocks
aggregate_data = function(.data) {
  .data$n_times = floor(.data$n_times / annual_subdivisions)
  every_twelfth = seq(annual_subdivisions, length(.data$ts), annual_subdivisions) # c(1, 13, 25, ...)
  .data$ts = .data$ts[every_twelfth]
  .data$eps = 1
  .data$cases = sapply(every_twelfth, function(mo) {
    ixs = seq_along(.data$cases)
    sum(.data$cases[ixs >= (mo - annual_subdivisions) & ixs < mo])
  })
  .data
}

my_simulate_data_list = function(...) {
  synth_data = simulate_data(...)
  
  synth_data_rds = readRDS(synth_data$datasets[1])
  file.remove(synth_data$datasets[1])
  return(synth_data_rds)
}

my_simulate_data = function(...) {
  synth_data = simulate_data(...)
  
  synth_data_rds = readRDS(synth_data$datasets[1])
  file.remove(synth_data$datasets[1])
  indx <- sapply(synth_data_rds, length)
  synth_df = lapply(synth_data_rds, function(x) {length(x) = max(indx); x}) %>%
    as.data.frame() %>%
    as_tibble() %>%
    drop_na()
  return(synth_df)
}

plot_labeller <- function(variable, value){
  var = case_match(variable, 
                   "seasonality_ratio" ~ "ε",
                   c("transmission_rates", "lambda") ~ "λ",
                   "tstar" ~ "t*",
                   "xi" ~ "ξ",
                   .default = as.character(variable))
  return(paste(var, "=", value))
}
