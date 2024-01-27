
c_posterior = "blue"

years = 365.25
annual_subdivisions = 12
n_traces = 1000 # maximum traces to plot

comparison_colors = c("lambda_nonseasonal_poisson" = "tomato",
                      "lambda_nonseasonal_negbin" = "orange",
                      "lambda_seasonal_poisson" = "navy",
                      "lambda_seasonal_negbin" = "steelblue",
                      "lambda_seasonal_negbin_ext" = "hotpink")

renamed_comparison_colors = c("Nonseasonal, poisson" = "tomato",
                              "Nonseasonal, negbin" = "orange",
                              "Seasonal, poisson" = "navy",
                              "Seasonal, negbin" = "steelblue",
                              "Seasonal, negbin_ext" = "hotpink")

rename_methods = names(comparison_colors)
names(rename_methods) = names(renamed_comparison_colors)

label_auto = function(x) {
  max_x = max(x, na.rm=T)
  if (max_x > 5e6) {
    paste(x/1e6, "M")
  } else if (max_x > 5e3) {
    paste(x/1e3, "k")
  } else {
    x
  }
}

label_auto2 = function(x) {
  case_when(x >= 1e6 ~ paste(x/1e6, "M"),
            x >= 1e3 ~ paste(x/1e3, "k"),
            TRUE ~ as.character(x))
}


# Use this to filter out chains that haven't converged
.calculate_rhat = function(.trace, .n_chains = n_chains) {
  if (is.null(.trace)) {
    return(NA)
  }
  chains = matrix(.trace, ncol = .n_chains)
  Rhat(chains)
}
calculate_rhat = function(x) {
  sapply(x, .calculate_rhat)
}

#' Aggregate time in twelve month blocks
aggregate_time = function(ts) {
  every_twelfth = seq(1, length(ts), annual_subdivisions) # c(1, 13, 25, ...)
  return(ts[every_twelfth])
}

aggregate_cases = function(cases) {
  every_twelfth = seq(1, length(cases), annual_subdivisions)
  cases = sapply(every_twelfth, function(mo) {
    ixs = seq_along(.data$cases)
    sum(.data$cases[ixs >= mo & ixs < (mo+annual_subdivisions)])
  })
  return(cases)
}

aggregate_eps = function(eps=NULL) {
  1
}

#' Aggregate in twelve month blocks
aggregate_data = function(.data) {
  .data$n_times = floor(.data$n_times / annual_subdivisions)
  every_twelfth = seq(1, length(.data$ts), annual_subdivisions) # c(1, 13, 25, ...)
  .data$ts = .data$ts[every_twelfth]
  .data$eps = 1
  .data$cases = sapply(every_twelfth, function(mo) {
    ixs = seq_along(.data$cases)
    sum(.data$cases[ixs >= mo & ixs < (mo+annual_subdivisions)])
  })
  .data
}

#' https://stackoverflow.com/a/53018594
with_timeout <- function(expr, cpu, elapsed){
  expr <- substitute(expr)
  envir <- parent.frame()
  setTimeLimit(cpu = cpu, elapsed = elapsed, transient = TRUE)
  on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
  eval(expr, envir = envir)
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
  if (variable == "seasonality_ratio") {
    var = "eps"
  } else if (variable == "transmission_rates") {
    var = "lambda"
  } else {
    var = "facet"
  }
  return(paste0(var, ": ", value))
}
