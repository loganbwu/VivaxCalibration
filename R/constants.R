
c_posterior = "blue"

years = 365.25
annual_subdivisions = 12
n_traces = 1000 # maximum traces to plot

comparison_colors = c("lambda_nonseasonal_poisson" = "tomato",
                      "lambda_nonseasonal" = "orange",
                      "lambda_seasonal_poisson" = "navy",
                      "lambda_seasonal" = "steelblue",
                      "lambda_seasonal_ext" = "hotpink")

# Use this to filter out chains that haven't converged
.calculate_rhat = function(.trace, .n_chains = n_chains) {
  if (is.null(.trace)) {
    return(NULL)
  }
  chains = matrix(.trace, ncol = .n_chains)
  Rhat(chains)
}
calculate_rhat = function(x) {
  sapply(x, .calculate_rhat)
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