#' Calaculate r-hat from a vector of samples
#'
#' Use this to filter out chains that haven't converged
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
