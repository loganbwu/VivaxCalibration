plot_labeller <- function(variable, value){
  var = case_match(variable, 
                   c("transmission_rates", "lambda") ~ "λ",
                   "phi" ~ "φ",
                   c("seasonality_ratio", "eps") ~ "ε",
                   "kappa" ~ "κ",
                   "phase" ~ "ψ",
                   "tstar" ~ "t*",
                   "xi" ~ "ξ",
                   .default = as.character(variable))
  return(paste(var, "=", value))
}

plot_labeller_novar <- function(variable, value){
  val = case_match(value, 
                   c("transmission_rates", "lambda") ~ "λ",
                   "phi" ~ "φ",
                   c("seasonality_ratio", "eps") ~ "ε",
                   "kappa" ~ "κ",
                   "phase" ~ "ψ",
                   "tstar" ~ "t*",
                   "xi" ~ "ξ",
                   .default = as.character(value))
  return(val)
}
