
c_posterior = "blue"

years = 365.25
annual_subdivisions = 12
n_traces = 100 # maximum traces to plot

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

names_params = c("lambda", "phi", "eps", "kappa", "phase", "tstar", "xi")
param_colors = brewer.pal(length(names_params), "Set2") %>% setNames(names_params)
