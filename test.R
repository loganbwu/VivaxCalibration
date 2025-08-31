install.packages("Rcpp", type="source", repos = "https://cran.rstudio.com/")
install.packages("rstan", type="source", repos = "https://cran.rstudio.com/")
library(rstan)

source(file = "R/models/temperate_v11.R")