R_functions = list.files("../R/functions", "\\.R$", recursive=TRUE, full.names=TRUE)
for (file in R_functions) { source(file) }
