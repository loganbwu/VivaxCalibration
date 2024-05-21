.plot_labeller <- function(variable, value){
  var = make_greek(variable)
  if (is.numeric(value)) {
    value = round(value, 6)
  }
  return(paste(var, "=", value))
}

.plot_labeller_years <- function(variable, value) {
  if (variable == "tstar") {
    x = plot_labeller(variable, paste(round(value/365.25, 2), "years"))
  } else {
    x = plot_labeller(variable, value)
  }
  return(x)
}

.plot_labeller_novar <- function(variable, value){
  return(make_greek(value))
}

plot_labeller = function(labels, multi_line=TRUE) {
  labels = lapply(names(labels), function(name) {
    .plot_labeller(name, labels[[name]])
  })
  if (multi_line) {
    labels
  }
  else {
    collapse_labels_lines(labels)
  }
}

plot_labeller_years = function(labels, multi_line=TRUE) {
  labels = lapply(names(labels), function(name) {
    .plot_labeller_years(name, labels[[name]])
  })
  if (multi_line) {
    labels
  }
  else {
    collapse_labels_lines(labels)
  }
}

plot_labeller_novar = function(labels, multi_line=TRUE) {
  labels = lapply(names(labels), function(name) {
    .plot_labeller_novar(name, labels[[name]])
  })
  if (multi_line) {
    labels
  }
  else {
    collapse_labels_lines(labels)
  }
}
