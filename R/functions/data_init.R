data_init = function(data, ...) {
  listargs = list(...)
  if (length(listargs) > 0 & is.list(listargs[[1]])) {
    listargs = listargs[[1]]
  }
  
  for (n in names(listargs)) {
    data[[n]] = listargs[[n]]
  }
  data
}
