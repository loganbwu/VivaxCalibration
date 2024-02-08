state_init = function(parameters, ...) {
  listargs = list(...)
  if (length(listargs) > 0 & is.list(listargs[[1]])) {
    listargs = listargs[[1]]
  }
  
  n_stages = parameters$n_dormant + 1
  names_compartments = c("S0", "I0",
                         paste0("Sl", seq_len(n_stages)),
                         paste0("Scl", seq_len(n_stages)),
                         paste0("Icl", seq_len(n_stages)))
  
  state = rep(0, length(names_compartments)) %>% setNames(names_compartments) %>% c(TrueShortIncubations=0, TrueLongIncubations=0, TrueRelapses=0)
  for (n in names(listargs)[names(listargs) %in% names(state)]) {
    state[n] = listargs[[n]]
  }
  state["S0"] = 1 - sum(state) + state["S0"]
  state
}
