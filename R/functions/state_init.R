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

#' Adapted for new clinical incidence model for use with temperate_4
state_init_v4 = function(parameters, ...) {
  listargs = list(...)
  if (length(listargs) > 0 & is.list(listargs[[1]])) {
    listargs = listargs[[1]]
  }
  
  n_stages = parameters$n_dormant + 1
  names_compartments = c("S0", "I0",
                         paste0("Sl", seq_len(n_stages)),
                         paste0("Scl", seq_len(n_stages)),
                         paste0("Icl", seq_len(n_stages)))
  
  state = rep(0, length(names_compartments)) %>%
    setNames(names_compartments) %>%
    c(ClinicalIncidence = 0)
  for (n in names(listargs)[names(listargs) %in% names(state)]) {
    state[n] = listargs[[n]]
  }
  state["S0"] = 1 - sum(state) + state["S0"]
  state
}

#' Adapted for new clinical incidence model. For use with version 7
state_init_2 = function(parameters, ...) {
  listargs = list(...)
  if (length(listargs) > 0 & is.list(listargs[[1]])) {
    listargs = listargs[[1]]
  }
  
  n_stages = parameters$n_dormant + 1
  names_compartments = c("S0", "I0",
                         paste0("Sl", seq_len(n_stages)),
                         paste0("Scl", seq_len(n_stages)),
                         paste0("Icl", seq_len(n_stages)))
  
  state = rep(0, length(names_compartments)) %>%
    setNames(names_compartments) %>%
    c(ClinicalPrimary = 0, ClinicalRelapse = 0)
  for (n in names(listargs)[names(listargs) %in% names(state)]) {
    state[n] = listargs[[n]]
  }
  state["S0"] = 1 - sum(state) + state["S0"]
  state
}

#' Adapted for new clinical incidence model including all primaries or relapses. For use with version 8
state_init_3 = function(parameters, ...) {
  listargs = list(...)
  if (length(listargs) > 0 & is.list(listargs[[1]])) {
    listargs = listargs[[1]]
  }
  
  n_stages = parameters$n_dormant + 1
  names_compartments = c("S0", "I0",
                         paste0("Sl", seq_len(n_stages)),
                         paste0("Scl", seq_len(n_stages)),
                         paste0("Icl", seq_len(n_stages)))
  
  state = rep(0, length(names_compartments)) %>%
    setNames(names_compartments) %>%
    c(AllPrimary = 0, AllRelapse = 0, ClinicalPrimary = 0, ClinicalRelapse = 0)
  for (n in names(listargs)[names(listargs) %in% names(state)]) {
    state[n] = listargs[[n]]
  }
  state["S0"] = 1 - sum(state) + state["S0"]
  state
}

