extract_summary = function(data, sim_out) {
  with(data, cbind(as.data.frame(
    summary(
      sim_out,
      pars = "y",
      probs = c(0.05, 0.5, 0.95)
    )$summary))) %>%
    rownames_to_column(var = "name") %>%
    mutate(index_i = name %>% str_extract("(?<=\\[)[0-9]+") %>% as.numeric(),
           time = data$ts[index_i],
           index_j = name %>% str_extract("(?<=,)[0-9]+") %>% as.numeric(),
           compartment = names(data$y0)[index_j],
           variable = name %>% str_remove("\\[.*")) %>%
    as_tibble() %>%
    select(time, variable=compartment, mean) %>%
    pivot_wider(names_from=variable, values_from=mean) %>%
    mutate(time = time,
           Infectious = rowSums(across(starts_with("I"))),
           Dormant = rowSums(across(matches(paste0("^Sl(", paste(seq_len(data$n_dormant), collapse="|"), ")$")))),
           Latent := !!rlang::sym(paste0("Sl", data$n_dormant+1)),
           Susceptible = S0,
           Prevalence = 1 - S0,
           # Below incidence rates are in units of people per day
           across(matches("True"), ~ (lead(.x) - .x) / as.numeric(lead(time) - time)),
           TrueCases = TrueShortIncubations + TrueLongIncubations + TrueRelapses,
           ReportedShortIncubations = TrueShortIncubations * data$alpha,
           ReportedLongIncubations = TrueLongIncubations * data$alpha,
           ReportedRelapses = TrueRelapses * (1-data$relapse_clinical_immunity) * data$alpha,
           ReportedCases = ReportedShortIncubations + ReportedLongIncubations + ReportedRelapses,
           # checksum = rowSums(across(-one_of("Cases"))),
           .keep = "used")
}


#' Adapted for updated clinical incidence model
#' Also split into primary cases and relapses
extract_summary_2 = function(data, sim_out) {
  with(data, cbind(as.data.frame(
    summary(
      sim_out,
      pars = "y",
      probs = c(0.05, 0.5, 0.95)
    )$summary))) %>%
    rownames_to_column(var = "name") %>%
    mutate(index_i = name %>% str_extract("(?<=\\[)[0-9]+") %>% as.numeric(),
           time = data$ts[index_i],
           index_j = name %>% str_extract("(?<=,)[0-9]+") %>% as.numeric(),
           compartment = names(data$y0)[index_j],
           variable = name %>% str_remove("\\[.*")) %>%
    as_tibble() %>%
    select(time, variable=compartment, mean) %>%
    pivot_wider(names_from=variable, values_from=mean) %>%
    mutate(time = time,
           Infectious = rowSums(across(starts_with("I"))),
           Dormant = rowSums(across(matches(paste0("^Sl(", paste(seq_len(data$n_dormant), collapse="|"), ")$")))),
           Latent := !!rlang::sym(paste0("Sl", data$n_dormant+1)),
           Susceptible = S0,
           Prevalence = 1 - S0,
           ClinicalPrimary = ClinicalPrimary - lag(ClinicalPrimary),
           ClinicalRelapse = ClinicalRelapse - lag(ClinicalRelapse),
           ClinicalIncidence = ClinicalPrimary + ClinicalRelapse,
           # checksum = rowSums(across(-one_of("Cases"))),
           .keep = "used")
}

