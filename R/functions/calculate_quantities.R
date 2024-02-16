calculate_quantities = function(x, data) {
  n_stages = names(x) %>%
    str_subset("Sl") %>%
    str_remove("Sl") %>%
    as.numeric() %>%
    max()
  x %>%
    mutate(time = time,
           Infectious = rowSums(across(starts_with("I"))),
           Dormant = rowSums(across(matches(paste0("^Sc?l(", paste(seq_len(data$n_dormant), collapse="|"), ")$")))),
           Latent = rowSums(across(matches(paste0("^Sc?l", data$n_dormant+1, "$")))),
           Susceptible = S0,
           Prevalence = 1 - S0,
           Total = rowSums(across(matches(paste0("^(S|I).{,5}$")))),
           # Below incidence rates are in units of people per day
           across(matches("True"), ~ (lead(.x) - .x) / as.numeric(lead(time) - time)),
           TrueCases = TrueShortIncubations + TrueLongIncubations + TrueRelapses,
           ReportedShortIncubations = TrueShortIncubations * data$alpha,
           ReportedLongIncubations = TrueLongIncubations * data$alpha,
           ReportedRelapses = TrueRelapses * (1-data$relapse_clinical_immunity) * data$alpha,
           ReportedCases = ReportedShortIncubations + ReportedLongIncubations + ReportedRelapses,
           .keep = "used")
}

calculate_quantities_2 = function(x, data) {
  n_stages = names(x) %>%
    str_subset("Sl") %>%
    str_remove("Sl") %>%
    as.numeric() %>%
    max()
  x %>%
    mutate(time = time,
           Infectious = rowSums(across(starts_with("I"))),
           Dormant = rowSums(across(matches(paste0("^Sc?l(", paste(seq_len(data$n_dormant), collapse="|"), ")$")))),
           Latent = rowSums(across(matches(paste0("^Sc?l", data$n_dormant+1, "$")))),
           Susceptible = S0,
           Prevalence = 1 - S0,
           Total = rowSums(across(matches(paste0("^(S|I).{,5}$")))),
           ClinicalPrimary = ClinicalPrimary - lag(ClinicalPrimary),
           ClinicalRelapse = ClinicalRelapse - lag(ClinicalRelapse),
           ClinicalIncidence = ClinicalPrimary + ClinicalRelapse,
           .keep = "used")
}