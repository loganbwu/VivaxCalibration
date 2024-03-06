extract_incidence = function(fit, data) {
  incidence = rstan::extract(fit, "incidence")[[1]]
  if (n_traces > dim(incidence)[1]) {
    sample_ix = sample(seq_len(dim(incidence)[1]), n_traces, replace=T)
  } else {
    sample_ix = seq_len(dim(incidence)[1])
  }
  incidence_sample = as_tibble(t(incidence[sample_ix,])) %>%
    mutate(j = row_number()) %>%
    pivot_longer(-j, names_to = "trace", values_to = "incidence") %>%
    drop_na(j) %>%
    group_by(j) %>%
    mutate(ts = data$ts[j],
           lower = quantile(incidence, 0.025, na.rm=T),
           upper = quantile(incidence, 0.975, na.rm=T),
           legend = "95% prediction interval")
}