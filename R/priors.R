# Define priors
duration = tribble(
  ~name, ~delta, ~n_dormant,
  "baseline", 1/162, 11,
  "extended", 1/240, 11
)

ascertainment = tribble(
  ~name, ~alpha_shape1, ~alpha_shape2,
  "baseline", 10, 90,
  "high", 90, 10
)

phenotype = tribble(
  ~name, ~p_silent_shape1, ~p_silent_shape2, ~p_long_shape1, ~p_long_shape2, ~p_RCI_shape1, ~p_RCI_shape2,
  "baseline", 15, 30, 22, 5, 1, 1,
  "no_immunity", 15, 30, 22, 5, 1, 50
)

region = tribble(
  ~name, ~N,
  "Dengzhou", 300000,
  "Guantang", 40000
)

# Combine them all and join the values
scenarios = expand_grid(
  duration = duration$name,
  ascertainment = ascertainment$name,
  phenotype = phenotype$name,
  region = region$name
) %>%
  left_join(duration, by=c("duration"="name")) %>%
  left_join(ascertainment, by=c("ascertainment"="name")) %>%
  left_join(phenotype, by=c("phenotype"="name")) %>%
  left_join(region, by=c("region"="name")) %>%
  mutate(n_changes = ((duration!="baseline") + (ascertainment!="baseline") + (phenotype!="baseline")),
         scenario = row_number(), .before = 0) %>%
  filter(n_changes <= 1) %>%
  mutate(name = paste(duration, "duration,", ascertainment, "ascertainment,", phenotype, "phenotype,", region))
