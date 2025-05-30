# Define priors
duration = tribble(
  ~name, ~delta, ~n_dormant,
  "baseline", 1/162, 11
)

ascertainment = tribble(
  ~name, ~alpha_shape1, ~alpha_shape2,
  "baseline", 10, 90
)

phenotype = tribble(
  ~name, ~p_silent_shape1, ~p_silent_shape2, ~p_long_shape1, ~p_long_shape2, ~p_RCI_shape1, ~p_RCI_shape2,
  "baseline", 15, 30, 22, 5, 1, 1,
  "no_immunity", 15, 30, 22, 5, 1, 999,
  "low_immunity", 15, 30, 22, 5, 333, 667,
  "high_immunity", 15, 30, 22, 5, 667, 333
)

# helper to visualise beta function
# xx = seq(0, 1, length.out=1000)
# yy = dbeta(xx, 667, 333)
# plot(xx, yy, 'l')

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
  # filter(n_changes <= 1) %>%
  mutate(name = paste(duration, "duration,", ascertainment, "ascertainment,", phenotype, "phenotype,", region),
         name_short = paste(duration, "duration,", ascertainment, "ascertainment,", phenotype, "phenotype"),
         name_shortest = case_when(duration == "extended" ~ "Extended dormancy",
                                   ascertainment == "high" ~ "High ascertainment",
                                   phenotype == "no_immunity" ~ "No relapse immunity",
                                   TRUE ~ "Baseline")) %>%
  mutate(name_short = name %>%
           str_replace_all("baseline [a-z]+", ".") %>%
           str_remove_all(", [A-Za-z]+$") %>%
           str_replace("\\., \\., \\.", "baseline") %>%
           str_remove(" phenotype") %>%
           str_replace("duration", "dormancy") %>%
           str_to_sentence(),
         name_shortest = name_short %>%
           str_remove_all("\\., ") %>%
           str_remove_all("\\.") %>%
           str_remove_all(", $") %>%
           str_replace("_", " ") %>%
           str_replace_all(", ", ",\n")
         ) %>%
  arrange(n_changes) %>%
  mutate_at(vars(matches("^name")), fct_inorder)
