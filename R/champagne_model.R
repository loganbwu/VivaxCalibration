library(pacman)
pacman::p_load("deSolve", "tidyverse", "viridis")

pv_champ_function <- function(t, state, parameters) { # with optional birth death
  with(as.list(c(state, parameters)), {
    # rate of change
    ds_0 <- -(1 - alpha * beta) * (lambda * (i_L + i_0) + delta) * s_0 + 
      (lambda * (i_0 + i_L) + delta) * alpha * beta * s_L + 
      alpha*beta*f*s_L + 
      gamma_L*s_L + 
      r*i_0
    
    di_0 <- -(lambda * (i_0 + i_L) + delta) * i_0 + 
      gamma_L*i_L - r*i_0
    
    ds_L <- -(1 - alpha*(1 - beta))*(lambda*(i_L + i_0) + delta + f)*s_L +
      alpha*(1-beta)*(lambda*(i_0+i_L) + delta)*s_0 - gamma_L*s_L + r*i_L
    
    di_L <- (1 - alpha)*(lambda*(i_L + i_0) + delta)*(s_0 + s_L) +
      (lambda*(i_L + i_0) + delta)*i_0 + (1 - alpha)*f*s_L - gamma_L*i_L - r*i_L
    
    # return the rate of change
    list(c(ds_0, di_0, ds_L, di_L))
  }) # end with(as.list ...
}

pv_champ_alpha <- 0.4 #prop of effective care
pv_champ_beta <- 0.4 #prop of radical cure
pv_champ_gamma_L <- 1 / 223 #liver stage clearance rate
pv_champ_delta <- 0.05 #prop of imported cases
pv_champ_lambda <- 0.04 #transmission rate
pv_champ_f <- 1 / 72 #relapse frequency
pv_champ_r <- 1 / 60 #blood stage clearance rate

pv_champ_parameters <- c(alpha = pv_champ_alpha, beta = pv_champ_beta, gamma_L = pv_champ_gamma_L, delta = pv_champ_delta, lambda = pv_champ_lambda, f = pv_champ_f, r = pv_champ_r)
state <- c(s_0 = 0.99, i_0 = 0.01, s_L = 0, i_L = 0)

pand_length <- 200

times <- seq(0, pand_length, by = 1)

pv_champ_outputs <- ode(y = state, times = times, func = pv_champ_function, parms = pv_champ_parameters) %>%
  as_tibble() %>%
  mutate(across(.fns = as.numeric))


pv_champ_long_outputs <- pivot_longer(pv_champ_outputs, -time, names_to = "state", values_to = "proportions")

pv_champ_plot <- ggplot(pv_champ_long_outputs, aes(x = time, y = proportions, colour = state)) +
  geom_line(size = 0.75) +
  theme_bw() +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since first infection") +
  ylab("Proportion in each compartment")
pv_champ_plot



# INFERENCE ------------

coeffs_fn <- function(alpha, beta, gamma_L, delta, lambda, f, r, I_eq, ...){
  c(
    (delta + gamma_L + r)*(delta*(1 - I_eq)*(delta + gamma_L + f) - r*I_eq*(delta + gamma_L + (1 - alpha*(1 - beta))*f/(1 - alpha))) + f*r*(delta + r)*I_eq,
    I_eq*(1 - I_eq)*((delta + gamma_L + f)*(delta + gamma_L + r) + delta*(2*delta + 2 * gamma_L + r + f)) - I_eq^2*r*(2*delta + 2*gamma_L + r + alpha*beta*f)/(1 - alpha),
    I_eq^2*(1 - I_eq)*(3*delta + 2*gamma_L + r + f) - I_eq^3*r/(1 - alpha),
    I_eq^3*(1 - I_eq)
  )
}

delta_fn <- function(alpha, beta, gamma_L, f, r, h) {
  
}



combinations <- expand_grid(
  alpha = c(0, 0.4, 0.8),
  beta = c(0, 0.4, 0.8),
  # p = seq(0, 1, length.out = 1001), # p_imported
  p = seq(0, 0.0001, length.out=10),
  f = pv_champ_f,
  gamma_L = pv_champ_gamma_L,
  r = pv_champ_r,
  incidence = 1:500
) %>%
  mutate(
    h = incidence/1000/365,
    I_eq = h*(1 - alpha)/r,
    delta = p*h/(1 - I_eq)
  ) %>%
  mutate(lambda = pmap(., coeffs_fn)) %>%
  mutate(lambda = map(lambda, polyroot),
         alpha=paste0("alpha:", alpha),
         beta=paste0("beta:", beta)
  ) %>%
  mutate(lambda = map(lambda, ~Re(.x[which(abs(Im(.x))<10^(-10))]))) %>%
  mutate(lambda = map_dbl(lambda, max)) %>%
  filter(lambda > 0)

ggplot(combinations, aes(x = p, y = incidence, fill = lambda)) +
  geom_raster() +
  facet_grid(beta ~ alpha, labeller = "label_parsed") +
  scale_fill_viridis() +
  theme_bw() +
  # scale_x_continuous(expand = expansion(0), limits = c(0, 1)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 500)) +
  labs(y="Incidence (per 1000 person-year)", x="Proportion of imported cases (p)", fill=expression(lambda))

# library(here)
# 
# ggsave(here("figures", "simulating_lambdas.pdf"), width = 6, height = 3)
