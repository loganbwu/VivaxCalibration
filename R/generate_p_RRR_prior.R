# Based on "The probability of Plasmodium vivax acute malarial illness following primary infection and relapse" - Ross et al.

# For Model variant B1^3 (Table 5), Ross et al. find that the 'ratio of probability for relapse compared to primary infection' is:
# shorter blood stage duration: 0.31 (mean; 0.09-0.68 for 95% credible interval)
# longer blood stage duration: 0.49 (mean; 0.13-0.93 for 95% credible interval)
# They use two parameter sets (shorter and longer) because of uncertainty in the duration of blood-stage infection.
# "The shorter estimate of 27 days was informed by a study estimating the mean durations of blood-stage infection in a cohort of children in endemic areas in Albinama, Papua New Guinea and in Thailand [7] (input parameter set 1). The longer duration comes from neurosyphilis patients deliberately infected as malaria therapy in the 1940s and 1950s. Among those infected for the first time with P. vivax via blood, a median of 76 days was estimated (input parameter set 2)." - suppporting information.
# Our parameter p_{RRR} is the relative reduction of relapse - 1 minus Ross's figure. Therefore:
# shorter blood stage duration: 0.69 (0.32-0.91), or
# longer blood stage duration: 0.51 (0.07-0.87).

# My model uses a beta distribution for the priors. We will use a simple search algorithm to derive two priors that match these results. Note that we have three data points but only two free parameters in the beta distribution so it may not be exact. Also, we do not know exactly what kind of credible intervals Ross used, but we assume equal-tailed.
# Findings: we match all three within a tolerance.

blood_stage = "longer"

if (blood_stage == "shorter") {
  ymean = 0.69
  ylower = 0.32
  yupper = 0.91
} else if (blood_stage == "longer") {
  ymean = 0.51
  ylower = 0.07
  yupper = 0.87
}

loss_fn = function(x) {
  shape1 = x[1]
  shape2 = x[2]
  lower = qbeta(0.025, shape1, shape2)
  upper = qbeta(0.975, shape1, shape2)
  .mean = shape1 / (shape1 + shape2)
  loss = (lower-ylower)^2 + (upper-yupper)^2 + 1*(.mean-ymean)^2
  return(loss)
}

result = optim(c(2, 2), loss_fn)
shape_optim = result$par

# Calculate mean
mean_optim = shape_optim[1] / (shape_optim[1] + shape_optim[2])
quantiles_optim = qbeta(c(0.025, 0.975), shape_optim[1], shape_optim[2])

# plot
xx = seq(0, 1, length.out=500)
yy = dbeta(xx, shape_optim[1], shape_optim[2])
plot(xx, yy, 'l')

# Findings - check against the comments on line 9-10; should be an APPROXIMATE match:
cat(paste0(blood_stage, " blood stage parameters: [", paste(round(shape_optim, 3), collapse=", "), "]\n", blood_stage, " blood stage p_RRR statistics: ", round(mean_optim, 3), " (", paste(round(quantiles_optim, 3), collapse="-"), " CI)"))
