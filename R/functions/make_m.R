#' Seasonal function defined by Griffin et al., 2010, "Reducing Plasmodium falciparum Malaria Transmission in Africa: A Model-Based Evaluation of Intervention Strategies" and then improved in 2015, "The Interaction between Seasonality and Pulsed Interventions against Malaria in Their Effects on the Reproduction Number".
#' @param eps Ratio of peak to trough
#' @param kappa How 'peaky' it is
#' @param m_0 Average value
make_m = function(eps, kappa, phase=0, m_0=1) {
  function(t) {
    m_0 * (eps + (1-eps)*pi/beta(0.5, kappa+0.5) *((1+sin(2*pi*(t-phase)/365.25))/2)^kappa)
  }
}
