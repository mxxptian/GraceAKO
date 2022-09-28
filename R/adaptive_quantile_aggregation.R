#' Function to apply adaptive version of the quantile aggregation method, Meinshausen et al. (2008)
#'
#' @param pvals a n_bootstrap by n_covariates matrix of the feature statistics, AGCs
#' @param gamma_min the pre-specified percentile value used for aggregation, range from 0 to 1

adaptive_quantile_aggregation = function(pvals, gamma_min=0.05){

  gammas = seq(gamma_min, 1.05, 0.05)
  list_Q = matrix(0,nrow = length(gammas), ncol = length(gammas))
  for (gamma in gammas) {
    list_Q = c(list_Q, fixed_quantile_aggregation(pvals, gamma))

  }

  return (min(1, (1 - log(gamma_min)) * min(list_Q)))

}
