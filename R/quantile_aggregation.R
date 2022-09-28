#' Function to apply the corresponding quantile aggregation function
#'
#' @param pvals a n_bootstrap by n_covariates matrix of the feature statistics, AGCs
#' @param gamma the pre-specified percentile value used for the fixed quantile aggregation, range from 0 to 1
#' @param gamma_min the pre-specified percentile value used for the adaptive quantile aggregation, range from 0 to 1
#' @param adaptive the statement to determine the quantile aggregation function
#'

quantile_aggregation = function(pvals, gamma=0.5, gamma_min=0.00001, adaptive=FALSE){
  if (adaptive == TRUE){
    return (adaptive_quantile_aggregation(pvals, gamma_min))
  }

  else{
    return (fixed_quantile_aggregation(pvals, gamma))
  }

}
