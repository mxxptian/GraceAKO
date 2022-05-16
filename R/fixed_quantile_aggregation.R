#' Function to apply quantile aggregation function based on Meinshausen et al (2008)
#'
#' @param pvals a n_bootstrap by n_covariates matrix of the feature statistics, AGCs
#' @param gamma the pre-specified percentile value used for aggregation, range from 0 to 1

fixed_quantile_aggregation = function(pvals, gamma = 0.3){


  converted_score = (1 / gamma) *  quantile(pvals, gamma)

  return (min(1, converted_score))
}
