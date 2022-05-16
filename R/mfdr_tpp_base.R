#' Function to compute the modified false discovery rate and true positive proportion for baseline model (e.g., Grace)
#'
#' @param q the pre-specified false discovery rate level
#' @param beta_true the true coefficient of graph-structured covariates
#' @param beta_est the estimated coefficient of graph-structured covariates
#'
#'

mfdr_tpp_base = function(q, beta_true, beta_est){
  true_beta = which(beta_true!=0)
  Ec = which(beta_true ==0)

  Sa_b = which(beta_est!=0)
  mFDR_b = length(intersect(Sa_b, Ec))/(length(Sa_b)+1/q)

  #compute TPP
  t_base = length(intersect(Sa_b, true_beta))
  TPP_base = t_base/max(1,length(true_beta))

  return(list(mFDR_baseline = mFDR_b, TPP_baseline = TPP_base))
}
