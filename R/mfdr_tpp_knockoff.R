#'  Function to compute the modified false discovery rate and true positive proportion for Grace-AKO
#'
#' @param q the pre-specified false discovery rate level
#' @param w a 1 by n_covariates vector of aggregated feature statistics
#' @param t a data-dependent threshold
#' @param beta_true the true coefficient of graph-structured covariates
#' @param beta_est the estimated coefficient of graph-structured covariates
#'

mfdr_tpp_knockoff = function(q, w, t, beta_true, beta_est){
  Sa = c()
  for(m in 1:length(w)){if(w[m]<=t){Sa = c(Sa, m)}}

  true_beta = which(beta_true!=0)
  Ec = which(beta_true == 0)
  mFDR_k = length(intersect(Sa, Ec))/(length(Sa)+1/q)

  #compute TPP with knockoffs

  tpp_k= length(intersect(Sa, true_beta))
  TPP_k = tpp_k/max(1,length(true_beta))

  return(list(mFDR_withknockoff = mFDR_k, TPP_withknockoff = TPP_k, l.sa = length(Sa)))
}
