#' Function of intermediary products derived from the conversion of feature statistics (single knockoffs) to aggregated feature statistics.
#'
#' @param test_score a vector of the feature statistics obatianed from single knockoff procedure
#' @param offset a modified value of the conversion procedure
#'
#' @return a vector of the intermediary products

empirical_pval = function(test_score, offset = 1){
  pvals = c()
  n_features = length(test_score)
  if (offset !=0 && offset!=1){
    return("'offset' must be either 0 or 1")
  }
  else{
    test_score_inv = -test_score
    for (i in 1:n_features){
      if (test_score[i] <= 0){
        pvals = c(pvals, 1)
      }
      else{

        pvals = c(pvals,(offset+sum(test_score_inv[i] >= test_score))/n_features)
      }
    }
  }
  return (pvals)
}

