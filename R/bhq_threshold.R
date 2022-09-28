#' Function for applying the Standard Benjamini-Hochberg for controlling False discovery rate
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param fdr the pre-specified false discovery rate level
#'
#' @return The threshold value corresponding to the Standard Benjamini-Hochberg procedure

bhq_threshold = function(pvals, fdr=0.1){

  n_features = length(pvals)
  pvals_sorted = sort(pvals)
  selected_index = 2 * n_features
  for (i in seq(n_features, 1, -1)){
    if (pvals_sorted[i] <= (fdr * i / n_features)){
      selected_index = i
      break
    }
  }

  if (selected_index <= n_features){
    return (pvals_sorted[selected_index])
  }

  else{
    return ('-1.0')
  }
}
