#' Function to apply the Benjamini-Hochberg-Yekutieli procedure for controlling FDR
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param reshaping_function the default value for reshaping function can be referred by Benjamini & Yekutieli (2001) otherwise, the reshaping function can be referred by Ramdas et al (2017)
#' @param fdr the pre-specified false discovery rate level
#' @return The threshold value corresponding tthe Benjamini-Hochberg-Yekutieli procedure


bhy_threshold = function(pvals,
                         reshaping_function = NULL, fdr = 0.1)
{

  n_features = length(pvals)
  p_vals_sorted = sort(pvals)
  selected_index = 2 * n_features


  if (is.null(reshaping_function)==TRUE)
  {
    temp = seq(n_features,1)
    sum_inverse = 0
    for (i in temp) {
      sum_inverse = sum_inverse + 1 / i
    }

    return (bhq_threshold(pvals, fdr / sum_inverse))
  }
  else{
    for (i in seq(n_features - 1, 0, -1)){
      if (p_vals_sorted[i] <= fdr * reshaping_function(i) / n_features){
        selected_index = i
        break
      }

    }

    if (selected_index <= n_features){
      return (p_vals_sorted[selected_index])
    }
    else{
      return ('-1.0')
    }

  }

}

