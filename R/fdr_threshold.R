#' Function to compute the corresponding threshold (BY/BH) for controlling FDR
#'
#' @param pvals a vector of the feature statistics, AGCs
#' @param fdr the pre-specified false discovery rate level
#' @param method the test procedure applied to compute the threshold. The defult method to compute threshold is the the Standard Benjamini-Hochberg procedure.
#' @param reshaping_function the default value for reshaping function can be referred by Benjamini & Yekutieli (2001) otherwise, the reshaping function can be referred by Ramdas et al (2017)
#'
#' @return The threshold value corresponding to the test procedure

fdr_threshold = function(pvals, fdr=0.1, method='bhq', reshaping_function=NULL){
  if (method == 'bhq'){
    # pvals_bhq = as.vector(unlist(pvals))
    return (bhq_threshold(pvals, fdr=fdr))
  }
  else{
    if(method == 'bhy'){
      # pvals_bhy = as.vector(unlist(pvals))
      return( bhy_threshold(
        pvals, fdr=fdr, reshaping_function=reshaping_function))
    }
    else{
      return('{} is not support FDR control method')
    }
  }
}
