#' Function to conduct Grace-AKO and Grace
#'
#' @import glmnet
#' @import lattice
#' @import dplyr
#' @import stats
#' @import quantile
#' @import Grace
#' @import snow
#' @import MonteCarlo
#' @import MASS
#' @import knockoff
#' @import corpcor
#' @import miscTools
#' @import pracma
#' @import ModelMetrics
#'
#' @param X standardized n (number of rows) by p (number of columns) design matrix.
#' @param Y centered n by 1 vector of the response variable.
#' @param L_1 p by p (weighted) processed matrix of a graph for Grace-AKO.
#' @param L p by p (weighted) adjacency matrix of a graph for Grace.
#' @param fdr the pre-specified false discovery rate level
#' @param lambda.L tuning parameters of the penalty weight matrix
#' @param lambda.1 tuning parameters of the L_1 penalty
#' @param lambda.2 tuning parameters of the ridge penalty
#' @param n_bootstraps the number of multiple knockoffs
#' @param gamma the pre-specified percentile value used for fixed aggregation
#' @param gamma_min the pre-specified percentile value used for adaptive aggregation

Grace_multi_knockoff = function(X,Y,L_1,L, fdr, lambda.L, lambda.1,
                                lambda.2, n_bootstraps = 25, gamma = 0.3, gamma_min = 0.1){
  A <- abs(L_1)

  ########################L estimated by data structure####################
  A_k1 = cbind(A, A)
  A_k2 = A_k1
  A_k = rbind(A_k1,A_k2)



  d_prec = apply(A_k,2,sum)


  L_prec = matrix(0, nrow = length(d_prec), ncol = length(d_prec))
  for (i in 1:length(d_prec)) {
    for (w in i:length(d_prec)) {
      if(i!=w){
        if(d_prec[i]!=0 & d_prec[w]!=0){

          L_prec[i,w] = -A_k[i,w]/sqrt(abs(d_prec[i]*d_prec[w]))
          L_prec[w,i] = L_prec[i,w]
        }
        else{
          L_prec[i,w]=0
          L_prec[w,i] = L_prec[i,w]
        }
      }

    }
  }
  diag(L_prec)[(d_prec!=0)] <- 1


  ###########################Grace baseline#######################
  #
  # select_var_coef_prec = grace.test.result_prec$beta
  #
  # Sa_grace_only = which(select_var_coef_prec != 0)

  # lasso = cv.glmnet(X,Y, alpha = 1)
  # ridge = cv.glmnet(X,Y,alpha = 0)

  cv_baseline = cvGrace(as.matrix(X), Y,
                        L,
                        lambda.L = lambda.L,
                        lambda.1 = lambda.1,
                        lambda.2 = lambda.2)
  #baseline
  grace_linear_out = grace(Y,
                           as.matrix(X),
                           L,
                           lambda.L = cv_baseline[1],
                           lambda.1  = cv_baseline[2],
                           lambda.2 = cv_baseline[3])

  beta_est_baseline = grace_linear_out$beta



  ###########################Grace with multi-knockoff#######################

  # n_bootstraps = 25
  pvals = matrix(0, ncol(X), n_bootstraps)
  #Multiple Knockoffs
  for (i in 1:n_bootstraps) {
    print(i)

    mu = apply(X,2,mean)


    covariance = tryCatch({suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,verbose=F)), nrow=ncol(X)))},
                          warning = function(w){}, error = function(e) {
                            stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
                          }, finally = {})

    gaussian_knockoffs = create.gaussian(X, mu, covariance)

    new_xg = cbind(Y, X, gaussian_knockoffs)


    cv_knockoff = cvGrace(as.matrix(new_xg[,-1]), new_xg[,1],  L_prec,
                          lambda.L = lambda.L,
                          lambda.1 = lambda.1,
                          lambda.2 = lambda.2)


    grace_linear_knockoff = grace(new_xg[,1], as.matrix(new_xg[,-1]), L_prec,
                                  lambda.L = cv_knockoff[1],
                                  lambda.1  = cv_knockoff[2],
                                  lambda.2 = cv_knockoff[3])


    beta_est_knockoff = grace_linear_knockoff$beta

    w = c()
    par = ncol(X)
    for (m in 1:ncol(X)) {
      w = c(w, abs(beta_est_knockoff[m])-abs(beta_est_knockoff[m+par]))
    }



    pvals[,i] = empirical_pval(w, offset = 0)

  }
  aggregated_pval = apply(pvals, 1, quantile_aggregation, gamma=gamma, gamma_min=gamma_min, adaptive=FALSE)

  threshold = fdr_threshold(aggregated_pval, fdr=fdr, method='bhq',
                            reshaping_function=NULL)



  return(list(t = threshold,
              w = aggregated_pval,
              beta_est_baseline = beta_est_baseline,
              beta_est_knockoff = beta_est_knockoff))
}

