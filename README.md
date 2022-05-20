# GraceAKO

`GraceAKO` is an R package to do variable selection with guarantee of false discovery rate control. This algorithm refers to the main ideas of Network-constraint regularization proposed in 2008 by Li and Li (https://academic.oup.com/bioinformatics/article/24/9/1175/206444) and Agregation of multiple knockoffs proposed by Tuan-Binh Nguyen in 2020 (https://arxiv.org/abs/2002.09269). The main characteristics of real data are finite sample and high dimensions. Grace-AKO incorporate the prior information, the pathways lying among the genes into Laplacian matrix and utilize the aggragation quantile function to compute the threshold value for multiple knockoffs. Grace-AKO can strictly control the FDR under the predefined level without much power loss. GraceAKO is developed to help to go through the algorithm and can also return the mediated results.

# Installation


You can install the development version of
`GraceAKO` from Github via the `devtools` package. I suppose using
the `remotes` package would work as well.

``` r
devtools::install_github("mxxptian/GraceAKO")
```

# Example


``` r
# set parameters
n = 200
g = 10
k = 10
rho = 0.7
p = g*(k+1)
c1 = c(5, rep(5/sqrt(10),10), -5, rep(-5/sqrt(10), 10), 3, rep(3/sqrt(10), 10), -3, rep(-3/sqrt(10), 10))
beta_true = c(c1, rep(0, p - length(c1)))

lambda.L = seq(0.1,2,0.5)
lambda.1 =  seq(110,150,10)
lambda.2 = seq(1,10,5)

fdr = 0.1
n_bootstraps = 25
gamma = 0.1

#generate data
l <- matrix(0,nrow = k+1,ncol = k+1)
l[1,] <- l[,1] <- rep(-1,k+1)
diag(l) <- 1
l[1,1] <- k

L <- matrix(0,nrow = g*(k+1), ncol=g*(k+1))
for (i in 1:g) {
  L[((i-1)*(k+1)+1):(i*(k+1)),((i-1)*(k+1)+1):(i*(k+1))] <- l
}

i = 1
X = matrix(0, nrow = n, ncol = g*(k+1))
for (i in 1:n) {
  TF <- rnorm(g)
  for (j in 1:g) {
    genes <- rnorm(k, rho*TF[j], sqrt(1-rho^2))
    X[i,((j-1)*(k+1)+1):(j*(k+1))] <- c(TF[j],genes)
  }
}

sigma.error = sqrt(sum(beta_true[which(beta_true!=0)]^2)/4)

Y <- c(X %*% beta_true + rnorm(n, sd = sigma.error))

L_1 = build_L1(L)

library(dplyr)
X = as.data.frame(X)
std.X = X %>% mutate_all(~scale((.)%>% as.vector))

Y = as.data.frame(Y)
std.Y = Y %>% mutate_all(~scale((.)%>% as.vector))

X = as.matrix(std.X)
Y1 = as.vector(std.Y)$Y
Y = as.vector(Y1)

#conduct Grace-AKO and Grace

result = Grace_multi_knockoff(X = X, Y = Y, L_1 = L_1, L = L, fdr = fdr, 
                              lambda.L = lambda.L, lambda.1 = lambda.1, 
                              lambda.2= lambda.2, n_bootstraps = n_bootstraps, gamma = gamma)
mfdr_tpp_base = mfdr_tpp_base(q=fdr, beta_true, result$beta_est_baseline)
mfdr_tpp_k = mfdr_tpp_knockoff(q = fdr, result$w, result$t, beta_true, result$beta_est_knockoff)
```
