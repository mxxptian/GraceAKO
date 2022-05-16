#' Funtion to construct the Laplacian matrix based on the graphical structure
#'
#' @param L p by p (weighted) adjacency matrix of a graph.
#'

build_L1 = function(L){
  diag(L) = 0
  for (i in 1:nrow(L)) {
    for (w in i:ncol(L)) {
      if(L[i,w]!=0){
        L[i,w] = 1
        L[w,i] = L[i,w]
      }
    }
  }
  return(L)
}
