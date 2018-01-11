HJ.Biplot <- function(X, dimension = 3, Scaling = 5, sup.rows = NULL, sup.cols = NULL, grouping=NULL){
  bip=PCA.Biplot(X, alpha = 2, dimension = dimension, Scaling = Scaling, sup.rows = sup.rows, sup.cols = sup.cols, grouping=grouping)
  bip$Type = "HJ"
  bip$Title="HJ Biplot"
  bip$call <- match.call()
  return(bip)
}
  