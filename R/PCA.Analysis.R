PCA.Analysis <- function(X, dimension = 3, Scaling = 5, ...) {
  acp=PCA.Biplot(X, alpha = 1, dimension = dimension, Scaling = Scaling, ...)
  unclass(acp)
  class(acp)="PCA.Analysis"
  acp$Title="Principal Components Analysis"
  return(acp)
}
