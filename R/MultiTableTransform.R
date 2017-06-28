MultiTableTransform <- function(X, InitTransform="Standardize columns") 
{
  ng = length(X) #Number of groups
  VarNames=NULL
  for (i in 1:ng){
    X[[i]] = TransformIni(X[[i]],InitTransform)}
  return(X)
}