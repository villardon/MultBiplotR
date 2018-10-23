# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

ContinuousProximities<- function(x, y=NULL, ysup=FALSE, transpose=FALSE, coef = "Pythagorean", r = 1) {
  
  distances = c("Pythagorean", "Taxonomic", "City", "Minkowski", "Divergence", "dif_sum", "Camberra", "Bray_Curtis", "Soergel", "Ware_Hedges", "Gower")
  if (is.numeric(coef)) coef = distances[coef]
  
  if (!is.null(y)){
    if (!(ncol(x)==ncol(y))) stop("Columns don't match")
  }
  
  Type="dissimilarity"
  
  if (is.null(y)) {Shape="Squared"}
  else{
    if (ysup) Shape="Squared"
    else Shape="Rectagular"}
  
  result= list()
  result$TypeData="Continuous"
  result$Type=Type
  result$Coefficient=coef
  result$Data=x
  result$SupData=y
  result$r=r
  
  if ((!is.null(y) & (!ysup)))
  result$Proximities=ContinuousDistances(x, y, coef = coef, r = r)
    else
  result$Proximities=SymmetricContinuousDistances(x,  coef = coef, r = r)

  result$SupProximities=NULL
  if (!is.null(y) & (ysup)) result$SupProximities=ContinuousDistances(x, y,   coef = coef, r = r)
  class(result)="proximities"
  return(result)
}
