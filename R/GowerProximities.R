# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Enero/2014

GowerProximities<- function(x, y=NULL, transformation=3) {
  
  transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
  if (is.numeric(transformation)) transformation=transformations[transformation]
  
  if (transformation==1) Type="similarity"
  else Type="dissimilarity"
  
  result= list()
  result$TypeData="Mixed"
  result$Type=Type
  result$Coefficient=NULL
  result$Transformation=transformation
  result$Data=x
  result$SupData=y
  result$Proximities=GowerSimilarities(x, y=NULL, transformation)
  result$SupProximities=NULL
  if (!is.null(y)) result$SupProximities=GowerSimilarities(x,y, transformation)
  class(result)="proximities"
  return(result)
}
