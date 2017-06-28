# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca


CategoricalProximities<- function(Data, SUP=NULL, coefficient="GOW" , transformation=3, ...) {
  
  nrow= dim(Data)[1]
  ncol= dim(Data)[2]
  factors=TRUE
  for (j in 1:ncol)
    factors = factors & is.factor(Data[[j]])
  if (!factors) stop("You must provide a data frame with categorical factors to calculate the similarity")
  
  
  coefficients = c("GOW", "ESK", "IOF", "OF", "GOO1", "GOO2", "GOO3", "GOO4", "GAM", "LIN", 
                   "AND", "SMI")
  
  if (is.numeric(coefficient)) coefficient=coefficients[coefficient]
  
  
  transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
  
  if (is.numeric(transformation)) transformation=transformations[transformation]
  
  if (transformation==1) Type="similarity"
  else Type="dissimilarity"
  
  x=as.matrix(ConvertFactors2Integers(Data))
  
  result= list()
  result$TypeData="Nominal"
  result$Type=Type
  result$Coefficient=coefficient
  result$Transformation=transformation
  result$Data=Data
  result$SupData=SUP
  result$Proximities=CategoricalDistances(x, y=NULL, coefficient=coefficient, transformation=transformation)
  result$SupProximities=NULL
  if (!is.null(SUP)) result$SupProximities=CategoricalDistances(x,y=SUP, coefficient, transformation)
  class(result)="proximities"
  return(result)
}
