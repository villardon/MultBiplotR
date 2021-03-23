# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2017
# Integer is treated as numeric unless otherwise is specified
# 

GowerProximities<- function(x, y=NULL, Binary=NULL, Classes=NULL, transformation=3, IntegerAsOrdinal=FALSE, BinCoef= "Simple_Matching", ContCoef="Gower", NomCoef="GOW", OrdCoef="GOW") {
  if (!is.data.frame(x)) stop("Main data is not organized as a data frame")
  
  NewX=AdaptDataFrame(x, Binary=Binary, IntegerAsOrdinal=IntegerAsOrdinal)
  
  if (is.null(y)) NewY=NewX
  else{
    if (!is.data.frame(y)) stop("Suplementary data is not organized as a data frame")
    NewY=AdaptDataFrame(y, Binary=Binary, IntegerAsOrdinal=IntegerAsOrdinal)
  }
  
  n = dim(NewX$X)[1]
  p = dim(NewX$X)[2]
  
  n1 = dim(NewY$X)[1]
  p1 = dim(NewY$X)[2]
  
  if (!(p==p1)) stop("Number of columns of the two matrices are not the same")
  
  
  transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
  if (is.numeric(transformation)) transformation=transformations[transformation]
  
  if (transformation==1) Type="similarity"
  else Type="dissimilarity"
  
  if ( (BinCoef== "Simple_Matching") & (ContCoef=="Gower") & (NomCoef=="GOW") & (OrdCoef=="GOW"))
  coefficient="Gower Similarity"
  else
    paste("Binary: ",BinCoef, ", Continuous: ", ContCoef, ", Nominal: ", NomCoef, ", Ordinal: ", OrdCoef)
  result= list()
  result$TypeData="Mixed"
  result$Type=Type
  result$Coefficient=coefficient
  result$Transformation=transformation
  result$Data=NewX$X
  result$SupData=NewY$X
  result$Types=NewX$Types
  result$Proximities=GowerSimilarities(NewX$X, y=NewY$X,  transformation=transformation, Classes=NewX$Types, BinCoef= BinCoef, ContCoef=ContCoef, NomCoef=NomCoef, OrdCoef=OrdCoef)
  rownames(result$Proximities)=rownames(x)
  colnames(result$Proximities)=rownames(x)
  result$SupProximities=NULL
  if (!is.null(y)) result$SupProximities=GowerSimilarities(x,y, transformation)
  class(result)="proximities"
  return(result)
}
