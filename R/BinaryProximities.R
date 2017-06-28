# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

BinaryProximities<- function(x, y=NULL, coefficient="Jaccard" , transformation=NULL, transpose=FALSE, ...) {
  
  if (is.data.frame(x)) x=Dataframe2BinaryMatrix(x, ...)
  
  if (transpose) x= t(x)
  
  if (!CheckBinaryMatrix(x)) stop("Input must be a binary matrix")
  if (!is.null(y)){
    if (is.data.frame(y)) y=Dataframe2BinaryMatrix(y, ...)
    if (!CheckBinaryMatrix(y)) stop("Input must be a binary matrix")
    if (transpose) y= t(y)
    if (!(ncol(x)==ncol(y))) stop("Columns don't match")
  }
  
  
  coefficients = c("Kulezynski", "Russell_and_Rao", "Jaccard", "Simple_Matching", "Anderberg", "Rogers_and_Tanimoto", "Sorensen_Dice_and_Czekanowski", 
                   "Sneath_and_Sokal", "Hamman", "Kulezynski2", "Anderberg2", "Ochiai", "S13", "Pearson_phi", "Yule", "Czekanowski", "Dice")
  
  if (is.numeric(coefficient)) coefficient=coefficients[coefficient]
  
  
  if (is.null(transformation))
    switch(coefficient, `Kulezynski` = {
      transformation=9
    },`Russell_and_Rao` = {
      transformation=3
    },`Jaccard` = {
      transformation=3
    },`Simple_Matching` = {
      transformation=3
    },`Anderberg` = {
      transformation=3
    },`Rogers_and_Tanimoto` = {
      transformation=3
    },`Sorensen_Dice_and_Czekanowski` = {
      transformation=3
    },`Sneath_and_Sokal` = {
      transformation=3
    },`Hamman` = {
      transformation=7
    },`Kulezynski2` = {
      transformation=3
    },`Anderberg2` = {
      transformation=3
    },`Ochiai` = {
      transformation=3
    },`S13` = {
      transformation=3
    },`Pearson_phi` = {
      transformation=7
    },`Yule` = {
      transformation=7
    },`Sorensen` = {
      transformation=3
    },`Dice` = {
      transformation=3
    })
  
  transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
  if (is.numeric(transformation)) transformation=transformations[transformation]
  
  if (transformation==1) Type="similarity"
  else Type="dissimilarity"
  
  if (is.null(y)) {Shape="Squared"}
  else{Shape="Rectagular"}
  
  result= list()
  result$TypeData="Binary"
  result$Type=Type
  result$Coefficient=coefficient
  result$Transformation=transformation
  result$Data=x
  result$SupData=y
  result$Proximities=BinaryDistances(x, y=NULL, coefficient, transformation)
  result$SupProximities=NULL
  if (!is.null(y)) result$SupProximities=BinaryDistances(x,y, coefficient, transformation)
  class(result)="proximities"
  return(result)
}
