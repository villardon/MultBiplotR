
# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

summary.Principal.Coordinates <- function(object, printdata=FALSE, printproximities=FALSE, printcoordinates=FALSE, printqualities=FALSE, ...){

  print("PRINCIPAL COORDINATES ANALYSIS")
  print(paste("Type of Data : ", object$TypeData))
  print(paste("Type of Proximity : ", object$Type))
  print(paste("Coefficient : ", object$Coefficient))
  print(paste("Transformation : ", object$Transformation))
  print("-----------")
  dims=dim(object$RowCoordinates)[2]
  
  if (printdata){ print("-----------")
                  print("DATA")
                  print(object$Data)}
  
  if (printproximities){ print("-----------")
                  print("PROXIMITIES")
                  print(object$Proximities)}
  
  print("-----------")
  print("Eigenvalues and explained Variance")
  result=matrix(0,dims,3)
  result[,1]=object$EigenValues[1:dims]
  result[,2]=object$Inertia[1:dims]
  result[,3]=cumsum(object$Inertia[1:dims])
  rownames(result)=colnames(object$RowCoordinates)
  colnames(result)=c("Eigenvalues", "Variance Explained", "Cummulative")
  print(result)
  
  if (printcoordinates){ print("-----------")
                         print("COORDINATES")
                         print(object$RowCoordinates)}
  
  if (printqualities){ print("-----------")
                         print("QUALITIES OF REPRESENTATION")
                         print(object$RowQualities)}
 
}