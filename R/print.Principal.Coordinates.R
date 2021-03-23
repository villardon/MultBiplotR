
# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

print.Principal.Coordinates <- function(x, ...){
  
  print("PRINCIPAL COORDINATES ANALYSIS")
  print(paste("Type of Data : ", x$TypeData))
  print(paste("Type of Proximity : ", x$Type))
  print(paste("Coefficient : ", x$Coefficient))
  print(paste("Transformation : ", x$Transformation))
  dims=dim(x$RowCoordinates)[2]
  
  print("-----------")
  print("Eigenvalues and explained Variance")
  result=matrix(0,dims,3)
  result[,1]=x$EigenValues[1:dims]
  result[,2]=x$Inertia[1:dims]
  result[,3]=cumsum(x$Inertia[1:dims])
  rownames(result)=colnames(x$RowCoordinates)
  colnames(result)=c("Eigenvalues", "Variance Explained", "Cummulative")
  print(result)
}