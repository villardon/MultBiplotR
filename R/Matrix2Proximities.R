#Converts a distance or a similarity matrix into an object of class "proximities" 
# suitable for Principal Coordinates or MDS (SMACOF)

Matrix2Proximities <- function(x, TypeData="User Provided",Type=c("dissimilarity", "similarity", "products"), Coefficient="None", Transformation="None", Data=NULL){
   if (class(x)=="dist") x=as.matrix(x)
  if (!is.matrix(x)) stop("The input is not a matrix")
  if (!(nrow(x)==ncol(x))) stop("The matrix must be squared")
  if (length(Type)>1) Type=Type[1]
  r=nrow(x)*(nrow(x))/2
  symmetrical=TRUE
  k=0
  while ((symmetrical) & (k<r))
    for (i in 1:(nrow(x)-1))
      for (j in (i+1):nrow(x)){
        k=k+1
        if (!(x[i,j]==x[j,i])) symmetrical=FALSE}
    if (!symmetrical) stop("The matrix must be symmetrical")
  if ((Type=="dissimilarity") & (sum(diag(x))>0)) stop("The diagonal elements are not zero")
  if ((Type=="product") & (sum(diag(x)>0))) stop("The diagonal elements must be positive")
  if ((Type=="products") & (sum(diag(x)>0))) stop("The diagonal elements must be positive")
      
      Res=list(TypeData=TypeData, Type=Type, Coefficient=Coefficient, Transformation=Transformation, Data=Data, Proximities=x)
      class(Res)="proximities"
  return(Res)
      
}