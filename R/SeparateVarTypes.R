
SeparateVarTypes <- function(X, TypeVar=NULL, TypeFit=NULL){
  if (!is.data.frame(X)) stop("X must be a data frame")
  VarTypes=c("c","b","n","o","f","a")
  FitTypes=c("wa","r","g","g1")

  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(TypeVar)){
    for (i in 1:p){
      if (class(X[[i]])=="numeric") TypeVar=c(TypeVar,"c")
      else
        if (class(X[[i]])=="factor") {
          if (length(levels(X[[i]]))==2) TypeVar=c(TypeVar,"b")
          else
            TypeVar=c(TypeVar,"n")}
      else 
        if (class(X[[i]])=="factor") TypeVar=c(TypeVar,"b")
    }
  }
  else
    if (length(TypeVar)==1) TypeVar=rep(TypeVar,p)
  
  if (length(TypeVar)!=p) stop("The vector with the Variable Types should have the same length as the number of variables")
  if (is.numeric(TypeVar)) TypeVar=VarTypes[TypeVar]
  pp=unique(TypeVar)
  if (!setequal(intersect(pp, VarTypes), pp)) {message=cat("Can not recognize some variable types: ",setdiff(pp,intersect(pp, VarTypes)))
                                               stop(message)}
  
  if (is.null(TypeFit)){
    for (i in 1:p){
      if (TypeVar[i]=="c") TypeFit=c(TypeFit,"r")
      if (TypeVar[i]=="b") TypeFit=c(TypeFit,"r")
      if (TypeVar[i]=="n") TypeFit=c(TypeFit,"r")
      if (TypeVar[i]=="o") TypeFit=c(TypeFit,"r")
      if (TypeVar[i]=="f") TypeFit=c(TypeFit,"wa")
      if (TypeVar[i]=="a") TypeFit=c(TypeFit,"g")
    }
  }
  else
    if (length(TypeFit)==1) TypeFit=rep(TypeFit,p)
  
  if (length(TypeFit)!=p) stop("The vector with the Variable Fittings should have the same length as the number of variables")
  if (is.numeric(TypeFit)) TypeFit=FitTypes[TypeFit]
  pp=unique(TypeFit)
  if (!setequal(intersect(pp, FitTypes), pp)) {message=cat("Can not recognize some fits: ",setdiff(pp,intersect(pp, FitTypes)))
                                               stop(message)}
  
  x=list("Continuous","Binary","Nominal","Ordinal","Frequency","Abundance")
  
  x$Continuous$Data=X[which(TypeVar=="c")]
  if (length(x$Continuous$Data)==0) x$Continuous=NULL
  else x$Continuous$Fit=TypeFit[which(TypeVar=="c")]
  
  x$Binary$Data=X[which(TypeVar=="b")]
  if (length(x$Binary$Data)==0) x$Binary=NULL
  else x$Binary$Fit=TypeFit[which(TypeVar=="b")]
  
  x$Nominal$Data=X[which(TypeVar=="n")]
  if (length(x$Nominal$Data)==0) x$Nominal=NULL
  else x$Nominal$Fit=TypeFit[which(TypeVar=="n")]
  
  x$Ordinal$Data=X[which(TypeVar=="o")]
  if (length(x$Ordinal$Data)==0) x$Ordinal=NULL
  else x$Ordinal$Fit=TypeFit[which(TypeVar=="o")]
  
  x$Frequency$Data=X[which(TypeVar=="f")]
  if (length(x$Frequency$Data)==0) x$Frequency=NULL
  else x$Frequency$Fit=TypeFit[which(TypeVar=="f")]
  
  x$Abundance$Data=X[which(TypeVar=="f")]
  if (length(x$Abundance$Data)==0) x$Abundance=NULL
  else x$Abundance$Fit=TypeFit[which(TypeVar=="f")]
  
  class(x)="SeparatedData"
  return(x)
}

