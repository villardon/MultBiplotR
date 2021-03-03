MultiTableStatistics <- function(X, dual=FALSE) 
{
  ng = length(X) #Number of groups
  StudyNames=names(X)
  if (!dual){
    x=NULL
    VarNames=NULL
    for (i in 1:ng){
      VarNames=c(VarNames,paste(colnames(X[[i]]),"_",StudyNames[i],sep=""))
      x=cbind(x, X[[i]])}
    colnames(x)=VarNames
    Statistics=list()
    Statistics$Non_Scaled_Data=x
    Statistics$Means = apply(x, 2, mean)
    Statistics$Medians = apply(x, 2, median)
    Statistics$Deviations = apply(x, 2, sd)
    Statistics$Minima = apply(x, 2, min)
    Statistics$Maxima = apply(x, 2, max)
    Statistics$P25 = apply(x, 2, quantile)[2, ]
    Statistics$P75 = apply(x, 2, quantile)[4, ]
    Statistics$GMean = mean(x)
    Statistics$nrows = dim(x[1])
    Statistics$ncols = dim(x[2])
  }
  
  if (dual){
    x=NULL
    IndNames=NULL
    for (i in 1:ng){
      IndNames=c(IndNames,paste(rownames(X[[i]]),"_",StudyNames[i],sep=""))
      x=rbind(x, X[[i]])
    }
    rownames(x)=IndNames
    Statistics=list()
    Statistics$Non_Scaled_Data=x
    Statistics$Means = apply(x, 2, mean)
    Statistics$Medians = apply(x, 2, median)
    Statistics$Deviations = apply(x, 2, sd)
    Statistics$Minima = apply(x, 2, min)
    Statistics$Maxima = apply(x, 2, max)
    Statistics$P25 = apply(x, 2, quantile)[2, ]
    Statistics$P75 = apply(x, 2, quantile)[4, ]
    Statistics$GMean = mean(x)
    Statistics$nrows = dim(x[1])
    Statistics$ncols = dim(x[2])
  }
  
  return(Statistics)
}

