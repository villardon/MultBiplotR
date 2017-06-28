# Takes a three-way list of matrices X and converts it to a two way matrix with
# the structure required by the Statis programs using a _ to separate variable and occassion
# or study. When whatlines is 1 the final matrix adds the rows of the three dimensional array, then the columns must
# be the same for all studies. When whatlines is 2 the columns are concatenated and then the number
# of rows must be the same for all studies.

Three2TwoWay <- function(X, whatlines=2){
  K=length(X)
  Studies=names(X)
  nr=matrix(0,K,1)
  nc=matrix(0,K,1)
  for (i in 1:K){
    nr[i]=dim(X[[i]])[1]
    nc[i]=dim(X[[i]])[2]
  }
  
  if ((whatlines == 2) & (sum(nr==mean(nr))!=K)) stop("The number of rows is not the same for all matrices")
  if ((whatlines == 1) & (sum(nc==mean(nc))!=K)) stop("The number of columns is not the same for all matrices")
  
  if (whatlines == 1){
    X2=X[[1]]
    labels=paste(rownames(X[[1]]), Studies[1], sep="_")
    for (i in 2:K){
      X2=rbind(X2,X[[i]])
      labels=c(labels, paste(rownames(X[[i]]), Studies[i], sep="_") )
    }
    rownames(X2)<-labels
  }
  
  if (whatlines == 2){
    X2=X[[1]]
    labels=paste(colnames(X[[1]]), Studies[1], sep="_")
    for (i in 2:K){
      X2=cbind(X2,X[[i]])
      labels=c(labels, paste(colnames(X[[i]]), Studies[i], sep="_") )
    }
    colnames(X2)<-labels
  }

  return(X2)
}


