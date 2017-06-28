BootstrapDistance <- function(D, W=diag(nrow(D)), nB=200, dimsol=2, ProcrustesRot=TRUE, method=c("Sampling", "Permutation")){
    if (length(method) > 1) method = method[1] 
    if (nrow(D)!=ncol(D)) stop("The distance matrix must be squared")
    if (sum(diag(D))>0) stop("That is not a distance matrix (no zeros in the diagonal) ")
    n=nrow(D)
    D <- as.matrix(D)
    B=-0.5*(diag(n)-matrix(1,n,n)/n) %*% D^2 %*% (diag(n)-matrix(1, n, n)/n) 
   
    SvdD <- svd(sqrt(W) %*% B %*% sqrt(W))
    PO <- solve(sqrt(W)) %*% SvdD$u%*%diag(sqrt(SvdD$d))
    PO <- PO[,1:dimsol]
    Bhat <- PO %*% t(PO)
    Dhat=EuclideanDistance(PO)
    Er <- D-Dhat
    Er <- as.dist(Er)
    Signos <- as.dist(matrix(1,n,n))
    
    Lambdas=matrix(0,n,nB)
    rownames(Lambdas)=rownames(D)
    Inertias=matrix(0,n,nB)
    Coordinates=list()
    Qualities=list()
    for (i in 1:n) {
      Coordinates[[i]]=matrix(0,dimsol,nB)
      rownames(Coordinates[[i]])=paste("Dim_",1:dimsol,sep="")
      colnames(Coordinates[[i]])=paste("Rep_",1:nB,sep="")
      Qualities[[i]]=matrix(0,dimsol,nB)
      rownames(Qualities[[i]])=paste("Dim_",1:dimsol,sep="")
      colnames(Qualities[[i]])=paste("Rep_",1:nB,sep="")
    }
    names(Coordinates)<-rownames(D)
    names(Qualities)<-rownames(D)
    
    for (i in 1: nB){
      if (method=="Sampling") ErE <- (sample(Er,length(Er),replace = TRUE))*Signos
      else{
      mm=sample.int(n*(n-1)/2, size = n*(n-1)/2, replace = FALSE, prob = NULL)
      ErE=Er[mm]*Signos}
      DE <- as.matrix(ErE) + Dhat 
      B=-0.5*(diag(n)-matrix(1,n,n)/n) %*% DE^2 %*% (diag(n)-matrix(1, n, n)/n) 
      SvdDE <- svd(sqrt(W) %*% B %*% sqrt(W))
      Lambdas[,i] <- SvdDE$d
      Inertias[,i] <- (SvdDE$d/sum(SvdDE$d))*100
      POE <- solve(sqrt(W)) %*%  SvdDE$u %*% diag(sqrt(SvdDE$d))
      SumCuadTot=apply(POE^2,1,sum)
      POE <- POE[,1:dimsol]
      rownames(POE) <- rownames(D)
      if (ProcrustesRot)
        POE=SimpleProcrustes(PO,POE)$Yrot
      colnames(POE) <- paste("Dim", 1:dimsol)
      RQI=((diag(1/SumCuadTot)) %*% POE[,1:dimsol]^2)*100
      rownames(RQI) <- rownames(D)
      for (j in 1:n)
        for (k in 1:dimsol){
          Coordinates[[j]][k,i]=POE[j,k]
          Qualities[[j]][k,i]=RQI[j,k]}
    }
    Info=paste("Bootstrap for PCoA based on ", method, "of distances")
    result=list(Info=Info, InitialDistance=D,Eigenvalues=Lambdas, Inertias=Inertias, Coordinates=Coordinates, Qualities=Qualities, NReplicates=nB)
    class(result)="PCoABootstrap"
    return(result)
  }


