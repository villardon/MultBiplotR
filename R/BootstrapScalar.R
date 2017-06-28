BootstrapScalar <- function(B, W=diag(nrow(B)), nB=200, dimsol=2, ProcrustesRot=TRUE, method=c("Sampling", "Permutation")){
  if (length(method) > 1) method = method[1] 
  if (nrow(B)!=ncol(B)) stop("The scalar pruducts matrix must be squared")
  n=nrow(B)
  B <- as.matrix(B)
  
  SvdD <- svd(sqrt(W) %*% B %*% sqrt(W))
  PO <- solve(sqrt(W)) %*% SvdD$u%*%diag(sqrt(SvdD$d))
  D<-EuclideanDistance(PO)
  PO <- PO[,1:dimsol]
  Bhat <- PO %*% t(PO)
  Er <- B-Bhat
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
    BE <- as.matrix(ErE) + Bhat - diag(diag(Bhat)) + diag(diag(B))
    SvdDE <- svd(sqrt(W) %*% BE %*% sqrt(W))
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
  Info=paste("Bootstrap for PCoA based on ", method, "of scalar products")
  result=list(Info=Info,InitialDistance=D,Eigenvalues=Lambdas, Inertias=Inertias, Coordinates=Coordinates, Qualities=Qualities, NReplicates=nB)
  class(result)="PCoABootstrap"
  return(result)
}

