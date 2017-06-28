BootstrapSmacof <- function(D, W=NULL, Model=c("Identity", "Ratio", "Interval", "Ordinal"), dimsol=2, maxiter=100, maxerror=0.000001, StandardizeDisparities=TRUE, ShowIter=TRUE,
                            nB=200, ProcrustesRot=TRUE, method=c("Sampling", "Permutation")){
  
  if (length(Model) > 1) Model = Model[1] 
  if (length(method) > 1) method = method[1] 
  if (nrow(D)!=ncol(D)) stop("The distance matrix must be squared")
  if (sum(diag(D))>0) stop("That is not a distance matrix (no zeros in the diagonal) ")
  
  
  D <- as.matrix(D)
  n=nrow(D)
  if (is.null(W)) W = matrix(1,n,n)
  
  Sol <- SMACOF(D, W=W, Model=Model, dimsol=dimsol, maxiter=maxiter, maxerror=maxerror, StandardizeDisparities=StandardizeDisparities, ShowIter=ShowIter)
  
  PO <- Sol$X
  PO <- PO[,1:dimsol]
  Dhat=Sol$Dhat
  Er <- D-Dhat
  Er <- as.dist(Er)
  Signos <- as.dist(matrix(1,n,n))
  
  RawStress=matrix(0,nB)
  stress1=matrix(0,nB)
  stress2=matrix(0,nB)
  sstress1=matrix(0,nB)
  sstress2=matrix(0,nB)
  Coordinates=list()
  Qualities=list()
  for (i in 1:n) {
    Coordinates[[i]]=matrix(0,dimsol,nB)
    rownames(Coordinates[[i]])=paste("Dim_",1:dimsol,sep="")
    colnames(Coordinates[[i]])=paste("Rep_",1:nB,sep="")
  }
  names(Coordinates)<-rownames(D)
  
  
  for (i in 1: nB){
    if (method=="Sampling") ErE <- (sample(Er,length(Er),replace = TRUE))*Signos
    else{
      mm=sample.int(n*(n-1)/2, size = n*(n-1)/2, replace = FALSE, prob = NULL)
      ErE=Er[mm]*Signos}
    DE <- as.matrix(ErE) + Dhat 
    
    Sol <- SMACOF(DE, W=W, X=Sol$X, Model=Model, dimsol=dimsol,maxiter=maxiter, maxerror=maxerror, StandardizeDisparities=StandardizeDisparities, ShowIter=ShowIter)
    
    RawStress[i] <- Sol$stress
    stress1[i] <- Sol$stress1
    stress2[i] <- Sol$stress2
    sstress1[i] <- Sol$sstress1
    sstress2[i] <- Sol$sstress2
    POE <- Sol$X
    rownames(POE) <- rownames(D)
    if (ProcrustesRot)
      POE=SimpleProcrustes(PO,POE)$Yrot
    colnames(POE) <- paste("Dim", 1:dimsol)
    
    for (j in 1:n)
      for (k in 1:dimsol){
        Coordinates[[j]][k,i]=POE[j,k]}
  }
  Info=paste("Bootstrap for SMACOF based on ", method, "of distances")
  result=list(Info=Info, InitialDistance=D,RawStress=RawStress, stress1=stress1, stress2=stress2, sstress1=sstress1, sstress2=sstress2, Coordinates=Coordinates, NReplicates=nB)
  class(result)="PCoABootstrap"
  return(result)
}


