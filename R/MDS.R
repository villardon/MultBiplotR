MDS <- function(Proximities, W=NULL, Model=c("Identity", "Ratio", "Interval", "Ordinal"), dimsol=2, maxiter=100, maxerror=0.000001, Bootstrap=FALSE,
               nB=200, ProcrustesRot=TRUE, BootstrapMethod=c("Sampling", "Permutation"), StandardizeDisparities=FALSE, ShowIter=FALSE){
  if (length(Model)>1) Model=Model[1]
  if (length(BootstrapMethod) > 1) BootstrapMethod = BootstrapMethod[1] 
  if (!(class(Proximities)=="proximities")) stop("You need a proximities matrix")
  Dimnames=paste("Dim", 1:dimsol)
  P=Proximities$Proximities
  n <- dim(P)[1]
  if (is.null(W)) W = matrix(1,n,n)
  sol=SMACOF(P, W=W, Model=Model, dimsol=dimsol, maxiter=maxiter, maxerror=maxerror, StandardizeDisparities=StandardizeDisparities, ShowIter=ShowIter)
  Proximities$Analysis="MDS"
  Proximities$Model=Model
  rownames(sol$X)=rownames(P)
  colnames(sol$X)=Dimnames
  Proximities$RowCoordinates=sol$X
  Proximities$RowQualities=NULL
  Proximities$Disparities=sol$Dh
  Proximities$RawStress=sol$stress
  Proximities$stress1=sol$stress1
  Proximities$stress2=sol$stress2
  Proximities$sstress1=sol$sstress1
  Proximities$sstress2=sol$sstress2
  Proximities$rsq=sol$rsq
  Proximities$rho=sol$rho
  Proximities$tau=sol$tau
  Proximities$Transformation=Model
  if (Bootstrap){
    Proximities$BootstrapInfo=BootstrapSmacof(P, W=W, nB=nB, dimsol=dimsol, Model=Model, maxiter=maxiter,  StandardizeDisparities=StandardizeDisparities, 
                                              maxerror=maxerror, ProcrustesRot=ProcrustesRot, method=BootstrapMethod, ShowIter=ShowIter)
  }
  class(Proximities) <- c("Principal.Coordinates","MDS")
  return(Proximities)
}