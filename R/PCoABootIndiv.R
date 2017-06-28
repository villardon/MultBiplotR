PCoABootIndiv <- function(PcSol,nrep=200,dimension=2){
  dist=PcSol$Proximities
  n=dim(PcSol$Proximities)[1]
  BootEigenVal=matrix(0,nrep,n)
  BootInertia=matrix(0,nrep,n)
  Config=array(0,c(n,dimension,nrep))
  for (k in 1:nrep){
    sample = sort(ceiling(n*runif(n)))
    common=intersect(1:n,sample)
    differ=setdiff(1:n,sample)
    Dist=list()
    Dist$Data=PcSol$Data[sample,]
    Dist$SupData=PcSol$Data[differ,]
    Dist$Proximities= PcSol$Proximities[sample,sample]
    Dist$SupProximities=PcSol$Proximities[differ,sample]
    class(Dist) <- "proximities"
    PC=PrincipalCoordinates(Dist)
    BootEigenVal[k,]=PC$EigenValues
    BootInertia[k,]=PC$Inertia
    conf=matrix(0,n,dimension)
    for (i in 1:n){
      conf[sample[i],]=PC$RowCoordinates[i,]
    }
    conf[differ,]=PC$SupRowCoordinates
    Config[,,k]=conf
  }
  PcSol$BootstrapEigen=BootEigenVal
  PcSol$BootstrapInertia=BootInertia
  PcSol$BootstrapConfigs=Config
  PcSol$BootstrapProcessing="None"
  return(PcSol)
}

                       