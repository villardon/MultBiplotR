# Takes an object from a PCoA for mexed types of data and constructs a biplot

MixedBiplot <- function(pco, tolerance = 1e-05, maxiterlogist = 100, penalization = 0.2, showiter = FALSE, IncludeConst=TRUE){
  
  Numericas=which(pco$Types=="numeric")
  if (length(Numericas)>0){
    X=pco$Data[, Numericas]
    n=dim(X)[1]
    p=dim(X)[2]
    pco=AddContVars2Biplot(pco,  X, dims=NULL, Scaling = 5, Fit=NULL)
    pco=c(pco, pco$ContSupVarsBiplot)
    pco$ContSupVarsBiplot=NULL
    pco$ColContributions = (pco$Structure)^2
    pco$RowContributions=pco$RowQualities
    
    sca = sum(pco$RowCoordinates^2)
    scb = sum(pco$ColCoordinates^2)
    sca = sca/n
    scb = scb/p
    scf = sqrt(sqrt(scb/sca))
    pco$RowCoordinates = pco$RowCoordinates * scf
    pco$ColCoordinates = pco$ColCoordinates/scf
    class(pco)="ContinuousBiplot"
    type=TRUE
  }
  
  Binarias=which(pco$Types=="binary")
  if (length(Binarias)>0){
    Y=as.matrix(pco$Data[, Binarias])
    pco=AddBinVars2Biplot(pco, Y, IncludeConst=IncludeConst, penalization=penalization, freq=NULL, tolerance = tolerance, maxiter = maxiterlogist)
  }
  
  Ordenadas=which(pco$Types=="ordered")
  if (length(Ordenadas)>0){
    Y=as.data.frame(pco$Data[, Ordenadas])
    colnames(Y)=colnames(pco$Data)[Ordenadas]
    pco=AddOrdVars2Biplot(pco, Y, tol = tolerance, maxiterlogist = maxiterlogist, penalization = penalization, showiter = showiter)
  }

  
  Factores=which(pco$Types=="factor")
  if (length(Factores)>0){
  Y=as.data.frame(pco$Data[, Factores])
  colnames(Y)=colnames(pco$Data)[Factores]
  class(model) = "NomSupVarsBiplot"
  
  }
  
  
return(pco)
}

# "TypeData"       "Type"           "Transformation" "Data"           "Proximities"    "Analysis"       "EigenValues"    "Inertia"        "RowCoordinates" "RowQualities"   "RawStress"      "stress1"       
# "stress2"        "sstress1"       "sstress2"       "rsq"            "rho"            "tau"     