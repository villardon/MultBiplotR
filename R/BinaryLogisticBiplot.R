BinaryLogisticBiplot <- function(x, dim = 2, compress = FALSE, init = "mca", 
                                 method = "EM", rotation = "none", 
                                 tol = 1e-04, maxiter = 100, penalization = 0.2, 
                                 similarity="Simple_Matching", ...) {
  ntot = dim(x)[1]
  
  if (is.null(rownames(x)))
    rownames(x) <- paste("R",1:ntot, sep="")
  
  if (compress) {
    NewTable <- ExtractTable(x)
    x = NewTable$Patterns
    freq = NewTable$Frequencies
  } else {
    freq = matrix(1, ntot, 1)
  }
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  if (init=="random"){
    a <- matrix(rnorm(n*dim), n, dim )
  }
  
  if (init=="mirt"){
    mod2 <- mirt(NewTable$Data,dim)
    a <- fscores(mod2,method = "EAP", full.scores = TRUE)
    if (compress) a=a[NewTable$Unique,]
  }
  
  if (init=="PCoA"){
    dis= BinaryProximities(x, coefficient=similarity)
    result=PrincipalCoordinates(dis, w = freq)
    a = result$RowCoordinates
  }
  
  if (init=="mca"){
    corr=CA(as.matrix(x),dim=dim)
    a=corr$RowCoordinates[,1:dim]
  }
  
  print("Fitting the model") 
  
  switch(method, EM = {
    LogBip = BinaryLogBiplotEM(x, freq, dimens = dim, aini=a,tol = tol, maxiter = maxiter, penalization = penalization, ...)
    if (compress) {LogBip$RowCoordinates = ExpandCoord(LogBip$RowCoordinates, NewTable)
    }
  }, Joint = {
    LogBip = BinaryLogBiplotJoint(x, freq, dimens = dim, ainit=a, tolerance = tol, maxiter = maxiter, penalization = penalization, ...)
  }, mirt = {
    LogBip = BinaryLogBiplotJoint(x, freq, dimens = dim, tolerance = tol, maxiter = maxiter, penalization = penalization, ...)
  }, External = {
    LogBip = BinaryLogBiplotJoint(x, freq, dimens = dim, tolerance = tol, maxiter = maxiter, penalization = penalization, ...)
  }, JointGD = {
    LogBip = BinaryLogBiplotGD(x, freq, dimens = dim, tolerance = tol, maxiter = maxiter, penalization = penalization, ...)
  },AlternatedGD = {
    LogBip = BinaryLogBiplotGD(x, freq, dimens = dim, tolerance = tol, maxiter = maxiter, penalization = penalization, ...)
  })
  
  LogBip$Type="Binary Logistic Biplot"
  LogBip$InitialConfig=init
  LogBip$Method=method
  LogBip$Rotation=rotation
  class(LogBip)="Binary.Logistic.Biplot"
  
  return(LogBip)
}



# method = c("Joint", "EM", "mirt", "external")
# rotation = c("none", "oblimin", "quartimin", "oblimax" ,"entropy",  "quartimax", "varimax",  "simplimax" ) see GPARotation
# init = c("mca", "random", "PCoA", "mirt")

