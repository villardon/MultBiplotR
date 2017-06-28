ScreePlot <-  function(bip, Type="bar"){
  if (is.null(bip$EigenValues)) stop("Scree plot is not available for this object")
  plot(bip$EigenValues, xlab="Axis", ylab="Eigenvalues", main=paste(bip$Title, "- Scree Plot"))
}