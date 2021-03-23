summary.PCA.Analysis <-function(object, latex=FALSE, ...){
  cat(" ###### Principal Components Analysis #######\n\n")
  cat("Transformation of the raw data:\n")
  print(object$Initial_Transformation)
  cat("\n Eigenvalues & Explained Variance (Inertia)\n")
  
  pp=cbind(object$EigenValues[1:object$Dimension], object$Inertia[1:object$Dimension], object$CumInertia[1:object$Dimension])
  
  colnames(pp)=c("Eigenvalue", "Exp. Var", "Cummulative")
  print(pp)
  
  if (object$Type=="FA"){
    cat(" Factor Analysis Loadings\n\n")
    print(object$ColCoordinates)
  }
  
  cat("\n\n STRUCTURE OF THE PRINCIPAL COMPONENTS\n")
  print(round(object$Structure, digits=3))
  
  if (latex){
    print(xtable(pp, caption="Explained Variance"))
    print(xtable(round(object$Structure, digits=3), caption="Correlations with the Principal Components"))
  }
  
}