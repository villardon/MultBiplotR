print.ContinuousBiplot <- function(x, ...){
  if (x$Type=="PCA")
    cat(" ###### Biplot for Principal Components Analysis #######\n")
  if (x$Type=="FA")
    cat(" ###### Biplot for Factor Analysis #######\n\n")
  if (x$Type=="HJ"){
    cat(" ###### HJ - Biplot #######\n\n")
  }
  cat("Call\n")
  print(x$call)
  cat("Transformation of the raw data:\n")
  print(x$Initial_Transformation)
  cat("Type of coordinates:\n")  
  if (x$alpha==2)
    tipo="Principal Normalization (Baricentric Scaling)"
  if (x$alpha==1)
  tipo="Row Principal Normalization (RMP-Biplot)"
  if (x$alpha==0)
    tipo="Column Principal Normalization (CMP-Biplot)"
  if (x$alpha==0.5)
    tipo="Symmetrical Normalization (SQRT - Biplot)"
  if ((x$alpha>0) & (x$alpha<1) & (x$alpha != 0.5))
    tipo=paste("Custom Normalization (Biplot con \alpha = ",gamma,")")
  print(tipo)
  cat("\n Eigenvalues & Explained Variance (Inertia)\n")
  
  if (x$Type!="LogFreqBiplot"){
  pp=cbind(x$EigenValues[1:x$Dimension], x$Inertia[1:x$Dimension], x$CumInertia[1:x$Dimension])
  colnames(pp)=c("Eigenvalue", "Exp. Var", "Cummulative")
  print(pp)
  }
  
}