print.Canonical.Biplot <- function(object, ...) {
  cat(" ###### 1-WAY CANONICAL/MANOVA Biplot #######\n\n")
  cat("Call\n")
  print(object$call)
  cat("________________________________________________\n\n")
  cat("Eigenvalues and explained variance\n")
  pp = cbind(1:length(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
  colnames(pp) = c("Axis", "Eigenvalue", "Explained Variance", "Cummulative")
  print(pp)
  cat("________________________________________________\n\n")
  cat(paste("Global contrast based on Wilks Lambda :", object$Wilksf, "\n"))
  cat(paste("p-value :", object$Wilksp, "\n"))
  
  if (!is.null(object$ContSupVarsBiplot))  cat("\nSome supplementary variables have been added (vector model)\n")
}