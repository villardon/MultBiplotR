print.Canonical.Biplot <- function(x, ...) {
  cat(" ###### 1-WAY CANONICAL/MANOVA Biplot #######\n\n")
  cat("Call\n")
  print(x$call)
  cat("________________________________________________\n\n")
  cat("Eigenvalues and explained variance\n")
  pp = cbind(1:length(x$EigenValues), x$EigenValues, x$Inertia, x$CumInertia)
  colnames(pp) = c("Axis", "Eigenvalue", "Explained Variance", "Cummulative")
  print(pp)
  cat("________________________________________________\n\n")
  cat(paste("Global contrast based on Wilks Lambda :", x$Wilksf, "\n"))
  cat(paste("p-value :", x$Wilksp, "\n"))
  
  if (!is.null(x$ContSupVarsBiplot))  cat("\nSome supplementary variables have been added (vector model)\n")
}