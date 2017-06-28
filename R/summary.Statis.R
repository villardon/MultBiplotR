summary.Statis <- function(object, ...) {
  
  cat("Call\n")
  print(object$call)
  cat("Correlations among occasions\n")
  print(object$RV)
  
  cat("Biplot for Statis\n")
  summary(object$Biplot)
  
}




print.Statis <- function(object, ...) {
  
  cat("Call\n")
  print(object$call)
  cat("Correlations among occasions\n")
  print(object$RV)
  cat("Biplot for Statis\n")
  print(object$Biplot)
}