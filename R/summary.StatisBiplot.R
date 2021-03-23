summary.StatisBiplot <- function(object, ...) {
  
  cat("Call\n")
  print(object$call)
  cat("Correlations among occasions\n")
  print(object$RV)
  
  cat("Biplot for Statis\n")
  summary(object$Biplot)
  
}




print.StatisBiplot <- function(x, ...) {
  
  cat("Call\n")
  print(x$call)
  cat("Correlations among occasions\n")
  print(x$RV)
  cat("Biplot for Statis\n")
  print(x$Biplot)
}