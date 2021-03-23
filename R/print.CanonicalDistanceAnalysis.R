print.CanonicalDistanceAnalysis <- function(x, ...){
  cat("Canonical Analysis of Distances\n")
  cat("Call\n")
  print(x$call)
  res=matrix(0,3,4)
  cat("\nAnalysis of Variance\n")
  res=c(x$TSS,  x$BSS, x$WSS, x$Fexp, x$pvalue)
  names(res)=c("Total SS", "Between SS", "Within SS", "F Exp.", "p-value")
  print(res)
  cat("\n\n")
}