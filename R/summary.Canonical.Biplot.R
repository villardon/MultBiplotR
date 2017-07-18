summary.Canonical.Biplot <- function(object, ...) {
	cat(" ###### 1-WAY CANONICAL/MANOVA Biplot #######\n\n")
	cat("________________________________________________\n\n")

	if (!is.null(object$LDA)){
	cat("LINEAR DISCRIMINANT ANALYSIS SUMMARY\n")
	print(object$LDA)
	cat("\n Percent of Correct Classifications\n")
	print(c(object$PercentCorrect, object$TotalPercentCorrect))
	cat("________________________________________________\n\n")}

	if (!is.null(object$MANOVA)){
	cat("MULTIVARIATE ANALYSIS OF VARIANCE\n\n")
	cat("Wilks\n")
	pp=summary(object$MANOVA, test="Wilks")
	print(pp)
	cat("\n Pillai\n")
	pp=summary(object$MANOVA, test="Pillai")
	print(pp)
	cat("\n Hotelling-Lawley\n")
	pp=summary(object$MANOVA, test="Hotelling-Lawley")
	print(pp)
	cat("\n Roy\n")
	pp=summary(object$MANOVA, test="Roy")
	print(pp)
	cat("________________________________________________\n\n")}

	cat(paste("Global contrast based on Wilks Lambda :", object$Wilksf, "\n"))
	cat(paste("p-value :", object$Wilksp, "\n"))
	cat("________________________________________________\n\n")

	cat("Analysis of Variance for each Variable\n")
	print(object$ANOVAS)
	cat("________________________________________________\n\n")
	cat("Eigenvalues and explained variance\n")
	pp = cbind(1:length(object$EigenValues), object$EigenValues, object$Inertia, object$CumInertia)
	colnames(pp) = c("Axis", "Eigenvalue", "Explained Variance", "Cummulative")
	print(pp)
	cat("________________________________________________\n\n")

	cat("Correlations between the canonical axis and the original variables\n")
	print(round(object$Structure_Correlations, digits=3))
	cat("________________________________________________\n\n")
	cat("Squared Correlations \n")
	print(round(object$Structure_Correlations^2, digits = 3))
	cat("________________________________________________\n\n")
	cat("Quality of representation of the group means \n")
	print(round(object$GroupContributions*100, digits=3))
	cat("________________________________________________\n\n")
	cat("Quality of representation of the group means - Cummulative\n")
	print(round(CumSum(object$GroupContributions)*100, digits=3))
	cat("________________________________________________\n\n")
	cat("Goodness of fit if the variables (to explain the group means) \n")
	print(round(object$ColContributions*100, digits=3))
	cat("________________________________________________\n\n")
	cat("Goodness of fit if the variables (to explain the group means) - Cummulative \n")
	print(round(CumSum(object$ColContributions*100), digits=3))
	cat("________________________________________________\n\n")
	cat("Goodness of fit if the variables (to explain the individuasls) - Cummulative \n")
	print(round(object$QLRVars*100, digits=3))
	cat("________________________________________________\n\n")
	
	if (!is.null(object$ContSupVarsBiplot)){  cat("\nSUPPLEMENTARY VARIABLES ADDED\n")
	  cat("\nCorrelations  of the Supplementary Variables with the Canonical Axis\n")
	  print(round(object$ContSupVarsBiplot$Structure, digits = 3))
	  cat("\nCorrelations  of the Supplementary Variables with the Canonical Axis\n")
	  print(round(object$ContSupVarsBiplot$Structure^2, digits = 3))
	  cat("\nGoodness of fit of the Supplementary Variables (R2)\n")
	  print(round(object$ContSupVarsBiplot$R2, digits = 3))
	  
	  }

}
