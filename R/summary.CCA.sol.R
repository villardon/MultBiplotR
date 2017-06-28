summary.CCA.sol <- function(object, ...) {
	dims=object$dimens
	cat(" \n\n###### CANONICAL CORRESPONDENCES ANALYSIS  #######\n\n")
	cat("________________________________________________\n\n")
cat("Explained Variance\n")
Exp=round(rbind(object$EigenValues[1:dims], object$Exp_from_relation[1:dims], object$Exp_from_species_const[1:dims], object$Exp_from_species_non_const[1:dims]), digits=4)
colnames(Exp)=colnames(object$Site_Scores)[1:dims]
rownames(Exp)=c("Eigenvalues", "Of Species-Environment Relation", "Of Species Data (constrained)", "Of Species Data (non-constrained)")
print(Exp)
cat("\nCummulative Explained Variance\n")
print(t(apply(Exp,1,cumsum))[2:4,])
}
