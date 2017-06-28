wa <- function(CCA.sol, transformed=FALSE){
	if (transformed) print(CCA.sol$Transformed_Weighted_Averages) else print(CCA.sol$Weighted_Averages)
}
