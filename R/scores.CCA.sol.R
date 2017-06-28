scores.CCA.sol <- function(CCA.sol){
	dims=CCA.sol$dimens
	scores=list()
	scores$sites=CCA.sol$Site_Scores[,1:dims]
	scores$species=CCA.sol$Species_Scores[,1:dims]
	scores$variables=CCA.sol$Env_Var_Scores[,1:dims]
	return(scores)
}