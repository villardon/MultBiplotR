CCA = function(P, Z, alpha = 1, dimens = 4) {
	# P abundance matrix, Z environmental variables
	# alpha : biplot decomposition,  
    # with alpha=1 the emphasis is on the sites
    # with alpha=0 the emphasis is on the species


	if (is.data.frame(P)) 
		P = as.matrix(P)
	if (is.data.frame(Z)) 
		Z = as.matrix(Z)

	# Setting the properties of data
	if (is.null(rownames(P))) 
		rownames(P) <- rownames(P, do.NULL = FALSE, prefix = "Si")
	if (is.null(colnames(P))) 
		colnames(P) <- colnames(P, do.NULL = FALSE, prefix = "Sp")

	SpeciesNames = colnames(P)
	SiteNames = rownames(P)

	if (is.null(rownames(Z))) 
		rownames(Z) <- SiteNames
	if (is.null(colnames(Z))) 
		colnames(Z) <- colnames(Z, do.NULL = FALSE, prefix = "Ev")
	EnvNames = colnames(Z)


	I = nrow(P)
	J = ncol(P)
	K = ncol(Z)
	
	CCA.sol = list()
	CCA.sol$call <- match.call()
	CCA.sol$nsites = I
	CCA.sol$nspecies = J
	CCA.sol$nvars = K
	CCA.sol$alpha = alpha
	CCA.sol$dimens = dimens

	N = sum(P)
	P = P/N

	r = apply(P, 1, sum)
	c = apply(P, 2, sum)
	PC = P - r %*% t(c)

	PC = sqrt(diag(1/r)) %*% PC %*% sqrt(diag(1/c))
	valsin = svd(PC)$d

	InTotal = sum(diag(t(PC) %*% PC))
	InerRows = apply(PC^2, 1, sum)
	InerCols = apply(PC^2, 2, sum)

	#Initial Weighted averages (before transformation of abundances)
	W1 = diag(1/c) %*% t(P) %*% Z
	rownames(W1) = SpeciesNames
    # Useful Parameters of the environmental variables
     CCA.sol$Minima_Z = apply(Z,2,min)
     CCA.sol$Maxima_Z = apply(Z,2,max)
     CCA.sol$Medians_Z = apply(Z,2,weighted.quantile, r,0.5 )
     CCA.sol$P25_Z = apply(Z,2,weighted.quantile, r,0.25 )
     CCA.sol$P75_Z = apply(Z,2,weighted.quantile, r,0.75 )
     
	# Centering and standardizing the environmental variables
	waZ = t(t(r) %*% Z)
	Z = (diag(I) - matrix(1, I, 1) %*% t(r)) %*% Z
	s11 = t(Z^2) %*% r
	Z = Z %*% diag(sqrt(as.vector(1/s11)))	
	rownames(Z) = SiteNames
	colnames(Z) = EnvNames
	# Initial Weighted averages (after transformation of abundances)
	W = diag(1/c) %*% t(P) %*% Z
	rownames(W) = SpeciesNames
	colnames(W) = EnvNames
	
    CCA.sol$Means_Z = waZ
	CCA.sol$Deviations_Z = sqrt(s11)

	S22 = t(Z) %*% diag(r) %*% Z
	sqrts22 = matrixsqrt(S22)
	insqrts22 = matrixsqrtinv(S22)

	WW = sqrt(diag(c)) %*% W %*% insqrts22
	dvs = svd(WW)

	rang = length(dvs$d)
	DimNames = "CCA1"
	for (i in 2:rang) DimNames = c(DimNames, paste("CCA", i, sep=""))

	# Principal coordinates for species and variables 
	A = diag(sqrt(1/c)) %*% dvs$u %*% diag(dvs$d)
	B = sqrts22 %*% dvs$v %*% diag(dvs$d)

	# Goodness of fit o f species and environmental variables
	# to explain the weighted means

	clresp = round((diag(1/apply(A^2, 1, sum)) %*% A^2) * 100, digits = 2)
	rownames(clresp) = SpeciesNames
	colnames(clresp) = DimNames
	clrvar = round((diag(1/apply(B^2, 1, sum)) %*% B^2) * 100, digits = 2)
	rownames(clrvar) = EnvNames
	colnames(clrvar) = DimNames

	clrespAFC = round((diag(1/InerCols) %*% (dvs$u %*% diag(dvs$d))^2) * 100, digits = 2)
	rownames(clrespAFC) = SpeciesNames
	colnames(clrespAFC) = DimNames

	# Explained Variance 
	# From the relation among species and environment
	varexp = (dvs$d^2/sum(dvs$d^2)) * 100
	names(varexp)=DimNames
	# Of the species (constrained)
	varexp2 = (dvs$d^2/InTotal) * 100
	names(varexp2)=DimNames
	# Of the species (non constrained)
	varexp3 = round((valsin^2/InTotal) * 100, digits = 2)
	names(varexp3)=DimNames

	A = diag(sqrt(1/c)) %*% dvs$u %*% diag(dvs$d)

	# Principal coordinates of the individuals as weighted averages
	x = diag(1/r) %*% P %*% diag(sqrt(1/c)) %*% dvs$u

	clrsitesAFC = round((diag(1/InerRows) %*% (diag(sqrt(r)) %*% x)^2) * 100, digits = 2)
	rownames(clrsitesAFC) = SiteNames
	colnames(clrsitesAFC) = DimNames
	clrsitesAFCacum = t(apply(clrsitesAFC, 1, cumsum))
	rownames(clrsitesAFCacum) = SiteNames
	colnames(clrsitesAFCacum) = DimNames

	# Principal coordinates of the individuals as linear combinations 
	x2 = Z %*% insqrts22 %*% dvs$v

	# Canonical coefficients
	cocan = insqrts22 %*% t(Z) %*% diag(r) %*% x
	interset = cor(x, Z)

	# Final cooordinates 
	A = diag(sqrt(1/c)) %*% dvs$u %*% diag(dvs$d^alpha)
	rownames(A) = SpeciesNames
	colnames(A) = DimNames
	B = sqrts22 %*% dvs$v %*% diag(dvs$d^(1 - alpha))
	rownames(B) = EnvNames
	colnames(B) = DimNames
	x = x %*% diag((dvs$d^(-1 * alpha)))
	rownames(x) = SiteNames
	colnames(x) = DimNames
	x2 = x2 %*% diag((dvs$d^(-1 * alpha)))
	rownames(x2) = SiteNames
	colnames(x2) = DimNames

	# Solution 
	

	CCA.sol$Weighted_Averages = W1
	CCA.sol$Transformed_Weighted_Averages = W
	CCA.sol$EigenValues=dvs$d^2
	CCA.sol$Site_Scores = x[,1:dimens]
	CCA.sol$LC_Site_Scores = x2[,1:dimens]
	CCA.sol$Species_Scores = A[,1:dimens]
	CCA.sol$Env_Var_Scores = B[,1:dimens]
	CCA.sol$Exp_from_relation= varexp[1:dimens]
	CCA.sol$Exp_from_species_const= varexp2[1:dimens]
	CCA.sol$Exp_from_species_non_const= varexp3[1:dimens]
	CCA.sol$Species_Quality=clresp[,1:dimens]
	CCA.sol$Species_Quality_AFC=clrespAFC[,1:dimens]
	CCA.sol$Variables_Quality= clrvar[,1:dimens]
    CCA.sol$Sites_Quality= clrsitesAFC[,1:dimens]
	class(CCA.sol)<- "CCA.sol"
	return(CCA.sol)
}
