CanonicalBiplot <- function(X, group, SUP = NULL, InitialTransform = 5, LDA=FALSE, MANOVA = FALSE) {

  cl <- match.call()

  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering",
                              "Column centering", "Standardize columns", "Row centering",
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(InitialTransform))
    InitialTransform = ContinuousDataTransform[InitialTransform]

	Bip = list() #Container for the solution
	Bip$call=cl
	# Setting the properties of data
	if (is.null(rownames(X)))
		rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "I")
	RowNames = rownames(X)
	if (is.null(colnames(X)))
		colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
	VarNames = colnames(X)

	Bip$Title = "Canonical/MANOVA Biplot"
	Bip$Type = "Canonical"
	Bip$Non_Scaled_Data = X
	Bip$Means = apply(X, 2, mean)
	Bip$Medians = apply(X, 2, median)
	Bip$Deviations = apply(X, 2, sd)
	Bip$Minima = apply(X, 2, min)
	Bip$Maxima = apply(X, 2, max)
	Bip$P25 = apply(X, 2, quantile)[2, ]
	Bip$P75 = apply(X, 2, quantile)[4, ]
	Bip$GMean = mean(as.matrix(X))
  Bip$Initial_Transformation=InitialTransform
	X = TransformIni(as.matrix(X), transform = InitialTransform) # Initial transformation
	rownames(X) <- RowNames
	if (is.factor(group)) {
		GroupNames = levels(group)
	}
	g = length(levels(group))
	n = dim(X)[1]
	m = dim(X)[2]
	r = min(c(g - 1, m))
	Bip$ncols=m
  Bip$nrows=n
  Bip$dim=r
  if (LDA) {Bip$LDA=lda(X,group)
	Bip$Predict=predict(Bip$LDA,X)$class
	Bip$ClassificationTable = table(group, Bip$Predict)
	Bip$PercentCorrect=diag(prop.table(Bip$ClassificationTable, 1))
	names(Bip$PercentCorrect)=GroupNames
	Bip$TotalPercentCorrect=sum(diag(prop.table(Bip$ClassificationTable)))
	names(Bip$TotalPercentCorrect)= "Total"}
  
  
  if (MANOVA) Bip$MANOVA = manova(X ~ group)

	DimNames = "Dim 1"
	for (i in 2:r) DimNames = c(DimNames, paste("Dim", i))
	Z = Factor2Binary(group) # Matrix of indicators
	ng = colSums(Z)
	S11 = t(Z) %*% Z
	Xb = solve(S11) %*% t(Z) %*% X
	B = t(Xb) %*% S11 %*% Xb
	S = t(X) %*% X - B
	Y = (S11^0.5) %*% Xb %*% matrixsqrtinv(S)
	SV = svd(Y)

	H = matrixsqrt(S) %*% SV$v[, 1:r] # Variable coordinates
	B = matrixsqrtinv(S) %*% SV$v[, 1:r] # Canonical Weigths
	J = Xb %*% B # Center Coordinates
	V = X %*% B # Individual Coordinates
	if (!is.null(SUP)) {
		VS = SUP %*% B
		rownames(VS)=rownames(SUP)
		colnames(VS)=DimNames
		# Bip$SupPredict=predict(Bip$LDA,SUP)$class
	}
	else {
		VS=NULL
		Bip$SupPredict=NULL}

	# Inertia, ANOVAs for each Canonical Variate and MANOVA
	sct = diag(t(V) %*% V)
	sce = diag(t(J) %*% S11 %*% J)
	scr = sct - sce
	fs = (sce/(g - 1))/(scr/(n - g))
	signif2 = df(fs, (g - 1), (n - g))

	vprop = SV$d[1:r]
	iner = (vprop^2/sum(vprop^2)) * 100
	acum = cumsum(iner)

	Bip$EigenValues = vprop
	Bip$Inertia = iner
	Bip$CumInertia = acum
	# colnames(Bip$EigenValues) <- c("Eigenvalue", "Explained Variance", "Cummulative")
	# rownames(Bip$EigenValues) <- DimNames

	lambda = vprop^2
	pill = 1/(1 + lambda)
	pillai = det(diag(pill))
	glh = g - 1
	gle = n - g
	t = ((glh^2 * m^2 - 4)/(m^2 + glh^2 - 5))^0.5
	w = gle + glh - 0.5 * (m + glh + 1)
	df1 = m * glh
	df2 = w * t - 0.5 * (m * glh - 2)
	Bip$Wilksf = ((1 - pillai^(1/t))/(pillai^(1/t))) * (df2/df1)
	Bip$Wilksp = 1 - pf(Bip$Wilksf, df1, df2)

	Bip$GroupContributions = diag(1/rowSums(J^2)) %*% J^2
	Bip$ColContributions = diag(1/rowSums(H^2)) %*% H^2

	Bip$ExplTotal = matrix(0, r, 1)
	Bip$RowContributions = matrix(0, n, r)
	Bip$QLRVars = matrix(0, m, r)

	SCT = sum(X^2)
	SCRows = rowSums(X^2)
	SCCols = colSums(X^2)

	for (j in 1:r) {
		Fitted = V[, 1:j] %*% t(H[, 1:j])
		residuals = X - Fitted
		Bip$ExplTotal[j] = 1 - sum(residuals^2)/SCT
		Bip$RowContributions[, j] = 1 - rowSums(residuals^2)/SCRows
		Bip$QLRVars[, j] = 1 - colSums(residuals^2)/SCCols
	}


	FitX = V %*% t(H)
	Resid = X - FitX

	SCR = sum(Resid^2)
	FIT = 1 - (SCR/SCT)

	sctotal = diag(t(X) %*% X)

	scdentro = diag(S)
	scentre = sctotal - scdentro
	fs = (scentre/glh)/(scdentro/gle)
	pval = 1 - pf(fs, glh, gle)

	Bip$ANOVAS = cbind(sctotal, scentre, scdentro, fs, pval)
	colnames(Bip$ANOVAS) <- c("Total", "Groups", "Error", "F", "p-val")

	falfau = qt(1 - (0.025), (n - g))
	falfab = qt(1 - (0.025/(g * m)), (n - g))
	falfam = sqrt(qf(1 - 0.05, m, (n - g - m + 1)) * (((n - g) * m)/(n - g - m + 1)))
	falfac = sqrt(qchisq(0.95, 2))

	Bip$UnivRad = falfau * diag(solve(sqrt(S11)))/sqrt(n - g)
	Bip$BonfRad = falfab * diag(solve(sqrt(S11)))/sqrt(n - g)
	Bip$MultRad = falfam * diag(solve(sqrt(S11)))/sqrt(n - g)
	Bip$ChisRad = falfac * diag(solve(sqrt(S11)))/sqrt(n - g)

	Bip$n = n
	Bip$p = m
	Bip$g = g
	Bip$X = X
	Bip$groups = group

	Bip$RowCoordinates = V
	rownames(Bip$RowCoordinates) = RowNames
	colnames(Bip$RowCoordinates) = DimNames
	Bip$Sup_Individual_Coord = VS
	Bip$ColCoordinates = H
	rownames(Bip$ColCoordinates) = VarNames
	colnames(Bip$ColCoordinates) = DimNames
	Bip$GroupCoordinates = J
	rownames(Bip$GroupCoordinates) = GroupNames
	colnames(Bip$GroupCoordinates) = DimNames
	Bip$Canonical_Weights = B
	rownames(Bip$Canonical_Weights) = VarNames
	colnames(Bip$Canonical_Weights) = DimNames
	Bip$Structure_Correlations = cor(X, V)
	rownames(Bip$Structure_Correlations) = VarNames
	colnames(Bip$Structure_Correlations) = DimNames
	rownames(Bip$GroupContributions) = GroupNames
	colnames(Bip$GroupContributions) = DimNames
	rownames(Bip$ColContributions) = VarNames
	colnames(Bip$ColContributions) = DimNames
	rownames(Bip$QLRVars) = VarNames
	colnames(Bip$QLRVars) = DimNames
	
	  NGroups=length(levels(group))
	  Bip$Clusters = group
	  Bip$ClusterNames = levels(group)
	
	palette(rainbow(NGroups))
	ClusterColors = palette()
	Bip$ClusterType="us"
	Bip$ClusterColors=ClusterColors
	

	class(Bip) <- "Canonical.Biplot"
	return(Bip)
}
