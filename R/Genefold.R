Genefold <- function(DATA, TRANSFORMATION = "ordinal", DIMENS = 2, LENGTHCONSTR = TRUE, LENGTHX = 1, VARCONSTR = TRUE, INITIAL = "rational", 
	MAXITER = 100, CONVERGENCE = 1e-05, HISTORY = FALSE, Constrained=FALSE) {
	n = dim(DATA)[1]
	m = dim(DATA)[2]
	eps = 1e-12

	if (VARCONSTR) 
		Umax = 1/sqrt(m * VarDisk(m))
	else Umax = 1/eps

	SOL = Initial(DATA, TRANSFORMATION, INITIAL, DIMENS, LENGTHX, Umax)
	SOL$Data = DATA
	SOL$D = DistUnfold(SOL$X, SOL$Y)
	L = Loss(SOL)

	conv = 0
	iter = 0
	if (HISTORY) 
		print(paste("Iteration: ", iter, "-   Loss: ", L))
	while ((conv == 0) & (MAXITER > 0)) {
		Linit = L
		SOL = UpdateX(SOL, LENGTHX, LENGTHCONSTR, eps)
		SOL$D = DistUnfold(SOL$X, SOL$Y)
		SOL = UpdateY(SOL)
		SOL$Y = SOL$Y - (matrix(1, m, m)/m) %*% SOL$Y
		SOL$D = DistUnfold(SOL$X, SOL$Y)

		if (TRANSFORMATION == "ordinal") 
			for (i in 1:n) SOL$GAMMA[i, ] = UpdateGAMMA(DATA[i, ], SOL$a[i] * SOL$D[i, ])

		SOL$a = Updatea(SOL$GAMMA, SOL$D, Umax)
		L = Loss(SOL)
		iter = iter + 1
		Diff = abs(Linit - L)

		if (HISTORY) 
			print(paste("Iteration: ", iter, "-   Loss: ", L, "   Diff: ", Diff))

		if (iter == MAXITER) 
			conv = 1
		if (Diff < CONVERGENCE) 
			conv = 1
	}
	
	if (Constrained) 
	  SOL$Analysis = "Constrained Unfolding - Genefold"
	else SOL$Analysis  = "Unfolding - Genefold"

	# Some goodness of fit indices and results
	C1IP = 0
	Rp = 0
	CiIP = matrix(0, n, 1)
	Rpi = matrix(0, n, 1)
	for (i in 1:n) {
		CiIP[i] = indexC1(DATA[i, ], SOL$D[i, ])
		C1IP = C1IP + CiIP[i]
		Rpi[i] = cor(DATA[i, ], SOL$D[i, ])
		Rp = Rp + Rpi[i]
	}

	vm = matrix(1, m, 1)
	vn = matrix(1, 1, n)
	v = matrix(1, m, m)
	Dm = SOL$D %*% v/m
	Dc = SOL$D - Dm

	if (TRANSFORMATION == "linear") {
		DATAm = DATA %*% v/m
		DATAc = DATA - DATAm
		a = ((DATAc * Dc) %*% vm)/((DATAc * DATAc) %*% vm)
		A = a %*% matrix(1, 1, m)
		SOL$GAMMA = DATA * A + (Dm - (A * DATAm))
		SOL$a = a
	} else for (i in 1:n) SOL$GAMMA[i, ] = lsqisotonic(DATA[i, ], SOL$D[i, ], vm)

	Dev = SOL$GAMMA - SOL$D
	S1 = ((Dev * Dev) %*% vm)/((SOL$D * SOL$D) %*% vm)
	S1 = sqrt(vn %*% S1/n)
	S2 = ((Dev * Dev) %*% vm)/((Dc * Dc) %*% vm)
	S2 = sqrt(vn %*% S2/n)
	C1IP = C1IP/n
	Rp = Rp/n

	SOL$Iterations = iter
    SOL$RecoveredOrders= CiIP
    SOL$MeanRecoveredOrders= C1IP
    SOL$Correlations=Rpi
    SOL$MeanCorrelations=Rp
    SOL$Loss=L
    SOL$Stress1=S1
    SOL$Stress2=S2
    rownames(SOL$Y)=colnames(DATA)
	class(SOL) <- "Genefold"
	return(SOL)
}




