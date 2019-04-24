
OrdinalLogisticFit <- function(y, x = NULL, penalization = 0.1, tol = 1e-04, maxiter = 200, show = FALSE) {
	n <- length(y)
	if (is.null(x)) {
		p = 0
	} else {
		if (!is.matrix(x)) {
			x = as.matrix(x)
		}
		p <- ncol(x)
	}
	J = max(y)
	npar = J + p - 1
	Y = matrix(0, n, J)
	for (i in 1:n) if (y[i] > 0) 
		Y[i, y[i]] = 1
	R = matrix(0, n, J)
	for (i in 1:n) {
		R[i, ] = cumsum(Y[i, ])
	}
	#Este if es mio porque vi que cuando la variable independiente tomaba solo dos categorias\n
	#el programa daba error, y por supuesto cuando tomaba solo una tambien. En este caso\n
#mostrar un mensaje de eliminar la variable\n
if (J > 2) {
		A = matrix(colSums(R[, 1:J - 1])/n, J - 1, 1)
	} else if (J == 2) {
		A = matrix(sum(R[, 1:J - 1])/n, J - 1, 1)
	} else {
		stop("There is a variable with the same value for all the items. Revise the data set.")
	}

	Beta = matrix(0, p, 1)
	B = rbind(A, Beta)
	err = 0.1
	iter = 0
	while ((err > tol) & (iter < maxiter)) {
		iter = iter + 1

		if (is.null(x)) 
			ETA = matrix(1, n, 1) %*% t(A)
		else ETA = matrix(1, n, 1) %*% t(A) - x %*% Beta %*% matrix(1, 1, (J - 1))

		PIA = exp(ETA)/(1 + exp(ETA))
		PIA = cbind(PIA, matrix(1, n, 1))
		PI = matrix(0, n, J)
		PI[, 1] = PIA[, 1]
		PI[, 2:J] = PIA[, 2:J] - PIA[, 1:(J - 1)]
		Rho = log(PIA[, 1:(J - 1)]/(PIA[, 2:J] - PIA[, 1:(J - 1)]))
		gRho = log(1 + exp(Rho))
		U = PIA[, 2:J]/(PIA[, 1:(J - 1)] * (PIA[, 2:J] - PIA[, 1:(J - 1)]))
		D = R[, 1:(J - 1)] - R[, 2:J] * (PIA[, 1:(J - 1)]/PIA[, 2:J])
		D2 = PIA[, 1:J] * (1 - PIA[, 1:J])
		P = array(0, c(n, J, npar))
		Q = array(0, c(n, J - 1, npar))
		GR = matrix(0, npar, 1)
		for (k in 1:npar) {
			if (k < J) 
				P[, k, k] = matrix(1, n, 1)
			else P[, 1:(J - 1), k] = (-1 * x[, k - (J - 1)]) %*% matrix(1, 1, (J - 1))
			Q[, , k] = P[, 1:(J - 1), k] * D2[, 1:(J - 1)] - P[, 2:J, k] * (PIA[, 1:(J - 1)]/PIA[, 2:J]) * D2[, 2:J]
			GR[k] = sum(D * U * Q[, , k])
		}
		# Expectation of the second derivative, (Fisher's scoring method is used)\n
		HESS = matrix(0, npar, npar)
		for (s in 1:npar) for (k in 1:npar) {
			HESS[s, k] = sum(U * Q[, , s] * Q[, , k])
		}
		HESS = HESS + 2 * diag(J + p - 1) * penalization
		Bnew = B + solve(HESS) %*% (GR - 2 * B * penalization)
		err = sum((Bnew - B)^2)
		B = Bnew
		A = matrix(B[1:(J - 1)], J - 1, 1)
		if (is.null(x))
		Beta = matrix(0, p, 1)
		else
		Beta = matrix(B[J:(J + p - 1)], p, 1)
		
		if (show) 
			print(c(iter, err))
	}

	L = sum(R[, 1:(J - 1)] * Rho - R[, 2:J] * gRho)

	Deviance = -2 * sum(L)

	model <- list()
	model$nobs = n
	model$J = J
	model$nvar = p
	model$Eta=ETA
	model$fitted.values = PI
	model$pred = matrix(max.col(PI), n, 1)
	model$Covariaces = solve(HESS)
	model$clasif = table(y, model$pred)
	model$PercentClasif = sum(y == model$pred)/n
	model$coefficients = Beta
	model$thresholds = A
	model$logLik = L
	model$penalization = penalization
	model$Deviance = Deviance
	model$iter = iter
	return(model)
}

