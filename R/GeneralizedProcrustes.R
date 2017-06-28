
GeneralizedProcrustes <- function(x, tolerance = 1e-05, maxiter = 100, Plot = FALSE) {

	# Calculating the matrix sizes
	n = dim(x)[1] 
	p = dim(x)[2]
	g = dim(x)[3]


	# Centering the initial matrices
	for (i in 1:g) x[, , i] = (diag(n) - ones(n)/n) %*% x[, , i]

	#Total sum of squares
	sct = sum(x^2)
	w = matrix(0, g, 1)
	for (i in 1:g) w[i] = sum(x[, , i]^2)
	w = sqrt(sum(w)/w)


	# Initialization of parameters (Rotation -identity- and scaling -1-)
	T = array(0, c(p, p, g))
	s = matrix(1, g, 1)
	for (i in 1:g) T[, , i] = diag(p)

	# Consensus matrix
	Transformed = array(0, c(n, p, g))
	for (j in 1:g) Transformed[, , j] = s[j] * x[, , j] %*% T[, , j]
	G = apply(Transformed, c(1, 2), sum)/g


	# Residual sum of squares
	Residuals = x
	rsq = 0
	for (j in 1:g) Residuals[, , j] = Transformed[, , j] - G
	rsq = sum(Residuals^2)

	rsqold = rsq

	# Iterative algorithm
	error = 1
	iter = 0
	history=NULL
	while ((error > tolerance) & (iter < maxiter)) {
		iter = iter + 1

		#This update is made as in Gower (1975)
		scg = sum(G^2)
		for (i in 1:g) s[i, 1] = (sct/(g * scg)) * (tr(t(G) %*% x[, , i] %*% T[, , i])/sum((x[, , i] %*% T[, , i])^2))

		for (j in 1:g) Transformed[, , j] = s[j] * x[, , j] %*% T[, , j]

		# Updating Rotations
		for (i in 1:g) {

			int = t(x[, , i]) %*% G
			udv = svd(int)
			T[, , i] = udv$v %*% t(udv$u)

			scg = sum(G^2)
			s[i, 1] = (sct/(g * scg)) * (tr(t(G) %*% x[, , i] %*% T[, , i])/sum((x[, , i] %*% T[, , i])^2))
			for (j in 1:g) Transformed[, , j] = s[j] * x[, , j] %*% T[, , j]
			G = apply(Transformed, c(1, 2), sum)/g

		}


		# Residual sum of squares
		for (j in 1:g) Residuals[, , j] = Transformed[, , j] - G
		rsq = sum(Residuals^2)

		error = ((rsqold - rsq)^2)/rsq^2
		history=rbind(history,c(iter, rsq, error))
		print(history[iter,])
		rsqold = rsq

	}

Proc=list()
Proc$History=history
Proc$X=x
Proc$RotatedX=Transformed
Proc$Scale=s
Proc$Rotations=T
Proc$Consensus=G

	if (Plot)  { dev.new()
		plot(G[, 1], G[, 2])
for (i in 1:g) plines(G[, 1:2], Transformed[, 1:2, i])}
	
	sct = sum(Transformed^2)
	scg = g * sum(G^2)
	gfit = 1 - (rsq/sct)
Proc$Fit=gfit
class(Proc)="GenProcustes"
return(Proc)
}

