Multiquad <- function(nnodes, dims) {
	Q = Hermquad(nnodes)
	I = patterns_eq(nnodes, dims)
	n = dim(I)[1]
	X2 = matrix(0, n, dims)
	A2 = matrix(1, n, 1)
	for (i in 1:n) for (j in 1:dims) {
		X2[i, j] = Q$X[I[i, j]]
		A2[i, 1] = A2[i, 1] * Q$W[I[i, j]]
	}
	Max = Q$W[1] * Q$W[round((nnodes + 1)/2)]/15
	keep = (A2 > Max)
	n2 = sum(keep)
	X = matrix(0, n2, dims)
	A = matrix(1, n2, 1)
	k = 0
	for (i in 1:n) if (keep[i]) {
		k = k + 1
		X[k, ] = X2[i, ]
		A[k, ] = A2[i, ]
	}
	QUAD = list()
	QUAD$X = X
	QUAD$A = A
	class(QUAD) <- "MultiGaussQuadrature"
	return(QUAD)
}
