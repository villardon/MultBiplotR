CumSum <- function(X, dimens = 1) {
	X = as.matrix(X)
	if (dimens == 2) 
		X = t(X)
	n = dim(X)[1]
	p = dim(X)[2]
	for (i in 1:n) X[i, ] = cumsum(X[i, ])
	if (dimens == 2) 
		X = t(X)
	return(X)
}