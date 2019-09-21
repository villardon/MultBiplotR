SimpleProcrustes <- function(X, Y, centre=FALSE) {
	# Z= s *
	nx = nrow(X)
	px = ncol(X)
	ny = nrow(Y)
	py = ncol(Y)

	if (nx != ny) 
		stop("Matrices should have the same number of rows: ", nrow(X), " or ", nrow(Y))

	if (px != py) {
		warning("Matrices should have the same number of columns: Adjusting to conform")
		add = abs(px - py)
		if (px < py) 
			X = cbind(X, matrix(0, nx, add))
		if (px > py) 
			Y = cbind(Y, matrix(0, ny, add))
		px = ncol(X)
		py = ncol(Y)
	}
	
	nx = nrow(X)
	px = ncol(X)
	ny = nrow(Y)
	py = ncol(Y)
	
	J = diag(nx) - matrix(1, nx,nx) / nx
    if (centre) {
    	X= J %*% X
    	Y= J %*% Y
    }
	Procrustes = list()
	Procrustes$xmean = apply(X,2,mean)
	Procrustes$X = X 
	Procrustes$Y = Y
	
	C = t(X) %*% J %*% Y
	SVD = svd(C)
	T = SVD$u %*% t(SVD$v)
	s = tr(C %*% T)/ tr(t(Y) %*% J %*% Y)
	t=apply((X-s*Y%*%T),2, sum) / ny
  Z= s*Y %*% T + matrix(1,nx, 1) %*% t
    
    Procrustes$Yrot=Z
    Procrustes$rotation=T
    Procrustes$translation=t
    Procrustes$scale=s
    Procrustes$rss = sum((X-Z)^2)
    Procrustes$fit = 1 - sum((X-Z)^2) / sum(X^2)
    Procrustes$correlations = cor(X,Z)
    class(Procrustes)="Procrustes"
    return(Procrustes)
}


tr <- function(x) (sum(diag(x)))