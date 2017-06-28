diagonal <- function(d) {
	n = length(d)
	D = diag(1, n, n)
	diag(D) <- d
	D
}

