ones <- function(n) {
	if (length(n) == 1) 
		res = matrix(1, n, n)
	else res = array(1, n)
	return(res)
}