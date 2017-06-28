zeros <- function(n) {
	if (length(n) == 1) 
		res = matrix(0, n, n)
	else res = array(0, n)
	return(res)
}
