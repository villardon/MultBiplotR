# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

Numeric2Binary <- function(y, name = "MyVar", cut = NULL) {
	n = length(y)
	if (is.null(cut)) 
		cut = median(y)
	Z = matrix(0, n, 1)
	for (i in 1:n) if (y[i] > cut) 
		Z[i] = 1
	colnames(Z) <- paste(name, " > ", round(cut, digits = 3), sep = "")
	return(Z)
}