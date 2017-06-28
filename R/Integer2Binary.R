# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

Integer2Binary <- function(y, name="My_Factor") {
	tt = table(y)
	y = as.integer(as.factor(y))
	ncat = max(y)
	n = length(y)
	Z = matrix(0, n, ncat)
	for (i in 1:n) Z[i, y[i]] = 1
	Nam = names(tt)
	for (j in 1:ncat) Nam[j] = paste(name, "-", names(tt)[j], sep = "")
	colnames(Z) <- Nam
	return(Z)
}