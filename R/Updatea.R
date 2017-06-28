#UPDATEA Update of the unfolding regression weights
#   UPDATEA computes the regression weights a_i under the restriction 
#   a_i <= c for all i and with c an upper bound.
#
#Using this file implies that you agree with the license (see license.txt)

#Author: K. Van Deun, Department of Psychology, Catholic University of
#Leuven (BELGIUM)

Updatea <- function(GAMMA, D, Lup) {
  n = dim(D)[1]
  m = dim(D)[2]
  
  v = matrix(1, m, m)
  vm = matrix(1, m, 1)
  GAMMAc = GAMMA - (GAMMA %*% v/m)
  Dc = D - (D %*% v/m)
  aup = ((GAMMAc * D) %*% vm)/((Dc * Dc) %*% vm)
  aup[which(aup < 0)] = 0
  aup[which(aup > Lup)] = Lup
  return(aup)
}
