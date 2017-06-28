# LOSS Calculate unfolding loss.
#   LOSS(a,GAMMA,D) Calculates the loss between the distances and optimal 
#   transformed data with the loss function given by 
#   n^(-1) \sum_i||GAMMA_i-a_i*D_i||²_J (_J denotes a centering of the 
#   deviation vector GAMMA_i-a_i*D_i).
#
# Using this file implies that you agree with the license (see license.txt)
# From K. Van Deun, Department of Psychology, Catholic University of
# Leuven (Belgium).

Loss <- function(SOL) {
  n = dim(SOL$GAMMA)[1]
  m = dim(SOL$GAMMA)[2]
  int = SOL$GAMMA - ((SOL$a %*% matrix(1, 1, m)) * SOL$D)
  v = matrix(1, m, m)/m
  int = int - (int %*% v)
  L = sum(sum(int * int))/n
  return(L)
}
