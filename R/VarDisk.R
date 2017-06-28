#  VARDISK Variance of distances on unit disk
#   VARDISK determines the distribution of the variance of the distance of 
#   a point in a disk to n other points in the disk using Monte Carlo
#   simulation. The function returns the average variance over 1000
#   samples.
#
# Using this file implies that you agree with the license (see license.txt)

# Author: K. Van Deun, Department of Psychology, Catholic University of
# Leuven (BELGIUM)

VarDisk <- function(n) {
  runs = 1000
  D = matrix(1, runs, n)
  for (i in 1:runs) {
    r_refp = sqrt(runif(1))
    theta_refp = 2 * pi * runif(1)
    refcoord = matrix(c(r_refp * cos(theta_refp), r_refp * sin(theta_refp)), 1, 2)
    r_ps = sqrt(runif(n))
    theta_ps = 2 * pi * runif(n)
    coord = cbind(r_ps * cos(theta_ps), r_ps * sin(theta_ps))
    D[i, ] = DistUnfold(refcoord, coord)
  }
  Dc = D - (D %*% matrix(1, n, n)/n)
  r = sum((Dc^2))/(n * runs)
  return(r)
}
