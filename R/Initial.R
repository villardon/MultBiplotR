#  Initial values for the genefold algorithm.
#   INIT calculates first the coordinates (X and Y); the transformed data 
#   GAMMA and regression weights a are calculated then conditional on the 
#   values for the coordinates. 
#   INITIAL determines whether the coordinates are random, semi-rational or 
#   rational (using an SVD of the data matrix -augmented for the column 
#  elements- with coordinates according to a preference sphere 
#   (see Van Deun, Heiser, and Delbeke(in press))
#
# Using this file implies that you agree with the license (see license.txt)

# Translation from : K. Van Deun, Department of Psychology, Catholic University of
# Leuven (Belgium)
Initial <- function(DATA, TRANSFORMATION, INITIAL = "rational", DIMENS = 2, LENGTHX, Lup) {
  n = dim(DATA)[1]
  m = dim(DATA)[2]
  z = matrix(0, 1, DIMENS)
  switch(INITIAL, rational = {
    ranking = matrix(0, n, m)
    for (i in 1:n) ranking[i, ] = rank(DATA[i, ])
    columnscore = sum(2:m)/(m - 1)
    columncoord = columnscore * matrix(1, m, m) + diag(as.vector(matrix(1, m, 1))) - diag(columnscore * as.vector(matrix(1, 
                                                                                                                         m, 1)))
    augmdata = rbind(ranking, columncoord)
    augmdata_c = augmdata - mean(1:m) * matrix(1, n + m, m)
    norm = DistUnfold(augmdata_c, matrix(0, 1, m))
    norm[which(norm == 0)] = 1
    prefsphere = (((norm)^(-1) * m) %*% matrix(1, 1, m)) * augmdata_c
    H = matrix(1, 1, n + m) %*% prefsphere
    H = matrix(1, n + m, 1) %*% H
    prefsphere = (H/(n + m)) - prefsphere
    # 2. Find coordinates using a principal coordinates approach 
    SOL = svd(prefsphere, DIMENS, DIMENS)
    coord = SOL$u %*% diag(SOL$d[1:DIMENS])^2
    
  }, `semi-rational` = {
    
    ranking = matrix(0, n, m)
    for (i in 1:n) ranking[i, ] = rank(DATA[i, ])
    columnscore = sum(2:m)/(m - 1)
    columncoord = columnscore * matrix(1, m, m) + diag(as.vector(matrix(1, m, 1))) - diag(columnscore * as.vector(matrix(1, 
                                                                                                                         m, 1)))
    augmdata = rbind(ranking, columncoord)
    augmdata_c = augmdata - mean(1:m) * matrix(1, n + m, m)
    norm = DistUnfold(augmdata_c, matrix(0, 1, m))
    norm[which(norm == 0)] = 1
    prefsphere = (((norm)^(-1) * m) %*% matrix(1, 1, m)) * augmdata_c
    H = matrix(1, 1, n + m) %*% prefsphere
    H = matrix(1, n + m, 1) %*% H
    prefsphere = (H/(n + m)) - prefsphere
    # 2. Find coordinates using a principal coordinates approach 
    SOL = svd(prefsphere, DIMENS, DIMENS)
    coord = SOL$u %*% diag(SOL$d[1:DIMENS]) ^2
    R = max(DistUnfold(coord, z))
    coord = coord + 0.5 * R * (matrix(runif((n + m) * DIMENS), n + m, DIMENS) - 0.5)
  }, random = {
    coord = matrix(runif((n + m) * DIMENS), n + m, DIMENS) - 0.5
  })
  #3. Set reference space (center Y;norm for row points maximal 1) and an
  # intermixed start
  
  c_coord = coord - ((matrix(1, n + m, m)/m) %*% coord[(n + 1):(n + m), ])
  R = max(DistUnfold(c_coord, z))
  norm_coord = c_coord/R
  X = norm_coord[1:n, ]
  Yinit = norm_coord[(n + 1):(n + m), ]
  Y = (median(DistUnfold(X, z))/median(DistUnfold(Yinit, z))) * Yinit
  
  #CALCULATION OF GAMMA AND a CONDITIONAL ON COORDINATES
  D = DistUnfold(X, Y)
  v = matrix(1, m, m)
  DATAc = DATA - (DATA %*% v/m)
  F = (DATAc * DATAc) %*% v
  GAMMA = DATA/(F^(0.5))
  a = Updatea(GAMMA, D, Lup)
  sol=list(X=X, Y=Y, GAMMA=GAMMA, a=a)
  return(sol)
}