UnfoldingSMACOF <- function(P, W = matrix(1, dim(P)[1], dim(P)[2]), Constrained = FALSE, 
                            ENV = NULL, model = "Ratio", condition = "Matrix", r = 2, maxiter = 100, 
                            tolerance = 1e-05, lambda = 1, omega = 0, X0, Y0){
  # Formulas from Heiser (Ecology book)
  n = dim(P)[1]
  m = dim(P)[2]
  DimNames = "Dim-1"
  for (i in 2:r) DimNames = c(DimNames, paste("Dim-", i, sep = ""))
  
  ColNames = colnames(P)
  RowNames = rownames(P)
  
  # Calculo de las distancias en la configuracion
  D = DistUnfold(X0, Y0)
  
  # Calculo de las disparidades optimas
  # [E,pstress]=cvpenalize(D,P,W,lambda,omega);
  
  DH = dhatsunfol(P, D, W, modelo = model, condicion = condition)
  cv2 = sd(as.vector(DH$Dh))/mean(DH$Dh)
  pstress = sum(sum((W * (D - DH$Dh))^2))^lambda * (1 + omega/cv2)
  
  Conf = UpdateUnfol(X0, Y0, W, DH$Dh)
  
  errorest = 1
  k = 0
  history = NULL
  while ((k <= maxiter) & (errorest > tolerance)) {
    k = k + 1
    if (Constrained) 
      Conf = UpdateUnfolConstr(X0, Y0, W, DH$Dh, ENV)
    else Conf = UpdateUnfol(X0, Y0, W, DH$Dh)
    D = DistUnfold(Conf$X, Conf$Y)
    #[E,pstress1]=cvpenalize(D,Dh,W,lambda,omega)
    DH = dhatsunfol(P, D, W, model, condition)
    cv2 = sd(as.vector(DH$Dh))/mean(DH$Dh)
    newpstress = sum(sum((W * (D - DH$Dh))^2))^lambda * (1 + omega/cv2)
    errorest = (pstress - newpstress)^2
    pstress = newpstress
    X0 = Conf$X
    Y0 = Conf$Y
    history = rbind(history, c(k, errorest))
    print(c(k, errorest))
  }
  
  rownames(Conf$X) = RowNames
  colnames(Conf$X) = DimNames
  rownames(Conf$Y) = ColNames
  colnames(Conf$Y) = DimNames
  
  dmean = sum(sum(D * W))/sum(sum(W))
  stress1 = sum(sum(((D - DH$Dh)^2) * W))/sum(sum(((D)^2) * W))
  stress2 = sum(sum(((D - DH$Dh)^2) * W))/sum(sum(((D - dmean)^2) * W))
  rs = cor(as.vector(DH$Dh * W), as.vector(D * W))
  rsq = rs^2
  
  sstress1 = sqrt(sum(sum(((D^2 - DH$Dh^2)^2) * W))/sum(sum(((D^2)^2) * W)))
  dmean2 = sum(sum((D^2) * W))/sum(sum(W))
  sstress2 = sqrt(sum(sum(((D^2 - DH$Dh^2)^2) * W))/sum(sum(((D^2 - dmean2)^2) * W)))
  rho = cor(as.vector(DH$Dh * W), as.vector(D * W), method = "spearman")
  tau = cor(as.vector(DH$Dh * W), as.vector(D * W), method = "kendall")
  
  fitcols = matrix(0, m, 1)
  for (i in 1:m) fitcols[i] = cor((DH$Dh * W)[, i], (D * W)[, i])^2
  rownames(fitcols) = ColNames
  colnames(fitcols) = "R-Squared"
  
  colnames(history) = c("Iteration", "Error")
  rownames(DH$Tol) = colnames(P)
  colnames(DH$Tol) = "Tolerance"
  if (Constrained) 
    Analysis = "Constrained Unfolding - SMACOF"
  else Analysis = "Unfolding - SMACOF"
  
  Unfold = list()
  Unfold$call <- match.call()
  Unfold$Analysis = Analysis
  Unfold$nsites = n
  Unfold$nspecies = m
  Unfold$nvars = NULL
  Unfold$alpha = NULL
  Unfold$dimens = r
  Unfold$Abundances = P
  Unfold$TransformedAbundances = P
  Unfold$Disparities = DH$Dh
  Unfold$Distances = D
  Unfold$Environment = ENV
  Unfold$TransfEnvironment = ENV
  Unfold$Weighted_Averages = NULL
  Unfold$Minima_Z = NULL
  Unfold$Maxima_Z = NULL
  Unfold$Medians_Z = NULL
  Unfold$P25_Z = NULL
  Unfold$P75_Z = NULL
  Unfold$Means = NULL
  Unfold$Deviations = NULL
  Unfold$IterationHistory = history
  Unfold$X = Conf$X
  Unfold$Y = Conf$Y
  Unfold$Tolerance = DH$Tol
  Unfold$RawStress = pstress
  Unfold$Stress1 = stress1
  Unfold$Stress2 = stress2
  Unfold$Sstress1 = sstress1
  Unfold$Sstress2 = sstress2
  Unfold$Pearson = rs
  Unfold$Spearman = rho
  Unfold$Kendall = tau
  Unfold$RSQ = rsq
  Unfold$ColumnsFit = fitcols
  Unfold$model = model
  Unfold$condition = condition
  class(Unfold) = "Unfolding"
  return(Unfold)
}
