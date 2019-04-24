OrdLogBipEM <- function (Data, freq=NULL, dim = 2, nnodes = 15, tol = 0.0001, maxiter = 100, maxiterlogist = 100,
          penalization = 0.2, show = FALSE, initial = 1, alfa = 1, Orthogonalize=TRUE, Varimax=TRUE, ...) 
{
  nrow= dim(Data)[1]
  ncol= dim(Data)[2]
  factors=TRUE
  for (j in 1:ncol)
    factors = factors & is.ordered(Data[[j]])
  if (!factors) stop("You must provide a data frame with ordered factors to calculate an Ordinal Logistic Biplot")
  
  CategoryNames=list()
  for (j in 1:ncol)
    CategoryNames[[j]]=levels(Data[[j]])
  
  initials = c("MCA", "MIRT", "RAND", "PCoAB", "PCoAO")
  if (is.numeric(initial)) {
    initial = initials[initial]
  }
  VarNames=colnames(Data)
  RowNames=rownames(Data)
  
  CategoryNames=list()
  for (j in 1:ncol)
    CategoryNames[[j]]=levels(Data[[j]])
  names(CategoryNames)=VarNames
  
  x=ConvertFactors2Integers(Data)
  Q = Multiquad(nnodes, dim)
  X = Q$X
  A = Q$A
  q = dim(X)[1]
  p = dim(x)[2]
  G = Dataframe2BinaryMatrix(Data)
  s = dim(G)[1]
  n = dim(G)[2]
  
  Ncats = apply(x,2,max)
  Maxcat = max(Ncats)
  par = list()
  par$coefficients = array(0, c(p, dim))
  par$thresholds = array(0, c(p, Maxcat - 1))
  par$fit = array(0, c(p, 9))
  dimnames(par$fit)[[1]] = dimnames(x)[[2]][1:dim(x)[2]]
  dimnames(par$fit)[[2]] = c("logLik", "Deviance", "df", "p-value", 
                             "PCC", "CoxSnell", "Macfaden", "Nagelkerke", "NullDeviance")
  if (show) {
    cat("Calculating initial coordinates -", initial, "\n")
  }
  if (initial == "MCA") {
    corr = CA(G, dim = dim, alpha = alfa)
    ability = corr$RowCoordinates[, 1:dim]
  }
  if (initial == "MIRT") {
    technical = list()
    technical$NCYCLES = maxiterlogist
    mod = mirt(x, model=dim, technical = technical, ...)
    ability = fscores(mod, full.scores = TRUE)[, 1:dim]
  }
  if (initial == "RAND") {
    ability = matrix(rnorm(nrow*dim), nrow, dim)
  }
  if (initial == "PCoAB") {
    D=BinaryProximities(G)
    PC=PrincipalCoordinates(D, dimension=dim)
    ability=PC$RowCoordinates
  }
  if (initial == "PCoAO") {
    
  }
  
  logLik = 0
  for (j in 1:p) {
    model = RidgeOrdinalLogistic(as.ordered(x[, j]), ability, tol = tol, maxiter = maxiterlogist, 
                       penalization = penalization, show = F)
    par$coefficients[j, ] = model$coefficients
    par$thresholds[j, 1:nrow(model$thresholds)] = model$thresholds
    par$fit[j, 1] = model$logLik
    par$fit[j, 2] = model$Deviance
    par$fit[j, 3] = model$df
    par$fit[j, 4] = model$pval
    par$fit[j, 5] = model$PercentClasif
    par$fit[j, 6] = model$CoxSnell
    par$fit[j, 7] = model$MacFaden
    par$fit[j, 8] = model$Nagelkerke
    par$fit[j, 9] = model$DevianceNull
    logLik = logLik + model$logLik
  }
  
  
  if (show) {
    cat("Calculating Iterations\n")
  }
  
  logLikold = logLik
  error = 1
  iter = 0
  
  while ((error > tol) & (iter < maxiter)) {
    iter = iter + 1
    if (show) cat("Iteration", iter, " - Calculating Abilities ")
    PT = EvalOrdlogist(X, par, Ncats)
    L = matrix(1, s, q)
    for (l in 1:s) for (k in 1:q) L[l, k] = prod(PT[k, ]^G[l,])
    Pl = L %*% A
    ability = matrix(0, s, dim)
    for (l in 1:s) for (j in 1:dim) {
      for (k in 1:q) ability[l, j] = ability[l, j] + X[k, 
                                                       j] * L[l, k] * A[k]
      ability[l, j] = ability[l, j]/Pl[l]
    }
    if (show) cat(" -> Fitting Variables")
    logLik = 0
    for (j in 1:p) {
      if (show) cat(" ", j)
      model = OrdinalLogisticFit(x[, j], ability, tol = tol, maxiter = maxiterlogist, 
                         penalization = penalization, show = F)
      par$coefficients[j, ] = model$coefficients
      par$thresholds[j, 1:nrow(model$thresholds)] = model$thresholds
      logLik = logLik + model$logLik
    }
    error = abs((logLikold - logLik)/logLik)
   
    logLikold = logLik
    if (show) {
      cat("\nIteration ", iter, "- Log-Lik:", logLik, 
                  " - Change:", error, "\n")
    }
  }
  
  if (Orthogonalize){
    if (show) cat("Orthogonalizing abilities\n")
    ability= OrthogonalizeScores(ability)
  }
  
  if (show) {
    cat("Fitting variables after convergence")
  }
  
  logLik = 0
  for (j in 1:p) {
    if (show) cat(" ", j)
    model = RidgeOrdinalLogistic(as.ordered(x[, j]), ability, tol = tol, maxiter = maxiterlogist, 
                                 penalization = penalization, show = F)
    par$coefficients[j, ] = model$coefficients
    par$thresholds[j, 1:nrow(model$thresholds)] = model$thresholds
    par$fit[j, 1] = model$logLik
    par$fit[j, 2] = model$Deviance
    par$fit[j, 3] = model$df
    par$fit[j, 4] = model$pval
    par$fit[j, 5] = model$PercentClasif
    par$fit[j, 6] = model$CoxSnell
    par$fit[j, 7] = model$MacFaden
    par$fit[j, 8] = model$Nagelkerke
    par$fit[j, 9] = model$DevianceNull
    logLik = logLik + model$logLik
  }
  
  rownames(par$coefficients) = VarNames
  colnames(par$coefficients) = paste("Dim_",1:dim,sep="")
  rownames(par$thresholds) = VarNames
  colnames(par$thresholds) = paste("C_",1:(Maxcat-1),sep="")
  
  d = sqrt(rowSums(par$coefficients^2) + 1)
  loadings = solve(diag(d)) %*% par$coefficients
  thresholds = solve(diag(d)) %*% par$thresholds
  r2 = rowSums(loadings^2)
  model = list()
  model$Title="Ordinal Logistic Biplot (EM)"
  if (maxiter>0) model$Type="Internal"
  else model$Type="External"
  if (maxiter>0) model$Fit="EM (Alternated Expectation-Maximization)"
  else model$Fit="Ordinal Logistic Regression fitted to a External Configuration"
  model$InitialConfiguration=initial
  model$Penalization=penalization
  model$NumberIterations=iter
  model$CategoryNames=CategoryNames
  
  model$RowCoordinates = ability
  rownames(model$RowCoordinates)=RowNames
  colnames(model$RowCoordinates)=paste("Dim_",1:dim,sep="")
  
  model$RowContributions=matrix(100/dim,nrow,dim)
  rownames(model$RowContributions)=RowNames
  colnames(model$RowContributions)=paste("Dim_",1:dim,sep="")

  model$ColumnParameters = par
  model$ColCoordinates = par$coefficients
  
  model$loadings = loadings
  rownames(model$loadings)=VarNames
  model$Communalities = matrix(r2, ncol,1)
  rownames(model$Communalities)=VarNames
  colnames(model$Communalities)="Communalities"
  model$ColContributions = loadings^2 
  rownames(model$ColContributions)=VarNames
  model$LogLikelihood = logLik
  model$Ncats = Ncats
  model$DevianceNull=sum(par$fit[,9])
  model$Deviance=sum(par$fit[,2])
  model$df=sum(par$fit[,3])
  model$Dif = (model$DevianceNull - model$Deviance)
  model$pval = 1 - pchisq(model$Dif, df = model$df)
  model$CoxSnell = 1 - exp(-1 * model$Dif/(nrow*ncol))
  model$Nagelkerke = model$CoxSnell/(1 - exp((model$DevianceNull/(-2)))^(2/(nrow*ncol)))
  model$MacFaden = 1 - (model$Deviance/model$DevianceNull)
  class(model) = "Ordinal.Logistic.Biplot"
  model$AIC=-2*model$LogLikelihood + 2 * model$df
  model$BIC=-2*model$LogLikelihood + model$df * log(nrow*ncol)
  return(model)
}





# Calculates the expected probabilities for a latent trait model
# with ordinal manifiest variables
EvalOrdlogist <- function(X, par, Ncats) {
  MaxCat = max(Ncats)
  dims = dim(par$coefficients)[2]
  nitems = dim(par$coefficients)[1]
  nnodos = dim(X)[1]
  Numcats = sum(Ncats)
  CumCats = cumsum(Ncats)
  eta = X %*% par$coefficients[1, ]
  ETA = matrix(1, nnodos, 1) %*% t(par$thresholds[1,1:(Ncats[1]-1)]) - eta %*% matrix(1, 1, (Ncats[1] - 1))	
  PIA = exp(ETA)/(1 + exp(ETA))
  PIA = cbind(PIA, matrix(1, nnodos, 1))
  PI = matrix(0, nnodos, Ncats[1])
  PI[, 1] = PIA[, 1]
  PI[, 2:Ncats[1]] = PIA[, 2:Ncats[1]] - PIA[, 1:(Ncats[1] - 1)]
  PT = PI
  for (j in 2:nitems) {
    eta = X %*% par$coefficients[j, ]
    ETA = matrix(1, nnodos, 1) %*% t(par$thresholds[j,1:(Ncats[j]-1)]) - eta %*% matrix(1, 1, (Ncats[j] - 1))			
    PIA = exp(ETA)/(1 + exp(ETA))
    PIA = cbind(PIA, matrix(1, nnodos, 1))
    PI = matrix(0, nnodos, Ncats[j])
    PI[, 1] = PIA[, 1]
    PI[, 2:Ncats[j]] = PIA[, 2:Ncats[j]] - PIA[, 1:(Ncats[j] - 1)]
    PT = cbind(PT, PI)
  }
  return(PT)
}
