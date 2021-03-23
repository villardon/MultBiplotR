
AddOrdVars2Biplot <- function(bip, Y, tol = 0.000001, maxiterlogist = 100, penalization = 0.2, 
                              showiter = TRUE, show=FALSE){
  
  nrow= dim(Y)[1]
  ncol= dim(Y)[2]
  
  CategoryNames=list()
  for (j in 1:ncol)
    CategoryNames[[j]]=levels(Y[[j]])
  
  Ncats = rep(0,ncol)
  for (j in 1:ncol) Ncats[j]= length(levels(Y[[j]]))
  Maxcat = max(Ncats)
  
  Y=ConvertFactors2Integers(Y)
  
  p = dim(Y)[2]
  dim=dim(bip$RowCoordinates)[2]
  VarNames=colnames(Y)
  par = list()
  par$coefficients = array(0, c(p, dim))
  par$thresholds = array(0, c(p, Maxcat - 1))
  par$fit = matrix(0, p, 9)
  rownames(par$fit) = colnames(Y)
  colnames(par$fit) = c("logLik", "Deviance", "df", "p-value", 
                        "PCC", "CoxSnell", "Macfaden", "Nagelkerke", "NullDeviance")
  ability=bip$RowCoordinates
  logLik = 0
  for (j in 1:p) {
    if (show) cat(" ", j)
    model = RidgeOrdinalLogistic(Y[, j], ability, tol = tol, maxiter = maxiterlogist, 
                                 penalization = penalization, show = FALSE)
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
  if (p>1){
    loadings = solve(diag(d)) %*% par$coefficients
    thresholds = solve(diag(d)) %*% par$thresholds}
  else
  {
    loadings = par$coefficients/d
    thresholds = par$thresholds/d}
  r2 = rowSums(loadings^2)
  model = list()
  model$Title="Ordinal Logistic Biplot (External)"
  
  model$Type="External"
  model$Fit="Ordinal Logistic Regression fitted to a External Configuration"
  model$Penalization=penalization
  model$CategoryNames=CategoryNames
  
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
  class(model) = "OrdSupVarsBiplot"
  model$AIC=-2*model$LogLikelihood + 2 * model$df
  model$BIC=-2*model$LogLikelihood + model$df * log(nrow*ncol)
  bip$OrdSupVarsBiplot=model
  return(bip)
  
}

