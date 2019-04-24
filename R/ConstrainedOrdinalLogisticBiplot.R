ConstrainedOrdinalLogisticBiplot <- function(Y, X, dim=2, Scaling=5,  tolerance = 1e-05, maxiter = 100, penalization = 0.1, show=FALSE){
  n <- nrow(Y)
  p <- ncol(Y)
  q <- ncol(X)
  r=min(p,q)
  
  VarNames=colnames(Y)
  RowNames=rownames(Y)
  Ncats=matrix(0,p,1)
  CategoryNames=list()
  for (j in 1:p){
    CategoryNames[[j]]=levels(Y[[j]])
    Ncats[j]= length(CategoryNames[[j]])
    }
  names(CategoryNames)=VarNames
  
  X=as.matrix(X)
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(Scaling)) 
    Scaling = ContinuousDataTransform[Scaling]
  
  
  Data = InitialTransform(X, transform = Scaling)
  X = Data$X
  if (Scaling=="Within groups standardization") Biplot$Deviations = Data$ColStdDevs
  
  Yh=matrix(0, n,p)
  Ajustes=list()
  for (i in 1:p){
    ajuste = RidgeOrdinalLogistic(Y[,i],X, tol = tolerance, maxiter = maxiter, penalization = penalization)
    Ajustes[[i]]=ajuste
    Yh[,i]=(ajuste$Eta[,1]-mean(ajuste$Eta[,1]))/sd(ajuste$Eta[,1])
  }
  
  names(Ajustes)=colnames(Y)
  
  SD = svd(Yh)
  EV = SD$d[1:r]^2
  Inertia = round((EV/sum(EV)) * 100, digits = 3)
  CumInertia = cumsum(Inertia)
  a = SD$u[,1:dim] %*% diag(SD$d[1:dim])
  freq = matrix(1, n, 1)

  Maxcat = max(Ncats)
  par = list()
  par$coefficients = array(0, c(p, dim))
  par$thresholds = array(0, c(p, Maxcat - 1))
  par$fit = array(0, c(p, 9))
  dimnames(par$fit)[[1]] = dimnames(Y)[[2]][1:dim(Y)[2]]
  dimnames(par$fit)[[2]] = c("logLik", "Deviance", "df", "p-value", 
                             "PCC", "CoxSnell", "Macfaden", "Nagelkerke", "NullDeviance")
  
  
  logLik = 0
  for (j in 1:p) {
    if (show) cat(" ", j)
    model = RidgeOrdinalLogistic(Y[, j], a, tol = tolerance, maxiter = maxiter, penalization = penalization, show = FALSE)
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
  model$Title="Constrained Ordinal Logistic Biplot"
  model$Type="Constrained"
  model$Fit="Ordinal Logistic Regression"
  model$Penalization=penalization
  model$NumberIterations=1
  model$CategoryNames=CategoryNames
  
  model$RowCoordinates = a
  rownames(model$RowCoordinates)=RowNames
  colnames(model$RowCoordinates)=paste("Dim_",1:dim,sep="")
  
  model$RowContributions=matrix(100/dim,n,dim)
  rownames(model$RowContributions)=RowNames
  colnames(model$RowContributions)=paste("Dim_",1:dim,sep="")
  
  model$ColumnParameters = par
  model$ColCoordinates = par$coefficients
  
  model$loadings = loadings
  rownames(model$loadings)=VarNames
  model$Communalities = matrix(r2, p,1)
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
  model$CoxSnell = 1 - exp(-1 * model$Dif/(n*p))
  model$Nagelkerke = model$CoxSnell/(1 - exp((model$DevianceNull/(-2)))^(2/(n*p)))
  model$MacFaden = 1 - (model$Deviance/model$DevianceNull)
  class(model) = "Ordinal.Logistic.Biplot"
  model$AIC=-2*model$LogLikelihood + 2 * model$df
  model$BIC=-2*model$LogLikelihood + model$df * log(n*p)
  model=AddContVars2Biplot(model,  X, dims=dim, Scaling = Scaling)
  return(model)
}

