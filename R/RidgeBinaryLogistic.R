RidgeBinaryLogistic <- function(y, X=NULL, data=NULL, freq = NULL, tolerance = 1e-05, maxiter = 100, 
                                penalization = 0.2, cte = FALSE, ref="first", bootstrap=FALSE, nmB=100,
                                RidgePlot=FALSE, MinLambda=0, MaxLambda=2, StepLambda=0.1) {
  # Fits a Binary Logistic Regression using a Ridge penalization
  
  if (class(y)=="formula"){
    X=model.matrix(y, data=data)
    y=model.frame(y, data=data)[1]
    cte=FALSE
  }
  
  name=colnames(y)

  if (is.factor(y[[1]])){
    if (length(levels(y[[1]]))!=2) stop("The number of levels of the response must be 2")
    name=names(y)
   y = Factor2Binary(y[[1]], Name=name)
   if (ref=="first") y=y[,2]
   else y=y[,1]
  }
  
  n=dim(X)[1]
  if (is.numeric(y)) {y= as.matrix(y)}
  if (is.matrix(X) & cte) X=cbind(rep(1,n),X)
  if (is.data.frame(X)) X=DataFrame2Matrix4Regression(X, Intercept=cte)
  if (is.vector(X)){
    X=matrix(X, length(X),1)
    if (cte) X=cbind(rep(1,length(X)),X)
  }
  
  
  X=as.matrix(X)
  
  n <- dim(X)[1]
  m <- dim(X)[2]
  
  if (is.null(freq)) 
    freq = matrix(1, n, 1)
  model = list()
  null = list()
  # Null Model
  null$beta = RidgeBinaryLogisticFit(y, matrix(1,n,1), freq, tolerance = tolerance, maxiter = maxiter, penalization = penalization)
  null$linterm = matrix(1,n,1) %*% null$beta
  null$fitted = exp(null$linterm)/(1 + exp(null$linterm))
  null$Deviance = -2 * sum(y * log(null$fitted)*freq + (1 - y) * log(1 - null$fitted)*freq)
  
  # Full Model
  model$y=y
  model$x=X
  model$Penalization=penalization
  model$beta = RidgeBinaryLogisticFit(y, X, freq, tolerance = tolerance, maxiter = maxiter, penalization = penalization)
  colnames(model$beta)="Beta"
  model$linterm = X %*% model$beta
  
  model$fitted = exp(model$linterm)/(1 + exp(model$linterm))
  model$residuals = y - model$fitted
  model$Prediction = as.numeric(model$fitted >= 0.5)
  
  v = (model$fitted * (1 - model$fitted))
  vv = diag(1, n, n)
  diag(vv) <- v
  Imod=diag(m)
  Imod[1,1]=0
  In = (t(X) %*% vv %*% X) + 2 * penalization * Imod
  model$Covariances = solve(In)
  
  model$Deviance = -2 * sum(y * log(model$fitted) * freq + (1 - y) * log(1 - model$fitted) * freq)
  model$NullDeviance=null$Deviance
  model$Dif=(null$Deviance - model$Deviance)
  model $df=length(model$beta)-length(null$beta)
  model $p=1-pchisq(model$Dif, df =  model$df)
  nn=sum(freq)
  model$CoxSnell=1-exp(-1*model$Dif/nn)
  model$Nagelkerke=model$CoxSnell/(1-exp((null$Deviance/(-2)))^(2/nn))
  model$MacFaden=1-(model$Deviance/null$Deviance)
  model$SSRes=sum((model$residuals^2)*freq)
  model$SSTot=sum((y^2)*freq)
  model$R2 = 1 - (sum((model$residuals^2)*freq)/sum((y^2)*freq))
  obs=NULL
  pred=NULL
  for (i in 1:n){
    obs=c(obs,rep(y[i],freq[i]))
    pred=c(pred, rep(model$Prediction[i],freq[i]))
  }
  model$Classification=table(obs,pred)
  model$PercentCorrect = sum(pred == obs)/nn
  
  
  if (RidgePlot){
    print("RidgePlot")
    lambdas=seq(MinLambda, MaxLambda, StepLambda)
    nlamb=length(lambdas)
    betaslambda=matrix(0, nlamb, m)
    for (j in 1:nlamb){
      betaslambda[j,] = RidgeBinaryLogisticFit(y, X, freq, tolerance = tolerance, maxiter = maxiter, penalization = lambdas[j])
    }
    plot(lambdas, betaslambda[,1], type="l", xlim=c(0,max(lambdas)), ylim=c(min(betaslambda),max(betaslambda)), panel.first=grid(), axes=F)
    for (j in 2:m)
      points(lambdas, betaslambda[,j], type="l", col=j)
    abline(v=penalization)
    abline(h=0)
    axis(1, pos=0)
    axis(2, pos=0)

  }
  
  if (bootstrap){
    model$bootstrap=TRUE
    betas=matrix(0, nmB, m)
    for (j in 1:nmB){
      betas[j,] = RidgeBinaryLogisticFit(y, X, freq, tolerance = tolerance, maxiter = maxiter, penalization = penalization)
    }
  }
  class(model) <- "RidgeBinaryLogistic"
  return(model)
}

