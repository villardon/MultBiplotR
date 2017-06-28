RidgeBinaryLogistic <- function(y, xd, freq = NULL, tolerance = 1e-05, maxiter = 100, penalization = 0.2, cte = TRUE) {
# Fits a Binary Logistic Regression

n=dim(xd)[1]
  if (is.numeric(y)) {y= as.matrix(y)}
  if (is.matrix(xd) & cte) xd=cbind(rep(1,n),xd)
  if (is.data.frame(xd)) xd=DataFrame2Matrix4Regression(xd, Intercept=cte)
  if (is.vector(xd)){
    xd=matrix(xd, length(xd),1)
    if (cte) xd=cbind(rep(1,length(xd)),xd)
  }

  xd=as.matrix(xd)

n <- dim(xd)[1]
  m <- dim(xd)[2]
	
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
  model$x=xd
  model$Penalization=penalization
	model$beta = RidgeBinaryLogisticFit(y, xd, freq, tolerance = tolerance, maxiter = maxiter, penalization = penalization)
	colnames(model$beta)="Beta"
	model$linterm = xd %*% model$beta

	model$fitted = exp(model$linterm)/(1 + exp(model$linterm))
	model$residuals = y - model$fitted
	model$Prediction = as.numeric(model$fitted >= 0.5)
	
	v = (model$fitted * (1 - model$fitted))
	vv = diag(1, n, n)
	diag(vv) <- v
	Imod=diag(m)
	Imod[1,1]=0
	In = (t(xd) %*% vv %*% xd) + 2 * penalization * Imod
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
  yy=NULL
  pred=NULL
  for (i in 1:n){
    yy=c(yy,rep(y[i],freq[i]))
    pred=c(pred, rep(model$Prediction[i],freq[i]))
  }
	model$Classification=table(yy,pred)
	model$PercentCorrect = sum(pred == yy)/nn
	class(model) <- "RidgeBinaryLogistic"
	return(model)
}

