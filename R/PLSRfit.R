PLSRfit <- function(Y, X, S=2, tolerance=0.000005, maxiter=100, show=FALSE){
  I1=dim(X)[1]
  J=dim(X)[2]
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  I=I1
  
  uini=svd(X)
  T=matrix(0, I1, S)
  U=matrix(0, I1, S)
  W=matrix(0, J, S)
  C=matrix(0, K, S)
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  
  # Initial Step
  t0=matrix(1, I,1)
  c0=t(Y) %*% t0/ sum(t0^2)
  
  Y=Y-t0%*%t(c0)
  xb=matrix(apply(X,2,mean), ncol=1)
  X=X-t0%*%t(xb)
  
  for (i in 1:S){
    error=1
    iter=0
    u=matrix(uini$u[,i],I1,1)
    while ((error>tolerance) & (iter<maxiter)){
      iter=iter+1
      w=(t(X) %*% u)/sum(u^2)
      w=w/sqrt(sum(w^2))
      t=X %*% w
      t=t/sum(t^2)
      c=t(Y) %*% t/ sum(t^2)
      newu= Y %*% c
      error=sum((u-newu)^2)
      u=newu
      if (show) print(c(i, iter, error))
    }
    T[,i]=t
    U[,i]=u
    W[,i]=w
    C[,i]=c
    p=t(X) %*% t/ sum(t^2)
    P[,i]=p
    q=t(Y) %*% u/ sum(u^2)
    Q[,i]=q
    X=X-t %*% t(p)
  }
  result=list(c0=c0, T=T, W=W, P=P, U=U, C=C, Q=Q)
  return(result)
}


PLSR <- function(Y, X, S=2,  InitTransform=5, grouping=NULL,  centerY=TRUE, scaleY=TRUE, tolerance=0.000005, maxiter=100, show=FALSE, CrossValidation=FALSE){
  
  if (is.data.frame(X)) X=as.matrix(X)
  
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(InitTransform)) 
    InitTransform = ContinuousDataTransform[InitTransform]
  
  result=list()
  I1=dim(X)[1]
  J=dim(X)[2]
  
  
  if (is.numeric(Y) & !is.matrix(Y)) {Y= matrix(Y, ncol=1)
  colnames(Y)="Response"}
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Component", 1:S)
  result$Method="PLSR"
  result$X=X
  result$Y=Y
  result$centerY=centerY
  result$scaleY=scaleY
  result$Initial_Transformation=InitTransform
  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1
  
  rownames(Y)<-rownames(X)
  
  
  if (centerY & !scaleY){
    Y=TransformIni(Y,transform=4)
  }
  if (centerY & scaleY){
    Y=TransformIni(Y,transform=5)
  }
  
  Data = InitialTransform(X, transform = InitTransform, grouping=grouping)
  X = Data$X
  if (InitTransform=="Within groups standardization") result$Deviations = Data$ColStdDevs
  
  result$ScaledX=X
  result$ScaledY=Y
  
  fit=PLSRfit(Y, X, S=S, tolerance=tolerance, maxiter=maxiter, show=show)
  
  rownames(fit$P)=xnames
  colnames(fit$P)=dimnames
  rownames(fit$Q)=ynames
  colnames(fit$Q)=dimnames
  rownames(fit$C)=ynames
  
  rownames(fit$T)=inames
  colnames(fit$T)=dimnames
  rownames(fit$U)=inames
  colnames(fit$U)=dimnames
  rownames(fit$W)=xnames
  colnames(fit$W)=dimnames
  rownames(fit$C)=ynames
  colnames(fit$C)=dimnames
  
  result$Intercept=fit$c0
  result$XScores=fit$T
  result$XWeights=fit$W
  result$XLoadings=fit$P
  result$YScores=fit$U
  result$YWeights=fit$C
  result$YLoadings=fit$Q
  result$RegParameters=fit$W%*% t(fit$C)
  
  t0=matrix(1, I,1)
  result$ExpectedY= t0%*%t(fit$c0) + fit$T %*% t(fit$C)
  result$R2=diag(t(result$ExpectedY)%*%result$ExpectedY)/diag(t(Y)%*%Y)
  result$XStructure=cor(result$X,fit$T)
  result$YStructure=cor(result$Y,fit$U)
  result$YXStructure=cor(result$Y,fit$T)
  
  if (CrossValidation){
    CrossR2=matrix(0,I,1)
    PsParameters=matrix(0,J,I)
    CrossParameters=matrix(0,J,I)
    t0=matrix(1, I-1,1)
    for (i in 1:I){
      Y1=matrix(Y[-i,], ncol=K)
      X1=X[-i,]
      fit=PLSRfit(Y1, X1, S=S, tolerance=tolerance, maxiter=maxiter, show=show)
      ExpectedY= t0%*%t(fit$c0) + fit$T %*% t(fit$C)
      CrossR2[i]=diag(t(ExpectedY)%*%ExpectedY)/diag(t(Y1)%*%Y1)
      CrossParameters[,i]=fit$W%*% t(fit$C)
      PsParameters[,i] = result$RegParameters +  (I-1) * ( result$RegParameters -fit$W%*% t(fit$C))
    }
    
    rownames(PsParameters)=xnames
    rownames(CrossParameters)=xnames
    result$Crossvalidation=CrossValidation
    result$CrossR2=CrossR2
    result$CrossParameters=CrossParameters
    result$PsParameters = PsParameters
    stdErr= apply(PsParameters, 1, sd)/sqrt(I)
    Z=result$RegParameters/stdErr
    pval=dnorm(Z)
    EI=result$RegParameters - 1.96 *stdErr
    ES=result$RegParameters + 1.96 *stdErr
    
    result$RegParameters=cbind(result$RegParameters, stdErr, Z , pval, EI, ES)
    colnames(result$RegParameters)=c("Beta", "Std. Error", "Z", "p-val", "CI : Lower", "CI : Upper")
  }
  
  class(result)="PLSR"
  return(result)
}

summary.PLSR <- function(plsr){
  
}


plot.PLSR <- function(plsr, ParameterBoxPlot=FALSE, ParameterCI=TRUE, Correlations=FALSE, ...){
  I=length(plsr$RegParameters[,1])
  VarLabels= rownames(plsr$RegParameters)
  if (plsr$CrossVlidation){
    if (ParameterBoxPlot){
      boxplot(t(plsr$CrossParameters), main="Box plot for the Jackkinfe pseudo-samples")
      abline(h=0, col="red")}
    
    if (ParameterCI){
      signif=plsr$RegParameters[,4]>0.05
      plot(1:I, plsr$RegParameters[,1], ylim=range(c(plsr$RegParameters[,5], 
      plsr$RegParameters[,6])), pch=19, xlab="Variables", ylab="Mean +/- CI", 
      main="Scatter plot with confidence intervals error", col=signif+1)
      arrows(1:I, plsr$RegParameters[,5], 1:I, plsr$RegParameters[,6], length=0.05, angle=90, code=3, col=signif+1)
      text(1:I, plsr$RegParameters[,6], labels = VarLabels, pos=3, col=signif+1)
      abline(h=0, col="red", lty=3)
     }
    
  }
}
