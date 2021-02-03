PLSR <- function(Y, X, S=2,  InitTransform=5, grouping=NULL,  centerY=TRUE, scaleY=TRUE, tolerance=0.000005, maxiter=100, show=FALSE, Validation=NULL, nB=500){
  
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
  
  if (Validation=="Cross"){
    CrossR2=matrix(0,I,1)
    PsParameters=matrix(0,J,I)
    CrossParameters=matrix(0,J,I)
    t0=matrix(1, I-1,1)
    for (i in 1:I){
      Y1=matrix(result$Y[-i,], ncol=K)
      X1=result$X[-i,]
      
      if (centerY & !scaleY){
        Y1=TransformIni(Y1,transform=4)
      }
      if (centerY & scaleY){
        Y1=TransformIni(Y1,transform=5)
      }
      Data = InitialTransform(X1, transform = InitTransform, grouping=grouping)
      X1 = Data$X
      fit=PLSRfit(Y1, X1, S=S, tolerance=tolerance, maxiter=maxiter, show=show)
      ExpectedY= t0%*%t(fit$c0) + fit$T %*% t(fit$C)
      CrossR2[i]=diag(t(ExpectedY)%*%ExpectedY)/diag(t(Y1)%*%Y1)
      CrossParameters[,i]=fit$W%*% t(fit$C)
      PsParameters[,i] = result$RegParameters +  (I-1) * ( result$RegParameters -fit$W%*% t(fit$C))
    }
    
    
    rownames(PsParameters)=xnames
    rownames(CrossParameters)=xnames
    
    result$Validation=Validation
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
  
  
  if (Validation=="Bootstrap"){
    CrossR2=matrix(0,nB,1)
    CrossParameters=matrix(0,J,nB)
    t0=matrix(1, I,1)
    for (i in 1:nB){
      muestra=sample.int(I, I, replace=TRUE)
      Y1=matrix(Y[muestra,], ncol=K)
      X1=X[muestra,]
      if (centerY & !scaleY){
        Y1=TransformIni(Y1,transform=4)
      }
      if (centerY & scaleY){
        Y1=TransformIni(Y1,transform=5)
      }
      Data = InitialTransform(X1, transform = InitTransform, grouping=grouping)
      X1 = Data$X
      fit=PLSRfit(Y1, X1, S=S, tolerance=tolerance, maxiter=maxiter, show=show)
      ExpectedY= t0%*%t(fit$c0) + fit$T %*% t(fit$C)
      CrossR2[i]=diag(t(ExpectedY)%*%ExpectedY)/diag(t(Y1)%*%Y1)
      CrossParameters[,i]=fit$W%*% t(fit$C)
    }
    
    if (!is.null(Validation)){
    rownames(CrossParameters)=xnames
    
    result$Validation=Validation
    result$CrossR2=CrossR2
    result$CrossParameters=CrossParameters
    stdErr= apply(CrossParameters, 1, sd)
    Z=result$RegParameters/stdErr
    pval=dnorm(Z)
    EI=result$RegParameters - 1.96 *stdErr
    ES=result$RegParameters + 1.96 *stdErr
    
    result$RegParameters=cbind(result$RegParameters, stdErr, Z , round(pval, digits=5), EI, ES)
    colnames(result$RegParameters)=c("Beta", "Std. Error", "Z", "p-val", "CI : Lower", "CI : Upper")
    }
    
  }
  
  
  class(result)="PLSR"
  return(result)
}
