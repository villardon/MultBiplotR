PLSR1BinFit <- function(Y, X, S=2, tolerance=0.000005, maxiter=100, show=FALSE, penalization=0, cte =TRUE, Algorithm=1){
  
  I1=dim(X)[1]
  J=dim(X)[2]
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  I=I2
  
  uini=svd(X)
  T=matrix(0, I1, S)
  U=matrix(0, I1, S)
  W=matrix(0, J, S)
  C=matrix(0, K, S)
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  
  freq=matrix(1,I,1)
  w=matrix(0,J,1)
  # El algotitmo 1 se basa en Bastien et al.
  if (Algorithm==1){ 
    for (j in 1:J){
      x=as.matrix(X[,j])
      colnames(x)=xnames[j]
      fit=RidgeBinaryLogistic(Y, x)
      w[j]=fit$beta[2]
    }
    w=w/sqrt(sum(w^2))
    t=X %*% w
    T[,1]=t
    W[,1]=w
    p=t(X) %*% t/ sum(t^2)
    P[,1]=p
    X1=X-t %*% t(p)
    
    for (i in 2:S){
      for (j in 1:J){
        x=as.matrix(cbind(T[,1:(i-1)],X[,j]))
        colnames(x)=c(dimnames[1:(i-1)] ,xnames[j])
        fit=RidgeBinaryLogistic(Y, x)
        w[j]=fit$beta[(i+1)]
      }
      w=w/sqrt(sum(w^2))
      t=X1 %*% w
      T[,i]=t
      W[,i]=w
      p=t(X) %*% t/ sum(t^2)
      P[,i]=p
      X1=X1-t %*% t(p)
    }
    fit=RidgeBinaryLogistic(Y, T, penalization=penalization, cte=cte)
    
  }
  
  # El algoritmo 2 
  if (Algorithm==2){
    # We have to take the constant into account
    if (cte){
    t0=matrix(1, nrow=I, ncol=1)
    fit=RidgeBinaryLogistic(Y, t, penalization=penalization, cte=FALSE)
    c0=fit$beta[1]
    xb=matrix(apply(X,2,mean), ncol=1)
    X=X-t0%*%t(xb)

    }
    
    for (i in 1:S){
      error=1
      iter=0
      u=runif(I)
      while ((error>tolerance) & (iter<maxiter)){
        iter=iter+1
        w=(t(X) %*% u)/sum(u^2)
        w=w/sqrt(sum(w^2))
        t=X %*% w
        fit=RidgeBinaryLogistic(Y, t, penalization=penalization, cte=cte)
        c=fit$beta[2]
        newu= fit$linterm
        error=sum((u-newu)^2)
        u=newu
        if (show) print(c(i, iter, error))
      }
      T[,i]=t
      #U[,i]=u
      W[,i]=w
      C[,i]=c
      p=t(X) %*% t/ sum(t^2)
      P[,i]=p
      X=X-t %*% t(p)
    }
    

    
    fit = list()
    null = list()
    # Null Model
    null$beta = c0
    null$linterm = matrix(1,I,1) %*% null$beta
    null$fitted = exp(null$linterm)/(1 + exp(null$linterm))
    null$Deviance = -2 * sum(Y * log(null$fitted) + (1 - Y) * log(1 - null$fitted))
    
    # Full Model
    fit$y=Y
    fit$x=T
    fit$Penalization=penalization
    fit$beta = matrix(c(c0, C), ncol=1)
    colnames(fit$beta)="Beta"
    T1=cbind(rep(1,I),T)
    fit$linterm =  T1 %*% fit$beta
    
    fit$fitted = exp(fit$linterm)/(1 + exp(fit$linterm))
    fit$residuals = Y - fit$fitted
    fit$Prediction = as.numeric(fit$fitted >= 0.5)
    
    v = (fit$fitted * (1 - fit$fitted))
    vv = diag(1, I, I)
    diag(vv) <- v
    Imod=diag(S+1)
    Imod[1,1]=0
    In = (t(T1) %*% vv %*% T1) + 2 * penalization * Imod
    fit$Covariances = solve(In)
    fit$Deviance = -2 * sum(Y * log(fit$fitted) * freq + (1 - Y) * log(1 - fit$fitted) * freq)
    fit$NullDeviance=null$Deviance
    fit$Dif=(null$Deviance - fit$Deviance)
    fit$df=length(fit$beta)-length(null$beta)
    fit$p=1-pchisq(fit$Dif, df =  fit$df)
    fit$CoxSnell=1-exp(-1*fit$Dif/I)
    fit$Nagelkerke=fit$CoxSnell/(1-exp((null$Deviance/(-2)))^(2/I))
    fit$MacFaden=1-(fit$Deviance/null$Deviance)
    fit$SSRes=sum((fit$residuals^2))
    fit$SSTot=sum((Y^2))
    fit$R2 = 1 - (sum((fit$residuals^2))/sum((Y^2)))
    fit$Classification=table(Y,fit$Prediction)
    fit$PercentCorrect = sum(fit$Prediction == Y)/I
    class(fit) <- "RidgeBinaryLogistic"
  }
  result=list(c0=c0, T=T, W=W, P=P, U=U, C=C, Q=Q, fit=fit)
  return(result)
}