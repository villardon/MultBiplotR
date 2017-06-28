VarimaxRot <- function (a){
  n=dim(a)[1]
  p=dim(a)[2]
  TT=diag(p)
  h=sqrt(apply(a^2, 1, sum))
  anorm=diag(1/h) %*% a
  simp=n*sum(anorm^4) - sum(apply(anorm^2,2,sum)^2)
  maxiter = 40
  iter=0
  err=simp

  while ((err>0.00001) & (iter<maxiter)){
    iter=iter+1
    for (k in 1:(p-1))
    for (l in (k+1):p){
      u=anorm[,k]^2-anorm[,l]^2
      vv= 2 * anorm[,k] * anorm[,l]
      aa=sum(u)
      bb=sum(vv)
      cc=sum(u^2-vv^2)
      dd=2*sum(u*vv)
      e=dd-2*aa*bb/n
      f=cc-(aa^2-bb^2)/n
      ang=atan2(e,f)/4
      T1=diag(p)
      T1[k,k]=cos(ang)
      T1[k,l]=-1*sin(ang)
      T1[l,k]=sin(ang)
      T1[l,l]=cos(ang)
      TT=TT %*% T1
    }
    anorm=diag(1/h) %*% a %*% TT
    simp2=n*sum(anorm^4) - sum(apply(anorm^2,2,sum)^2)
    err = abs(simp2-simp)
    simp=simp2
  }
  b= a%*% TT
  colnames(TT)=colnames(a)
  rownames(TT)=colnames(a)
  dimnames(b)=dimnames(a)
  res=list(loadings=b, rotmat=TT)
  
  return(res)
}

