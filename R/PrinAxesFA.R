PrinAxesFA <- function(R, dimsol=3, method=1, tol=0.0001,  MaxIter=50){
  n=dim(R)[1]
  methods=c("A1", "HSC", "MC")
  if (is.numeric(method)) method=methods[method]
  
  switch(method, A1 = {
    ico=rep(1,n)
  }, HSC={
    ico =apply(abs(R*(1-diag(n))), 1, max)
  }, MC={
    ico=rep(1,n)
    for (i in 1:n){
      rc=matrix(R[i, -i], n-1, 1)
      rx=R[-i, -i]
      ico[i]=t(rc) %*% solve(rx) %*% rc
    }
  })
  
  R=R*(1-diag(n))+diag(ico)
  sol=svd(R)
  A=sol$v[,1:dimsol]%*% diag(sqrt(sol$d[1:dimsol]))
  Comunal=apply(A^2, 1, sum)
  error=sum((ico-Comunal)^2)
  ico=Comunal
  iter=0
  
  while ((error > tol) & (iter < MaxIter)){
    iter=iter+1
    R=R*(1-diag(n))+diag(ico)
    sol=svd(R)
    A=sol$v[,1:dimsol]%*% diag(sqrt(sol$d[1:dimsol]))
    Comunal=apply(A^2, 1, sum)
    error=sum((ico-Comunal)^2)
    ico=Comunal 
  }
  rownames(A)=rownames(R)
  colnames(A)=paste("Factor",1:dimsol, sep="_")
  eigenvalues=apply(A^2, 2, sum)
  names(eigenvalues)=paste("Factor",1:dimsol, sep="_")
  return(list(A=A, Eigenvalues=eigenvalues, SCT= n))
}

