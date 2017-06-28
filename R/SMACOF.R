SMACOF <- function(P, X=NULL, W=NULL, Model=c("Identity", "Ratio", "Interval", "Ordinal"), dimsol=2, maxiter=100, maxerror=0.000001, StandardizeDisparities=TRUE, ShowIter=FALSE){
  if (length(Model)>1) Model=Model[1]
  n <- dim(P)[1]
  if (is.null(W)) W = matrix(1,n,n)
  W=W-diag(diag(W))
  V=-1*W+diag(apply(W,2,sum))
  W=as.dist(W)
  # Initial solution (Principal coordinates)
  if (is.null(X)){
  b <- -0.5 * (diag(n) - matrix(1, n, n)/n) %*% P^2 %*% (diag(n) - matrix(1, n, n)/n)
  solut <- svd(b)
  X <- solut$u %*% diag(sqrt(solut$d))
  X=X[,1:dimsol]
  }
  k=0

  Z=X 

  D=as.dist(EuclideanDistance(X))
  Dh=Dhats(P,D,W,Model, Standardize=StandardizeDisparities)
  stress=sum(((D-Dh)^2)*W)
  error=1
  while ((k <= maxiter) & (error > maxerror)) {
    k=k+1
    BZ=as.matrix(-1*W*Dh/D)
    BZ[which(is.na(BZ))]<-0
    bzii=-1*apply(BZ,2,sum)
    BZ=BZ+diag(bzii)
    X=ginv(V)%*%BZ%*%Z
    Z=X
    D=as.dist(EuclideanDistance(Z))
    Dh=Dhats(P,D,W,Model, Standardize=StandardizeDisparities)
    newstress=sum(((D-Dh)^2)*W)
    error = abs(stress-newstress)
    if (ShowIter) print(c(k, error))
    stress=newstress
  }

  scresid=sum((((D-Dh)^2)*W))
  scdis=sum((((D)^2)*W))
  dmean=sum(D*W)/sum(W)
  scdesv=sum(((D-dmean)^2)*W)

  stress1=sqrt(scresid/scdis)
  stress2=sqrt(scresid/scdesv)
  sstress1=sqrt(sum(((D^2-Dh^2)^2)*W)/sum(((D^2)^2)*W))
  dmean2=sum((D^2)*W)/sum(W)
  sstress2=sqrt(sum(((D^2-Dh^2)^2)*W)/sum(((D^2-dmean2)^2)*W))
  rsq=cor(Dh*W,D*W)^2
  rho=cor(Dh*W,D*W,method="spearman")
  tau=cor(Dh*W,D*W,method="kendall")
  res=list(X=X,D=D,Dh=Dh,stress=stress, stress1=stress1, stress2=stress2, sstress1=sstress1, sstress2=sstress2, rsq=rsq, rho=rho, tau =tau )
  class(res)="SMACOF"
  return(res)
}

