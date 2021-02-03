# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2014

PrincipalCoordinates <- function (Proximities, w = NULL, dimension = 2, tolerance=0.0001, Bootstrap=FALSE, BootstrapType=c("Distances", "Products"), nB=200, ProcrustesRot=TRUE, BootstrapMethod=c("Sampling", "Permutation")) {
  if (length(BootstrapType) > 1) BootstrapType = BootstrapType[1] 
  if (length(BootstrapMethod) > 1) BootstrapMethod = BootstrapMethod[1] 
  r=dimension
  if (!(class(Proximities)=="proximities")) stop("You need a proximities matrix")
  Dimnames=paste("Dim", 1:r)
  
  if (is.null(Proximities$Proximities))
    dis=Proximities$D
  else 
    dis=Proximities$Proximities
  n <- dim(dis)[1]

  if (is.null(w)) w = matrix(1, 1, n)/n
  else w= matrix(w/sum(w), nrow=1)
  
  if (length(w)!=n) stop("The vector of weigths is not correct")
  
  # Scalar Product matrix
  b <- -0.5 * (diag(n) - matrix(1, n, 1) %*% w) %*% dis^2 %*% (diag(n) - matrix(1, n, 1) %*% w)
  
  Dw = diag(as.vector(w))
  
  B=sqrt(Dw) %*% b %*% sqrt(Dw)

  solut <- svd(B)
  
  Inertia = (solut$d/sum(solut$d)) * 100
  g <- diag(as.vector(1/sqrt(w))) %*% solut$u %*% diag(sqrt(solut$d))
  rownames(g)=rownames(dis)
  ra=sum(as.numeric(solut$d>tolerance))
  st <- apply(g^2, 1, sum)
  qlr <- diag(1/st) %*% (g^2)
  qlr=round(qlr[, 1:r]*100, digits=2)
  rownames(qlr)=rownames(dis)
  colnames(qlr)=Dimnames
  cumqlr=t(apply(qlr,1,cumsum))
  Proximities$Analysis="Principal Coordinates"
  Proximities$EigenValues = solut$d
  Proximities$Inertia = Inertia
  Proximities$RowCoordinates = g[,1:r]
  colnames(Proximities$RowCoordinates)=Dimnames
  Proximities$RowQualities = qlr
  
  if (!is.null(Proximities$SupProximities)) {
    t=dim(Proximities$SupProximities)
    SupCoord=0.5* ((matrix(1,t,1) %*% matrix(diag(b),1,n)) - Proximities$SupProximities^2) %*% g %*% diag(solut$d^(-1))
    SupCoord=SupCoord[,1:ra]
    st <- apply(SupCoord^2, 1, sum)
    qlr <- diag(1/st) %*% (SupCoord^2)
    qlr=round(qlr[, 1:r]*100, digits=2)
    rownames(qlr)=rownames(Proximities$SupProximities)
    colnames(qlr)=Dimnames
    Proximities$SupRowCoordinates = SupCoord[,1:r]
    colnames(Proximities$SupRowCoordinates)=Dimnames
    Proximities$SupRowQualities = qlr
  }
  
  Dh= as.dist(dis)
  D=as.dist(EuclideanDistance( Proximities$RowCoordinates))
  stress=sum(((D-Dh)^2))
  scresid=sum((((D-Dh)^2)))
  scdis=sum((((D)^2)))
  dmean=sum(D)/length(D)
  scdesv=sum(((D-dmean)^2))
  Proximities$RawStress=stress
  Proximities$stress1=sqrt(scresid/scdis)
  Proximities$stress2=sqrt(scresid/scdesv)
  Proximities$sstress1=sqrt(sum(((D^2-Dh^2)^2))/sum(((D^2)^2)))
  dmean2=sum((D^2))/length(D)
  Proximities$sstress2=sqrt(sum(((D^2-Dh^2)^2))/sum(((D^2-dmean2)^2)))
  Proximities$rsq=cor(Dh,D)^2
  Proximities$rho=cor(Dh,D,method="spearman")
  Proximities$tau=cor(Dh,D,method="kendall")
  
  if (Bootstrap){
    if (BootstrapType=="Distances") Proximities$BootstrapInfo=BootstrapDistance(dis, W=diag(nrow(dis)), nB=nB, dimsol=dimension, ProcrustesRot=ProcrustesRot, method=BootstrapMethod)
    if (BootstrapType=="Products") Proximities$BootstrapInfo=BootstrapScalar(b, W=diag(nrow(dis)), nB=nB, dimsol=dimension, ProcrustesRot=ProcrustesRot, method=BootstrapMethod)
  }
  class(Proximities) <- "Principal.Coordinates"
  return(Proximities)
}