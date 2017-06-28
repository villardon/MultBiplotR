# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

WeightedPCoA <- function(Proximities, weigths = matrix(1,dim(Proximities$Proximities)[1],1), dimension = 2, tolerance=0.0001) {
  r=dimension
  if (!(class(Proximities)=="proximities")) stop("You need a proximities matrix")
  Dimnames=NULL
  for (j in 1:r) Dimnames=c(Dimnames,paste("Dim",j))
  dis=Proximities$Proximities
  n <- dim(dis)[1]
  weigths=weigths/sum(weigths)
  H= diag(n) - matrix(1, n, 1)%*%t(weigths)
  b <- -0.5 * H %*% dis^2 %*% H
  solut <- svd(b)
  Inertia = (solut$d/sum(solut$d)) * 100
  g <- solut$u %*% diag(sqrt(solut$d))
  rownames(g)=rownames(dis)
  ra=sum(as.numeric(solut$d>tolerance))
  st <- apply(g^2, 1, sum)
  qlr <- diag(1/st) %*% (g^2)
  qlr=round(qlr[, 1:r]*100, digits=2)
  rownames(qlr)=rownames(dis)
  colnames(qlr)=Dimnames
  cumqlr=t(apply(qlr,1,cumsum))
  Proximities$EigenValues = solut$d
  Proximities$Inertia = Inertia
  Proximities$RowCoordinates = g[,1:r]
  colnames(Proximities$RowCoordinates)=Dimnames
  Proximities$RowQualities = qlr
  class(Proximities) <- "Principal.Coordinates"
  return(Proximities)
  
}

