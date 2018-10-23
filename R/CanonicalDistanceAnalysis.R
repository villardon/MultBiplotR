CanonicalDistanceAnalysis <- function(Prox, group, dimens=2, Nsamples=1000, PCoA="Standard", ProjectInd=TRUE){
  cl <- match.call()
  if (!is.factor(group)) stop("The grouping variable must be a factor")
  if (class(Prox)!="proximities") stop("The D argument must be an object with Proximities")
  D=Prox$Proximities
  print(PCoA)
  PCoAs= c("Standard", "Weighted", "WPCA")
  if (is.numeric(PCoA)) PcoA=PCoAs(PCoA)
  Result=list()

  Result$call=cl
  Result$Type="CDA"
  Result$Distances=D
  Result$Groups=group
  GroupNames = levels(group)
  g = length(levels(group))
  n = dim(D)[1]

  G = Factor2Binary(group) # Matrix of indicators
  ng=diag(t(G) %*% G)

  D=0.5*D^2
  F=diag(1/ng) %*% (t(G) %*% D %*% G) %*% diag(1/ng)
  f=matrix(diag(F),g,1)
  DE = (2*F - f %*% matrix(1, 1, g) - t(f %*% matrix(1, 1, g)))
  DE=0.5*DE

  TSS=sum(D)/n
  BSS=(t(ng) %*% DE %*% ng)/n
  WSS=TSS-BSS
  Fexp=(BSS/(g-1))/(WSS/(n-g))
  Result$TSS= TSS
  Result$BSS = BSS
  Result$WSS = WSS
  Result$glt=n-1
  Result$glb=g-1
  Result$glw=n-g
  Result$Fexp=Fexp
  SamplingDist=rep(0,Nsamples)
  for (i in 1:Nsamples){
    mysample=sample(1:n)
    DS=D[mysample,mysample]
    FS=diag(1/ng) %*% (t(G) %*% DS %*% G) %*% diag(1/ng)
    f=matrix(diag(FS),g,1)
    DES = 0.5 * (2*FS - f %*% matrix(1, 1, g) - t(f %*% matrix(1, 1, g)))
    TSS=sum(D)/n
    BSS=(t(ng) %*% DES %*% ng)/n
    WSS=TSS-BSS
    SamplingDist[i]=(BSS/(g-1))/(WSS/(n-g))
  }

  # Result$SamplingDist=SamplingDist
  Result$pvalue=sum((SamplingDist - Fexp)>0)/Nsamples
  Result$Nsamples= Nsamples

  switch(PCoA, Standard = {
    H=(diag(g) - matrix(1, g, g)/g)
    B =  H %*% DE %*% H
  },Weighted = {
    H=(diag(g) - matrix(1, g, 1) %*% matrix(ng, 1, g)/n)
    B = H %*% DE %*% H
  },WPCA = {
    # Not finished (Have to be revised before it can be used)
    H=(diag(n) - matrix(1, n, n)/n)
    B =  H %*% D %*% H
  })

  solut <- svd(B)
  Result$ExplainedVariance = (solut$d/sum(solut$d)) * 100
  Y = solut$u %*% diag(sqrt(solut$d))
  d0=apply(Y^2,1, sum)
  rownames(Y)=GroupNames
  st <- apply(Y^2, 1, sum)
  Result$MeanCoordinates=Y[,1:dimens]
  colnames(Result$MeanCoordinates)=paste("Dim",1:dimens)
  qlr <- diag(1/st) %*% (Result$MeanCoordinates^2)
  Result$Qualities=round(qlr[, 1:dimens]*100, digits=2)
  rownames(Result$Qualities)=GroupNames
  colnames(Result$Qualities)=paste("Dim",1:dim(Result$Qualities)[2])
  Result$CummulativeQualities=t(apply(Result$Qualities,1,cumsum))

 if (ProjectInd){
   x=as.matrix(Prox$Data)
   Means= diag(1/ng) %*% t(G) %*% x
   Di=ContinuousProximities(Means, y=x, coef = Prox$Coefficient, r = Prox$r)$Proximities^2
   Y = Y[,1:dimens]
   Yi=-0.5 * solve(t(Y)%*%Y) %*% t(Y) %*% H %*% t(Di-matrix(1,n,1) %*% matrix(d0,1,g))
   Yi=t(Yi)
   rownames(Yi)=rownames(D)
   Result$RowCoordinates=Yi
 }

  Result=AddCluster2Biplot(Result, ClusterType="us", Groups=group )
  class(Result)="CanonicalDistanceAnalysis"
  return(Result)
}

