PCA.Bootstrap <- function(X, dimens=2, Scaling = "Standardize columns", B=1000, type="np"){
  # X is supposed to be transformed previously
  # If data are not standardized only non parametric Bootstrap should be used
  # For the moment, only reflection of the eigenvectors is considered
  # Some kind of rotation (Procrustes) may be needed 
  # type = c("np", "pa", "spper", "spres")
  # np: Non Parametric
  # pa: parametric (resuduals are suppose to have a normal distribution)
  # spper: Semi-parametric Residuals are permutated
  # spres: Semi-parametric Residuals are resampled
  
  if (is.data.frame(X)) X=as.matrix(X)
  n= dim(X)[1]
  p= dim(X)[2]
  
  rnames=rownames(X)
  cnames=colnames(X)
  dimnames=paste("Dim", 1:dimens)
  
  res=list()
  res$Type=type
  res$InitTransform=Scaling
  res$InitialData=X  
  XT=InitialTransform(X, transform = Scaling)$X
  res$TransformedData=XT
  res$Type=type

  # Initial Calculation of the Principal Components
  acp=svd(XT, nu=dimens, nv=dimens)
  rownames(acp$u)=rnames
  colnames(acp$u)=dimnames
  rownames(acp$v)=cnames
  colnames(acp$v)=dimnames
  V=acp$v
  res$InitialSVD=acp
  res$InitScores=XT %*% acp$v
  res$InitCorr= cor(XT, res$InitScores)
  # Initialization of the matrices of bootstrap estimates
  res$Samples=matrix(0, B, n)
  res$EigVal=matrix(0, B, min(n,p))
  res$Inertia=matrix(0, B, min(n,p))
  res$Us=array(0, dim=c(n, dimens, B))
  res$Vs=array(0, dim=c(p, dimens, B))
  res$As=array(0, dim=c(n, dimens, B))
  res$Scores=array(0, dim=c(n, dimens, B))
  res$Bs=array(0, dim=c(p, dimens, B))
  res$Struct=array(0, dim=c(p, dimens, B))
  
  if (type=="np"){
    # Resampling individuals
    for (i in 1:B){
      samp=sample.int(n, size=n, replace=TRUE)
      res$Samples[i,]=samp
      XB=X[samp,]
      XB=InitialTransform(XB, transform = Scaling)$X
      acpB=svd(XB, nu=dimens, nv=dimens)
      signs=sign(diag(t(V) %*% acpB$v))
      acpB$v= acpB$v * matrix(1, p,1 ) %*% matrix(signs, ncol=dimens)
      acpB$u= acpB$u * matrix(1, n,1 ) %*% matrix(signs, ncol=dimens)
      rownames(acpB$u)=rnames
      colnames(acpB$u)=dimnames
      rownames(acpB$v)=cnames
      colnames(acpB$v)=dimnames
      
      res$EigVal[i,]=acpB$d
      res$Inertia[i,]=100*(acpB$d^2)/sum(acpB$d^2)
      res$Us[,,i]=acpB$u
      res$Vs[,,i]=acpB$v
      res$As[,,i]= XB %*% acpB$v
      res$Bs[,,i]= t(XB) %*% acpB$u
      res$Struct[,,i]=cor(XB, res$As[,,i])
      res$Scores[,,i]= XT %*% acpB$v
    }
  }
  class(res)="PCA.Bootstrap"
  return(res)
}

plot.PCA.Bootstrap <- function(pcaboot, Eigenvalues=TRUE, Inertia=FALSE, EigenVectors=TRUE, Structure=TRUE, Squared=TRUE, Scores=TRUE, ColorInd="black", TypeScores="ch"){
  require(Hmisc)
  rnames=rownames(pcaboot$InitialData)
  cnames=colnames(pcaboot$InitialData)
  if (Eigenvalues){
    boxplot(pcaboot$EigVal^2, main="Box-Plots for the Eigenvalues", xlab="Components", ylab="Eigenvalues")
    points(1:length(pcaboot$InitialSVD$d),pcaboot$InitialSVD$d^2, cex=1.5, pch=4, col="red")
    MedEig=apply(pcaboot$EigVal^2, 2,mean)
    ICPerc=apply(pcaboot$EigVal^2, 2,quantile, c(0.025, 0.975))
    expr = errbar(1:length(MedEig), MedEig, ICPerc[1,], ICPerc[2,], add=FALSE, pch=16, xlab="Components", ylab="Eigenvalues")
    points(1:length(MedEig),pcaboot$InitialSVD$d^2, cex=1.5, pch=4, col="red")
    title(main="5% CI for the Eigenvalues based on Percentiles")
    sdev=apply(pcaboot$EigVal^2, 2,sd)
    expr = errbar(1:length(MedEig), MedEig, MedEig-1.96*sdev,  MedEig+1.96*sdev, add=FALSE, pch=16, xlab="Components", ylab="Eigenvalues")
    points(1:length(MedEig),pcaboot$InitialSVD$d^2, cex=1.5, pch=4, col="red")
    title(main="5% CI for the Eigenvalues based on Moments")
    }
  
  
  dimens=dim(pcaboot$Vs)[2]
  
  if (Inertia){
    boxplot(pcaboot$Inertia, main="Bootstrap plot: Percent of Explained Variance")}
  
  if (EigenVectors){
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Vs[,i,]))
      maxim=apply(bb, 2,max)
      boxplot(bb, main=paste("PC",i))
      points(1:length(pcaboot$InitialSVD$v[,i]),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(MedEig), maxim, cnames, pos=3, cex=0.7, col="blue")
      abline(0,0, lty="dashed", col="red")
    }
    mtext("Box-Plots for the Eigenvectors", outer = TRUE, cex = 1.5)
    par(op)
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Vs[,i,]))
      media=apply(bb, 2,mean)
      ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
      expr = errbar(1:length(media), media, ICPerc[1,], ICPerc[2,], add=FALSE, pch=16, xlab="Variables", ylab="Coefficients")
      points(1:length(media),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(media), ICPerc[2,], cnames, pos=3, cex=0.8, col="blue")
      title(main=paste("PC",i))
      abline(0,0, lty="dashed", col="red")
    }
    mtext("5% CI for Eigenvectors based on Percentiles", outer = TRUE, cex = 1.5)
    par(op)
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Vs[,i,]))
      media=apply(bb, 2,mean)
      sdev=apply(bb, 2,sd)
      expr = errbar(1:length(media), media, media-1.96*sdev, media+1.96*sdev, add=FALSE, pch=16, xlab="Variables", ylab="Coefficients")
      points(1:length(media),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(media), media+1.96*sdev, cnames, pos=3, cex=0.8, col="blue")
      title(main=paste("PC",i))
      abline(0,0, lty="dashed", col="red")
    }
    mtext("5% CI for Eigenvectors based on Moments", outer = TRUE, cex = 1.5)
    par(op)

  }
  

  if (Structure){
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]))
      maxim=apply(bb, 2,max)
      boxplot(bb, main=paste("PC",i), ylim=c(-1, 1))
      # points(1:length(pcaboot$InitialSVD$v[,i]),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(MedEig), maxim, cnames, pos=3, cex=0.7, col="blue")
      abline(0,0, lty="dashed", col="red")
    }
    mtext("Box-Plots for the Correlations", outer = TRUE, cex = 1.5)
    par(op)
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]))
      media=apply(bb, 2,mean)
      ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
      expr = errbar(1:length(media), media, ICPerc[1,], ICPerc[2,], add=FALSE, pch=16, xlab="Variables", ylab="Correlations", ylim=c(-1, 1))
      #points(1:length(MedEig),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(media), ICPerc[2,], cnames, pos=3, cex=0.8, col="blue")
      title(main=paste("PC",i))
      abline(0,0, lty="dashed", col="red")
    }
    mtext("5% CI for Correlations based on Percentiles", outer = TRUE, cex = 1.5)
    par(op)
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]))
      media=apply(bb, 2,mean)
      sdev=apply(bb, 2,sd)
      expr = errbar(1:length(media), media, media-1.96*sdev, media+1.96*sdev, add=FALSE, pch=16, xlab="Variables", ylab="Correlations", ylim=c(-1, 1))
      #points(1:length(media),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(media), media+1.96*sdev, cnames, pos=3, cex=0.8, col="blue")
      title(main=paste("PC",i))
      abline(0,0, lty="dashed", col="red")
    }
    mtext("5% CI for Correlations based on Moments", outer = TRUE, cex = 1.5)
    par(op)
    
  }
  
  if (Squared){
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,])^2)
      maxim=apply(bb, 2,max)
      boxplot(bb, main=paste("PC",i), xlab="Variables", ylab="Squared Correlations", ylim=c(0, 1))
      # points(1:length(pcaboot$InitialSVD$v[,i]),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(MedEig), maxim, cnames, pos=3, cex=0.7, col="blue")
      abline(0,0, lty="dashed", col="red")
    }
    mtext("Box-Plots for the Squared Correlations", outer = TRUE, cex = 1.5)
    par(op)
    
    op <- par(mfrow=c(dimens,1),oma = c(0, 0, 2, 0)) 
    for (i in 1:dimens){
      bb=t(as.matrix(pcaboot$Struct[,i,]^2))
      media=apply(bb, 2,mean)
      ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
      expr = errbar(1:length(media), media, ICPerc[1,], ICPerc[2,], add=FALSE, pch=16, xlab="Variables", ylab="Squared Correlations", ylim=c(0, 1))
      #points(1:length(MedEig),pcaboot$InitialSVD$v[,i], cex=1.5, pch=4)
      text(1:length(media), ICPerc[2,], cnames, pos=3, cex=0.8, col="blue")
      title(main=paste("PC",i))
      abline(0,0, lty="dashed", col="red")
    }
    mtext("5% CI for Squared Correlations based on Percentiles", outer = TRUE, cex = 1.5)
    par(op)
  }
  
  if (Scores){
    bb1=t(as.matrix(pcaboot$Scores[,1,]))
    media1=apply(bb1, 2,mean)
    ICPerc1=apply(bb1, 2,quantile, c(0.025, 0.975))
    bb2=t(as.matrix(pcaboot$Scores[,2,]))
    media2=apply(bb2, 2,mean)
    ICPerc2=apply(bb2, 2,quantile, c(0.025, 0.975))
    plot(media1, media2, col=ColorInd, cex=0.7, pch=16, xlab = "PC1", ylab="PC2")
    text(media1, media2, rnames, pos=3, col=ColorInd, cex=0.5)
    points(pcaboot$InitScores[,1], pcaboot$InitScores[,2], col=ColorInd, cex=0.7, pch=5)
    text(pcaboot$InitScores[,1], pcaboot$InitScores[,2], paste("i",rnames, sep=""), pos=3, col=ColorInd, cex=0.5)
    title(main="Row Scores")
    n=dim(pcaboot$InitialData)[1]
    
    if (TypeScores=="el")
    for (i in 1:n){
      dat=t(as.matrix(pcaboot$Scores[i,,]))
      E=ConcEllipse(dat, confidence=0.95)
      plot(E, col=ColorInd[i])
    }
    
    if (TypeScores=="ch")
    for (i in 1:n){
      dat=t(as.matrix(pcaboot$Scores[i,,]))
      F=Fraction(dat, confidence=0.95)
      plot(F, col=ColorInd[i])
    }
    
    if (TypeScores=="st")
      for (i in 1:n){
        dat=t(as.matrix(pcaboot$Scores[i,,]))
        F=Fraction(dat, confidence=0.95)
        plot(F, col=ColorInd[i], type="st")
      }
    
  }
  
}

summary.PCA.Bootstrap <- function(pcaboot){
  
  cat(" ###### Bootstrap for Principal Components Analysis #######\n\n")
  rnames=rownames(pcaboot$InitialData)
  cnames=colnames(pcaboot$InitialData)
  cat("Transformation of the raw data:\n")
  print(pcaboot$InitTransform)
  
  cat("\n\nEigenvalues\n")
  MedEig=apply(pcaboot$EigVal^2, 2,mean)
  ICPerc=apply(pcaboot$EigVal^2, 2,quantile, c(0.025, 0.975))
  sdev=apply(pcaboot$EigVal^2, 2,sd)
  m1=cbind(round(pcaboot$InitialSVD$d^2,3), round(MedEig, 3), round(t(ICPerc), 3), round(MedEig-1.96*sdev, 3), round(MedEig+1.96*sdev, 3))
  rownames(m1)=paste("Dim", 1:length(MedEig))
  colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
  print(m1)
  
  cat("\n\nAccounted Variance\n")
  MedEig=apply(pcaboot$Inertia, 2,mean)
  ICPerc=apply(pcaboot$Inertia, 2,quantile, c(0.025, 0.975))
  sdev=apply(pcaboot$Inertia, 2,sd)
  Inertia=100*pcaboot$InitialSVD$d^2/sum(pcaboot$InitialSVD$d^2)
  m1=cbind(round(Inertia,3), round(MedEig, 3), round(t(ICPerc), 3), round(MedEig-1.96*sdev, 3), round(MedEig+1.96*sdev, 3))
  rownames(m1)=paste("Dim", 1:length(MedEig))
  colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
  print(m1)
  
  cat("\n\nEigenvector Coefficients\n")
  dimens=dim(pcaboot$InitialSVD$u)[2]
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(pcaboot$Vs[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(pcaboot$InitialSVD$v[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
  
  cat("\n\nCorrelations with the components\n")

  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(pcaboot$Struct[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(pcaboot$InitCorr[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
    
  
  cat("\n\nSquared Correlations with the components\n")
 
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(pcaboot$Struct[,i,]^2))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))

    m1=cbind(round(pcaboot$InitCorr[,i]^2,3), round(media, 3), round(t(ICPerc), 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5")
    print(m1)
  }
  
  cat("\n\nRow Scores\n")
  
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(pcaboot$Scores[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(pcaboot$InitScores[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
  
}


