
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