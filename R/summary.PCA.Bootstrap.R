summary.PCA.Bootstrap <- function(object, ...){
  
  cat(" ###### Bootstrap for Principal Components Analysis #######\n\n")
  rnames=rownames(object$InitialData)
  cnames=colnames(object$InitialData)
  cat("Transformation of the raw data:\n")
  print(object$InitTransform)
  
  cat("\n\nEigenvalues\n")
  MedEig=apply(object$EigVal^2, 2,mean)
  ICPerc=apply(object$EigVal^2, 2,quantile, c(0.025, 0.975))
  sdev=apply(object$EigVal^2, 2,sd)
  m1=cbind(round(object$InitialSVD$d^2,3), round(MedEig, 3), round(t(ICPerc), 3), round(MedEig-1.96*sdev, 3), round(MedEig+1.96*sdev, 3))
  rownames(m1)=paste("Dim", 1:length(MedEig))
  colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
  print(m1)
  
  cat("\n\nAccounted Variance\n")
  MedEig=apply(object$Inertia, 2,mean)
  ICPerc=apply(object$Inertia, 2,quantile, c(0.025, 0.975))
  sdev=apply(object$Inertia, 2,sd)
  Inertia=100*object$InitialSVD$d^2/sum(object$InitialSVD$d^2)
  m1=cbind(round(Inertia,3), round(MedEig, 3), round(t(ICPerc), 3), round(MedEig-1.96*sdev, 3), round(MedEig+1.96*sdev, 3))
  rownames(m1)=paste("Dim", 1:length(MedEig))
  colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
  print(m1)
  
  cat("\n\nEigenvector Coefficients\n")
  dimens=dim(object$InitialSVD$u)[2]
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(object$Vs[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(object$InitialSVD$v[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
  
  cat("\n\nCorrelations with the components\n")
  
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(object$Struct[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(object$InitCorr[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
  
  
  cat("\n\nSquared Correlations with the components\n")
  
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(object$Struct[,i,]^2))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    
    m1=cbind(round(object$InitCorr[,i]^2,3), round(media, 3), round(t(ICPerc), 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5")
    print(m1)
  }
  
  cat("\n\nRow Scores\n")
  
  for (i in 1:dimens){
    cat("\nPrincipal Component :",i,"\n")
    bb=t(as.matrix(object$Scores[,i,]))
    media=apply(bb, 2,mean)
    ICPerc=apply(bb, 2,quantile, c(0.025, 0.975))
    sdev=apply(bb, 2,sd)
    m1=cbind(round(object$InitScores[,i],3), round(media, 3), round(t(ICPerc), 3), round(media-1.96*sdev, 3), round(media+1.96*sdev, 3))
    colnames(m1)=c("Initial", "Bootstrap Mean", "CI- P2.5", "CI- P97.5", "CI- M EI", "CI- M ES")
    print(m1)
  }
  
}
