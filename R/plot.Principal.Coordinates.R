# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Marzo/2021

plot.Principal.Coordinates <- function(x, A1 = 1, A2 = 2, LabelRows=TRUE, WhatRows = NULL, RowCex=1, RowPch=16, Title="",
                                       RowLabels = NULL, RowColors = NULL, ColColors=NULL, ColLabels=NULL, SizeQualInd = FALSE, SmartLabels = TRUE, 
                                       ColorQualInd = FALSE, ColorQual="black", PlotSup=TRUE, Bootstrap=FALSE, 
                                       BootstrapPlot=c("Ellipse", "CovexHull", "Star"), margin=0, 
                                       PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, CexClustCenters=1, LegendClust=TRUE,
                                       ClustCenters = FALSE, UseClusterColors = TRUE, ShowAxis=FALSE, PlotBinaryMeans=FALSE,
                                       MinIncidence=0, ShowBox=FALSE,  ColorSupContVars=NULL, ColorSupBinVars=NULL, ColorSupOrdVars=NULL,
                                       TypeScale = "Complete", SupMode="s", PlotSupVars=FALSE, ...){
  
  if (is.null(ColColors)) 
    ColorVar = "black"
  if (is.null(ColLabels)) 
    ColLabels=colnames(x$Data)
  
  if (length(BootstrapPlot)>1) BootstrapPlot = BootstrapPlot[1]
  
  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }
  
  a = x$RowCoordinates[,c(A1,A2)]
  n = dim(a)[1]
  if (is.null(WhatRows)) WhatRows = matrix(1, n, 1)
  if (is.null(RowColors)) RowColors = matrix(1, n, 1)
  if (is.null(RowLabels)) RowLabels =rownames(a)
  
  #op=par(mai=c(0,0,0.5,0))
  #op=par(mar=c(1, 1, 1, 1) + 0.1)
  
  xmin = min(a[, 1])
  xmax = max(a[, 1])
  ymin = min(a[, 2])
  ymax = max(a[, 2])
  xrang=abs(xmax-xmin)
  yrang=abs(ymax-ymin)
 
  if (x$Analysis == "Principal Coordinates"){
  xlab=paste("(Dim", A1,"(",round(x$Inertia[A1], digits=1),"%)")
  ylab=paste("(Dim",A2,"(",round(x$Inertia[A2], digits=1),"%)")}
  else{
    xlab=paste("Dim", A1)
    ylab=paste("Dim", A2)
  }
  
  P = rbind(a, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))
  plot(P[, 1], P[, 2], cex = 0, asp = 1, xaxt = xaxt, yaxt = yaxt, xlab="",ylab="", bty="n", ...)
  
  # plot(a[, 1], a[, 2], cex = 0, asp = 1,  xaxt = xaxt, yaxt = yaxt, xlab=xlab,ylab=ylab, bty="n", ...)
  
  if (x$Analysis == "Principal Coordinates"){
    title(main = paste("Principal Coordinates", Title, sep="   "), omi = c(0, 0, 0, 0))}
  else{
    title(main = paste("MDS", Title, sep="   "), omi = c(0, 0, 0, 0))
  }
  
  if (ShowBox) rect(xmin, ymin, xmax, ymax)
  
  if (PlotClus) {
    RowColors=PlotBiplotClusters(a, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, CexClustCenters=CexClustCenters, Legend=LegendClust, ...)
    if (x$ClusterType=="gm"){
      ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
      if (UseClusterColors) RowColors = ColorInd2
      PchInd=rep(16,n)
    }
  }

  qlrrows = x$RowQualities[, A1] + x$RowQualities[, A2]
  if (ColorQualInd==TRUE) RowColors = cscale(qlrrows, seq_gradient_pal("white", ColorQual))
  
  neje = dim(a)[2]
  points(a[, 1], a[, 2], col = RowColors, cex=WhatRows*RowCex, pch=RowPch, asp=1, ...)
  
  if (PlotSupVars) 
    plot.Supplementary.Variables(x, F1=A1, F2=A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=SupMode, TypeScale=TypeScale, ColorSupContVars=ColorSupContVars, ColorSupBinVars=ColorSupBinVars, ColorSupOrdVars=ColorSupOrdVars )
  
  
  if (LabelRows)  
      if (SmartLabels) 
        textsmart(cbind(a[, 1], a[, 2]), RowLabels, CexPoints = WhatRows*RowCex, ColorPoints = RowColors)
  else text(a[, 1], a[, 2], labels = RowLabels, cex=WhatRows*RowCex , col = RowColors)
  
  if ((PlotBinaryMeans) & (x$TypeData=="Binary") & (!is.null(x$Data))){
    tot=apply(x$Data, 2, sum)
    cc= diag(1/tot) %*% t(x$Data) %*% a
    cc=cc[which(tot>MinIncidence),]
    rownames(cc)=ColLabels[which(tot>MinIncidence)]
    points(cc[, 1], cc[, 2], col = ColColors, pch=16)
    textsmart(cbind(cc[, 1], cc[, 2]), ColorPoints=ColColors, colnames(x$Data), CexPoints=0.7)
  }
  
  if((PlotSup) & !is.null(x$SupRowCoordinates)){
    b = x$SupRowCoordinates[,c(A1,A2)]
    points(b[, 1], b[, 2], col = 2)
    if (LabelRows)  
      if (SmartLabels) 
        textsmart(cbind(b[, 1], b[, 2]), CexPoints = 1, ColorPoints = 2)
    else text(b[, 1], b[, 2],labels=rownames(b), col = 2)
  }

  if ((Bootstrap==TRUE) & (!is.null(x$BootstrapInfo)))
  plot.PCoABootstrap(x$BootstrapInfo, BootstrapPlot=BootstrapPlot , confidence=0.95, Colors=RowColors)
  # par(op)
}
