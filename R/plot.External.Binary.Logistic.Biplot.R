plot.External.Binary.Logistic.Biplot <- function(x, F1 = 1, F2 = 2, ShowAxis=FALSE, margin=0.1, PlotVars = TRUE, PlotInd = TRUE, WhatRows = NULL, 
                                                 WhatCols = NULL, LabelRows=TRUE, LabelCols=TRUE, RowLabels = NULL, ColLabels = NULL, 
                                                 RowColors = NULL, ColColors = NULL, Mode = "s", TickLength= 0.01,
                                                 RowCex = 0.8, ColCex = 0.8, SmartLabels = FALSE, MinQualityRows = 0, MinQualityCols = 0, 
                                                 dp = 0, PredPoints=0, SizeQualRows = FALSE, ShowBox=FALSE, SizeQualCols = FALSE, ColorQualRows = FALSE,
                                                 ColorQualCols = FALSE, PchRows = NULL, PchCols = NULL, PlotClus = FALSE, TypeClus = "ch", ClustConf=1,  
                                                 Significant=FALSE, alpha=0.05, Bonferroni=FALSE, PlotSupVars = TRUE, ...) {
  
  a = x$RowCoordinates[,c(F1,F2)]
  n = dim(a)[1]
  neje = dim(a)[2]
  
  p = dim(x$ColumnParameters)[1]
  
  # Determining what rows to plot
  if (is.null(WhatRows)) 
    WhatRows = matrix(1, n, 1)
  else
    if (!CheckBinaryVector(WhatRows)) {
      AllRows = matrix(0, n, 1)
      AllRows[WhatRows]=1
      WhatRows=AllRows
    }

  WhatRows=as.logical(WhatRows)
  
  RowQualities=(x$RowQualities[,F1] +x$RowQualities[,F2])/100

  # WhatRows=WhatRows & (RowQualities>MinQualityRows)
  
  # Determining what columns to plot
  if (is.null(WhatCols)) 
    WhatCols = matrix(1, p, 1)
  else
    if (!CheckBinaryVector(WhatCols)){
      AllCols = matrix(0, p, 1)
      AllCols[WhatCols]=1
      WhatCols=AllCols
    }
  WhatCols=as.logical(WhatCols)     
  
  if (Significant){
    if (Bonferroni)
      WhatCols= WhatCols & (x$VarInfo$Bonferroni<alpha)
    else
      WhatCols= WhatCols & (x$VarInfo$pvalues<alpha)
  }
  

  #WhatCols=WhatCols & (x$VarInfo$Nagelkerke>MinQualityCols)
  
  if (is.null(RowColors)) 
    RowColors = matrix("blue", n, 1)
  if (length(RowColors)==1)
    RowColors = rep(RowColors, n)
  
  if (is.null(ColColors)) 
    ColColors = matrix("black", p, 1)
  
  if (is.null(RowLabels)) 
    RowLabels = rownames(x$RowCoordinates)
  
  if (is.null(ColLabels)) 
    ColLabels = rownames(x$ColumnParameters)
  
  # Determining sizes and colors of the points
  if (length(RowCex == 1)) 
    RowCex = rep(RowCex, n)
  
  if (length(ColCex == 1)) 
    ColCex = rep(ColCex, p)
  
  qlrcols = x$VarInfo$Nagelkerke
  qlrrows = x$RowQualities[, F1] + x$RowQualities[, F2]
  
  if (SizeQualRows) 
    RowCex = cscale(qlrrows, rescale_pal())
  if (SizeQualCols) 
    ColCex = cscale(qlrcols, rescale_pal())
  
  if (ColorQualRows) 
    RowColors = cscale(qlrrows, seq_gradient_pal("white", "red"))
  if (ColorQualCols) 
    ColColors = cscale(qlrcols, seq_gradient_pal("white", "blue"))
  
  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }
  
  if ((margin < 0) | (margin > 0.3)) 
    margin = 0
  
  xmin = min(a[, 1])
  xmax = max(a[, 1])
  ymin = min(a[, 2])
  ymax = max(a[, 2])
  P = rbind(a, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))
  
  
  
  
  op=par(mai=c(0,0,0.5,0))
  op=par(mar=c(1, 1, 1, 1) + 0.1)
  

  plot(P[, 1], P[, 2], asp=1, xaxt = xaxt, yaxt = yaxt, cex=0, bty="n", xlab=paste("Dimension",F1), ylab=paste("Dimension",F2), main="External Logistic Biplot", ...)
  if (ShowBox) rect(xmin, ymin, xmax, ymax)
  if (PlotClus) {
    RowColors=PlotBiplotClusters(a, x$Clusters, TypeClus = "ch", ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = TRUE, ClustConf=ClustConf, ...)
  }

  if (PlotInd){
    for (i in 1:n)
      if (WhatRows[i]){
        points(a[i, 1], a[i, 2], col = RowColors[i], cex=RowCex[i], pch = PchRows, ...)
        if (LabelRows)
        text(a[i, 1], a[i, 2], labels = RowLabels[i], col = RowColors[i], pos=1, cex=RowCex[i], ...)
      }
  }
  

  if (PlotVars){
    for (i in 1:p)
      if (WhatCols[i]){
        PlotBinaryVar(b0=x$ColumnParameters[i,1], bi1=x$ColumnParameters[i,F1+1], bi2=x$ColumnParameters[i,F2+1], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                      mode=Mode, Color = ColColors[i], label=ColLabels[i], tl=TickLength, CexPoint=ColCex[i])
      }
  }
  
  for (idp in dp){
    if ((idp > 0) & (idp < (p + 1))) 
      if (WhatCols[idp]){
        g = matrix(c(x$ColumnParameters[idp,F1+1], x$ColumnParameters[idp,F2+1]),2,1)
        nn = (t(g) %*% g)
        scal <- (a %*% g)/nn[1, 1]
        Dscal <- diag(as.vector(scal))
        Fpr <- Dscal %*% matrix(rep(1, nrow(a)), ncol = 1) %*% t(g)
        nrFpr <- nrow(Fpr)
        
        dlines(a , Fpr, color=ColColors[idp])
      }
    }
  
  for (idp in PredPoints){
    if ((idp > 0) & (idp < (n + 1))){
      for (j in 1:p){
        if (WhatCols[j]){
          g = matrix(c(x$ColumnParameters[j,F1+1], x$ColumnParameters[j,F2+1]),2,1)
          nn = (t(g) %*% g)
          scal <- (a[idp,] %*% g)/nn[1, 1]
          
          Fpr <- scal %*% t(g)
          nrFpr <- nrow(Fpr)
          dlines(matrix(a[idp,],1,2) , Fpr, color=ColColors[j])
        }
      }
    }
  }
  
  if (PlotSupVars) 
    plot.Supplementary.Variables(x, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=Mode)
  
par(op)
}



