# Plots a biplot for continuous data

plot3d.ContinuousBiplot <- function(x, A1 = 1, A2 = 2, A3=3, ShowAxis = TRUE, margin = 0, PlotVars = TRUE, PlotInd = TRUE, WhatInds=NULL, WhatVars=NULL, LabelVars = TRUE, 
                                  LabelInd = TRUE, IndLabels=NULL, VarLabels=NULL, mode = "a", CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL, LabelPos=1,
                                  SmartLabels = FALSE, MinQualityInds = 0, MinQualityVars = 0.3, dp = 0, PredPoints=0, PlotAxis = FALSE, TypeScale = "Complete", 
                                  ValuesScale = "Original", SizeQualInd = FALSE, SizeQualVars = FALSE, ColorQualInd = FALSE, ColorQualVars = FALSE, PchInd = NULL, 
                                  PchVar = NULL, PlotClus = FALSE, TypeClus="ch", ClustConf=1, ClustCenters=FALSE, UseClusterColors=TRUE, PlotSupVars=TRUE,  ...) {
  
  modes=c("p", "a", "b", "h", "ah", "s")
  if (is.numeric(mode)) 
    mode = modes[mode]
  TypeScales=c("Complete", "StdDev", "BoxPlot")
  if (is.numeric(TypeScale)) 
    TypeScale = TypeScales[TypeScale]
  ValuesScales=c("Original", "Transformed")
  if (is.numeric(ValuesScale)) 
    ValuesScale = ValuesScales[ValuesScale]
  
  # Obtaining coordinates and qualities for the representation
  A = x$RowCoordinates[, c(A1, A2, A3)]
  B = x$ColCoordinates[, c(A1, A2, A3)]
  n = dim(A)[1]
  p = dim(B)[1]
  
  if (is.null(IndLabels)) IndLabels=rownames(A)
  if (is.null(VarLabels)) VarLabels=rownames(B)
  
  # Determining what rows to plot
  if (is.null(WhatInds)) 
    WhatInds = matrix(1, n, 1)
  else
    if (!CheckBinaryVector(WhatInds)) {
      AllRows = matrix(0, n, 1)
      AllRows[WhatInds]=1
      WhatInds=AllRows
    }
  WhatInds=as.logical(WhatInds)
  
  # Determining what columns to plot
  if (is.null(WhatVars)) 
    WhatVars = matrix(1, p, 1)
  else
    if (!CheckBinaryVector(WhatVars)){
      AllCols = matrix(0, p, 1)
      AllCols[WhatVars]=1
      WhatVars=AllCols
    }
  WhatVars=as.logical(WhatVars)   
  
  qlrcols = x$ColContributions[, A1] + x$ColContributions[, A2] + x$ColContributions[, A3]
  qlrrows = x$RowContributions[, A1] + x$RowContributions[, A2] + x$RowContributions[, A3]
  
  WhatInds = WhatInds & (qlrrows>MinQualityInds*100)
  WhatVars = WhatVars & (qlrcols>MinQualityVars*100)
  
  # Determining sizes and colors of the points
  if (is.null(CexInd)) 
    CexInd = rep(0.8, n)
  else if (length(CexInd == 1)) 
    CexInd = rep(CexInd, p)
  
  if (is.null(CexVar)) 
    CexVar = rep(0.8, p)
  else if (length(CexVar == 1)) 
    CexVar = rep(CexVar, p)
  
  if (is.null(PchInd)) 
    PchInd = rep(1,n)
  else if (length(PchInd == 1)) 
    PchInd = rep(PchInd, p)
  
  if (is.null(PchVar)) 
    PchVar = rep(16, p)
  else if (length(PchVar == 1)) 
    PchVar = rep(PchVar, p)
  
  if (is.null(ColorInd)) 
    if (is.null(x$ColorInd)) 
      ColorInd = rep("red",n)
  else ColorInd = x$ColorInd
  if (length(ColorInd)==1) ColorInd = rep(ColorInd,n)
  
  if (is.null(ColorVar)) 
    ColorVar = rep("black", p)
  if (SizeQualInd) 
    CexInd = cscale(qlrrows, rescale_pal())
  if (SizeQualVars) 
    CexVar = cscale(qlrcols, rescale_pal())
  
  if (ColorQualInd) 
    ColorInd = cscale(qlrrows, seq_gradient_pal("white", "red"))
  if (ColorQualVars) 
    ColorVar = cscale(qlrcols, seq_gradient_pal("white", "blue"))
  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }
  
  if ((margin < 0) | (margin > 0.3)) 
    margin = 0
  
  if (PlotVars & PlotInd){
    P = rbind(A, B)
  }
  else
    if (PlotVars) P=B
  else P=A
  
  xmin = min(P[, 1])
  xmax = max(P[, 1])
  ymin = min(P[, 2])
  ymax = max(P[, 2])
  zmin = min(P[, 3])
  zmax = max(P[, 3])
  xrang=abs(xmax-xmin)
  yrang=abs(ymax-ymin)
  zrang=abs(zmax-zmin)
  fact=0.03
  xmin=xmin-fact*xrang
  xmax=xmax+fact*xrang
  ymin=ymin-fact*yrang
  ymax=ymax+fact*yrang
  zmin=zmin-fact*zrang
  zmax=zmax+fact*zrang
  
  if (xmax <0 ) xmax=xmax*(-1)
  
#   P = rbind(P, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
#   P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))

  plot3d(0,0,0, cex = 0, asp = 1, xlab = paste("Dimensions", A1,"-",A2,"-",A3), ylab = "", main = x$Title, box=TRUE, cex.axis=0.5, ...)
  title(main = x$Title, omi = c(0, 0, 0, 0))

#   if (PlotClus) {
#     ColorInd=PlotBiplotClusters(A, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, ...)
#     if (x$ClusterType=="gm"){
#       ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
#       if (UseClusterColors) ColorInd = ColorInd2
#       PchInd=rep(16,n)
#       
#     }
#   }
  
  if (PlotInd) 
    for (i in 1:n)
      if (WhatInds[i]){
        points3d(A[i, 1], A[i, 2], A[i, 3], col = ColorInd[i], cex=CexInd[i], pch = PchInd[i], ...)
      }
  
  if (LabelInd)
    for (i in 1:n)
      if (WhatInds[i]) text3d(A[i, 1], A[i, 2], A[i, 3], IndLabels[i], cex = CexInd[i], col = ColorInd[i], pos = LabelPos, ...)
  

if (PlotVars) {
  HH=matrix(0,2*p,3)
  for (i in 1:p)
    HH[2*i,]=B[i,]
  
  segments3d(HH, cex = CexVar, pch = PchVar, col = ColorVar)
  text3d(B[, 1], B[, 2], B[, 3], rownames(B), adj = 0, col = ColorVar, cex = CexVar)
}

M=matrix(c(xmax, 0,0, 0, ymax, 0, 0, 0, zmax), 3,3)
xlabel=paste("A",A1)
ylabel=paste("A",A2)
zlabel=paste("A",A3)
text3d(M, texts=c(xlabel, ylabel, zlabel))


#   if (PlotVars) {
#     if (mode=="s")
#       Scales = GetBiplotScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)
#     
#     for (j in 1:p) if (WhatVars[j]) 
#       VarBiplot(B[j, 1], B[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = VarLabels[j], mode = mode, CexPoint = CexVar[j], Color = ColorVar[j], 
#                 ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j], ...)
#     
#     for (idp in dp)
#       if ((idp > 0) & (idp < (p + 1))) {
#         g = B[idp, ]
#         nn = (t(g) %*% g)
#         scal <- (A %*% g)/nn[1, 1]
#         Dscal <- diag(as.vector(scal))
#         Fpr <- Dscal %*% matrix(rep(1, nrow(A)), ncol = 1) %*% t(g)
#         nrFpr <- nrow(Fpr)
#         dlines(A, Fpr, color=ColorVar[idp])
#       }
#     
#     for (idp in PredPoints)
#       if ((idp > 0) & (idp < (n + 1)))
#         for (j in 1:p){
#           g = B[j, ]
#           nn = (t(g) %*% g)
#           scal <- (A[idp,] %*% g)/nn[1, 1]
#           Fpr <- scal %*% t(g)
#           nrFpr <- nrow(Fpr)
#           dlines(matrix(A[idp,],1,2) , Fpr, color=ColorVar[j])
#         }
#   }
#   
#   if (PlotSupVars) 
#     plot.Supplementary.Variables(x, F1=A1, F2=A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode, TypeScale=TypeScale)

if (ShowAxis) abclines3d(0, 0, 0, a = diag(3), col = "gray")  

}


