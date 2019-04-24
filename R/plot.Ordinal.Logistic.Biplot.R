plot.Ordinal.Logistic.Biplot <- function(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0, PlotVars = TRUE, PlotInd = TRUE, LabelVars = TRUE, LabelInd = TRUE, mode = "a", 
                                         CexInd = NULL, CexVar = NULL, ColorInd = NULL, ColorVar = NULL, SmartLabels = TRUE, MinQualityVars = 0, dp = 0, PredPoints=0, PlotAxis = FALSE, TypeScale = "Complete", 
                                         ValuesScale = "Original", SizeQualInd = FALSE, SizeQualVars = FALSE, ColorQualInd = FALSE, ColorQualVars = FALSE, PchInd = NULL, PchVar = NULL, 
                                         PlotClus = FALSE, TypeClus="ch", ClustConf=1, ClustCenters=FALSE, UseClusterColors=TRUE, ClustLegend=TRUE, ClustLegendPos="topright", TextVarPos =1, 
                                         PlotSupVars=FALSE, ...){
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
  A = x$RowCoordinates[, c(A1, A2)]
  B = x$ColCoordinates[, c(A1, A2)]
  thresholds=x$ColumnParameters$thresholds
  # B=diag(1/sqrt(apply(B^2,1,sum))) %*% B
  rownames(B)=rownames(x$ColCoordinates)
  n = dim(A)[1]
  p = dim(x$ColumnParameters$coefficients)[1]
  
  qlrcols = x$ColContributions[, A1] + x$ColContributions[, A2]
  qlrrows = x$RowContributions[, A1] + x$RowContributions[, A2]
  
  # Determining sizes and colors of the points
  if (is.null(CexInd)) 
    CexInd = 0.5
  
  if (is.null(CexVar)) 
    CexVar = rep(0.8, p)
  else if (length(CexVar == 1)) 
    CexVar = rep(CexVar, p)
  
  if (is.null(PchInd)) 
    PchInd = 1
  if (is.null(PchVar)) 
    PchVar = rep(16, p)
  
  if (is.null(ColorInd)) 
    if (is.null(x$ColorInd)) 
      ColorInd = rep("red",n)
  
  if (length(ColorInd)==1) ColorInd = rep(ColorInd,n)
  
  if (is.null(ColorVar)) 
    ColorVar = rep("black", p)
  if (length(ColorVar)==1) ColorVar = rep(ColorVar,p)

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
  
  if (margin < 0)
    margin = 0
  
  P=A
  
  xmin = min(P[, 1])
  xmax = max(P[, 1])
  ymin = min(P[, 2])
  ymax = max(P[, 2])
  if (xmax <0 ) xmax=xmax*(-1)
  
  P = rbind(P, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))
  plot(P[, 1], P[, 2], cex = 0, asp = 1, xlab = paste("Axis", A1), ylab = paste("Axis", A2), xaxt = xaxt, yaxt = yaxt, ...)
  title(main = x$Title, omi = c(0, 0, 0, 0))
  
  if (PlotClus) {
    ColorInd=PlotBiplotClusters(A, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, Legend=ClustLegend, LegendPos=ClustLegendPos,   ...)
    if (x$ClusterType=="gm"){
      ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
      if (UseClusterColors) ColorInd = ColorInd2
      PchInd=16
    }
  }
  
  if (PlotInd) 
    points(A[, 1], A[, 2], cex = CexInd, col = ColorInd, pch = PchInd, ...)
  
  if (LabelInd) 
    if (SmartLabels) 
      textsmart(cbind(A[, 1], A[, 2]), CexPoints = CexInd, ColorPoints = ColorInd, ...)
  else text(A[, 1], A[, 2], rownames(A), cex = CexInd, col = ColorInd, pos = 1, ...)
  
  if (PlotVars) {

    for (j in 1:p) if (qlrcols[j] > MinQualityVars) {
      
      OrdVarBiplot(B[j, 1], B[j, 2], thresholds[j,1:(x$Ncats[j]-1)], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = rownames(B)[j], mode = mode, CexPoint = CexVar[j], Color = ColorVar[j], 
                  PchPoint = PchVar[j], textpos=TextVarPos,   ...)}

    for (idp in dp)
      if ((idp > 0) & (idp < (p + 1))) {
        g = B[idp, ]
        nn = (t(g) %*% g)
        scal <- (A %*% g)/nn[1, 1]
        Dscal <- diag(as.vector(scal))
        Fpr <- Dscal %*% matrix(rep(1, nrow(A)), ncol = 1) %*% t(g)
        nrFpr <- nrow(Fpr)
        
        dlines(A, Fpr, color=ColorVar[idp])
      }
    
    for (idp in PredPoints)
      if ((idp > 0) & (idp < (n + 1)))
        for (j in 1:p){
          g = B[j, ]
          nn = (t(g) %*% g)
          scal <- (A[idp,] %*% g)/nn[1, 1]
          Fpr <- scal %*% t(g)
          nrFpr <- nrow(Fpr)
          dlines(matrix(A[idp,],1,2) , Fpr, color=ColorVar[j])
        }
  }
  
  if (PlotSupVars) 
    plot.Supplementary.Variables(x, F1=A1, F2=A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode)
  
  
}