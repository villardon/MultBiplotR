plot.CanonicalDistanceAnalysis<- function(x, A1 = 1, A2 = 2, ScaleGraph = TRUE,  ShowAxis=FALSE, ShowAxes=FALSE, LabelAxes=TRUE, margin=0.1, PlotAxis = FALSE, ShowBox=TRUE,
                                          PlotGroups = TRUE, LabelGroups=TRUE, CexGroup=1.5, PchGroup=16, ColorGroup = NULL, voronoi = TRUE, VoronoiColor="black",
                                          PlotInd = TRUE, LabelInd=TRUE, CexInd = 0.8, PchInd = 3, ColorInd = NULL,  WhatInds=NULL, IndLabels=NULL,
                                          PlotVars = TRUE, LabelVar=TRUE, CexVar = NULL, PchVar = NULL, ColorVar=NULL, WhatVars=NULL, VarLabels=NULL,
                                          mode="a", TypeScale = "Complete", ValuesScale = "Original", SmartLabels=TRUE,
                                          AddLegend =TRUE, LegendPos="topright", PlotCircle = TRUE, ConvexHulls = FALSE, TypeCircle = "M",
                                          MinQualityVars = 0, dpg = 0, dpi=0,  PredPoints=0, PlotClus = TRUE, TypeClus = "ch", ClustConf = 1, CexClustCenters=1,
                                          ClustCenters = FALSE, UseClusterColors = TRUE, ...){


  if (is.null(x$ClusterType))  x=AddCluster2Biplot(x, ClusterType="us", Groups=x$Groups)

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
  n = dim(A)[1]
  M = x$MeanCoordinates[, c(A1, A2)]
  g = dim(M)[1]

  if (is.null(rownames(A))) rownames(A)=paste("i", 1:n, sep="")

  if (is.null(IndLabels)) IndLabels=rownames(A)

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

  # Determining sizes and colors of the points
  if (is.null(CexInd))
    CexInd = rep(0.8, n)
  else if (length(CexInd == 1))
    CexInd = rep(CexInd, n)

  if (is.null(PchInd))
    PchInd = rep(3,n)
  else if (length(PchInd == 1))
    PchInd = rep(PchInd, n)

  if (is.null(ColorInd))
    ColorInd = as.numeric(x$Groups)
  else if (length(ColorInd == 1))
    ColorInd = rep(ColorInd, n)


  if (is.null(CexGroup))
    CexGroup = rep(1.5, g)
  else if (length(CexGroup == 1))
    CexGroup = rep(CexGroup, g)

  if (is.null(PchGroup))
    PchGroup = rep(16,g)
  else if (length(PchGroup == 1))
    PchGroup = rep(PchGroup, g)

  if (is.null(ColorGroup))
    ColorGroup = 1:g


  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }

  P=M

  xmin = min(P[, 1])
  xmax = max(P[, 1])
  ymin = min(P[, 2])
  ymax = max(P[, 2])
  if (xmax <0 ) xmax=xmax*(-1)

  P = rbind(P, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))
  XLabel=paste("Dimension", A1, " (", round(x$ExplainedVariance[A1], digits=2),"%)", sep="" )
  YLabel=paste("Dimension", A2, " (", round(x$ExplainedVariance[A2], digits=2),"%)", sep="" )

  if (x$Type == "PERMANOVA")
    Main="PERMANOVA : Graphical representation"
  else
    Main=paste("Canonical Distance Analysis (p-value =", round(x$pvalue, digits=5),")")
  plot(P[, 1], P[, 2], cex = 0, asp = 1, xlab = XLabel, ylab = YLabel, xaxt = xaxt, yaxt = yaxt, main=Main, axes=ShowAxes,  ...)

  if (PlotClus) {
    RowColors=PlotBiplotClusters(A, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, CexClustCenters=CexClustCenters, ...)
    if (x$ClusterType=="gm"){
      ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
      if (UseClusterColors) RowColors = ColorInd2
      PchInd=rep(16,n)
    }
  }

  if (PlotInd)
    points(A[, 1], A[, 2], cex = CexInd, col = RowColors, pch = PchInd, ...)

  if (LabelInd)
    if (SmartLabels)
      textsmart(cbind(A[, 1], A[, 2]), CexPoints = CexInd, ColorPoints = RowColors, ...)
  else text(A[, 1], A[, 2], rownames(A), cex = CexInd, col = ColorInd, pos = 1, ...)

  if (PlotGroups)
    points(M[, 1], M[, 2], cex = CexGroup, col = ColorGroup, pch = PchGroup, ...)

  if (LabelGroups)
  if (SmartLabels)
    textsmart(cbind(M[, 1], M[, 2]), CexPoints = CexGroup, ColorPoints = ColorGroup, ...)
  else text(M[, 1], M[, 2], rownames(M), cex = CexGroup, col = ColorGroup, pos = 1, ...)

  if (ConvexHulls) {
    lev = levels(x$Groups)
    for (i in 1:nlevels(x$Groups)) {
      XP = A[which(x$Groups == lev[i]), ]
      XP = cbind(XP[, 1], XP[, 2])
      hpts <- chull(XP)
      hpts <- c(hpts, hpts[1])
      lines(XP[hpts, ], col = ColorGroup[i])
    }
  }


}
