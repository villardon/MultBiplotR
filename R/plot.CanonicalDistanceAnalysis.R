plot.CanonicalDistanceAnalysis<- function(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE, PlotVars = TRUE, PlotInd = TRUE, 
                                          LabelInd=TRUE, CexGroup=1, PchGroup=16, margin=0.1, AddLegend=TRUE, ShowAxes=FALSE, LabelAxes=TRUE, LabelGroups=TRUE,
                                          PlotCircle = TRUE, ConvexHulls = FALSE, TypeCircle = "M", ColorGroups = NULL, ColorVars = NULL, LegendPos="topright",
                                          ColorInd = NULL, voronoi = TRUE, mode="a", TypeScale = "Complete", ValuesScale = "Original",
                                          MinQualityVars = 0, dpg = 0, dpi=0,  PredPoints=0, PlotAxis = FALSE, CexInd = NULL, CexVar = NULL,
                                          PchInd = NULL, PchVar = NULL, ColorVar=NULL, ShowAxis=FALSE, VoronoiColor="black", ShowBox=TRUE, ...){
  
  
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
    PchInd = rep(1,n)
  else if (length(PchInd == 1)) 
    PchInd = rep(PchInd, p)
}