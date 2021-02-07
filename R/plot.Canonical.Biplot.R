plot.Canonical.Biplot <- function(x, A1 = 1, A2 = 2, ScaleGraph = TRUE, PlotGroups = TRUE, PlotVars = TRUE, PlotInd = TRUE, WhatInds=NULL,
                                  WhatVars=NULL, WhatGroups=NULL, IndLabels=NULL, VarLabels=NULL, GroupLabels=NULL, AbbreviateLabels=FALSE,
                                  LabelInd=TRUE, LabelVars = TRUE, CexGroup=1, PchGroup=16, margin=0.1, AddLegend=FALSE, ShowAxes=FALSE, LabelAxes=FALSE, LabelGroups=TRUE,
                                  PlotCircle = TRUE, ConvexHulls = FALSE, TypeCircle = "M", ColorGroups = NULL, ColorVars = NULL, LegendPos="topright",
                                  ColorInd = NULL, voronoi = TRUE, mode="a", TypeScale = "Complete", ValuesScale = "Original",
                                  MinQualityVars = 0, dpg = 0, dpi=0, dp=0,  PredPoints=0, PlotAxis = FALSE, CexInd = NULL, CexVar = NULL,
                                  PchInd = NULL, PchVar = NULL, ColorVar=NULL, ShowAxis=FALSE, VoronoiColor="black", ShowBox=FALSE, ShowTitle=TRUE,
                                  PlotClus = FALSE, TypeClus = "ch", ClustConf = 1, ClustCenters = FALSE,  UseClusterColors = TRUE, CexClustCenters=1, ...) {


  modes=c("p", "a", "b", "h", "ah", "s")
  if (is.numeric(mode))
    mode = modes[mode]
  TypeScales=c("Complete", "StdDev", "BoxPlot")
  if (is.numeric(TypeScale))
    TypeScale = TypeScales[TypeScale]
  ValuesScales=c("Original", "Transformed")
  if (is.numeric(ValuesScale))
    ValuesScale = ValuesScales[ValuesScale]

  if (is.null(ColorGroups))
    ColorGroups = x$ClusterColors

  if (is.null(ColorInd))
    ColorInd = x$ClusterColors[as.integer(x$groups)]


  if (LabelAxes){
    xlabel = paste("CV", A1, "(", round(x$Inertia[A1], digits=2),"%)")
    ylabel = paste("CV", A2, "(", round(x$Inertia[A2], digits=2),"%)")
    Title="Canonical/MANOVA Biplot"}
  else{
    xlabel = ""
    ylabel = ""
     
    Title=paste(x$Title, " / ", A1, "-", A2, " (", round(x$Inertia[A1]+x$Inertia[A2], digits=2),"%)")
  }


  scf = 1
  if (ScaleGraph) {
    sca = sum(x$GroupCoordinates^2)
    scb = sum(x$ColCoordinates^2)
    sca = sca/x$g
    scb = scb/x$p
    scf = sqrt(sqrt(scb/sca))
  }

  J = x$GroupCoordinates[, c(A1, A2)] * scf
  H = x$ColCoordinates[, c(A1, A2)]/scf
  V = x$RowCoordinates[, c(A1, A2)] * scf

  g = dim(J)[1]
  p = dim(H)[1]
  n = dim(V)[1]

  if (is.null(GroupLabels)) GroupLabels=rownames(J)
  if (is.null(IndLabels)) IndLabels=rownames(V)
  if (is.null(VarLabels)) VarLabels=rownames(H)

  if (AbbreviateLabels){
    GroupLabels=abbreviate(GroupLabels, minlength = 5L)
    IndLabels=abbreviate(IndLabels, minlength = 5L)
    VarLabels=abbreviate(VarLabels, minlength = 5L)
  }

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

  # Determining what groups to plot
  if (is.null(WhatGroups))
    WhatGroups = matrix(1, p, 1)
  else
    if (!CheckBinaryVector(WhatGroups)){
      AllCols = matrix(0, p, 1)
      AllCols[WhatGroups]=1
      WhatGroups=AllCols
    }
  WhatVars=as.logical(WhatVars)

  if (is.null(CexInd))
    CexInd = 0.8

  if (is.null(CexVar))
    CexVar = rep(0.8, x$p)
  else if (length(CexVar == 1))
    CexVar = rep(CexVar, p)

  if (is.null(PchInd))
    PchInd = 1
  if (is.null(PchVar))
    PchVar = rep(16, p)

  if (is.null(ColorVar))
    ColorVar = rep("black", p)

  if (!LabelVars) VarLabels=rep(" ", p)

  if (!is.null(x$Sup_Individual_Coord)){VS=x$Sup_Individual_Coord*scf}

  if (ShowAxis) {
    xaxt = "s"
    yaxt = "s"
  } else {
    xaxt = "n"
    yaxt = "n"
  }

  if (PlotInd | PlotClus)
  PP = rbind(J, H, V)
  else
  PP = rbind(J, H)

  xmin = min(PP[, 1])
  xmax = max(PP[, 1])
  ymin = min(PP[, 2])
  ymax = max(PP[, 2])
  xrang=abs(xmax-xmin)
  yrang=abs(ymax-ymin)
  fact=0.03
  xmin=xmin-fact*xrang
  xmax=xmax+fact*xrang
  ymin=ymin-fact*yrang
  ymax=ymax+fact*yrang

  if (xmax <0 ) xmax=xmax*(-1)

  PP = rbind(PP, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
  PP = rbind(PP, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))

  plot(PP[, 1], PP[, 2], cex = 0, asp = 1, main = Title, xlab = xlabel, ylab = ylabel, xaxt = xaxt, yaxt = yaxt, axes=ShowAxes, ...)
  #op=par(mai=c(0,0,0.5,0))
  #op=par(mar=c(1, 1, 1, 1) + 0.1)

  if (ShowBox) rect(xmin, ymin, xmax, ymax)

  if (PlotClus) {
    ColorInd=PlotBiplotClusters(V, x$Clusters, TypeClus = TypeClus, ClusterColors = x$ClusterColors, ClusterNames=x$ClusterNames, centers = ClustCenters, ClustConf=ClustConf, CexClustCenters=CexClustCenters, ...)
    if (x$ClusterType=="gm"){
      ColorInd2=rgb((x$P %*% t(col2rgb(x$ClusterColors)))/255)
      if (UseClusterColors) ColorInd = ColorInd2
      PchInd=rep(16,n)
    }
  }

  if (PlotInd) {
    points(V[, 1], V[, 2], col = ColorInd, pch=PchInd, cex=CexInd)
    if (LabelInd)
      text(V[, 1], V[, 2], IndLabels, adj = 0, col = ColorInd, cex = 0.8)
    if (!is.null(x$Sup_Individual_Coord)){
      points(VS[, A1], VS[, A2], cex = 0.5, col = "Black")
      text(VS[, A1], VS[, A2], paste(rownames(VS),"*",sep=""), adj = 0, col = "Black", cex = 0.8)
    }
  }

  if (PlotGroups) {
    points(J[, 1], J[, 2], cex = CexGroup, pch = PchGroup, col = ColorGroups)
    if (LabelGroups)
      text(J[, 1], J[, 2], GroupLabels, adj = 0, col = ColorGroups, cex = CexGroup)
  }

  if (AddLegend) legend(LegendPos, rownames(J), fill = ColorGroups)

  if (PlotVars) {
    qlrcols = x$ColContributions[, A1] + x$ColContributions[, A2]
    if (mode=="s")
      Scales = GetBiplotScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)

    for (j in 1:p) if (qlrcols[j] > MinQualityVars)
      VarBiplot(H[j, 1], H[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = VarLabels[j], mode = mode, CexPoint = CexVar[j], Color = ColorVar[j],
                ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j], ...)

    for (idp in dpg)
      if ((idp > 0) & (idp < (p + 1))) {
        g = H[idp, ]
        nn = (t(g) %*% g)
        scal <- (J %*% g)/nn[1, 1]
        Dscal <- diag(as.vector(scal))
        Fpr <- Dscal %*% matrix(rep(1, nrow(J)), ncol = 1) %*% t(g)
        nrFpr <- nrow(Fpr)

        dlines(J, Fpr, color=ColorGroups)
      }

    if (PlotInd){
      for (idp in dp)
        if ((idp > 0) & (idp < (p + 1))) {
          g = H[idp, ]
          nn = (t(g) %*% g)
          scal <- (V %*% g)/nn[1, 1]
          Dscal <- diag(as.vector(scal))
          Fpr <- Dscal %*% matrix(rep(1, nrow(V)), ncol = 1) %*% t(g)
          nrFpr <- nrow(Fpr)

          dlines(V, Fpr, color=ColorInd)
        }
    }

    for (idp in PredPoints)
      if ((idp > 0) & (idp < (n + 1)))
        for (j in 1:p){
          g = H[j, ]
          nn = (t(g) %*% g)
          scal <- (J[idp,] %*% g)/nn[1, 1]
          Fpr <- scal %*% t(g)
          nrFpr <- nrow(Fpr)
          dlines(matrix(J[idp,],1,2) , Fpr, color=ColorVar[j])
        }

  }


  if (PlotCircle) {
    switch(TypeCircle, U = {
      radius = x$UnivRad
    }, B = {
      radius = x$BonfRad
    }, M = {
      radius = x$MultRad
    }, C = {
      radius = x$ChisRad
    })
    radius=radius*scf;

    for (i in 1:x$g) {
      Circle(radius[i], c(J[i, 1], J[i, 2]), col = ColorGroups[i])}
  }

  if (ConvexHulls) {
    lev = levels(x$groups)
    for (i in 1:nlevels(x$groups)) {
      XP = V[which(x$groups == lev[i]), ]
      XP = cbind(XP[, 1], XP[, 2])
      hpts <- chull(XP)
      hpts <- c(hpts, hpts[1])
      lines(XP[hpts, ], col = ColorGroups[i])
    }
  }

  if (voronoi) {
    tv = deldir(J[, 1], J[, 2], rw = c(xmin, xmax, ymin, ymax))
    plot(tv, add = TRUE, wlines = "tess", xlim = c(xmin, xmax), ylim = c(ymin, ymax))
  }

  if (!is.null(x$ContSupVarsBiplot))
    plot(x$ContSupVarsBiplot, F1=A1, F2=A2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode, TypeScale=TypeScale)
  #par(op)
}

