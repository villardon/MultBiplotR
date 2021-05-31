plot3dCanonicalBiplot <- function(Bip, A1 = 1, A2 = 2, A3=3,  ScaleGraph = TRUE, PlotGroups = TRUE, PlotVars = TRUE, PlotInd = TRUE, 
                                  LabelInd=TRUE, CexGroup=1, PchGroup=16 , margin=0.1, AddLegend=FALSE, ShowAxes=FALSE, LabelAxes=TRUE, LabelGroups=TRUE,
                                  PlotCircle = TRUE, ConvexHulls = FALSE, TypeCircle = "M", ColorGroups = NULL, ColorVars = NULL, LegendPos="topright",
                                  ColorInd = NULL, mode="a", TypeScale = "Complete", ValuesScale = "Original",
                                  MinQualityVars = 0, dpg = 0, dpi=0,  PredPoints=0, PlotAxis = FALSE, CexInd = NULL, CexVar = NULL,
                                  PchInd = NULL, PchVar = NULL, ColorVar=NULL, ShowAxis=TRUE, ColorAxis="gray", ...) {
  
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
    ColorGroups = 1:nlevels(Bip$groups)
  
  if (is.null(ColorInd)) 
    ColorInd = ColorGroups[as.integer(Bip$groups)]
  
  if (LabelAxes){
    xlabel = paste("Axis", A1)
    ylabel = paste("Axis", A2)
    zlabel = paste("Axis", A3)
    Title="Canonical/MANOVA Biplot"}
  else{
    xlabel = ""
    ylabel = ""
    zlabel = ""
    Title=paste("Canonical Biplot", " / ", A1, "-", A2, "-", A3, " (", round(Bip$Eigenvalues[A1,2]+Bip$Eigenvalues[A2,2]+Bip$Eigenvalues[A3,2], digits=2),"%)")
  }
  
  
  scf = 1
  if (ScaleGraph) {
    sca = sum(Bip$GroupCoordinates^2)
    scb = sum(Bip$ColCoordinates^2)
    sca = sca/Bip$g
    scb = scb/Bip$p
    scf = sqrt(sqrt(scb/sca))
  }
  
  J = Bip$GroupCoordinates[, c(A1, A2, A3)] * scf
  H = Bip$ColCoordinates[, c(A1, A2, A3)]/scf
  V = Bip$RowCoordinates[, c(A1, A2, A3)] * scf
  
  g = dim(J)[1]
  p = dim(H)[1]
  n = dim(V)[1]
  
  
  if (is.null(CexInd)) 
    CexInd = 0.8
  
  if (is.null(CexVar)) 
    CexVar = 0.8
  else if (length(CexVar) == 1) 
    CexVar = rep(CexVar, p)
  
  if (is.null(PchInd)) 
    PchInd = 1
  
  if (is.null(PchVar)) 
    PchVar = rep(16, p)
  
  if (is.null(ColorVar)) 
    ColorVar = rep("black", p)
  
  if (!is.null(Bip$Sup_Individual_Coord)){VS=Bip$Sup_Individual_Coord*scf}
  
  PP = rbind(J, H, V)
  xmin = min(PP[, 1])
  xmax = max(PP[, 1])
  ymin = min(PP[, 2])
  ymax = max(PP[, 2])
  zmin = min(PP[, 3])
  zmax = max(PP[, 3])
  
  #   xrang=abs(xmax-xmin)
  #   yrang=abs(ymax-ymin)
  #   zrang=abs(zmax-zmin)
  #   fact=0.03
  #   xmin=xmin-fact*xrang
  #   xmax=xmax+fact*xrang
  #   ymin=ymin-fact*yrang
  #   ymax=ymax+fact*yrang
  #   zmin=zmin-fact*zrang
  #   zmax=zmax+fact*zrang
  #   
  #   if (xmax <0 ) xmax=xmax*(-1)
  Title="Canonical/MANOVA Biplot"
  
  open3d()
  plot3d(0, 0, 0, cex = 0, asp = 1, main = Title, box=FALSE, cex.axis=0.5, axes=FALSE, ...)
  
  if (PlotInd) {
    points3d(V[, A1], V[, A2], V[, A3], col = ColorInd, pch=PchInd, cex=CexInd, ...)
    if (LabelInd)
      text3d(V[, A1], V[, A2], V[, A3], rownames(V), adj = 0, col = ColorInd, cex = 0.8)
    if (!is.null(Bip$Sup_Individual_Coord)){
      points(VS[, A1], VS[, A2], VS[, A3], cex = 0.5, col = "Black")
      text(VS[, A1], VS[, A2], VS[, A3], paste(rownames(VS),"*",sep=""), adj = 0, col = "Black", cex = 0.8)
    }
  }
  
  if (PlotGroups) {
    points3d(J[, A1], J[, A2], J[, A3], cex = CexGroup, pch = PchGroup, col = ColorGroups)
    if (LabelGroups)
      text3d(J[, A1], J[, A2], J[, A3], rownames(J), adj = 0, col = ColorGroups, cex = 1.5)
  }
  
  if (AddLegend) legend3d(LegendPos, rownames(J), fill = ColorGroups)
  
  if (PlotVars) {
    p=dim(H)[1]
    HH=matrix(0,2*p,3)
    for (i in 1:p)
      HH[2*i,]=H[i,c(A1,A2,A3)]
    segments3d(HH, cex = CexVar, pch = PchVar, col = ColorVar)
    text3d(H[, A1], H[, A2], H[, A3], rownames(H), adj = 0, col = ColorVar, cex = CexVar)
  }
  
  M=matrix(c(xmax, 0,0, 0, ymax, 0, 0, 0, zmax), 3,3)
  text3d(M, texts=c(xlabel, ylabel, zlabel), col=ColorAxis)
  
  # rgl.bbox(color = c("#333377", "white"), emission = "#333377", 
  #          specular = "#3333FF", shininess = 5, alpha = 0.8 )
  # bbox3d(color = c("#333377", "black"), emission = "#333377", 
  #        specular = "#3333FF", shininess = 5, alpha = 0.8)
  # #   
  #   if (PlotVars) {
  #     qlrcols = Bip$QLRVarMeans[, A1] + Bip$QLRVarMeans[, A2]+ Bip$QLRVarMeans[, A3]
  #     if (mode=="s")
  #       Scales = GetBiplotScales(Bip, TypeScale = TypeScale, ValuesScale = ValuesScale)
  #     
  #     for (j in 1:p) if (qlrcols[j] > MinQualityVars) 
  #       VarBiplot(H[j, 1], H[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = rownames(H)[j], mode = mode, CexPoint = CexVar[j], Color = ColorVar[j], 
  #                 ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j], ...)
  #     
  #     for (idp in dpg)
  #       if ((idp > 0) & (idp < (p + 1))) {
  #         g = H[idp, ]
  #         nn = (t(g) %*% g)
  #         scal <- (J %*% g)/nn[1, 1]
  #         Dscal <- diag(as.vector(scal))
  #         Fpr <- Dscal %*% matrix(rep(1, nrow(J)), ncol = 1) %*% t(g)
  #         nrFpr <- nrow(Fpr)
  #         
  #         dlines(J, Fpr, color=ColorGroups)
  #       }
  #     
  #     for (idp in dpi)
  #       if ((idp > 0) & (idp < (p + 1))) {
  #         g = H[idp, ]
  #         nn = (t(g) %*% g)
  #         scal <- (V %*% g)/nn[1, 1]
  #         Dscal <- diag(as.vector(scal))
  #         Fpr <- Dscal %*% matrix(rep(1, nrow(V)), ncol = 1) %*% t(g)
  #         nrFpr <- nrow(Fpr)
  #         
  #         dlines(V, Fpr, color=ColorInd)
  #       }
  #     
  #     for (idp in PredPoints)
  #       if ((idp > 0) & (idp < (n + 1)))
  #         for (j in 1:p){
  #           g = H[j, ]
  #           nn = (t(g) %*% g)
  #           scal <- (J[idp,] %*% g)/nn[1, 1]
  #           Fpr <- scal %*% t(g)
  #           nrFpr <- nrow(Fpr)
  #           dlines(matrix(J[idp,],1,2) , Fpr, color=ColorVar[j])
  #         }
  #     
  #   }
  #   
  
  if (PlotCircle) {
    switch(TypeCircle, U = {
      radius = Bip$UnivRad
    }, B = {
      radius = Bip$BonfRad
    }, M = {
      radius = Bip$MultRad
    }, C = {
      radius = Bip$ChisRad
    }) 
    radius=radius*scf;
    for (i in 1:Bip$g) spheres3d(J[i, A1], J[i, A2], J[i, A3], radius[i], col = ColorGroups[i])
  }
  
  if (ConvexHulls) {
    lev = levels(Bip$groups)
    for (i in 1:nlevels(Bip$groups)) {
      XP = V[which(Bip$groups == lev[i]), ]
      XP = cbind(XP[, A1], XP[, A2], XP[, A3])
      ts.surf <- t(convhulln(XP))
      rgl.triangles(XP[ts.surf,1],XP[ts.surf,2],XP[ts.surf,3],col=ColorGroups[i],alpha=.5)
    }
  }
  if (ShowAxis) abclines3d(0, 0, 0, a = diag(3), col = ColorAxis, linewiidth=2)
  
}

