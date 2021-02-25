plot.Unfolding <- function(x, A1 = 1, A2 = 2, ShowAxis = FALSE, margin = 0.1, PlotSites = TRUE, PlotSpecies = TRUE, PlotEnv = TRUE, LabelSites = TRUE, 
                           LabelSpecies = TRUE, LabelEnv = TRUE, SpeciesQuality = FALSE, MinQualityVars = 0, dp = 0, PlotAxis = FALSE, TypeScale = "Complete", ValuesScale = "Original", 
                           mode = "h", CexSites = NULL, CexSpecies = NULL, CexVar = NULL, ColorSites = NULL, ColorSpecies = NULL, ColorVar = NULL, PchSites = NULL, PchSpecies = NULL, 
                           PchVar = NULL, SizeQualSites = FALSE, SizeQualSpecies = FALSE, SizeQualVars = FALSE, ColorQualSites = FALSE, ColorQualSpecies = FALSE, ColorQualVars = FALSE, 
                           SmartLabels = FALSE, PlotTol = FALSE, ...) {
  dims = x$dimens
  A = x$Y[, c(A1, A2)]
  
  X = x$X[, c(A1, A2)]
  
  SpeciesNames = rownames(A)
  
  SiteNames = rownames(X)
  
  QualitySpecies = x$ColumnsFit
  
  QualitySites = (x$Sites_Quality[, A1] + x$Sites_Quality[, A2])/100

  
  if (is.null(x$B)) {
    B = x$Env_Var_Scores[, c(A1, A2)]
    EnvNames = rownames(B)
    QualityVars = x$Var_Fit
    
    if (!is.null(x$nvars)) 
      if (is.null(CexVar)) 
        CexVar = rep(0.8, x$nvars)
    else if (length(CexVar == 1)) 
      CexVar = rep(CexVar, x$nvars)
    
    if (is.null(PchVar)) 
      PchVar = rep(16, x$nvars)
    
    if (!is.null(x$nvars)) 
      if (is.null(ColorVar)) 
        ColorVar = rep("blue", x$nvars)
    if (SizeQualVars) 
      CexVar = cscale(QualityVars, rescale_pal())
    if (ColorQualVars) 
      ColorVar = cscale(QualityVars, seq_gradient_pal("white", "darkblue"))
  }
    
    if (is.null(CexSites)) 
      CexSites = 0.8
    if (is.null(CexSpecies)) 
      CexSpecies = 0.8
    
    if (is.null(PchSites)) 
      PchSites = 3
    if (is.null(PchSpecies)) 
      PchSpecies = 1
    
    
    if (is.null(ColorSites)) 
      ColorSites = "black"
    if (is.null(ColorSpecies)) 
      ColorSpecies = "red"
    
    if (SizeQualSites) 
      CexSites = cscale(QualitySites, rescale_pal())
    if (SizeQualSpecies) 
      CexSpecies = cscale(QualitySpecies, rescale_pal())
    
    if (ColorQualSites) 
      ColorSites = cscale(QualitySites, seq_gradient_pal("white", "black"))
    if (ColorQualSpecies) 
      ColorSpecies = cscale(QualitySpecies, seq_gradient_pal("white", "red"))

    
    if (ShowAxis) {
      xaxt = "s"
      yaxt = "s"
    } else {
      xaxt = "n"
      yaxt = "n"
    }
    
    if ((margin < 0) | (margin > 0.3)) 
      margin = 0
    
    
    P = rbind(A, X)
    xmin = min(P[, 1])
    xmax = max(P[, 1])
    ymin = min(P[, 2])
    ymax = max(P[, 2])
    P = rbind(P, c(xmin - (xmax - xmin) * margin, ymin - (ymax - ymin) * margin))
    P = rbind(P, c(xmax + (xmax - xmin) * margin, ymax + (ymax - ymin) * margin))
    xlabel = paste("Dimension", A1)
    ylabel = paste("Dimension", A2)
    plot(P[, A1], P[, A2], cex = 0, asp = 1, main = x$Analysis, xlab = xlabel, ylab = ylabel, , xaxt = xaxt, yaxt = yaxt, ...)
    
    if (PlotSpecies) 
      points(A[, A1], A[, A2], cex = CexSpecies, col = ColorSpecies, pch = PchSpecies, ...)
    if (LabelSpecies) 
      if (SmartLabels) 
        textsmart(cbind(A[, A1], A[, A2]), CexPoints = CexSpecies, ColorPoints = ColorSpecies)
    else text(A[, A1], A[, A2], rownames(A), cex = CexSpecies, col = ColorSpecies, pos = 1)
    
    if (PlotSites) 
      points(X[, A1], X[, A2], cex = CexSites, col = ColorSites, pch = PchSites)
    if (LabelSites) 
      if (SmartLabels) 
        textsmart(cbind(X[, A1], X[, A2]), CexPoints = CexSites, ColorPoints = ColorSites)
    else text(X[, A1], X[, A2], rownames(X), cex = CexSites, col = ColorSites, pos = 1)
    
    
    if (!is.null(x$nvars)) 
      if (PlotEnv) {
        Scales = GetCCAScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)
        for (j in 1:x$nvars) if (QualityVars[j] > MinQualityVars) 
          VarBiplot(B[j, 1], B[j, 2], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = rownames(B)[j], mode = mode, CexPoint = CexVar[j], 
                    Color = ColorVar[j], ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j])
        
        
        if ((dp > 0) & (dp < (x$nvars + 1))) {
          g = B[dp, ]
          nn = (t(g) %*% g)
          scal <- (A %*% g)/nn[1, 1]
          Dscal <- diag(as.vector(scal))
          Fpr <- Dscal %*% matrix(rep(1, nrow(A)), ncol = 1) %*% t(g)
          nrFpr <- nrow(Fpr)
          dlines(A, Fpr)
        }
      }
    
    if (PlotTol) 
      for (i in 1:x$nspecies) Circle(x$Tolerance[i], origin = c(A[i, 1], A[i, 2]), col = "red")

}