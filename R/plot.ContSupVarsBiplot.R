plot.ContSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                   ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, PchVar=1, ColorVar=NULL, ...){
 
  if (is.null(ColorVar)) ColorVar=3
  B=x$ColCoordinates[,c(F1,F2)]
  b0=x$b0
  VarLabels=rownames(x$ColCoordinates)

  if (mode=="s")
    Scales = GetBiplotScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)
  
  p=dim(x$ColCoordinates)[1]
  if (is.numeric(B)) B=matrix(B, nrow=p)
  for (j in 1:p) 
    VarBiplot(B[j, 1], B[j, 2], b0=b0[j], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = VarLabels[j], mode = mode, 
              ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j], tl=0.001, Color = ColorVar, ...)
  
}