
plot.OrdSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                  ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ColorVar=NULL, ...){
  p=dim(x$ColumnParameters)[1]
  if (is.null(ColorVar)) ColorVar=rep("green", p)
  nr=dim(x$ColCoordinates)[1]
  B = matrix(x$ColCoordinates[, c(F1, F2)], nrow=nr)
  rownames(B)=rownames(x$ColCoordinates)
  thresholds= matrix(x$ColumnParameters$thresholds, nrow=nr)
  
  p = dim(x$ColumnParameters$coefficients)[1]
  for (j in 1:p)
    OrdVarBiplot(B[j, 1], B[j, 2], thresholds[j,1:(x$Ncats[j]-1)], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode, label=rownames(B)[j], Color=ColorVar[j])
}

