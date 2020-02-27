plot.BinSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                  ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ColorVar=NULL, ...){
  p=dim(x$ColumnParameters)[1]
  if (is.null(ColorVar)) ColorVar=rep("black", p)
  ColLabels=rownames(x$VarInfo)
  for (i in 1:p)
    PlotBinaryVar(b0=x$ColumnParameters[i,1], bi1=x$ColumnParameters[i,F1+1], bi2=x$ColumnParameters[i,F2+1], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                  mode=mode, label=ColLabels[i], Color=ColorVar[i])
}
