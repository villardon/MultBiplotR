plot.Supplementary.Variables <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                         ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, 
                                         ColorSupContVars="black", ColorSupBinVars="blue", ColorSupOrdVars="red", 
                                         ModeSupContVars="a", ModeSupBinVars="a", 
                                         ModeSupOrdVars="a", ...){
   if (!is.null(x$ContSupVarsBiplot))
    plot(x$ContSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=ModeSupContVars, TypeScale=TypeScale, ColorVar=ColorSupContVars) 
  
  if (!is.null(x$BinSupVarsBiplot))
    plot(x$BinSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=ModeSupBinVars, ColorVar=ColorSupBinVars)

  
  if (!is.null(x$OrdSupVarsBiplot))
    plot(x$OrdSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=ModeSupOrdVars, ColorVar=ColorSupOrdVars) 
  
}
