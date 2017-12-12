plot.Supplementary.Variables <- function(bip, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                         ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ...){
  if (!is.null(bip$ContSupVarsBiplot))
    plot(bip$ContSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode, TypeScale=TypeScale) 
  
  if (!is.null(bip$BinSupVarsBiplot))
    plot(bip$BinSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode) 
  
  if (!is.null(bip$OrdSupVarsBiplot))
    plot(bip$OrdSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode) 
  
}
