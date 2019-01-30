plot.StatisBiplot <- function(x, A1 = 1, A2 = 2, PlotType="Biplot", PlotRowTraj=FALSE, PlotVarTraj=FALSE, LabelTraj='Begining',
                        VarColorType="Biplot", VarColors=NULL, VarLabels=NULL, RowColors=NULL,
                        TableColors=NULL, RowRandomColors=FALSE, ...) {
  PlotTypes=c("Biplot", "Correlations", "Contributions", "InterStructure")
  if (is.numeric(PlotType)) PlotType=PlotTypes[PlotType]
  
  VarColorTypes=c("Biplot", "ByTable", "ByVar")
  if (is.numeric(VarColorType)) VarColorType=VarColorTypes[VarColorType]
  ColorVar=VarColors
  
  LabelTrajectories=c("Begining", "End", "None")
  if (is.numeric(LabelTraj)) LabelTraj=LabelTrajectories[LabelTraj]
  
  nind=dim(x$Biplot$RowCoordinates)[1]
  if (RowRandomColors) RowColors=colors()[sample.int(657,nind)]
  
  if (x$SameVar){
    nvars =x$NTables * x$NVars
    if (VarColorType=="ByTable"){
      if (is.null(TableColors)) TableColors=1:x$NTables
      ColorVar=rep(TableColors[1], x$NVars)
      for (i in 2:x$NTables) ColorVar=c(ColorVar, rep(TableColors[i], x$NVars))
      VarLabels=rep(x$VarLabels, x$NTables)
    }
    
    if (VarColorType=="ByVar"){
      if (is.null(VarColors)) VarColors=1:x$NVars
      ColorVar=rep(VarColors, x$NTables)
      VarLabels=rep(x$TableLabels[1], x$NVars)
      for (i in 2:x$NTables) VarLabels = c(VarLabels, rep(x$TableLabels[i], x$NVars))
     }
  }
  
  if (PlotType=="Biplot")
  plot(x$Biplot, A1=A1, A2=A2, ColorVar=ColorVar, VarLabels=VarLabels, ColorInd=RowColors, ...)
  
  if (PlotType=="Correlations")
  CorrelationCircle(x$Biplot, Colors=ColorVar, Labs=VarLabels, cex=0.7)
  
  if (PlotType=="Contributions")
  ColContributionPlot(x$Biplot, Colors=ColorVar, Labs=VarLabels, cex=0.7, MinQuality = 0.7)
  
  if (x$SameVar){
    if (VarColorType=="ByTable")
      legend("topright",legend=x$TableLabels, fill=TableColors, cex=0.7)
    if (VarColorType=="ByVar")
      legend("topright",legend=x$VarLabels, fill=VarColors, cex=0.7)
  }
  
  if (!x$SameVar) PlotVarTraj=FALSE
  
  print(PlotRowTraj)
  if (PlotRowTraj){
    print("The projection of trajectories will be available soon")
    TrajNames=names(x$TrajInd)
    for (i in 1:length(x$TrajInd))
      points(x$TrajInd[[i]][,A1], x$TrajInd[[i]][,A2], type="l", col=RowColors[i])
    
    if (LabelTraj=='Begining') text(x$TrajInd[[i]][1,A1], x$TrajInd[[i]][1,A2], label=TrajNames[i], col=RowColors[i])
  }
  
  if (PlotVarTraj){
    print("The projection of trajectories will be available soon")
  }
  
}
