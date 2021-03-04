plot.StatisBiplot <- function(x, A1 = 1, A2 = 2, PlotType="Biplot", PlotRowTraj=FALSE, PlotVarTraj=FALSE, LabelTraj='Begining',
                              VarColorType="ByVar", VarColors=NULL, VarLabels=NULL, RowColors=NULL,
                              TableColors=NULL, RowRandomColors=FALSE, TypeTraj="line", ...) {
  PlotTypes=c("Biplot", "Correlations", "Contributions", "InterStructure")
  if (is.numeric(PlotType)) PlotType=PlotTypes[PlotType]
  
  VarColorTypes=c("Biplot", "ByTable", "ByVar")
  if (is.numeric(VarColorType)) VarColorType=VarColorTypes[VarColorType]
  ColorVar=VarColors
  
  if (is.null(TableColors)) TableColors=1:x$NTables
  if (is.null(VarColors)) VarColors=1:x$NVars[1]
  
  LabelTrajectories=c("Begining", "End", "None")
  if (is.numeric(LabelTraj)) LabelTraj=LabelTrajectories[LabelTraj]
  
  nind=dim(x$Biplot$RowCoordinates)[1]
  
  noc=x$NTables
  
  if (RowRandomColors) RowColors=colors()[sample.int(657,nind)]
  
  if (is.null(RowColors)) RowColors=rep("black", nind)
  
  if (x$SameVar){
    nvars =x$NTables * x$NVars
    if (VarColorType=="ByTable"){
      ColorVar=rep(TableColors[1], x$NVars)
      for (i in 2:x$NTables) ColorVar=c(ColorVar, rep(TableColors[i], x$NVars))
      VarLabels=rep(x$VarLabels, x$NTables)
    }
    
    if (VarColorType=="ByVar"){
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
    ColContributionPlot(x$Biplot, Colors=ColorVar, Labs=VarLabels, cex=0.7, MinQuality = 0, AddSigns2Labs=FALSE)
  
  
  if (PlotType=="InterStructure"){
    par(mfrow=c(2,1))
    plot(x$InterStructure[,1],x$InterStructure[,2], xlim=c(0,1.1), asp=1, cex=0.1, col=TableColors, 
         main="InterStructure", xlab="Dim 1", ylab="Dim 2")
    arrows(0,0, x$InterStructure[,1],x$InterStructure[,2], col=TableColors)
    text(x$InterStructure[,1],x$InterStructure[,2],labels = x$TableLabels, pos=4, col=TableColors)
    
    barplot(x$InerInter, main="Explained Variance (Interstructure)")
    par(mfrow=c(1,1))
  }
  
  if (PlotType!="InterStructure"){
    if (x$SameVar){
      if (VarColorType=="ByTable")
        legend("topright",legend=x$TableLabels, fill=TableColors, cex=0.7)
      if (VarColorType=="ByVar")
        legend("topright",legend=x$VarLabels, fill=VarColors, cex=0.7)
    }
  }
  
  if (!x$SameVar) PlotVarTraj=FALSE
  
  
  if (PlotRowTraj){
    PlotTrajectory(x$TrajInd, Centers=x$Biplot$RowCoordinates, A1=A1, A2=A2, TypeTraj=TypeTraj, ColorTraj=RowColors, LabelTraj=LabelTraj)
  }
  
  if (PlotVarTraj){
    VarColors=1:x$NVars[1]
    PlotTrajectory(x$TrajVar, Centers=NULL, A1=A1, A2=A2, TypeTraj="line", ColorTraj=VarColors, LabelTraj=LabelTraj)
  }
  
}
