plot.Statis <- function(x, A1 = 1, A2 = 2, PlotRowTraj=FALSE, PlotVarTraj=FALSE, LabelTraj='begining',...) {
  plot(x$Biplot, A1=A1, A2=A2, ...)
  
  if (!x$SameVar) PlotVarTraj=FALSE

  if (PlotRowTraj){
    print("The projection of trajectories will be available soon")
  }
  
  }
