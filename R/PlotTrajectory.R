PlotTrajectory <- function(Traj, Centers=NULL, A1=1, A2=2, TypeTraj="line", ColorTraj=NULL, LabelTraj="Begining"){
  ntraj=length(Traj)
  noc=dim(Traj[[1]])[1]
  TrajNames=names(Traj)
  if (TypeTraj=="line"){
    for (i in 1:ntraj){
      points(Traj[[i]][,A1], Traj[[i]][,A2], type="l", col=ColorTraj[i])
      if (LabelTraj=='Begining') text(Traj[[i]][1,A1], Traj[[i]][1,A2], label=TrajNames[i], col=ColorTraj[i])
      if (LabelTraj=='End') text(Traj[[i]][noc,A1], Traj[[i]][noc,A2], label=TrajNames[i], col=ColorTraj[i])
    }
  }
  
  if (TypeTraj=="star"){
    for (i in 1:ntraj)
      for (j in 1:noc)
        segments(Centers[i,A1], Centers[i,A2], Traj[[i]][j,A1], Traj[[i]][j,A2], col= ColorTraj[i])
  }
}
