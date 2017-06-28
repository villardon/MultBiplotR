# Defines the points for a trajectory using two factors
# A separate line will be plotted for each category of the MainFactor
# The order of the points in the trajectory is given by the second factor
# The value is an object of class "trajectory" that contains the factors
# labels and colors to plot the trajectory

DefineTrajectory <- function(MainFactor, SecondaryFactor=NULL){
  
  n=length(MainFactor)
  if (is.factor(MainFactor)){
    MainFactor=droplevels(MainFactor)
    TrajNames=levels(MainFactor)
  }
  else {
    MainFactor=as.factor(MainFactor)
    MainFactor=droplevels(MainFactor)
    TrajNames=levels(MainFactor)
  }
  
  if (!is.null(SecondaryFactor))
    if (is.factor(SecondaryFactor)){
      AB <- fac.combine(list(MainFactor,SecondaryFactor), combine.levels=TRUE)
      SecondaryFactor=droplevels(SecondaryFactor)
      SecondaryFactor=as.integer(SecondaryFactor)
    }
  else{
    SecondaryFactor=1:n}
  
  Trajectories=list()
  for (i in 1:length(TrajNames)){
    Trajectories[[i]]=SecondaryFactor[which(MainFactor==TrajNames[i])]
  }
  
names(Trajectories)=TrajNames
}