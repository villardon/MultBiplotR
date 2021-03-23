BasicDescription <- function(X, groups=NULL, SortByGroups=FALSE, na.rm=FALSE, Intervals=TRUE){
  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(groups)) {
    groups=as.factor(rep(1,n))
    levels(groups)="Complete Sample"}

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  Levels=levels(groups)

  if (!Intervals) {ncols=8
  statnames=c("n", "Mean", "Std. Dev", "Median", "min", "max", "range", "S.E.")}
  else {
    ncols=10
    statnames=c("n", "Mean", "Std. Dev", "Median", "min", "max", "range", "S.E.", "CI-left", "CI-right")}

  Descriptives=list()

  if (SortByGroups==TRUE){
    for (i in 1:g){
      XX = as.matrix(X[which(groups == Levels[i]), ])
      Stats=matrix(0, p, ncols)
      colnames(Stats)=statnames
      Stats[,1]=apply(XX,2,Valid)
      Stats[,2]=apply(XX,2,mean)
      Stats[,3]=apply(XX,2,sd)
      Stats[,4]=apply(XX,2,median)
      Stats[,5]=apply(XX,2,min)
      Stats[,6]=apply(XX,2,max)
      Stats[,7]=Stats[,6]-Stats[,5]
      Stats[,8]= Stats[,3]/sqrt(Stats[,1])
      if (Intervals){
        Stats[,9:10]=t(apply(XX,2,ConfidenceInterval))
      }

      rownames(Stats)=colnames(X)
      Descriptives[[i]]=Stats
    }
    names(Descriptives) = Levels
  }


  if (SortByGroups==FALSE){
    for (i in 1:p){
      XX = as.numeric(X[,i])
      Stats = matrix(0, g, ncols)
      colnames(Stats)=statnames
      Stats[,1]=by(XX,groups,Valid)
      Stats[,2]=by(XX,groups,mean)
      Stats[,3]=by(XX,groups,sd)
      Stats[,4]=by(XX,groups,median)
      Stats[,5]=by(XX,groups,min)
      Stats[,6]=by(XX,groups,max)
      Stats[,7]=Stats[,6]-Stats[,5]
      Stats[,8]= Stats[,3]/sqrt(Stats[,1])

      if (Intervals){
        ppp=by(XX,groups,ConfidenceInterval)
        interv=matrix(0,g,2)
        for (j in 1:g)
          interv[j,]=ppp[j][[1]]
        Stats[,9:10]=interv
      }
      rownames(Stats)=Levels
      Descriptives[[i]]=Stats
    }
    names(Descriptives)= colnames(X)
  }
  return(Descriptives)
}

Valid <- function(v){
  length(Which)
}
