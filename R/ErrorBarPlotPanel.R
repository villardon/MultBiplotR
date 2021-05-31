ErrorBarPlotPanel <- function(X, groups=NULL, nrows=NULL, panel=TRUE, GroupsTogether=TRUE, 
                              Confidence=0.95, p.adjust.method="None", UseANOVA=FALSE, Colors="blue",
                              Title="Error Bar Plot",  sort=TRUE, ...){

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (is.vector(X)) X=matrix(X, ncol=1)

  separated=!panel
  
  k=0
  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(groups)) {
    groups=as.factor(rep(1,n))
    levels(groups)="Complete Sample"}

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  if (length(Colors)==1) Colors=rep(Colors,g)
  Levels=levels(groups)
  varnames=colnames(X)

  switch(p.adjust.method,
         None={
           Confidence=Confidence
         },
         Sidak={
           Confidence=Confidence^(1/g)
         },
         Bonferroni={
           Confidence=1 - ((1-Confidence)/(g*(g-1)/2))
         })

  if (GroupsTogether){

    if (is.null(nrows))
      nrows=round(sqrt(p))

    ncols=ceiling(p/nrows)
   print(c(p, nrows, ncols))
    if (separated==FALSE)
      op=par(mfrow=c(nrows, ncols) ,oma = c(0, 0, 2, 0))

    for (j in 1:p){
      if (separated==TRUE) dev.new()

      XX = as.numeric(X[,j])
      Means=tapply(XX,groups,mean)
      met=mean(XX)
      if (UseANOVA){
        Desv=sqrt(sum((tapply(XX, groups, sd )^2)*(table(groups)-1))/(sum(table(groups))-g))
        df=n-g
      }
      else{
        Desv=NULL
        df=n-1}
      
      int=tapply(XX,groups,ConfidenceInterval, Confidence=Confidence, Desv=Desv, df = df)
      interv=matrix(0, g, 2)
      for (i in 1:g){
        interv[i,]=int[i][[1]]
      }
      
      if (sort){
        oi=sort(Means, index.return = TRUE)$ix
        Levels=Levels[oi]
        Means=Means[oi]
        interv=interv[oi,]
        Colors=Colors[oi]
      }
       Hmisc::errbar(Levels,Means,interv[,1],interv[,2], xlab=varnames[j], errbar.col=Colors, main=Title,  ...)
      title(varnames[j])
    }

    if (separated==FALSE)
      par(op)
  }
  else{
    if (is.null(nrows))
      nrows=round(sqrt(g))

    ncols=ceiling(g/nrows)

    if (separated==FALSE)
      op=par(mfrow=c(nrows, ncols))

    for (i in 1:g){
      XX = as.matrix(X[which(groups == Levels[i]), ])
      x=numeric()
      grupvar=numeric()
      for (j in 1:p){
        x=c(x, as.numeric(XX[,j]))
        grupvar=c(grupvar, rep(j, length(XX[,j])))
      }
      grupvar=as.factor(grupvar)
      levels(grupvar)=varnames

      me=by(x,grupvar,mean)
      int=by(x,grupvar,ConfidenceInterval)
      Means=matrix(0, p, 1)
      interv=matrix(0, p, 2)

      for (j in 1:p){
        Means[j]=me[j][[1]]
        interv[j,]=int[j][[1]]
      }
      if (separated==TRUE) dev.new()
      errbar(as.factor(varnames),Means,interv[,1],interv[,2], main=Levels[i], errbar.col=Colors, ...)
      title(Levels[i])
    }
  }
  # mtext("Error Bar Panel", outer = TRUE, cex = 1.5)
}
