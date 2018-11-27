PostHocGraph <- function(pvalmatrix, version=1, color=c("red3", "red","gray97"), nrows=NULL, panel=TRUE, title="Post-Hoc p-values", ...){
  nvar=dim(pvalmatrix)[1]
  ncomp=dim(pvalmatrix)[2]
  ng=(1+sqrt(1+8*ncomp))/2
  varnames=rownames(pvalmatrix)
  compnames=colnames(pvalmatrix)
  groupnames=strsplit(compnames[1], ":")[[1]][1]
  for (i in 1:(ng-1))
    groupnames=c(groupnames, strsplit(compnames[i], ":")[[1]][2])


  if (is.null(nrows))
    nrows=round(sqrt(nvar))
  ncols=ceiling(nvar/nrows)


  if (version ==3){
    require(lattice)
    op=par(mfrow=c(nrows, ncols))
  }
  else require(gplots)
  P=list()

  for (l in 1:nvar){
    P[[l]]=matrix(1, ng, ng)
    k=0
    for (i in 1:(ng-1))
      for (j in (i+1):ng){
        k=k+1
        P[[l]][i,j]=pvalmatrix[l,k]
        P[[l]][j,i]=P[[l]][i,j]
      }
    rownames(P[[l]])=groupnames
    colnames(P[[l]])=groupnames
  }

  for (l in 1:nvar){
   title2=paste(title,":", varnames[l])
    if (version==1)
      heatmap(P[[l]], Rowv=NA, scale="none",symm=TRUE, revC=FALSE, breaks=c(0, 0.01, 0.05, 1), col=color, margins=c(7,7), main=title2, ...)

   if (version==2){
      #dev.new()
      heatmap.2(P[[l]], Rowv=FALSE, scale="none",symm=TRUE, dendrogram="none", revC=FALSE, breaks=c(0, 0.01, 0.05, 1), col=color, key=FALSE, density.info="none",
                cellnote=round(P[[l]], digits=3), notecol="black", sepcolor="black", trace="none", margins=c(5,5), main=title2, lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lhei=c(0.3, 4, 0.3), ...)}

   if (version==3) levelplot(P[[l]])

  }
  names(P)=varnames
save(P, file="~/Dropbox/0 Daniela/P.rda")

  if (version ==3){
    par(op)
  }
}
