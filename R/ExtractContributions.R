# Extract the contributions of a biplot summing for the components listed in Axes.
ExtractContributions <- function(bip, Type=3, Ordered=TRUE, Axes=c(1,2), MinContribution=0){
  Types=c("Rows", "Columns", "Both")
  if (is.numeric(Type)) Type=Types[Type]
  result=list()
  result$What=Type
  ColName="Axes-"
  for (i in 1 : length(Axes)) ColName =paste(ColName, Axes[i],"-")
  if (Ordered) result$Ordered="yes"
  else result$Ordered="no"
  if ((Type=="Rows") | (Type=="Both")){
    RowCont=bip$RowContributions[,Axes]
    RowCont=apply(RowCont,1,sum)
    if (Ordered) RowCont=sort(RowCont,decreasing=TRUE)
    RowCont=RowCont[RowCont>MinContribution]
    Nombres=names(RowCont)
    result$RowContributions=matrix(RowCont,length(RowCont),1)
    rownames(result$RowContributions)=Nombres
    colnames(result$RowContributions)=ColName
  }
  
  if ((Type=="Columns") | (Type=="Both")){
    ColCont=bip$ColContributions[,Axes]
    ColCont=apply(ColCont,1,sum)
    if (Ordered) ColCont=sort(ColCont,decreasing=TRUE)
    ColCont=ColCont[ColCont>MinContribution]
    Nombres=names(ColCont)
    result$ColContributions=matrix(ColCont, length(ColCont) , 1)
    rownames(result$ColContributions)=Nombres
      colnames(result$ColContributions)=ColName
  }
  return(result)
}