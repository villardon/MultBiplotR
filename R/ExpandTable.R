# This function expands the patterns of a table of patterns and frequencies

ExpandTable <- function(table){
  if (!(class(table)=="TableFrequencies")) stop("You must provide a Table with Frequencies to extend")
  ncomp=dim(table$Patterns)[1]
  p=dim(table$Patterns)[2]
  nexp=sum(table$Frequencies)
  x=matrix(0,nexp, p)
  for (i in 1:ncomp){
    indices=as.integer(table$EqualRows[[i]])
    for (j in 1:length(indices))
      x[indices[j],]=table$Patterns[i,]
  }	
  colnames(x) <- colnames(table$Patterns)
  x[table$Order,]=x
  rownames(x)=table$OriginalNames
  return(x)
}




ExpandCoord <- function(Coord, table){
  if (!(class(table)=="TableFrequencies")) stop("You must provide a Table with Frequencies to extend")
  ncomp=dim(Coord)[1]
  nexp=sum(table$Frequencies)
  x=matrix(0,nexp, ncomp)
  for (i in 1:ncomp){
    indices=as.integer(table$EqualRows[[i]])
    for (j in 1:length(indices))
      x[indices[j],]=Coord[i,]
  }	
  colnames(x) <- colnames(Coord)
  x[table$Order,]=x
  rownames(x)=table$OriginalNames
  return(x)
}
