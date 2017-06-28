
ExtractIndices <- function(table){
  if (!(class(table)=="TableFrequencies")) stop("You must provide a Table with Frequencies to extend")
  ncomp=dim(table$Patterns)[1]
  indices=integer()
  for (i in 1:ncomp){
    indices[i]=table$EqualRows[[i]][1]
  }	
  return(indices)
}
