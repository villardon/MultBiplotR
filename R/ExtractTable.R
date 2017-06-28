# Extracs the unique patterns and the frequencies of a data matrix

ExtractTable <- function(x) {
  
  result = list()
  n = dim(x)[1]
  p = dim(x)[2]
  
  if (is.null(rownames(x)))
    Names <- paste("R",1:n, sep="")
  else
    Names=rownames(x)
  
  rownames(x)=1:n
  x=x[do.call(order, as.data.frame(x)), ]
  Orden=as.integer(rownames(x))
  rownames(x)=Names[Orden]

  x=as.matrix(x)
  result$Data=x
  suma=sum(x)
  
  if (is.na(suma)){
    maximo=max(x, na.rm = TRUE)
    x[which(is.na(x))]=maximo +1
  }
  
  i=1
  patterns=NULL
  frequencies=NULL
  RowNames=NULL
  EqualRows = list()
  l=0
  while (i < n) {
    patterns = rbind(patterns, x[i, ])
    RowNames=rbind(RowNames, Names[i])
    k = 0
    freq = 0
    WhatRows = NULL
    while (((i+k)<n) & (sum(x[i, ] == x[i + k, ]) == p) ){
      freq = freq + 1
      WhatRows = c(WhatRows, i+k)
      k = k + 1
    }
    frequencies=rbind(frequencies,freq)
    l=l+1
    EqualRows[[l]]=WhatRows
    i=i+k
  }
  
  if (sum(x[i, ] == patterns[l, ]) == p) {
    frequencies[l] = frequencies[l] + 1
    EqualRows[[l]] = c(EqualRows[[l]], i)
  } else {
    patterns = rbind(patterns, x[i, ])
    frequencies=rbind(frequencies, 1)
    EqualRows[[l + 1]] = i
    RowNames=c(RowNames,Names[i])
  }
  
  if (is.na(suma)){
    patterns[which(patterns==(maximo+1))]=NA
  }
  
  rownames(patterns)<-RowNames
  rownames(frequencies) <- RowNames
  
  result$Order=Orden
  result$OriginalNames=Names
  result$Patterns = patterns
  result$Frequencies = frequencies
  result$EqualRows = EqualRows
  result$Unique=rep(0, length(EqualRows))
  for (i in 1:length(EqualRows)) result$Unique[i]=EqualRows[[i]][1]
  
  class(result) <- "TableFrequencies"
  return(result)
  
}


CompressCoord <- function(Coord, Table){
  Coord=Coord[Table$Unique,]
  
}
