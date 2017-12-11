AdaptDataFrame<- function(x, Binary=NULL, IntegerAsOrdinal=FALSE) {
  
  Names=colnames(x)
  Classes=as.character(sapply(x,class))
  
  Integer=which(Classes=="integer")
  if (length(Integer)>0){
    for (i in 1:length(Integer))
      if (IntegerAsOrdinal)
        x[[Integer[i]]]=as.ordered(x[[Integer[i]]])
      else
        x[[Integer[i]]]=as.numeric(x[[Integer[i]]])}
  
  if (length(Binary)>0){
    for (i in 1:length(Binary)){
      if (CheckBinaryVector(x[[Binary[i]]])){
        Classes[[Binary[i]]]="binary"
        class(x[[Binary[i]]])="numeric"
        }
      else
      stop("Some binary variable is not correctly declared")
    }
  }
  
  Factors=which(Classes=="factor")
  if (length(Factors)>0){
    for (i in 1:length(Factors)){
      if (length(levels(x[[Factors[i]]]))==2){
        x[[Factors[i]]]=Factor2Binary(x[[Factors[i]]], Name=Names[Factors[i]])[,2]
        class(x[[Factors[i]]])="numeric"
      Classes[[Factors[i]]]="binary"}
    }
  }
  
  Ordered=which(Classes=="c(\"ordered\", \"factor\")")
  if (length(Ordered)>0){
    for (i in 1:length(Ordered))
      Classes[[Ordered[i]]]="ordered"
  }
  
  return(list(X=x, Types=Classes))
  
}

