ConvertFactors2Integers <- function(x){
  if (is.data.frame(x)){
    p=dim(x)[2]
    for (i in 1:p)
      if (is.factor(x[,i]))
        x[,i]=as.numeric(x[,i])}
  return(x)
}