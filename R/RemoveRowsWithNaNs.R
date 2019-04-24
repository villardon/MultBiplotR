RemoveRowsWithNaNs <- function(x, cols=NULL) {
  p=dim(x)[2]
  if (is.null(cols))
    cols=1:p
  sumas=apply(x[, cols],1,sumisna)
  good=which(sumas==0)
  newx=x[good,]
  return(newx)
}

sumisna <- function(x){
  return(sum(is.na(x)))
}

