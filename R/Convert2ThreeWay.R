# Converts a two-dimensional matrix into a list where
# each cell is the two dimensional data matrix  
# for an occasion or group

Convert2ThreeWay <- function(x, groups, columns = FALSE, RowNames=NULL) {
  if (is.data.frame(x)) 
    x = as.matrix(x)
  
  n = dim(x)[1]
  p = dim(x)[2]
  g = length(levels(groups))
  n2=length(groups)
  
  if ((columns==FALSE) & (n != n2)) stop("Number of rows must be equal to the number of elements in groups")
  if ((columns==TRUE) & (p != n2)) stop("Number of columns must be equal to the number of elements in groups")
  
  
  levellab = levels(groups)
  
  ColumnNames = colnames(x)
  if (is.null(ColumnNames)) {
    for (i in 1:p) ColumnNames = c(ColumnNames, paste("C", i, sep = ""))
    colnames(x) = ColumnNames
  }
  
  
  if (is.null(RowNames)){
    RowNames = rownames(x)
    if (is.null(RowNames)) {
      for (i in 1:n) RowNames = c(RowNames, paste("R", i, sep = ""))
      rownames(x) = RowNames
    }
  }

  
  if (columns) x=t(x)
  
  X = list()
  Sizes = matrix(0,g, 1)
  for (i in 1:g) {
    X[[i]] = x[which(groups == levellab[i]), ]
    rownames(X[[i]])=RowNames[which(groups == levellab[i])]
    Sizes[i] = dim(X[[i]])[1]
    if (columns) X[[i]]=t(X[[i]])
  }
  names(X) <- levellab
  return(X)
}
