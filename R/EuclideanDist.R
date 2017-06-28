EuclideanDistance <- function (x) 
{
  n=nrow(x)
  dis=matrix(0,n,n)
  for (i in 1:(n-1)) for (j in (i+1):n) {  
      dis[i, j] = sqrt(sum((x[i,] - x[j,])^2))
      dis[j, i]=dis[i, j]}
  return(dis)
}