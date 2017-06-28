OrthogonalizeScores <- function(scores){
  dimens=dim(scores)[2]
  n=dim(scores)[1]
  for (i in 1:dimens)
    scores[,i]=scores[,i]-mean(scores[,i])
  cov=t(scores) %*% scores
  dec=svd(cov)
  newscores = sqrt(n) * scores %*% dec$v %*% solve(diag(sqrt(dec$d)))
  return(newscores)
}