wcor=function (d1,d2, w = rep(1, nrow(d1))/nrow(d1)) 
{
  s <- sum(w)
  m1 <- sum(d1 * w)/s
  m2 <- sum(d2 * w)/s
  (sum(d1 * d2 * w)/s - m1 * m2)/sqrt((sum(d1^2 * w)/s - m1^2) * (sum(d2^2 * w)/s - m2^2))
}
