
Games_Howell <- function(data, group)	{
  OK <- complete.cases(data, group)
  data <- data[OK]
  group <- factor(group[OK])
  n <- tapply(data, group, length)
  a <- length(n)
  phi.e <- sum(n)-a
  Mean <- tapply(data, group, mean)
  Variance <- tapply(data, group, var)
    t.df <- combn(a, 2, function(ij) {
      t <- abs(diff(Mean[ij]))/sqrt(sum(Variance[ij]/n[ij]))
      df <- sum(Variance[ij]/n[ij])^2/sum((Variance[ij]/n[ij])^2/(n[ij]-1))
      return(c(t, df))} )
    t <- t.df[1,]
    df <- t.df[2,]
    p <- ptukey(t*sqrt(2), a, df, lower.tail=FALSE)
    Games.Howell <- cbind(t, df, p)
    rownames(Games.Howell) <- combn(levels(group), 2, paste, collapse=":")
    return(Games.Howell=Games.Howell)
}
