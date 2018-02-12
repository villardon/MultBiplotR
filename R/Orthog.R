Orthog<-function (y, order, recenter = TRUE, rescale = TRUE, adjnames = TRUE) 
{
  sd <- function(x, na.rm = FALSE) {
    if (is.matrix(x)) 
      apply(x, 2, sd, na.rm = na.rm)
    else if (is.vector(x)) 
      sqrt(var(x, na.rm = na.rm))
    else if (is.data.frame(x)) 
      sapply(x, sd, na.rm = na.rm)
    else sqrt(var(as.vector(x), na.rm = na.rm))
  }
  n <- nrow(y)
  if (missing(order)) 
    order <- 1:ncol(y)
  y <- y[, order]
  p <- ncol(y)
  if (is.data.frame(y)) {
    numeric <- unlist(lapply(y, is.numeric))
    if (!all(numeric)) 
      stop("all columns of y must be numeric")
  }
  ybar <- colMeans(y)
  ysd <- sd(y)
  z <- scale(y, center = TRUE, scale = FALSE)
  z <- qr.Q(qr(z))
  zsd <- sd(z)
  if (rescale) 
    z <- z %*% diag(ysd/zsd)
  if (recenter) 
    z <- z + matrix(rep(ybar, times = n), ncol = p, byrow = TRUE)
  rownames(z) <- rownames(y, do.NULL = FALSE)
  colnames(z) <- colnames(y, do.NULL = FALSE)
  if (adjnames) {
    for (j in 2:p) {
      colnames(z)[j] <- paste(colnames(z)[j], ".", sep = "", 
                              paste(1:(j - 1), collapse = ""))
    }
  }
  z
}
