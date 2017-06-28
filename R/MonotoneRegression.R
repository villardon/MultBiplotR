# Performs weighted isotonic regression using
# the non-negative weights in W.
MonotoneRegression <- function(x, y, w = NULL) {
  if (is.null(w)) 
    w = matrix(1, length(x))
  n = length(x)
  xyord = cbind(x, y)
  ord = do.call(order, as.data.frame(cbind(x, y)))
  xyord = xyord[ord, ]
  iord = 1:n
  iord[ord] = 1:n
  
  # Initialize fitted values to the given values.
  yhat = xyord[, 2]
  block = 1:n
  w = w[ord]
  
  # Merge zero-weight points with preceding pos-weighted point (or
  # with the following pos-weighted point if at start).
  
  posWgts = (w > 0)
  if (any(!posWgts)) {
    idx = cumsum(posWgts)
    idx[(idx == 0)] = 1
    w = w[posWgts]
    yhat = yhat[posWgts]
    block = idx[block]
  }
  # If all blocks are monotonic, then we're done.
  diffs = diff(yhat)
  monotonic = all(diffs >= 0)
  
  while (!monotonic) {
    
    
    # Otherwise, merge blocks of non-increasing fitted values, and set the
    # fitted value within each block equal to a constant, the weighted mean
    # of values in that block.
    
    idx = cumsum(c(1, as.numeric(diffs > 0)))
    nblocks = max(idx)
    sumyhat = 1:nblocks
    wnew = 1:nblocks
    for (i in 1:nblocks) {
      sumyhat[i] = sum(w[idx == i] * yhat[idx == i])
      wnew[i] = sum(w[idx == i])
    }
    w = wnew
    yhat = sumyhat/w
    block = idx[block]
    diffs = diff(yhat)
    monotonic = all(diffs >= 0)
    
  }
  
  # Broadcast merged blocks out to original points, and put back in the
  # original order and shape.
  yhat = yhat[block]
  yhat = yhat[iord]
  return(yhat)
}
