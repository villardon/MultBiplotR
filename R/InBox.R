InBox <- function(x, y, xmin, xmax, ymin, ymax) {
  result = FALSE
  if ((x >= xmin) & (x <= xmax) & (y >= ymin) & (y <= ymax)) 
    result = TRUE
  return(result)
}
