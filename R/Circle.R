Circle <- function (radius = 1, origin = c(0, 0), col=1, ...) 
{
  t <- seq(-pi, pi, by = 0.01)
  a <- origin[1]
  b <- origin[2]
  r <- radius
  x <- a + r * cos(t)
  y <- b + r * sin(t)
  points(x, y, type = "l",col=col, ...)
}