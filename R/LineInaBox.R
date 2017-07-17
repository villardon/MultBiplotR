LineInaBox <- function(a=0, b=1, xmin=0, xmax=1, ymin=0, ymax=1){
  # b is the slope of the line
  ini=NULL
  final=NULL
  b1=1
  b2=b
  x1 = xmin
  y1 = a + b * xmin
  if ((y1 >= ymin) & (y1 <= ymax)){
    if ((x1 * b1 + y1 * b2) <= 0)
      ini = c(x1, y1)
  else
    final = c(x1, y1)}
  
  x1 = xmax
  y1 = a + b * xmax
  if ((y1 >= ymin) & (y1 <= ymax)){
    if ((x1 * b1 + y1 * b2) <= 0)
      ini = c(x1, y1)
   else
    final = c(x1, y1)}
  
  
  x1 = (ymin-a)/b
  y1 = ymin
  if ((x1 >= xmin) & (x1 <= xmax)){
    if ((x1 * b1 + y1 * b2) <= 0)
      ini = c(x1, y1)
  else
    final = c(x1, y1)}
  
  x1 = (ymax-a)/b
  y1 = ymax
  if ((x1 >= xmin) & (x1 <= xmax)){
    if ((x1 * b1 + y1 * b2) <= 0)
      ini = c(x1, y1)
  else
    final = c(x1, y1)}
  return(list(Ini=ini, Final=final))
  
}

