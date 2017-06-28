NiceNumber <- function(x=6, round=TRUE){
  expon = floor(log10(x))
  f = x/(10^expon)
  if (round){ 
    if (f<1.5) nf=1
    else if (f < 3) nf=2
    else if (f < 7) nf=5
    else nf=10}
  else{
    if (f<1) nf=1
    else if (f < 2) nf=2
    else if (f < 5) nf=5
    else nf=10
  }
  result=nf*(10^expon)
  return(result)
}