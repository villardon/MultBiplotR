PrettyTicks <- function(min=-3, max=3, ntick=5){
  range=NiceNumber(max-min, FALSE)
  d=NiceNumber(range/(ntick-1), TRUE)
  graphmin=floor(min/d)*d
  graphmax=ceiling(max/d)*d
  nfrac=max(c(-1*floor(log10(d)), 0))
  ticks=round(seq(from=graphmin, to=graphmax+0.5*d, by = d), digits=nfrac)
  result=list(ticks=ticks, labels=as.character(ticks))
  return(result)
}