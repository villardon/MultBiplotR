# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca

OrdinalDistances <- function(x, y=NULL, coefficient= "GOW", transformation="sqrt(1-S)") {
  # if (!is.matrix(x)) stop("Input must be a matrix")
  
  coefficients = c("GOW", "ESK", "IOF", "OF", "GOO1", "GOO2", "GOO3", "GOO4", "GAM", "LIN", 
                   "AND", "SMI")
  if (is.numeric(coefficient)) coefficient=coefficients[coefficient]
  
  
  nx=dim(x)[1]
  px=dim(x)[2]
  if (is.null(y)) {y=x
  Supl=FALSE}
  ny=dim(y)[1]
  py=dim(y)[2]
  
  
  ncat=rep(0,px)
  ranks=rep(0,px)
  for (i in 1:px){
    ncat[i]=length(levels(x[[i]]))
    ranks[i]= max(as.numeric(c(x[[i]],y[[i]])))-min(as.numeric(c(x[[i]],y[[i]])))
  }
  
  sim=matrix(0, nx,ny)
  for (i in 1:nx)
    for (j in 1:ny){
      sim[i,j]=1-abs(as.numeric(x[i,])-as.numeric(y[j,]))/ranks
    }
  sim = sim/px




transformations= c("Identity", "1-S", "sqrt(1-S)", "-log(s)", "1/S-1", "sqrt(2(1-S))", "1-(S+1)/2", "1-abs(S)", "1/(S+1)")
switch(transformation, `Identity` = {
  dis=sim
}, `Identity` = {
  dis=sim
}, `1-S` = {
  dis=1-sim
}, `sqrt(1-S)` = {
  dis = sqrt(1 - sim)
}, `-log(s)` = {
  dis=-1*log(sim)
}, `1/S-1` = {
  dis=1/sim -1 
}, `sqrt(2(1-S))` = {
  dis== sqrt(2*(1 - sim))
}, `1-(S+1)/2` = {
  dis=1-(sim+1)/2
}, `1-abs(S)` = {
  dis=1-abs(sim)
}, `1/(S+1)` = {
  dis=1/(sim)+1
})

rownames(dis)=rownames(x)
colnames(dis)=rownames(x)
return(dis)
}

