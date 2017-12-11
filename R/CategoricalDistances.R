# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca

CategoricalDistances <- function(x, y=NULL, coefficient= "GOW", transformation="sqrt(1-S)") {
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
  
  for (i in 1:px)
  ncat[i]=length(levels(x[[i]]))
  
  
  switch(coefficient, GOW = {
    sim=matrix(0, nx,ny)
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=sum(x[i,]==y[j,])
      }
    sim = sim/px
  }, ESK = {
    

    sim=matrix(0, nx,ny)
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+1
        else sim[i,j]=sim[i,j] + (ncat[k]^2/(ncat[k]^2+2))
      }
    sim = sim/px
    
  }, IOF = {
    
    Freq=list()
    for (j in 1:px)
      Freq[[j]]=table(x[,j])
    
    sim=matrix(0, nx,ny)
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+1
        else sim[i,j]=sim[i,j] + (1/(1+log10(ncat[k]/Freq[[k]][x[i,k]]) * log10(ncat[k]/Freq[[k]][y[i,k]])))
      }
    sim = sim/px
    
  }, OF = {
    
    Freq=list()
    for (j in 1:px)
      Freq[[j]]=table(x[,j])
    sim=matrix(0, nx,ny)
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+1
        else sim[i,j]=sim[i,j] + (1/(1+log10(Freq[[k]][x[i,k]]) * log10(Freq[[k]][y[i,k]])))
      }
    sim = sim/px
    
  }, GOO1 = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+ ( 1- sum(Prob[[k]][which(Prob[[k]]<=Prob[[k]][x[i,k]])]^2))
      }
    
    sim = sim/px
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
  }, GOO2 = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+ ( 1- sum(Prob[[k]][which(Prob[[k]]>=Prob[[k]][x[i,k]])]^2))
      }
    
    sim = sim/px
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
  }, GOO3 = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+ ( 1- Prob[[k]][x[i,k]]^2)
      }
    
    sim = sim/px
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
  }, GOO4 = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]+ (Prob[[k]][x[i,k]]^2)
      }
    
    sim = sim/px
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
  }, GAM = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    w=sum(1/ncat)
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j] - w * ((Prob[[k]][x[i,k]]* log2(Prob[[k]][x[i,k]])) + ((1-Prob[[k]][x[i,k]])*log2(1-Prob[[k]][x[i,k]])))
      }
    
    sim = sim/(px*w)
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
    
  }, LIN = {
    
    Prob=list()
    for (j in 1:px)
      Prob[[j]]=table(x[,j])/nx
    
    sim=matrix(0, nx,ny)
    
    # Not sure I'm doing this right
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px){
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j]- 2*log(Prob[[k]][x[i,k]])
          else sim[i,j]=sim[i,j]- 2*log(Prob[[k]][x[i,k]]+Prob[[k]][y[i,k]])
        }
      }
    sim = sim/px
    print(sim)
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
    
  }, AND = {
    # Not ready yet
  }, SMI = {
    
    Freq=list()
    for (j in 1:px)
      Freq[[j]]=table(x[,j])
    
    w=sum(1/ncat)
    sim=matrix(0, nx,ny)
    
    for (i in 1:nx)
      for (j in 1:ny){
        sim[i,j]=0
        for (k in 1:px)
          if (x[i,k]==y[j,k]) sim[i,j]=sim[i,j] - w * ((Prob[[k]][x[i,k]]* log2(Prob[[k]][x[i,k]])) + ((1-Prob[[k]][x[i,k]])*log2(1-Prob[[k]][x[i,k]])))
      }
    
    sim = sim/(px*w)
    if (!Supl) sim=sim-diag(diag(sim))+diag(rep(1,nx))
    
  })
  
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

