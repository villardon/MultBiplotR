# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

GowerSimilarities <- function(x, y=NULL, Classes=NULL, transformation="sqrt(1-S)", BinCoef= "Simple_Matching", ContCoef="Gower", NomCoef="GOW", OrdCoef="GOW") {
  
  n = dim(x)[1]
  p = dim(x)[2]
  n1 = dim(y)[1]
  
  sim=matrix(0,n,n1)
  
  ng=0
  
  Numericas = which(Classes=="numeric")
  if (length(Numericas)>0){
    ng=ng+1
  NumX=as.matrix(x[, Numericas])
  NumY=as.matrix(y[, Numericas])
  NumDis=ContinuousDistances(NumX, NumY,  coef = ContCoef)
  if (ContCoef != "Gower") NumDis=NumDis/max(NumDis)
  else
  NumSim=(length(Numericas) - NumDis)/length(Numericas)
  sim=sim+NumSim}
  
  Binarias=which(Classes=="binary")
  if (length(Binarias)>0){
    ng=ng+1
  BinX=as.matrix(x[, Binarias])
  BinY=as.matrix(y[, Binarias])
  BinSim=BinaryDistances(BinX, BinY, coefficient= "Simple_Matching", transformation="Identity")
  sim=sim+BinSim}
  

  Factores=which(Classes=="factor")
  if (length(Factores)>0){
    ng=ng+1
  CatX=as.data.frame(x[, Factores])
  CatY=as.data.frame(y[, Factores])
  CatSim <- CategoricalDistances(CatX, y=CatY, coefficient= NomCoef, transformation="Identity")
  sim=sim+CatSim}
  
  
  Ordinales=which(Classes=="ordered")
  if (length(Factores)>0){
    ng=ng+1
  OrdX=as.data.frame(x[, Ordinales])
  OrdY=as.data.frame(y[, Ordinales])
  OrdSim <- OrdinalDistances(CatX, y=CatY, coefficient= OrdCoef, transformation="Identity")
  sim=sim+OrdSim}
  
  sim=sim/ng
  
  
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
  
  return(dis)
  
}




