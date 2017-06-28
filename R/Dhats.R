Dhats <- function(P,D,W, Model=c("Identity", "Ratio", "Interval", "Ordinal"), Standardize=TRUE){
  if (length(Model)>1) Model=Model[1]
  P=as.dist(P)
  D=as.dist(D)
  W=as.dist(W)

  if (Model =="Identity"){
    Dh= P}
  
  if (Model =="Ratio"){
    b=sum(P*W*D) / sum((P^2)*W)
    Dh=P*b
  }
  
  if (Model =="Interval"){
    xb=sum(P*W)/sum(W)
    yb=sum(D*W)/sum(W)
    b=sum((P-xb)*W*(D-yb)) / sqrt(sum(((P-xb)^2)*W)*sum(((D-yb)^2)*W))
    a=yb-b*xb
    Dh=P*b
  }
  
  if (Model =="Ordinal"){
    yhat=MonotoneRegression(P,D)
    no=nrow(as.matrix(D))
    Dh=yhat*as.dist(matrix(1,no,no))
  }
  
  if (Standardize)
  Dh=((Dh*W)/sqrt(sum((Dh^2)*W)))*sqrt(sum(W))
  
  return(Dh)
}