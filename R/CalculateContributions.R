CalculateContributions <- function(X, A, B){
  I=dim(X)[1]
  J=dim(X)[2]
  I1=dim(A)[1]
  S=dim(A)[2]
  J1=dim(B)[1]
  S1=dim(B)[2]
  if (!I==I1) stop("Number of rows of the data and coordinates must be the same")
  if (!J==J1) stop("Number of columns of the data and coordinates must be the same")
  if (!S==S1) stop("Number of columns of both coordinate matrices must be the same")
  CumRowContributions=matrix(0,I,S)
  CumColContributions=matrix(0,J,S)
  CumFit=matrix(0,S,1)
  for (i in 1:S){
    Xesp= A[,1:i] %*% t(B[,1:i])
    CumFit[i]=sum(Xesp^2)/sum(X^2)
    CumRowContributions[,i] = apply(Xesp^2,1,sum)/apply(X^2,1,sum)
    CumColContributions[,i] = apply(Xesp^2,2,sum)/apply(X^2,2,sum)
  }
  
  RowContributions=CumRowContributions
  ColContributions=CumColContributions
  Fit=CumFit
  
  for (i in 2:S){
    Fit[i]=CumFit[i]-CumFit[i-1]
    RowContributions[,i] = CumRowContributions[,i] - CumRowContributions[,(i-1)] 
    ColContributions[,i] = CumColContributions[,i] - CumColContributions[,(i-1)]
  }
  
  RowContributions[which(RowContributions<0)]=0
  ColContributions[which(RowContributions<0)]=0
  
  RowContributions=round(RowContributions,digits=4)*100
  ColContributions=round(ColContributions,digits=4)*100
  
  rownames(RowContributions)=rownames(X)
  colnames(RowContributions)=colnames(A)
  rownames(ColContributions)=colnames(X)
  colnames(ColContributions)=colnames(A)
  colnames(Fit)="Percentage"
  rownames(Fit)=colnames(A)
  
  Structure=round(cor(X,A),digits=4)
    result=list(Fit=Fit, RowContributions=RowContributions, ColContributions=ColContributions, Structure=Structure)
  return(result)
}