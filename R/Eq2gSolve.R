Eq2gSolve <- function(a,b,c){
  solns<-matrix(rep(NA,length(a)))
  solns<-cbind(solns,solns)
  colnames(solns)<-c("soln 1","soln 2")
  if(a==0 && b!=0){
    solns[1,1]= (-1)*c/b
    solns[1,2]= (-1)*c/b;
  }else if(a==0 && b==0){
    print("Coefficients a and b of the polinomial are zero")
  }else if((b^2 -4*a*c) < 0){
    print("Polinomial has complex roots")
  }else{
    solns[1,1]<-((-1)*b + sqrt(b^2 - 4*a*c))/(2*a)
    solns[1,2]<-((-1)*b - sqrt(b^2 - 4*a*c))/(2*a)
  }
  solns
}