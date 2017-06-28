"dlines" <-
  function(SetA,SetB,lin="dotted",color="black", ...) {
    np<-nrow(SetA)
    if (length(color)==1) color = rep(color, np)
    for(i in 1:np) lines(rbind(SetA[i,1:2],SetB[i,1:2]),lty=lin,col=color[i], ...) 
    return(NULL)
  }