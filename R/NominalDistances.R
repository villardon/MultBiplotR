#X<-pruebanominal1[,-5]
#la linea siguiente es para correr la funcion
#S<-DistNom(X,method=5,similarity=TRUE,diag=TRUE,upper=TRUE)

#programita para la distancia overlap primera aproximacion
NominalDistances<-function(X,method=1,diag=FALSE,upper=FALSE,similarity=TRUE){
  if(class(X)!="data.frame"&class(X)!="matrix"){
    stop("X is not a valid object")
  }
  
  for (i in 1: ncol(X)){
    if(class(X[,i])!="factor")
      stop("X is not a factor")  
  }
  
  METHODS = c("Overlap", "Eskin", "IOF", "OF","Goodall3","Lin")
  if (all((1:6) != method)) {
    cat("1 = Overlap \n")
    cat("2 = Eskin \n")
    cat("3 = IOF \n")
    cat("4 = OF \n")
    cat("5 = Goodall3\n")
    cat("6 = Lin\n")
    cat("Select an integer (1-6): ")
    method <- as.integer(readLines(n = 1))
  }
  if (all((1:6) != method)) 
    (stop("Non convenient method number"))  
  
  A<-function(X1){
    counts<-length(unique(X1))
    return(counts)
  }
  nf<-apply(X,2,A)
  #nf numero de atributos de cada categoria
  w=1/ncol(X)
  #w pesos para algunas de las funciones
  
  #fk la frecuencia del atributo
  #ver como la voy a contar a la frecuencia del atributo 
  #para que sea conveniente
  #is the table function 
  #frecuencies
  frec<-list()
  prob<-list()
  N<-c()
  prob2<-list()
  Logprob<-list()
  for (i in 1:ncol(X)){
    frec[[i]]<-table(X[,i]) 
    N[i]<-sum(table(X[,i]))
    prob[[i]]<-frec[[i]]/N[i]
    Logprob[[i]]<-log(prob[[1]])
    prob2[[i]]<-(frec[[i]]*(frec[[i]]-1))/(N[i]*(N[i]-1))
  }
  
  #esto puede servir para devolver la lista
  #esto es para la similaridad IOF
  Logar<-lapply(frec,log)
  
  SumLogar<-lapply(Logar,prod)
  
  Z1<-function(Z){
    Z<-Z+1
    return(Z)
  }
  SumLogar1<-lapply(SumLogar,Z1)
  Z2<-function(A){
    
    A<-N[1]/A
    return(A)
  }
  
  frecinvn<-lapply(frec,Z2)
  Logar<-lapply(frecinvn,log)
  ProdLogar<-lapply(Logar,prod)
  ProdLogar1<-lapply(ProdLogar,Z1)
  sumprob<-lapply(prob,sum)
  
  Logprob2<-lapply(prob2,log)
  
  
  
  D<-matrix(,nrow=nrow(X),ncol=nrow(X))
  
  if(method==1){
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            a[k]<-w*1
          }
          else{
            a[k]<-w*0
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    D<-as.dist(D)
  }
  
  #programita para la distancia Eskin primera aproximacion
  
  else if (method==2){
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            a[k]<-w*1
          }
          else{
            a[k]<-w*((nf[k]^2)/((nf[k]^2)+2))
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    
    D<-as.dist(D)
  }
  
  else if (method==3){
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            a[k]<-w*1
          }
          else{
            a[k]<-w*(1/SumLogar1[[k]][1])
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    
    D<-as.dist(D)
  }
  
  else if (method==4){
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            a[k]<-w*1
          }
          else{
            a[k]<-w*(1/ProdLogar1[[k]][1])
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    
    D<-as.dist(D)
  }
  
  
  else if (method==5){
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            Z<-as.list(prob2[[k]])
            
            a[k]<-w*(1-Z[[as.character(X[i,k])]])
          }
          else{
            a[k]<-0
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    
    D<-as.dist(D)
  }
  #Lin
  else if (method==6){
    #lin
    w1<-1/sum(unlist(Logprob))
    
    for (i in 1:(nrow(X)-1)){
      for(j in (i+1):nrow(X)){
        a<-c()
        for(k in 1:ncol(X)){
          if(as.character(X[i,k])==as.character(X[j,k])){
            Z<-as.list(prob[[k]])
            
            a[k]<-w1*2*(log(Z[[as.character(X[i,k])]]))
          }
          else{
            Z<-as.list(prob[[k]])
            a[k]<-w1*2*(log(Z[[as.character(X[i,k])]]+Z[[as.character(X[j,k])]]))
          }
        }
        
        D[j,i]<-sum(a)
      }
    }
    D<-as.dist(D)
  }
  #tengo que revisar bien estos algoritmos pero creo que esta bien
  #con esto ya puedo hacer la funcion (de manera mas eficiente)
  
  Unos<-matrix(rep(1,nrow(X)*nrow(X)),ncol=nrow(X),nrow=nrow(X))
  Unos<-as.dist(Unos)
  D1=sqrt(1-D)
  if(similarity==TRUE){
    attr(D, "Size") <- nrow(X)
    attr(D, "Labels") <- rownames(X)
    attr(D, "Diag") <- diag
    attr(D, "Upper") <- upper
    attr(D, "method") <- METHODS[method]
    attr(D, "call") <- match.call()
    class(D) <- "dist"
    message("D is a similarity matrix. If you want an distance matrix, should write as argument in the function 'similarity = FALSE'")
    return(D)
  }
  if(similarity==FALSE){
    attr(D1, "Size") <- nrow(X)
    attr(D1, "Labels") <- rownames(X)
    attr(D1, "Diag") <- diag
    attr(D1, "Upper") <- upper
    attr(D1, "method") <- METHODS[method]
    attr(D1, "call") <- match.call()
    class(D1) <- "dist"
    message("D is an distance matrix")
    return(D1)
  }
}