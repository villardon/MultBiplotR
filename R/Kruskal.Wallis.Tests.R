Kruskal.Wallis.Tests <- function(X, groups, posthoc="none", alternative="two.sided", digits=4){

  n=dim(X)[1]
  p=dim(X)[2]

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  Levels=levels(groups)
  varnames=colnames(X)
  txt=capture.output(cat("\n\n\n***** KRUSKAL-WALLIS WITH POST-HOC COMPARISOSNS *****\n\n"))

  Comparisons=combn(levels(groups), 2, paste, collapse=":")

  Summary.Kruskal=matrix(0, p, 3)

  Summary.posthoc=list()
  for (j in 1:length(posthoc)) {
    Summary.posthoc[[j]]=matrix(0, p , g*(g-1)/2)
    rownames(Summary.posthoc[[j]])=varnames
    colnames(Summary.posthoc[[j]])=Comparisons
  }
  names(Summary.posthoc)=posthoc

  for (i in 1:p){
    txt=c(txt, capture.output(cat("*************************************************************************\n\n")))
    txt=c(txt, capture.output(cat("******* Variable : ", varnames[i], "-------------\n\n")))
    av=kruskal.test(X[,i], g=groups)

    Summary.Kruskal[i,1] = round(av$statistic, digits=digits)
    Summary.Kruskal[i,2] = round(av$parameter, digits=digits)
    Summary.Kruskal[i,3] = round(av$p.value, digits=digits)
    txt=c(txt, capture.output(print(av)))

    for (j in 1:length(posthoc)){
        ph=dunn.test(X[,i], g=groups, method=posthoc[j], kw=FALSE)
        txt=c(txt, capture.output(ph=dunn.test(X[,i], g=groups, method=posthoc[j], kw=FALSE)))
        Summary.posthoc[[j]][i,]=round(ph$P.adjusted, digits=digits)
      }
    txt=c(txt, capture.output(cat("*************************************************************************\n\n")))
  }
  txt=c(txt, capture.output(cat("\n\n SUMMARY OF RESULTS ****************************************************\n\n")))
  rownames(Summary.Kruskal)<- varnames
  colnames(Summary.Kruskal) <- c("Chi-squared", "df", "p.value")
  txt=c(txt, capture.output(cat("Kruskal-Wallis ****************************************************\n\n")))
  txt=c(txt, capture.output(print(Summary.Kruskal)))
  txt=c(txt, capture.output(cat("\n\n POST-HOC COMPARISONS - Dunn Test *************************************************\n\n")))
  testnames=names(Summary.posthoc)
  for (i in 1:length(Summary.posthoc)){
  txt=c(txt, capture.output(print(testnames[i])))
  txt=c(txt, capture.output(print(t(Summary.posthoc[[i]]))))
  }
  txt=c(txt, capture.output(cat("*************************************************************************\n\n")))


  return(txt)
}

