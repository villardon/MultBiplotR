PlotBiplotClusters <- function(A, Groups = ones(c(nrow(A), 1)), TypeClus = "st",
                               ClusterColors = NULL, ClusterNames = NULL, centers =
                                 TRUE, ClustConf = 1, Legend = FALSE, LegendPos =
                                 "topright", CexClustCenters=1, ...){
  
  TypeClusters=c("ch", "el", "st", "oc")
  if (is.numeric(TypeClus)) TypeClus=TypeClusters[TypeClus]
  g = length(levels(Groups))
  
  if (is.null(ClusterColors)) {
    if (g > 1) {
      palette(rainbow(g))
      ClusterColors = palette()
    } else ClusterColors = "blue"
  }
  
  if (is.null(ClusterNames)) {
    ClusterNames=levels(Groups)}
  
  levellab = levels(Groups)
  Sizes = zeros(c(g, 1))
  Means = zeros(c(g, 2))
  X = list()
  for (i in 1:g) {
    X[[i]] = A[which(Groups == levellab[i]), ]
    Sizes[i] = length(which(Groups == levellab[i]))
    if (is.matrix(X[[i]]))
    Means[i, ] = apply(X[[i]], 2, mean)
    else
    Means[i, ] = X[[i]]
  }
  
  
  ColorInd=matrix(0,dim(A)[1],1)
  for (i in 1:(dim(A)[1]))
    ColorInd[i]=ClusterColors[as.integer(Groups[i])]
  
  if (centers){
    Centers=matrix(0,g,2)
    Centers[,1]=tapply(A[,1],Groups,mean)
    Centers[,2]=tapply(A[,2],Groups,mean)
    points(Centers[,1], Centers[,2], cex=1.5, pch=16, col=ClusterColors)
    text(Centers[,1], Centers[,2], labels=ClusterNames, col =ClusterColors, cex=CexClustCenters, pos =4)
  }
  
  for (i in 1:g)
    if (Sizes[i]>2){
      if (TypeClus == "el"){
        elip=ConcEllipse(X[[i]], confidence=ClustConf)
        plot(elip, col=ClusterColors[i])
      }
      else{
        fr=Fraction(X[[i]], confidence=ClustConf)
        plot(fr, type=TypeClus, col=ClusterColors[i])
      }
    }
  
  LegendPositions= c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")
  if (is.numeric(LegendPos)) LegendPos=LegendPositions[LegendPos]
  if (Legend){
    legend(x=LegendPos, ClusterNames , col = ClusterColors,text.col = "black", pch = 16, cex=0.7)
  }
  
  return(ColorInd)
  
}

