AddCluster2Biplot <- function(Bip, NGroups=3, ClusterType="hi", Groups=NULL, Original=FALSE, ...){
  types = c("Hierarchical", "K-Means", "Gaussian Mixture", "User Privided")

  ClusterTypes=c("hi", "km", "gm", "us")
  if (is.numeric(ClusterType))
    ClusterType = ClusterTypes[ClusterType]

  Bip$ClusterType=ClusterType

  if (ClusterType=="hi"){
    if (NGroups<2) stop("You need at least 2 groups")

    if (Original) {
      if (class(Bip)=="ContinuousBiplot")
        distances = dist(Bip$Scaled_Data)
      if (class(Bip)=="External.Binary.Logistic.Biplot")
        distances=as.dist(Bip$Proximities)}
    else distances = dist(Bip$RowCoordinates)

    hc=hclust(distances, ...)
    Bip$Clusters = as.factor(cutree(hc,NGroups))
    Bip$ClusterNames="Cluster 1"
    for (k in 2:NGroups)
      Bip$ClusterNames =c(Bip$ClusterNames, paste("Cluster",k))
    hcd=as.dendrogram(hc)
    Bip$ClusterObject=hc

  }

  if (ClusterType=="km"){

    if (Original)
      hc=kmeans(Bip$Scaled_Data, centers=NGroups, ...)
      else
      hc=kmeans(Bip$RowCoordinates, centers=NGroups, ...)

    Bip$ClusterObject=hc
    Bip$ClusterNames="Cluster 1"
    for (k in 2:NGroups)
      Bip$ClusterNames =c(Bip$ClusterNames, paste("Cluster",k))
    Bip$Clusters = as.factor(hc$cluster)
  }

  if (ClusterType=="gm"){
    if (Original)
      hc=MGC(Bip$Scaled_Data, NG=NGroups, ...)
    else
      hc=MGC(Bip$RowCoordinates, NG=NGroups, ...)

    Bip$Clusters = as.factor(hc$Classification)
    Bip$P=hc$P
  }

  if (ClusterType=="us"){
    if (!is.factor(Groups)) stop("Groups have to be defined as a factor")
    NGroups=length(levels(Groups))
    Bip$Clusters = Groups
    Bip$ClusterNames = levels(Groups)
  }

  palette(rainbow(NGroups))
  ClusterColors = palette()

  Bip$ClusterColors=ClusterColors

  if (ClusterType=="hi"){
    clusMember = cutree(hc, NGroups)
    colLab <- function(n) {
      if (is.leaf(n)) {
        a <- attributes(n)
        labCol <- ClusterColors[clusMember[which(names(clusMember) == a$label)]]
        attr(n, "nodePar") <- c(a$nodePar, lab.col = labCol)
      }
    n
    }
    Bip$Dendrogram = dendrapply(hcd, colLab)
  }
  return(Bip)
}



