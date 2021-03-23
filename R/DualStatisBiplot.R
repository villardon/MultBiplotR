DualStatisBiplot <- function(X, InitTransform = "Standardize columns", dimens=2, SameInd=FALSE) {
  mycall=match.call()
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center", "Normalized residuals from independence")
  if (is.numeric(InitialTransform)) 
    InitialTransform = ContinuousDataTransform[InitialTransform]
  
  ng = length(X) #Number of groups
  StudyNames=names(X)
  nri = matrix(0, ng, 1) #Number of rows in each group
  nci = matrix(0, ng, 1) #Number of cols in each group
  for (i in 1:ng) {
    nri[i] = dim(X[[i]])[1]
    nci[i] = dim(X[[i]])[2]
  }
  nc = nci[1]
  if (sum(nci == nc) < ng) 
    stop("The number of columns (variables) must be the same in all ocassions")
  
  if (SameInd){
    nr = nri[1]
    if (sum(nri == nr) < ng) 
      stop("The number of individuals can not be the same in all ocassions (use SameInd=FALSE)")
  }
  #  Extracting the names of the occasions
  OccNames=names(X)
  if (is.null(OccNames)) {
    for (i in 1:ng) OccNames= c(OccNames, paste("Occasion_",i,sep=""))
  }
  
  # Initial transformation of data and calculation of statistics for the biplot
  
  BiplotStatis=MultiTableStatistics(X, dual=TRUE)
  BiplotStatis$call <- mycall
  BiplotStatis$Type="DualStatis"
  BiplotStatis$Title="Biplot induced by Dual-Statis"
  BiplotStatis$Initial_Transformation=InitTransform
  BiplotStatis$nrows=sum(nri)
  BiplotStatis$ncols=nc
  StatisRes=list()
  StatisRes$Title="Dual STATIS-ACT Biplot"
  StatisRes$Type="Dual STATIS-ACT"
  StatisRes$NTables=ng
  StatisRes$NCols=nc
  if (SameInd)
    StatisRes$NRows=nr
  else
    StatisRes$NRows=nri
  
  StatisRes$VarLabels=colnames(X[[1]])
  StatisRes$TableLabels=names(X)
  if (SameInd) StatisRes$RowLabels=rownames(X[[1]])
  
  X=MultiTableTransform(X, InitTransform = InitTransform, dual=TRUE)
  #Calculation of the objects
  Wt = list()
  for (i in 1:ng) {
    Wt[[i]] = t(X[[i]]) %*% X[[i]]
  }
  
  # Calculation of the scalar products
  P = matrix(0, ng, ng)
  for (i in 1:ng) for (j in 1:ng) {
    P[i, j] = sum(diag(Wt[[i]] %*% t(Wt[[j]])))
    P[j, i] = P[i, j]
  }
  rownames(P)=StudyNames
  colnames(P)=StudyNames
  
  #RV: correlations among the occasions
  RV = solve(diag(sqrt(diag(P)))) %*% P %*% solve(diag(sqrt(diag(P))))
  rownames(RV)=StudyNames
  colnames(RV)=StudyNames
  
  # This information will be common to all the techniques, then a single
  # function to plot is necassary
  StatisRes$Data = X
  StatisRes$Objects = Wt
  StatisRes$ScalarProducts = P
  StatisRes$RV = RV
  Inter = svd(RV)
  StatisRes$EigInter = Inter$d
  names(StatisRes$EigInter)=paste("Dim",1:length(StudyNames))
  
  StatisRes$InerInter=(Inter$d/sum(Inter$d))*100
  names(StatisRes$InerInter)=paste("Dim",1:length(StudyNames))
  StatisRes$InterStructure = -1 * Inter$u %*% diag(sqrt(Inter$d))
  rownames(StatisRes$InterStructure)=StudyNames
  colnames(StatisRes$InterStructure)=paste("Dim", 1:ng)
  # Weigths for the Compromise
  StatisRes$Weights = abs(Inter$u[, 1])/sum(abs(Inter$u[, 1]))
  
  # Compromise
  W=matrix(0,nc,nc)
  for (i in 1:ng)
    W=W + StatisRes$Weights[i]*Wt[[i]]
  StatisRes$Compromise=W
  
  # Euclidean Configuration of the compromise
  Intra = svd(W)
  r=sum(Intra$d>0.00001)
  
  BiplotStatis$EigenValues=Intra$d[1:r]
  BiplotStatis$Inertia=(Intra$d[1:r]/sum(diag(Intra$d[1:r])))*100
  BiplotStatis$CumInertia = cumsum(BiplotStatis$Inertia)
  BiplotStatis$alpha=1
  BiplotStatis$Dimension=dimens
  # Compromise-Consensus Coordinates and contributions
  BiplotStatis$ColCoordinates = Intra$u[,1:dimens] %*% diag(sqrt(Intra$d[1:dimens]))
  sf=apply((Intra$u[,1:r] %*% diag(sqrt(Intra$d[1:r])))^2,1,sum)
  BiplotStatis$ColContributions=(diag(1/sf) %*% BiplotStatis$ColCoordinates^2)*100
  rownames(BiplotStatis$ColCoordinates)=colnames(X[[1]])
  colnames(BiplotStatis$ColCoordinates)=paste("Dim", 1:dimens)
  rownames(BiplotStatis$ColContributions)=colnames(X[[1]])
  colnames(BiplotStatis$ColContributions)=paste("Dim", 1:dimens)
# ----------------Voy por aquí ---------
  # Trajectories for individuals (old kind)
  trajin=list()
  for (j in 1:ng)
    trajin[[j]] = Wt[[j]] %*% Intra$u[,1:r] %*% diag(sqrt(1/Intra$d[1:r]))
  
  StatisRes$TrajInd=list()
  for (i in 1:nr){
    Traj=NULL
    for (j in 1:ng)
      Traj=rbind(Traj , trajin[[j]][i,1:dimens])
    rownames(Traj)=StudyNames
    colnames(Traj)=paste("Dim", 1:dimens)
    StatisRes$TrajInd[[i]]=Traj
  }
  names(StatisRes$TrajInd)=rownames(X[[1]])
  
  # Squared cosines of the individual objects projected onto the compromise
  contbt=list()
  for (j in 1:ng){
    sf=apply(trajin[[j]]^2,1,sum)
    contbt[[j]] = (diag(1/sf) %*% trajin[[j]]^2)*100
  }
  
  StatisRes$ContribInd=list()
  for (i in 1:nr){
    Traj=NULL
    for (j in 1:ng)
      Traj=rbind(Traj , contbt[[j]][i,1:dimens])
    rownames(Traj)=StudyNames
    colnames(Traj)=paste("Dim", 1:dimens)
    StatisRes$ContribInd[[i]]=Traj
  }
  names(StatisRes$ContribInd)=rownames(X[[1]])
  
  # Coordinates and contributions of the variables on the biplot
  # First we calculate the principal coordinates in order to calculate the contributions
  trajvar=list()
  for (j in 1:ng)
    trajvar[[j]]=t(X[[j]]) %*% Intra$u[,1:r]
  # Coordinates of the variables on the biplot
  contvar=list()
  for (j in 1:ng){
    sf=apply(trajvar[[j]]^2,1,sum)
    contvar[[j]] = (diag(1/sf) %*% trajvar[[j]]^2)*100
    rownames(contvar[[j]])=paste(colnames(X[[j]]), "-", StudyNames[j], sep="")
    
  }
  
  # Standard coodinates to represent on the biplot (RMP-Biplot)
  # Standard coordinates are needed to represent the biplot axes with graded scales
  for (j in 1:ng){
    trajvar[[j]]=t(X[[j]]) %*% Intra$u[,1:r] %*% diag(1/Intra$d[1:r])
    rownames(trajvar[[j]])=paste(rownames(trajvar[[j]]), "-", StudyNames[j], sep="")
  }
  
  
  # This is just a reorganization of the column coordinates in order to represent them
  # as a biplot using the standard procedure
  
  BiplotStatis$ColCoordinates=NULL
  for (j in 1:ng)
    BiplotStatis$ColCoordinates=rbind(BiplotStatis$ColCoordinates, trajvar[[j]][,1:dimens])
  
  colnames(BiplotStatis$ColCoordinates)<- paste("Dim ",1:dimens)
  
  BiplotStatis$ColContributions=NULL
  for (j in 1:ng)
    BiplotStatis$ColContributions=rbind(BiplotStatis$ColContributions, contvar[[j]][,1:dimens])
  colnames(BiplotStatis$ColContributions)<- paste("Dim ",1:dimens)
  
  # Rescaling the biplot for better visualization. Used to be optional
  # but now I do for every biplot. 
  # The row coordinates are multiplied for a constant and the coumn coordinates divided
  # by the same constant. The scalar (inner product) in which biplot interpretation
  # is based does not change, but the visuatization is much better.
  # In a biplot with scales for each varibale, the rescaling does not affect 
  # the calculation of the scale marks.
  # The rescaling is performed with all the dimensiond retained in the final solution
  
  sca = sum(BiplotStatis$RowCoordinates^2)
  scb = sum(BiplotStatis$ColCoordinates^2)
  p=sum(nci)
  sca = sca/nr
  scb = scb/p
  scf = sqrt(sqrt(scb/sca))
  BiplotStatis$RowCoordinates = BiplotStatis$RowCoordinates * scf
  BiplotStatis$ColCoordinates = BiplotStatis$ColCoordinates/scf
  BiplotStatis$Scale_Factor=scf
  
  for (i in 1:nr){
    StatisRes$TrajInd[[i]]=StatisRes$TrajInd[[i]]*scf
  }
  
  BiplotStatis$ClusterType="us"
  BiplotStatis$Clusters = as.factor(matrix(1,nrow(BiplotStatis$RowContributions), 1))
  BiplotStatis$ClusterColors="blue"
  BiplotStatis$ClusterNames="Cluster 1"
  
  class(BiplotStatis) ="ContinuousBiplot"
  # If the variables are the same for all the occasions or studies
  # additional trajectories for the variables can be constructed
  # selecting conveniently the var coordinates. This are original
  # in our work, together with the simultaneous representation
  # of rows and columns.
  StatisRes$Biplot=BiplotStatis
  StatisRes$SameIndiv = TRUE
  StatisRes$SameInd = SameInd
  if (SameInd){
    p=nci[1]
    StatisRes$TrajVar=list()
    for (i in 1:p){
      Traj=NULL
      for (j in 1:ng)
        Traj=rbind(Traj , trajvar[[j]][i,1:dimens]/scf)
      rownames(Traj)=StudyNames
      colnames(Traj)=paste("Dim", 1:dimens)
      StatisRes$TrajVar[[i]]=Traj}
    names(StatisRes$TrajVar)=colnames(X[[1]])
  }
  
  # Calculation of the Inertia of each study (or occasion) accounted for the STATIS solution.
  # Comparing this with the inertia accounted by the biplot for each separate study
  # we have an index of the goodness of the consensus (or compromise) in explaining the 
  # study.
  StatisRes$Biplot$Structure = cor(StatisRes$Biplot$Non_Scaled_Data, StatisRes$Biplot$RowCoordinates)
  
  class(StatisRes) = "StatisBiplot"
  return(StatisRes)
}
