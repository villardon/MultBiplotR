
plot.NomSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                  ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ...){
  x$penalization
  for (i in 1:p){
    nameVar = "gears"
    nominalVar = x$Y[[1]]
    estimRows = bb$RowCoordinates
    library(NominalLogisticBiplot)
    plotNominalVariable(nameVar,nominalVar,estimRows,planex = 1,planey = 2,
                        xi=xmin,xu=xmax,yi=ymin,yu=imax,CexVar=0.7,ColorVar="magenta",PchVar=0.7,
                        addToPlot=TRUE,QuitNotPredicted=TRUE,ShowResults=TRUE,
                        linesVoronoi=TRUE,LabelVar=TRUE,tol = 1e-04, maxiter = 100,
                        penalization = 0.2)
  }
}
