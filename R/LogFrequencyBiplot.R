
LogFrequencyBiplot <- function(x, Scaling=2, logoffset=1, freqoffset=logoffset,  ...){
   if (is.data.frame(x)) x=as.matrix(x)
   x=log(x+logoffset)
   W=x+freqoffset
   ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", "Column centering", "Standardize columns", "Row centering", 
                               "Standardize rows", "Divide by the column means and center", "Normalized residuals from independence")
   if (is.numeric(Scaling)) 
     Scaling = ContinuousDataTransform[Scaling]
   if (Scaling=="Double centering") stop("Double Centering is not compatible with log-frequencies")
   if (Scaling=="Divide by the column means and center") stop("Divide by the column means and center is not compatible with log-frequencies")
   if (Scaling=="Normalized residuals from independence") stop("Normalized residuals from independence are not compatible with log-frequencies")
   Data = InitialTransform(x, transform = Scaling)
   X = Data$X
   
   NewBiplot=CrissCross(X,w=W, ...)
   NewBiplot$Non_Scaled_Data=x
   NewBiplot$Data=X
   NewBiplot$Means=apply(x,2,mean)
   NewBiplot$Medians=apply(x,2,median)
   NewBiplot$Deviations=apply(x,2,sd)
   NewBiplot$Minima=apply(x,2,min)
   NewBiplot$Maxima=apply(x,2,max)
   NewBiplot$Gmean=mean(x)
   NewBiplot$alpha=0.5
   NewBiplot$Initial_Transformation = Scaling
   
   NewBiplot$Type="LogFreqBiplot"
  return(NewBiplot)
  
}

