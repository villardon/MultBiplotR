PlotOrdinalResponses <- function(olb, A1=1, A2=2, inf = -12, sup = 12, Legend=TRUE, WhatVars=NULL){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  A = olb$RowCoordinates[, c(A1, A2)]
  B = olb$ColumnParameters$coefficients[, c(A1, A2)]
  names=rownames(B)
  thresholds=olb$ColumnParameters$thresholds
  
  n = dim(A)[1]
  p = dim(olb$ColumnParameters$coefficients)[1]
  if (is.null(WhatVars)) WhatVars=1:p
  nf = ceiling(sqrt(length(WhatVars)))
  nc = ceiling(length(WhatVars)/nf)
  olb$Communalities=round(olb$Communalities, digits=3)
  op <- par(mfrow=c(nf,nc))
  for (j in WhatVars)
  OrCoor=OrdVarCoordinates(tr=thresholds[j,1:(olb$Ncats[j]-1)], c(B[j, 1], B[j, 2]), 
                           inf = inf, sup = sup, plotresponse=T, label=names[j], labx=paste("Comunality =", olb$Communalities[j]), 
                           catnames=olb$CategoryNames[[j]], Legend=Legend)
}