GetBiplotScales <- function(Biplot, nticks = 4, TypeScale = "Complete", ValuesScale = "Original") {
  # TypeScale = c("Complete", "StdDev", "BoxPlot")
  # ValuesScale = c("Original", "Transformed")
  # OptimalScales
  p = Biplot$ncols
  Ticks = list()
  Labels = list()
  
  if ((Biplot$Type=="LogFreqBiplot") & (ValuesScale == "Original")){
    Minima=exp(Biplot$Minima)
    Maxima=exp(Biplot$Maxima)
  }
  else
  {
    Minima=Biplot$Minima
    Maxima=Biplot$Maxima
  }
  
  for (j in 1:p) {
    switch(Biplot$Initial_Transformation, `Raw Data` = {
      if (ValuesScale == "Original"){
        Labels[[j]] = cbreaks(c(Minima[j], Maxima[j]), pretty_breaks(nticks))$breaks
        if (Biplot$Type=="LogFreqBiplot")
        {Ticks[[j]] = log(Labels[[j]]) }
        else
        { Ticks[[j]] = Labels[[j]]}
      }else{
        Labels[[j]] = cbreaks(c(Minima[j], Maxima[j]), pretty_breaks(nticks))$breaks
        Ticks[[j]] = Labels[[j]]
      }
      
    }, `Substract the global mean` = {
      if (ValuesScale == "Original"){
        Labels[[j]] = cbreaks(c(Minima[j], Maxima[j]), pretty_breaks(nticks))$breaks
        if (Biplot$Type=="LogFreqBiplot")
        {Ticks[[j]] = log(Labels[[j]])-Biplot$Gmean }
        else
        { Ticks[[j]] = Labels[[j]]-Biplot$Gmean}}
      else{
        Labels[[j]] = cbreaks(c(Minima[j], Maxima[j]), pretty_breaks(nticks))$breaks
        Ticks[[j]] = Labels[[j]]-Biplot$Gmean
      }
      
    }, `Double centering` = {
      
      Minima=min(Biplot$Scaled_Data[,j])
      Maxima=max(Biplot$Scaled_Data[,j])
      Labels[[j]] = cbreaks(c(Minima, Maxima), pretty_breaks(nticks))$breaks
      Ticks[[j]] = Labels[[j]]
      
    }, `Column centering` = {
      
      switch(TypeScale, Complete = {
        if (ValuesScale == "Original") {
          Labels[[j]] = cbreaks(c(Biplot$Minima[j], Biplot$Maxima[j]), pretty_breaks(nticks))$breaks
          Ticks[[j]] = (Labels[[j]] - Biplot$Means[j])
        } else {
          Ticks[[j]] = Labels[[j]] = cbreaks(c(Biplot$Minima[j] - Biplot$Means[j], Biplot$Maxima[j] - Biplot$Means[j]), pretty_breaks(nticks))$breaks
          Labels[[j]] = round(Ticks[[j]], digits = 2)
        }
      }, StdDev = {
        if (ValuesScale == "Original") {
          Labels[[j]] = round(c(Biplot$Means[j] - 3 * Biplot$Deviations[j], Biplot$Means[j] - 2 * Biplot$Deviations[j], Biplot$Means[j] - Biplot$Deviations[j], 
                                Biplot$Means[j] + Biplot$Deviations[j], Biplot$Means[j] + 2 * Biplot$Deviations[j], Biplot$Means[j] + 3 * Biplot$Deviations[j]), digits = 2)
          Ticks[[j]] = Labels[[j]] + Biplot$Means[j]
        } else {
          Ticks[[j]] = round(c(-3 * Biplot$Deviations[j], -2 * Biplot$Deviations[j], -1 * Biplot$Deviations[j], Biplot$Deviations[j], 2 * Biplot$Deviations[j], 
                               3 * Biplot$Deviations[j]), digits = 2)
          Labels[[j]] = Ticks[[j]]
        }
      }, BoxPlot = {
        if (ValuesScale == "Original") {
          Labels[[j]] = c(Biplot$Minima[j], Biplot$P25[j], Biplot$Median[j], Biplot$P75[j], Biplot$Maxima[j])
          Ticks[[j]] = Labels[[j]] - Biplot$Means[j]
        } else {
          Ticks[[j]] = c(Biplot$Minima[j], Biplot$P25[j], Biplot$Median[j], Biplot$P75[j], Biplot$Maxima[j])
          Ticks[[j]] = round((Ticks[[j]] - Biplot$Means[j]), digits = 2)
          Labels[[j]] = Ticks[[j]]
        }
      })
      
    }, `Standardize columns` = {
      
      switch(TypeScale, Complete = {
        if (ValuesScale == "Original") {
          Labels[[j]] = cbreaks(c(Biplot$Minima[j], Biplot$Maxima[j]), pretty_breaks(nticks))$breaks
          Ticks[[j]] = (Labels[[j]] - Biplot$Means[j])/Biplot$Deviations[j]
        } else {
          Ticks[[j]] = c(-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3)
          Labels[[j]] = round(Ticks[[j]], digits = 2)
        }
      }, StdDev = {
        if (ValuesScale == "Original") {
          Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
          Labels[[j]] = round(Ticks[[j]] * Biplot$Deviations[j] + Biplot$Means[j], digits = 1)
        } else {
          Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
          Labels[[j]] = round(Ticks[[j]], digits = 2)
        }
      }, BoxPlot = {
        if (ValuesScale == "Original") {
          Labels[[j]] = c(Biplot$Minima[j], Biplot$P25[j], Biplot$Median[j], Biplot$P75[j], Biplot$Maxima[j])
          Ticks[[j]] = (Labels[[j]] - Biplot$Means[j])/Biplot$Deviations[j]
        } else {
          Ticks[[j]] = c(Biplot$Minima[j], Biplot$P25[j], Biplot$Median[j], Biplot$P75[j], Biplot$Maxima[j])
          Ticks[[j]] = round((Ticks[[j]] - Biplot$Means[j])/Biplot$Deviations[j], digits = 2)
          Labels[[j]] = Ticks[[j]]
        }
      })
      
    }, `Divide by the column means and center` = {
      Minima=min(Biplot$Scaled_Data[,j])
      Maxima=max(Biplot$Scaled_Data[,j])
      Labels[[j]] = cbreaks(c(Minima, Maxima), pretty_breaks(nticks))$breaks
      Ticks[[j]] = Labels[[j]]
      
    }, `Normalized residuals from independence` = {
      switch(TypeScale, Complete = {
        if (ValuesScale == "Original") {
          Ticks[[j]] = ""
        } else {
          Labels[[j]] = ""
        }
      }, StdDev = {
        if (ValuesScale == "Original") {
          Ticks[[j]] = ""
        } else {
          Labels[[j]] = ""
        }
      }, BoxPlot = {
        if (ValuesScale == "Original") {
          Ticks[[j]] = ""
        } else {
          Labels[[j]] = ""
        }
        
      })
    })
    # 	if (Biplot$Type=="LogFreqBiplot") Labels[[j]]=exp(Labels[[j]])
  }
  
  return(list(Ticks = Ticks, Labels = Labels))
}


