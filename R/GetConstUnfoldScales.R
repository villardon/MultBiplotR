GetConstUnfoldScales <- function(Unfold, nticks = 4, TypeScale = "Complete", ValuesScale = "Original") {
  p = dim(Unfold$Environment)[2]
  Ticks = list()
  Labels = list()
  for (j in 1:p) {
    switch(TypeScale, Complete = {
      if (ValuesScale == "Original") {
        Labels[[j]] = cbreaks(c(Unfold$Minima_Z[j], Unfold$Maxima_Z[j]), pretty_breaks(nticks))$breaks
        Ticks[[j]] = (Labels[[j]] - Unfold$Means_Z[j])/Unfold$Deviations_Z[j]
      } else {
        Ticks[[j]] = c(-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3)
        Labels[[j]] = round(Ticks[[j]], digits = 2)
      }
    }, StdDev = {
      if (ValuesScale == "Original") {
        Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
        Labels[[j]] = round(Ticks[[j]] * Unfold$Deviations_Z[j] + Unfold$Means_Z[j], digits = 1)
      } else {
        Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
        Labels[[j]] = round(Ticks[[j]], digits = 2)
      }
    }, BoxPlot = {
      if (ValuesScale == "Original") {
        Labels[[j]] = c(Unfold$Minima_Z[j], Unfold$P25_Z[j], Unfold$Median_Z[j], Unfold$P75_Z[j], Unfold$Maxima_Z[j])
        Ticks[[j]] = (Labels[[j]] - Unfold$Means[j])/Unfold$Deviations[j]
      } else {
        Ticks[[j]] = c(Unfold$Minima[j], Unfold$P25[j], Unfold$Median[j], Unfold$P75[j], Unfold$Maxima[j])
        Ticks[[j]] = round((Ticks[[j]] - Unfold$Means[j])/Unfold$Deviations[j], digits = 2)
        Labels[[j]] = Ticks[[j]]
      }
    })
  }
  return(list(Ticks = Ticks, Labels = Labels))
}