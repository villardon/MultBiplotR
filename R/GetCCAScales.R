GetCCAScales <- function(CCA, nticks = 7, TypeScale = "Complete", ValuesScale = "Original") {
	p = CCA$nvars
	Ticks = list()
	Labels = list()

	for (j in 1:p) {
		switch(TypeScale, Complete = {
			if (ValuesScale == "Original") {
				Labels[[j]] = cbreaks(c(CCA$Minima_Z[j], CCA$Maxima_Z[j]), pretty_breaks(nticks))$breaks
				Ticks[[j]] = (Labels[[j]] - CCA$Means_Z[j])/CCA$Deviations_Z[j]
			} else {
				Ticks[[j]] = c(-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3)
				Labels[[j]] = round(Ticks[[j]], digits = 2)
			}
		}, StdDev = {
			if (ValuesScale == "Original") {
				Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
				Labels[[j]] = round(Ticks[[j]] * CCA$Deviations_Z[j] + CCA$Means_Z[j], digits = 1)
			} else {
				Ticks[[j]] = c(-3, -2, -1, 1, 2, 3)
				Labels[[j]] = round(Ticks[[j]], digits = 2)
			}
		}, BoxPlot = {
			if (ValuesScale == "Original") {
				Labels[[j]] = c(CCA$Minima_Z[j], CCA$P25_Z[j], CCA$Median_Z[j], CCA$P75_Z[j], CCA$Maxima_Z[j])
				Ticks[[j]] = (Labels[[j]] - CCA$Means[j])/CCA$Deviations[j]
			} else {
				Ticks[[j]] = c(CCA$Minima[j], CCA$P25[j], CCA$Median[j], CCA$P75[j], CCA$Maxima[j])
				Ticks[[j]] = round((Ticks[[j]] - CCA$Means[j])/CCA$Deviations[j], digits = 2)
				Labels[[j]] = Ticks[[j]]
			}
		})
	}
return(list(Ticks=Ticks, Labels=Labels))
}
