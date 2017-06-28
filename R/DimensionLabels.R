DimensionLabels <- function(dimens, Root = "Dim") {
	DimNames = paste(Root, "1")
	for (i in 2:dimens) DimNames = c(DimNames, paste(Root, i))
	return(DimNames)
}