shepard.Unfolding <- function(unf) {
  plot(unf$Distances, unf$Disparities)
  plot(unf$TransformedAbundances, unf$Distances)
  plot(unf$Abundances, unf$Distances)
}