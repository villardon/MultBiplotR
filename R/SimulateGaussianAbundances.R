SimulateGaussianAbundances <- function(nsites = 100, nspecies = 15, dimens = 2, maximum = 10, tol = 1) {
  
  SiteNames = "S1"
  for (i in 2:nsites) SiteNames = c(SiteNames, paste("S", i, sep = ""))
  SpeciesNames = "Sp1"
  for (i in 2:nspecies) SpeciesNames = c(SpeciesNames, paste("Sp", i, sep = ""))
  
  X = matrix(rnorm(nsites * dimens), nsites, dimens)
  U = matrix(rnorm(nspecies * dimens), nspecies, dimens)
  A = rbind(X, U)
  plot(A[, 1], A[, 2], cex = 0, asp = 1)
  points(X[, 1], X[, 2], cex = 0.2)
  text(X[, 1], X[, 2], SiteNames)
  points(U[, 1], U[, 2], col = "red", pch = 16, , cex = 0.2)
  text(U[, 1], U[, 2], SpeciesNames, col = "red")
  
  for (i in 1:nspecies) Circle(tol, origin = c(U[i, 1], U[i, 2]), color = "red")
  AB = round(maximum * exp(-0.5 * (1/tol) * DistUnfold(X, U)^2), digits = 1)
  
  rownames(AB) = SiteNames
  colnames(AB) = SpeciesNames
  dev.new()
  unf = Unfolding(AB, model = "Ratio", condition=2, weight=2)
  plot(unf, PlotTol = T)
  Pr = SimpleProcrustes(X, unf$X)
  dev.new()
  plot(Pr)
  return(AB)
}
