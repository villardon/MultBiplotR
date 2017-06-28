UpdateMultBiplot <- function(){
  update.packages(c("scales", "geometry", "deldir", "rgl", "mirt", "GPARotation", "MASS", "kde2d", "lattice", "splom", "dae"))
  install.packages("http://biplot.usal.es/classicalbiplot/multbiplot-in-r/multbiplotr.gz", repos = NULL, type="source")
}
