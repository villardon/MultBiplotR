UpdateMultBiplot <- function(dependencies=TRUE, upgrade_dependencies = TRUE){
  if (!require(devtools))   
    install.packages("devtools")
  library(devtools)
  install_github("villardon/MultBiplotR", dependencies = dependencies, upgrade_dependencies=upgrade_dependencies )
}
