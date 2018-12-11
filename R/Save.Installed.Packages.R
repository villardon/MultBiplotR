Save.Installed.Packages <- function(){
  nombres=rownames(installed.packages())
  n=length(nombres)
  packs=paste("c('", nombres[1], "'", sep="")
  for (i in 2:n)
    packs=paste(packs, nombres[i], sep="', '")
  packs=paste(packs, "')", sep="")
  file=paste("InstalledPacks_",paste(version$major, version$minor), ".rda", sep="")
  save(nombres, file=file)
  return(nombres)
}


