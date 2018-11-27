Save.Installed.Packages <- function(){
  nombres=rownames(installed.packages())
  n=length(nombres)
  packs=paste("c('", nombres[1], "'", sep="")
  for (i in 2:n)
    packs=paste(packs, nombres[i], sep="', '")
  packs=paste(packs, "')", sep="")
  file=paste("InstalledPacks",paste(version$major, version$minor), sep="")
  save(packs, file=file)
}


