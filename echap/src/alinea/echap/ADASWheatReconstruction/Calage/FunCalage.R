#
#
#       Fonctions generiques utiles pour le calage
#
#
#
#Calcule un intervalle de confiance, NA si pas de donnees ou pas de variance
#
ConfInt <- function(x) {
  if (length(na.omit(x))> 1 & length(unique(na.omit(x))) > 1){
    t <- t.test(x)$conf.int
    t[2]-t[1]
  }
  else NA
}
#
#plot point + intervalle de confiance
#
PointIc <- function(x,y,ci,col=1,...) {
  points(x,y,col=col,...)
  segments(x,y-ci/2,x,y+ci/2,col=col)
}
#
#fonction generique renvoyant un data frame avec resultats de by
#bycols = choix des colones (par nom) de dat qu'il faut grouper
#FUN doit renvoyer une matrice ou un dataframe avec des noms de colones ou un vecteur avec des noms
#FUN recoit en entree la soous partie mat  du tableau et ...
#
tableby <- function(dat,bycols,FUN,...) {
  bylist <- sapply(bycols,function(col) dat[,col],simplify=F)
  data.frame(do.call("rbind",by(dat,bylist,function(mat) {
    res=FUN(mat,...)
    if (!is.null(dim(res))) 
      noms <- colnames(res)
    else {
      noms <- names(res)
      dim(res) <- c(1,length(res))
    }
    bymat <- mat[seq(nrow(res)),bycols]
    if (nrow(mat) < nrow(res)) for (i in 2:nrow(res))
      bymat[i,] <- bymat[1,]
    res <- cbind(bymat,res)
    colnames(res) <- c(bycols,noms)
    rownames(res) <- seq(nrow(res))
    res
    })))
}
#
#lecturefichier excell
#
readExcell <- function(file,onglet="Feuil1",asis=FALSE) {
res <- NULL
if (require(RODBC)) {
ch <- odbcConnectExcel(file)
if (asis)
  res <- sqlQuery(ch,paste("select * from [",onglet,"]",sep=""))
else
  res <- sqlQuery(ch,paste("select * from [",onglet,"$]",sep=""))
odbcClose(ch)}
res
}
#
#conversion des dates excell en num de jour dans un data frame
#
formatJ <- function(dat) {
  coldates <- grep("date",colnames(dat))
  if (length(coldates) > 1)
    dat[,coldates] <- apply(format(dat[,coldates],"%j"),2,as.numeric)
  else if (length(coldates) > 0)
    dat[,coldates] <- as.numeric(format(dat[,coldates],"%j"))
  dat
}
