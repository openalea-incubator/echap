#
#
#           Compilations donnees formes des feuilles (curvature, leaf width profil)
#
# scan compilation
# with Nflig added
#
write.csv(scanleafdb, 'Compil_scan.csv', row.names=FALSE)
#
# Courbures
#
# homogenise data, split curvdb into leaf meta-info and xy interpolated data
#
leafcols <- c('label','Source', 'N', 'var', 'rank','relative_ranktop')
ptcols <- c(leafcols, 'XY', paste('Pt',1:18,sep=''))
df <- do.call('rbind', sapply(names(curvdb), function(g) {
  curv <- curvdb[[g]]
  curv$label = g
  if (!'relative_ranktop'%in%colnames(curv))
    curv[['relative_ranktop']] <- as.numeric(NA)
  curv[,ptcols]}, simplify=FALSE))
# 
plants <- split(df,list(df$Source,df$N),drop=TRUE)
iplants <- seq(plants)
#
leaves <- do.call('rbind',mapply(function(pl,iplant) {
  mat=pl[pl$rank > 0 & pl$XY > 0,leafcols]
  mat$iplant=iplant
  mat},plants,iplants,SIMPLIFY=FALSE))
xyl <- do.call('rbind',mapply(curv2xy,plants,iplants,SIMPLIFY=FALSE))
# add inerv
leaves$inerv <- seq(nrow(leaves))
#
# add Plant notations (nff, Nflig) and HS,SSI from phend
# add Lb, Ab, Wb from scan data
# estimate ranktop + missing relative_ranktop (+ keep track of TT,HSest,nffest)
#
leaves <- add_dimphen(leaves, phend, scanleafdb, TTlin, TTem, TToM)
#
# estimate age / HS
# 
leaves <- do.call('rbind',lapply(split(leaves,leaves$Source), function(x) {
  x$dHS <- x$HS -  mean(x$HS,na.rm=TRUE)
  x$dSSI <-  x$SSI - mean(x$SSI,na.rm=TRUE)
  x}))
leaves$dHS <- ifelse(is.na(leaves$dHS),0,leaves$dHS)
leaves$dSSI <- ifelse(is.na(leaves$dSSI),0,leaves$dSSI)
leaves$dagep <- ifelse(abs(leaves$dHS) > 0,leaves$dHS, leaves$dSSI) 
leaves$age <- leaves$HSest + leaves$dagep - leaves$rank + 1 # TO DO use HS/SSI instead
# use age to re-estimate HS
leaves$HS <- leaves$age + leaves$rank + 1
#
# Export xydb for echap/adel pipelines
#
cols <- c('variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','Nflig','HS','inerv','x','y','hins','side')
#
xyTremie <- merge(xyl,leaves)[,c('variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','Nflig','HS','inerv','xr','yr','hins','side')]
xyTremie$x <- xyTremie$xr
xyTremie$y <- xyTremie$yr
write.csv(xyTremie[,cols],'xydb_Boigneville_Tremie12_Tremie13.csv', row.names=FALSE)
#
#
xyMercia <- merge(xyl_2010[,c('inerv','xr','yr','hins','side')],leaves_2010)
xyMercia$x <- xyMercia$xr
xyMercia$y <- xyMercia$yr
xyMercia$variety <- xyMercia$label
xyMercia$variety_code <- as.numeric(as.factor(xyMercia$variety))
xyMercia$harvest <- as.numeric(as.factor(xyMercia$Source))
write.csv(xyMercia[,cols],'xydb_Grignon2010.csv', row.names=FALSE)
#
# Export data for analysis wth tino/corine data pipeline/ median leaf shape
#
matphi <- xyl[,c('iplant','rank','phiP')]
matphi <- matphi[!duplicated(matphi),]
bl_e <- merge(leaves,matphi)[, c('inerv','label','Source','plant','rank','ranktop','age','lmax','wmax','A','Agr','mass','stat','phiP')]
write.csv(bl_e, 'leafshape_metainfo_Boigneville_2012_2013.csv',row.names=FALSE)
#
mati <- leaves[,c('iplant','rank','inerv')]
nervs_e <- merge(xyl,mati)[,c('inerv','phiS','s','x','y','xprot','yprot','xrot','yrot','xr','yr','hins','side')]
write.csv(nervs_e, 'leafshape_xydata_Boigneville_2012_2013.csv',row.names=FALSE)
#
# Process/export imageJ data
#
ijplants <- split(ij_curvdb,list(ij_curvdb$variety,ij_curvdb$Source,ij_curvdb$N),drop=TRUE)
iplants <- seq(ijplants)
#
#check
par(mfrow=c(5,5),mar=c(4,2,1,1))
mapply(function(pl,lab) plot(pl$x,pl$y,col=pl$serie+1,pch=16,xlab=lab),ijplants[1:25], names(ijplants)[1:25])
mapply(function(pl,lab) plot(pl$x,pl$y,col=pl$serie+1,pch=16,xlab=lab),ijplants[26:50], names(ijplants)[26:50])
mapply(function(pl,lab) plot(pl$x,pl$y,col=pl$serie+1,pch=16,xlab=lab),ijplants[51:63], names(ijplants)[51:63])
#
# split leaf info / xy data
# compute dHS for growing plants
leafcols <- c('variety', 'Source', 'N', 'serie', 'ntop_lig', 'green_leaves')
ijleaves <- do.call('rbind',mapply(function(pl,iplant) {
  mat=pl[pl$serie > 0, leafcols]
  mat <- mat[!duplicated(mat$serie),]
  mat$iplant <- iplant
  mat$dHS <- 0
  if (min(mat$ntop_lig) <= 0) {
    lig <- pl[pl$ntop_lig == 1,]
    Llig <- sum(sqrt(diff(lig$x)^2 + diff(lig$y)^2))
    last <- pl[pl$ntop_lig == 0,]
    Llast <- sum(sqrt(diff(last$x)^2 + diff(last$y)^2))
    mat$dHS <- Llast / Llig
  }
  mat},ijplants,iplants,SIMPLIFY=FALSE))
# add inerv
ijleaves$inerv <- seq(nrow(ijleaves))
# extract xy
ijxyl <- do.call('rbind',mapply(ijcurv2xy,ijplants,iplants,SIMPLIFY=FALSE))
# add inerv to xy and phi_P to leaves
mati <- ijleaves[,c('iplant','serie','inerv')]
ijxyl <- merge(ijxyl,mati)
matphi <- ijxyl[,c('iplant','serie','phiP')]
matphi <- matphi[!duplicated(matphi),]
ijleaves <- merge(ijleaves,matphi)
#
# compute HS,age, assume rank, ranktop, N2
#
#add HSest, nffest (genotypic estimates) and leaf meta info columns
ijleaves <- ijphen(ijleaves, TTlin, TTem, TToM)
#
# assume N2
ijleaves$N2 <- ifelse(ijleaves$HSest < ijleaves$nffest | ijleaves$dHS != 0, floor(ijleaves$HSest), round(ijleaves$nffest))
# rank
ijleaves$rank <- ijleaves$N2 + 1 - ijleaves$ntop_lig
ijleaves$ranktop <- ijleaves$nffest - ijleaves$rank + 1
#HS
# use opposite of delta green leaves as a proxy for dHS
dGL <- unsplit(lapply(split(ijleaves, ijleaves$variety),function(x) {
  moy <- mean(x$green_leaves,na.rm=TRUE)
  ifelse(is.na(x$green_leaves),0, x$green_leaves-moy)}),ijleaves$variety)
dHSm <- ijleaves$HSest - sapply(ijleaves$HSest,floor)
ijleaves$HS <- ijleaves$N2 + ifelse(ijleaves$dHS <=0, dHSm, ijleaves$dHS) - dGL
#age
ijleaves$age <- ijleaves$HS - ijleaves$rank
# Export data for analysis wth tino/corine data pipeline/ median leaf shape
#
bl_ij <- ijleaves[, c('inerv','label','Source','plant','rank','ranktop','age','lmax','wmax','A','Agr','mass','stat','phiP')]
write.csv(bl_ij, 'leafshape_metainfo_Boigneville_2011.csv',row.names=FALSE)
#
nervs_ij <- ijxyl[,c('inerv','phiS','s','x','y','xprot','yprot','xrot','yrot','xr','yr','hins','side')]
write.csv(nervs_ij, 'leafshape_xydata_Boigneville_2011.csv',row.names=FALSE)
