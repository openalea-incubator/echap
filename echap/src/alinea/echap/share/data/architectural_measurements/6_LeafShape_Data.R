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
