#
#
#           Compilations donnees formes des feuilles (curvature, leaf width profil)
#
# scan compilation
# with Nflig added
#
write.csv(scanleafdb, 'Compil_scan.csv', row.names=FALSE)

# Courbures
#
# add Plant notations (nff, Nflig) and HS,SSI from phend
# add Lb, Ab, Wb from scan data
# estimate HSest
xydb <- curvdb
scand <- scanleafdb[scanleafdb$id_Axe=='MB',]
scand <- split(scand,scand$variety)
#
for (g in names(xydb)) {
  xydb[[g]] <- merge(xydb[[g]],phend[[g]][,c('Source', 'N', 'nff','Nflig','HS','SSI')], all.x=TRUE)
  xydb[[g]] <- merge(xydb[[g]],scand[[g]][,c('Source','N','rank','stat','lmax','wmax','A_bl','A_bl_green')], all.x=TRUE)
  if (!'relative_ranktop'%in%colnames(xydb[[g]]))
    xydb[[g]][['relative_ranktop']] <- as.numeric(NA)
  #
  dates <- sapply(xydb[[g]]$Source, function(s) format(as.Date(strsplit(s,split='_')[[1]][3], '%d%m%y'),'%d/%m/%Y'))
  TT <- TTlin[[g]]$TT[match(dates,TTlin[[g]]$Date)]
  xydb[[g]]$TT <- TT
  TTflag <- aggregate(TTem[[g]], list(TTem[[g]]$nff),max)
  dTTflag <- diff(TTflag$TTlig) /  diff(TTflag$nff)
  nffM <- mean(sapply(split(TTem[[g]],TTem[[g]]$N),function(x) x$nff[1]))
  xydb[[g]]$HSest <- slopeM[[g]] * (TT + ifelse(is.na(xydb[[g]]$nff), 0, dTTflag * (xydb[[g]]$nff - nffM)) - TToM[[g]])
  xydb[[g]]$nffest <- ifelse(is.na(xydb[[g]]$nff),nffM, xydb[[g]]$nff)
  # ranktop estimation
  xydb[[g]]$ranktop <-  xydb[[g]]$nffest -  xydb[[g]]$rank + 1
  # relative_ranktop
  hs <- ifelse(!is.na(xydb[[g]]$HS), xydb[[g]]$HS, pmin(xydb[[g]]$nffest, xydb[[g]]$HSest))
  xydb[[g]]$relative_ranktop <- ifelse(!is.na(xydb[[g]]$relative_ranktop), xydb[[g]]$relative_ranktop, hs - xydb[[g]]$rank + 1)
  # restore stem
  xydb[[g]][xydb[[g]]$Organ == 0, c('ranktop', 'relative_ranktop')] <- 0
}
#
# Create xydb format
#
leafcols <- c('label','Source', 'N', 'var', 'nff', 'Nflig', 'HS','nffest', 'HSest','rank','ranktop', 'relative_ranktop', 'stat','lmax','wmax','A_bl', 'A_bl_green')
ptcols <- c(leafcols, 'XY', paste('Pt',1:18,sep=''))
xydb <- do.call('rbind', sapply(names(xydb), function(g) {dim <- xydb[[g]]; dim$label = g; dim[,ptcols]}, simplify=FALSE))
# add inerv and separate leaf_info from points
plants <- split(xydb,list(xydb$Source,xydb$N),drop=TRUE)
iplants <- seq(plants)
leaves <- do.call('rbind',mapply(function(pl,iplant) {mat=pl[pl$rank > 0 & pl$XY > 0,leafcols];mat$iplant=iplant;mat},plants,iplants,SIMPLIFY=FALSE))
xyl <- do.call('rbind',mapply(curv2xy,plants,iplants,SIMPLIFY=FALSE))
# add cols
leaves$inerv <- seq(nrow(leaves))
leaves$age <- leaves$HSest - leaves$rank + 1 # TO DO use HS/SSI instead
#
leaves$variety <- leaves$label
leaves$variety_code <- as.numeric(as.factor(leaves$variety))
leaves$harvest <- as.numeric(as.factor(leaves$Source))
leaves$plant <- leaves$N
leaves$A <- leaves$A_bl
leaves$Agr <- leaves$A_bl_green
leaves$nmax <- leaves$nffest
leaves$N2 <- leaves$Nflig
leaves$mass <- NA
#
xyl$phiS <- pi / 2
xyl$xprot <- xyl$x
xyl$yprot <- xyl$y
xyl$xrot <- xyl$x
xyl$yrot <- xyl$y
xyl$xr <- xyl$x
xyl$yr <- xyl$y
#
# Export xydb
#
xyTremie <- merge(xyl,leaves)[,c('variety','variety_code','harvest','plant','rank','ranktop','relative_ranktop','HS','inerv','x','y','hins','side')]
write.csv(xyTremie,'xydb_Boigneville_Tremie12_Tremie13.csv', row.names=FALSE)
#
# Export data for analysis wth tino/corine data / median shape
#
matphi <- xyl[,c('iplant','rank','phiP')]
matphi <- matphi[!duplicated(matphi),]
blade <- merge(leaves,matphi)[, c('inerv','label','Source','plant','rank','ranktop','age','lmax','wmax','A','Agr','mass','stat','phiP')]
write.csv(blade, 'export_blade_echap.csv',row.names=FALSE)
#
mati <- leaves[,c('iplant','rank','inerv')]
nervs <- merge(xyl,mati)[,c('inerv','phiS','s','x','y','xprot','yprot','xrot','yrot','xr','yr')]
write.csv(nervs, 'export_nervs_echap.csv',row.names=FALSE)
