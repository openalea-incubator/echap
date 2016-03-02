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
# Export for analysis wth tino/corine data
#
cols <- c('label','Source', 'N', 'var', 'nff', 'Nflig', 'HS','nffest', 'HSest','rank','ranktop', 'relative_ranktop', 'stat','lmax','wmax','A_bl', 'A_bl_green', 'XY', paste('Pt',1:18,sep=''))
write.csv(do.call('rbind', sapply(names(xydb), function(g) {dim <- xydb[[g]]; dim$label = g; dim[,cols]}, simplify=FALSE)), 'Compil_curvature.csv',row.names=FALSE)
