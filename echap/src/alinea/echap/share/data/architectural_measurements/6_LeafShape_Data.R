#
#
#           Compilations donnees formes des feuilles (curvature, leaf width profil)
#
#
# Courbures
#
xydb <- curvdb
# add Plant notations (nff, Nflig) from plant notations
# add Lb, Ab, Wb from scan data
# add HS from phend
# estimate HSest
for (g in names(xydb)) {
  xydb[[g]] <- merge(xydb[[g]],msdb[[g]][,c('Source', 'N', 'nff', 'Nflig')], all.x=TRUE)
  xydb[[g]] <- merge(xydb[[g]],phend[[g]][,c('Source', 'N', 'HS')], all.x=TRUE)
  cols <- c('Source', 'N', 'rank','relative_ranktop','stat','lmax','wmax','A_bl','A_bl_green')
  xydb[[g]] <- merge(xydb[[g]],scandim[[g]][,colnames(scandim[[g]])%in%cols], all.x=TRUE)
  for (c in cols)
    if (!c%in%colnames(xydb[[g]]))
      xydb[[g]][[c]] <- as.numeric(NA)
  dates <- sapply(xydb[[g]]$Source, function(s) format(as.Date(strsplit(s,split='_')[[1]][3], '%d%m%y'),'%d/%m/%Y'))
  TT <- TTlin[[g]]$TT[match(dates,TTlin[[g]]$Date)]
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
# Export
#
cols <- c('label','Source', 'N', 'var', 'nff', 'Nflig', 'HS', 'nffest', 'HSest','rank','ranktop', 'relative_ranktop', 'stat','lmax','wmax','A_bl', 'A_bl_green', 'XY', paste('Pt',1:18,sep=''))
write.csv(do.call('rbind', sapply(names(xydb), function(g) {dim <- xydb[[g]]; dim$label = g; dim[,cols]}, simplify=FALSE)), 'Compil_curvature.csv',row.names=FALSE)
