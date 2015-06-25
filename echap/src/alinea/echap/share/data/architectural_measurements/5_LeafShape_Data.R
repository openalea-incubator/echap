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
  cols <- c('Source', 'N', 'rank','stat','lmax','wmax','A_bl')
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
}
#
# Export
#
cols <- c('label','Source', 'N', 'var', 'nff', 'Nflig', 'HS', 'HSest','rank','stat','lmax','wmax','A_bl', paste('Pt',1:18,sep=''))
write.csv(do.call('rbind', sapply(names(xydb), function(g) {dim <- xydb[[g]]; dim$label = g; dim[,cols]}, simplify=FALSE)), 'Compil_curvature.csv',row.names=FALSE)
