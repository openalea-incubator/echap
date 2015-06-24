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
  dates <- sapply(xydb[[g]]$Source, function(s) format(as.Date(strsplit(s,split='_')[[1]][3], '%d%m%y'),'%d/%m/%Y'))
  TT <- TTlin[[g]]$TT[match(dates,TTlin[[g]]$Date)]
  xydb[[g]]$HSest <- slopeM[[g]] * (TT - TToM[[g]])
}

  #to do :add missing columns to homogenise data and bind
