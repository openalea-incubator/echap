#
#                      Reading input files
#
#
genos <- c('Mercia', 'Rht3', 'Tremie12', 'Tremie13')
#
prefix <- c('MerciaRht3', 'MerciaRht3','Tremie12', 'Tremie13')
names(prefix) <- genos
#
varname <- c('Mercia','Rht3','Tremie','Tremie')
names(varname) <- genos
#
# TT linear since sowing
#
TTlin <- sapply(genos, function(g) read.csv2(paste(prefix[g], '_TTlin_sowing.csv', sep='')), simplify=FALSE)
#
# dynamic notation tagged plants treated archi
#
tagged <- sapply(genos, function(g) readTagged(paste(prefix[g],'_suivis_plantes_baguees.txt', sep=''), TTlin[[g]], varname[g]), simplify=FALSE)
#
# dynamic notation (nec/symptom) tagged plant treated symptom
#
symtagged <- sapply(c('Tremie12', 'Tremie13'), function(g) readSymTagged(paste('../',g,'_treated_symptom_tagged.csv',sep=''), TTlin[[g]]), simplify=FALSE)
#
# destructive samplings
#
# misses : - scans Mercia/Rht3 09/06/2011 (tagged plants): images non analysees
#          - silhouettes Mercia/Rht3
#          - 2 scans Tremie 13  (Benjamin)
#
# scans
#
scanned <- list(Tremie12 = paste('sampled_plants', c('090312','020412','110412', '090512'),sep='_'),
                Tremie13 = paste('sampled_plants', c('220413','030513'), sep='_'))
#
scandb <- sapply(names(scanned), function(g) readScanned(prefix[g], scanned[[g]]),simplify=FALSE)
#
# scan compilation
#
comp <- do.call('rbind',mapply(function(x,name) {
  cols <- c('prelevement', 'N','id_Axe', 'rank','stat','lmax','wmax','A_bl', 'A_bl_green')
  x <- x[,cols[cols%in%colnames(x)]]
  if (!'lmax'%in% colnames(x))
    x['lmax'] <- NA
  if (!'wmax'%in% colnames(x))
    x['wmax'] <- NA
  x['variety'] = name
  x
  }, scandb, names(scandb),SIMPLIFY=FALSE))
#
write.csv(comp, 'Compil_scan.csv', row.names=FALSE)
# silhouettes
#
curvature <- list(Tremie12 = paste('sampled_plants', c('090312','110412', '090512', '120612'),sep='_'),
                  Tremie13 = paste('sampled_plants', c('290413'),sep='_'))
#
curvdb <- sapply(names(curvature), function(g) readCurv(prefix[g], curvature[[g]]), simplify=FALSE)
#
# notations (sheath, internode, diameter, stage)
#
notations <- list(Mercia = c('tagged_plants_010611','tagged_plants_090611'),
                  Rht3 = c('tagged_plants_010611','tagged_plants_090611'),
                  Tremie12 = c(scanned[['Tremie12']][-2],
                               'sampled_plants_120612',
                               'tagged_plants_120712'),
                  Tremie13 = scanned[['Tremie13']])
#
notdb <- sapply(names(notations), function(g) readNotations(prefix[g], notations[[g]], varname[g]), simplify=FALSE)
# special reader for notations ssi Tremie13 02/04/2013
dat <- readTagged('Tremie13_notations_sampled_plants_020413.txt', TTlin$Tremie13, 'Tremie')
notdb$Tremie13$sampled_plants_020413 <- dat[,-grep('^Lg', colnames(dat))]
