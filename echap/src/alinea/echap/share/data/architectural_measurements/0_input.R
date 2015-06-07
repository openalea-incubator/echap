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
# dynamic notation tagged plants
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
#          - 2 scans Tremie 13  (Benjamin)
#
# scans
#
scanned <- list(Tremie12 = paste('scanned_plants', c('090312','020412','110412', '090512'),sep='_'),
                Tremie13 = paste('scanned_plants', c('220413','030513'), sep='_'))
#
scandb <- sapply(names(scanned), function(g) readScanned(prefix[g], scanned[[g]]),simplify=FALSE)
#
# notations (sheath, internode, diameter, stage)
#
notations <- list(Mercia = c('tagged_plants_090611'),
                  Rht3 = c('tagged_plants_090611'),
                  Tremie12 = c(scanned[['Tremie12']][-2],
                               'silhouette_plants_120612',
                               'tagged_plants_120712'),
                  Tremie13 = scanned[['Tremie13']])
#
notdb <- sapply(names(notations), function(g) readNotations(prefix[g], notations[[g]], varname[g]), simplify=FALSE)
#
# notations ssi Tremie13 02/04/2013
dat <- readTagged('Tremie13_ssi_020413.txt', TTlin$Tremie13, 'Tremie')
notdb$Tremie13$ssi_sample_020413 <- dat[,-grep('^Lg', colnames(dat))]
