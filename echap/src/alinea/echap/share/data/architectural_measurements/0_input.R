#
#                      Reading input files
#
source('preprocessing_tools.R')
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
# silhouettes
#
# For Tremie 13 on 29/04/2013, ranks were extimated a posteriori.
# Original data are given for relative rank top only.
# Estimation was done manually knowing that after HSfit, plants of tremie 13 were at stage 9.4 for nff=11, and stage 9.8 for nff=12, ie, las ligulated leaf was always 9.
#Leaf 9 was guessed after manual inspection of images (Christian, February 2016)
#
curvature <- list(Tremie12 = paste('sampled_plants', c('090312','110412', '090512', '120612'),sep='_'),
                  Tremie13 = paste('sampled_plants', c('290413'),sep='_'))
#
curvdb <- sapply(names(curvature), function(g) readCurv(prefix[g], curvature[[g]]), simplify=FALSE)
#
# Mercia Rht3 data acquired with ImageJ
#
ij_curvature <- list(Mercia = paste('sampled_plants', c('270411','010611'),sep='_'),
                  Rht3 = paste('sampled_plants', c('270411','010611'),sep='_'))
ij_curvdb <- do.call('rbind',sapply(names(ij_curvature), function(g) readCurv_ij(g, ij_curvature[[g]]), simplify=FALSE))
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
#
# Plant notations (subset of notations)
#
plantdb <- lapply(notdb, function(x) plant_notations(x))
#
# leaf counts notation (subset of plant notation, completed for sampling with no notation using inspection of images)
#
leafcdb <- lapply(plantdb, function(x) x[,c('Source','N', 'nff', 'Nflig','Nfvis','stem_half_area')])
# for scan data Tremie 12 of sampled plant april 2, nflig = last measured leaf, hypothethise one growing leaf
data <- scandb$Tremie12[scandb$Tremie12$Source=='sampled_plants_020412',]
nfl <- do.call('rbind', lapply(split(data, data$N), function(x) data.frame(Source=x$Source[1],N=x$N[1],nff=NA, Nflig=max(x$rank), Nfvis=1, stem_half_area=NA)))
leafcdb$Tremie12 <- rbind(leafcdb$Tremie12, nfl)
# for curvature Tremie 13 on 29 april, nflig=9 by convention (see input.R), hypothethise one growing leaf
data <- curvdb$Tremie13[curvdb$Tremie13$Source=='sampled_plants_290413',]
nfl <- do.call('rbind', lapply(split(data, data$N), function(x) data.frame(Source=x$Source[1],N=x$N[1],nff=NA, Nflig=9, Nfvis=1, stem_half_area=NA)))
leafcdb$Tremie13 <- rbind(leafcdb$Tremie13, nfl)
#
# Extract 'other than leaf width profile' data from scan, homogenise and mix with leaf counts
#
scanleafdb <- do.call('rbind',mapply(function(x,name) {
  cols <- c('Source','prelevement', 'rep','N','id_Axe', 'rank','stat','lmax','wmax','A_bl', 'A_bl_green')
  x <- x[,cols[cols%in%colnames(x)]]
  if (!'lmax'%in% colnames(x))
    x['lmax'] <- NA
  if (!'wmax'%in% colnames(x))
    x['wmax'] <- NA
  x['variety'] = name
  nfl <- leafcdb[[name]]
  nfl$id_Axe='MB'
  x <- merge(x,nfl, all.x=TRUE)
  x
  }, scandb, names(scandb),SIMPLIFY=FALSE))
#
# import curvature / leaf area data from Corinne/Tino experiment 2010
#
leaves_2010 <-  read.csv('leafshape_metainfo_Grignon2010.csv', sep=',', dec='.')
xyl_2010 <-  read.csv('leafshape_xydata_Grignon2010.csv', sep=',', dec='.')
