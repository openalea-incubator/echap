#
#
#              Generation of synthetic HS_Data for reconstruction, using leaf blade profiles per nff
#
source('1_HS_plantes_baguees.R')
source('2_Blade_length_profiles.R')
#
# tagged plants : regenerate HS data with leaf profiles
# -----------------------------------------------------
#
dimPMf <- lapply(LbM, function(Lb) do.call('rbind',lapply(names(Lb), function(nff) data.frame(nff=as.numeric(nff), rank = Lb[[nff]]$ranks, x=Lb[[nff]]$L))))
dimPMMf <- lapply(LbMM, function(Lb) {Lb$x=Lb$L;Lb$rank=Lb$ranks;Lb})
#
dimPfitf <- mapply(function(dim,name) do.call('rbind', lapply(split(dim, dim$N), function(mat) final_length(mat,dimPMf[[name]], dimPMMf[[name]]))), dimP, names(dimP), SIMPLIFY=FALSE)
# HS,SSI per plant
phendbf <- sapply(genos, function(g) pheno(tagged[[g]], dimPfitf[[g]]), simplify=FALSE)
#
# ajout estimation SSI terrain 27/04 Mercia/Rht3
#
SSIobsMercia <- data.frame(Date='27/04/2011', N=c(8,11), SSI=c(7.3, 8.2))
SSIobsRht3 <- data.frame(Date='27/04/2011', N=c(1,4,5,6,8,11), SSI=c(7.4, 6.3,6.2,7.1,7.4,6.2))
phendbf$Mercia$SSI[match(paste(SSIobsMercia$Date,SSIobsMercia$N),paste(as.character(phendbf$Mercia$Date), phendbf$Mercia$N))] <- SSIobsMercia$SSI
phendbf$Rht3$SSI[match(paste(SSIobsRht3$Date,SSIobsRht3$N),paste(as.character(phendbf$Rht3$Date), phendbf$Rht3$N))] <- SSIobsRht3$SSI
phendbf$Mercia$GL = phendbf$Mercia$HS - phendbf$Mercia$SSI
phendbf$Rht3$GL = phendbf$Rht3$HS - phendbf$Rht3$SSI
#check
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(phendbf, function(p) {
  plot(c(0,2500),c(0,13),type='n')
  lapply(split(p,p$N), function(x) points(x$TT,x$HS,col=x$nff,pch=16))
  lapply(split(p,p$N), function(x) points(x$TT,x$GL,col=x$nff,pch=16))
})
#
# Export
#
write.csv(do.call('rbind', sapply(genos, function(g) {phen <- phendbf[[g]]; phen$label = g; phen[,c('Date', 'label', 'Trt', 'Rep', 'N', 'Axis', 'Nflig', 'Nfvis', 'nff', 'TT', 'HS', 'SSI', 'GL')]}, simplify=FALSE)), 'Compil_Pheno_treated_archi_tagged.csv',row.names=FALSE)
#
# HS/SSI data from detructive samplings
# -------------------------------------
#
phend <- NULL
#scan samples Tremie 12
dat <- scandb$Tremie12[scandb$Tremie12$id_Axe=='MB',1:9]
not <- notdb$Tremie12[grep('scanned',names(notdb$Tremie12))]
nfl <- do.call('rbind',lapply(not, function(x) data.frame(prelevement=x$Date, plant=x$N,id_Axe=x$axe,Nflig=x$Nflig, Nfvis=x$Nfvis)))
nfl <- na.omit(nfl[nfl$id_Axe=='MB',])
dat <- merge(dat,nfl)
phen <- do.call('rbind',lapply(split(dat,list(dat$prelevement,dat$plant), drop=TRUE), function(x) pheno_scan(x, LbMM$Tremie12)))
phen <- merge(phen, TTlin$Tremie12)
phend$Tremie12 <- phen
#scan samples Tremie13
dat <- scandb$Tremie13[scandb$Tremie13$id_Axe=='MB',]
dat$Nfvis <- 1#force HS computing (all sampling occured before flag leaf)
dat$A_bl_green <- dat$A_bl * dat$pcent_green / 100
phen <- do.call('rbind',lapply(split(dat,list(dat$prelevement,dat$plant), drop=TRUE), function(x) pheno_scan(x, LbMM$Tremie13)))
phen <- merge(phen, TTlin$Tremie13)
phend$Tremie13 <- phen
# silhouette data Tremie12 12/06/12
dat <- notdb$Tremie12$silhouette_plants_120612[,1:9]
dat <- dat[dat$axe=='MB',]
phen <- data.frame(Date=dat$Date, N=dat$N, nff=dat$Nflig, HS=dat$Nflig, SSI=dat$Nflig - dat$Nfvert, GL=dat$Nfvert)
phen <- merge(phen, TTlin$Tremie12)
phend$Tremie12 <- rbind(phend$Tremie12, phen)
# ssi sample Tremie13 02/04/2013
dat <- notdb$Tremie13$ssi_sample_020413
phen <- pheno_ssi(dat)[,c('Date', 'N','nff','HS','SSI','GL','TT')]
phend$Tremie13 <- rbind(phend$Tremie13, phen)
#
#check
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(names(phendbf), function(g) {
  p <- phendbf[[g]]
  plot(c(0,2500),c(0,14),type='n')
  lapply(split(p,p$N), function(x) points(x$TT,x$HS,col=x$nff,pch=16))
  lapply(split(p,p$N), function(x) points(x$TT,x$SSI,col=x$nff,pch=16,cex=0.7))
  lapply(split(p,p$N), function(x) points(x$TT,x$GL,col=x$nff,pch=16))
  if (g %in% names(phend)){
    phen <- phend[[g]]
    coul <- ifelse(is.na(phen$nff),8,nff)
    symb <- ifelse(is.na(phen$nff),16,1)
    points(phen$TT, phen$HS,pch=symb,col=coul)
    points(phen$TT, phen$SSI, pch=symb, col=coul,cex=0.7)
     points(phen$TT, phen$GL, pch=symb, col=coul)
  }
})
#
# Export
#
write.csv(do.call('rbind', sapply(names(phend), function(g) {phen <- phend[[g]]; phen$label = g; phen[,c('Date', 'label', 'N', 'nff', 'TT', 'HS', 'SSI', 'GL')]}, simplify=FALSE)), 'Compil_Pheno_treated_archi_sampled.csv',row.names=FALSE)
