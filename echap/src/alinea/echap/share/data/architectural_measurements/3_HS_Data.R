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
