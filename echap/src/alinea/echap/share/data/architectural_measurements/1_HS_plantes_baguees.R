#
#          Preprocessing 1: generation of HS data per plant on tagged plants
#
source('preprocessing_tools.R')
source('0_input.R')
#
# compute mature leaf lengths per leaf per plant
#
dimP <- lapply(tagged, dimTagged)
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(dimP, function(dim) {
  plot(c(0,15),c(0,30),type='n')
  lapply(split(dim,dim$N), function(x) points(x$rank,x$L,col=x$N,pch=16,type='b'))
})
#
# add interpolated-scaled length to yield more HS data
#
dimPM <- lapply(dimP, function(x) aggregate(x$L,list(nff=x$nff, rank=x$rank),mean, na.rm=TRUE))
dimPMM <- lapply(dimP, function(x) aggregate(x$L,list(rank=x$rank),mean, na.rm=TRUE))
#
dimPfit <- mapply(function(dim,name) do.call('rbind', lapply(split(dim, dim$N), function(mat) final_length(mat,dimPM[[name]], dimPMM[[name]]))), dimP, names(dimP), SIMPLIFY=FALSE)
#
# HS,SSI per plant
#
phendb <- sapply(genos, function(g) pheno(tagged[[g]], dimPfit[[g]]), simplify=FALSE)
#
# ajout estimation SSI terrain 27/04 Mercia/Rht3
#
SSIobsMercia <- data.frame(Date='27/04/2011', N=c(8,11), SSI=c(7.3, 8.2))
SSIobsRht3 <- data.frame(Date='27/04/2011', N=c(1,4,5,6,8,11), SSI=c(7.4, 6.3,6.2,7.1,7.4,6.2))
phendb$Mercia$SSI[match(paste(SSIobsMercia$Date,SSIobsMercia$N),paste(as.character(phendb$Mercia$Date), phendb$Mercia$N))] <- SSIobsMercia$SSI
phendb$Rht3$SSI[match(paste(SSIobsRht3$Date,SSIobsRht3$N),paste(as.character(phendb$Rht3$Date), phendb$Rht3$N))] <- SSIobsRht3$SSI
phendb$Mercia$GL = phendb$Mercia$HS - phendb$Mercia$SSI
phendb$Rht3$GL = phendb$Rht3$HS - phendb$Rht3$SSI
#
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(phendb, function(p) {
  plot(c(0,2500),c(0,13),type='n')
  lapply(split(p,p$N), function(x) points(x$TT,x$HS,col=x$nff,pch=16))
  lapply(split(p,p$N), function(x) points(x$TT,x$GL,col=x$nff,pch=16))
})
