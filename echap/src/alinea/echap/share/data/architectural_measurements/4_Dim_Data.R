#
#
#           Compilations donees dimensions
#
#
# extraction raw dimensions per rank from notations files
#
Lsheath <- lapply(notdb, function(x) dim_notations(x,'sheath_length', 'sheath'))
Linternode <- lapply(notdb, function(x) dim_notations(x,'internode_length', 'internode'))
Hcol <- lapply(notdb, function(x) dim_notations(x,'Hcol_', 'col'))
Lbnot <- lapply(notdb, function(x) dim_notations(x,'blade_length', 'blade'))
#
#inspect
view_notdim(Lsheath)
view_notdim(Linternode)
view_notdim(Hcol,c(0,80))
view_notdim(Lbnot)
#
# extraction hcol,lmax from silhouettes
#
hclbc <- lapply(curvdb, HcLbCurv)
# compare Hcol to Hcol from notations for Tremie12
hcnot <- Hcol$Tremie12
Lbnot <- Lbnot$Tremie12
new <- hclbc$Tremie12
new$Source <- sapply(new$Source,function(x) sub('silhouette','scanned',x))
comph <- merge(hcnot,new)
}

#
# Dimensions tagged plants
#
# initiatlise with blade length from dynamic notations
dimtdb <- lapply(tagged, function(x) {df <- dimTagged(x);df$Lb <- df$L;df[,-grep('L$',colnames(df))]})
# sheath, internode, col from final sampling (always occcured after full expansion)
for (g in genos) {
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Lsheath[[g]], 'Ls')
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Linternode[[g]], 'Li')
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Hcol[[g]], 'Hc')
}
# check
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(dimtdb, function(dim) {
  plot(c(0,15),c(0,70),type='n')
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Lb,col=x$N,pch=16,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Ls,col=x$N,pch=1,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Li,col=x$N,pch=16,cex=0.7,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Hc,col=x$N,pch=1,cex=0.7,type='b'))
})
#
# Export
#
write.csv(do.call('rbind', sapply(genos, function(g) {dim <- dimtdb[[g]]; dim$label = g; dim}, simplify=FALSE)), 'Compil_Dim_treated_archi_tagged.csv',row.names=FALSE)
#
# Dimension data from other detructive samplings
# ----------------------------------------------
#
# scanned leaves
scans12 <- scandb$Tremie12[scandb$Tremie12$id_Axe=='MB',1:10]
scans12$lmax <- ifelse(scans12$stat < 3, scans12$lmax,NA)
scans12$wmax <- ifelse(scans12$stat < 2, scans12$wmax,NA)
scans12$A_bl<- ifelse(scans12$stat < 2, scans12$A_bl,NA)
phen12 <- phend$Tremie12[grep('scanned',phend$Tremie12$Source),]
dat <- merge(scans12,phen12)
# add dimensions from notations

dat <- dat[dat$rank <=dat$Nflig,]

Lb12 <- data.frame(Source=dat$Source, N=dat$N,nff=dat$nff,organ='blade',rank=dat$rank, L=dat$lmax)
Lb12 <- Lb12[!is.na(Lb12$L),]

#
# plant level variables
# ToDo : add area/cumul length par plant from scanned data
plant_dim <- lapply(notdb, function(x) plant_notations(x)) 
