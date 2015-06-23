#
#
#           Compilations donees dimensions
#
#
# extraction raw dimensions per rank from notations files
#
Lsheath <- lapply(notdb, function(x) dim_notations(x,'sheath_length', 'sheath'))
Linternode <- lapply(notdb, function(x) dim_notations(x,'internode_length', 'internode'))
#Hcol are directe notation oof collar heights
Hcol <- lapply(notdb, function(x) dim_notations(x,'Hcol_', 'col'))
#Hins are internode + sheath lengths from notation of internode+sheath dimensions
Hins <- lapply(notdb, function(x) dim_notations(x,'Hins_', 'col'))
Lbnot <- lapply(notdb, function(x) dim_notations(x,'blade_length', 'blade'))
# extraction plant level notations main stem
msdim <- lapply(notdb, function(x) plant_notations(x))
#
# Visual check notations
view_dim(Lsheath)
view_dim(Linternode)
view_dim(Hcol,ylim=c(0,80))
view_dim(Hins,ylim=c(0,80))
view_dim(Lbnot)
#
# extraction hcol,Lb from silhouttes
# add rank to Tremie13 data
TT <- TTlin$Tremie13$TT[TTlin$Tremie13$Date=='29/04/2013']
HS <- (TT - TToM$Tremie13) * slopeM$Tremie13
curvdb$Tremie13$rank <- HS - curvdb$Tremie13$ranktop + 1
hclbc <- lapply(curvdb, HcLbCurv)
# add nff if known
nff <- lapply(msdim,function(x) x[,c('Source','N','nff')])
for (g in names(hclbc))
  hclbc[[g]] <- merge(hclbc[[g]], nff[[g]], all.x=TRUE)
#
view_dim(hclbc,'Hcol',ylim=c(0,80))
view_dim(hclbc,'Lb')
#
# extraction of main stem dimension data from scandb
scandim <- lapply(scandb, function(x) {
  cols <- c('Source', 'N','id_Axe', 'rank','stat','lmax','wmax','A_bl')
  x[x$id_Axe=='MB',cols[cols%in%colnames(x)]]})
#
# Dimensions tagged plants
#
# initiatlise with blade length from dynamic notations
dimtdb <- lapply(tagged, function(x) {df <- dimTagged(x);df$Lb <- df$L;df[,-grep('L$',colnames(df))]})
# sheath, internode, col from notations on final sampling (always occcured after full expansion)
for (g in genos) {
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Lsheath[[g]], 'Ls')
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Linternode[[g]], 'Li')
  dimtdb[[g]] <- add_dim(dimtdb[[g]], Hcol[[g]], 'Hc')
  for (w in c('Lsheath','Linternode','Hcol')) {
    data <- get(w)
    if (length(grep('tagged', as.character(data[[g]]$Source))) > 0)
      data[[g]] <- data[[g]][-grep('tagged', as.character(data[[g]]$Source)),]
    if (!is.null(data[[g]]))
      if (nrow(data[[g]]) <= 0)
        data[[g]] <- NULL
    assign(w,data)
  }
}
# check
par(mfrow=c(2,2),mar=c(4,4,1,1))
lapply(dimtdb, function(dim) {
  plot(c(0,15),c(0,50),type='n')
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Lb,col=x$N,pch=16,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Ls,col=x$N,pch=1,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Li,col=x$N,pch=16,cex=0.7,type='b'))
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Hc/2,col=x$N,pch=1,cex=0.7,type='b'))
})
#
# Export
#
write.csv(do.call('rbind', sapply(genos, function(g) {dim <- dimtdb[[g]]; dim$label = g; dim}, simplify=FALSE)), 'Compil_Dim_treated_archi_tagged.csv',row.names=FALSE)
#
# Dimension data from other detructive samplings
# ----------------------------------------------
#
# Consolidation of Tremie12 data by crossing of sources
#
# scanned blade length versus notation blade length
comp <- merge(Lbnot$Tremie12,scandim$Tremie12)[,c('Source','N','rank','L','lmax')]
plot(comp$lmax,comp$L,xlim=c(0,30), ylim=c(0,30))
head(comp[order(-abs(comp$L-comp$lmax)),],20)
#conc : scan data okay for blade data, Notation are simply redundant : no merge
# scanned blade lengths versus silhouette
comp <- merge(hclbc$Tremie12,scandim$Tremie12)[,c('Source','N','rank','Lb','lmax')]
plot(comp$lmax,comp$Lb,xlim=c(0,30), ylim=c(0,30))
head(comp[order(-abs(comp$Lb-comp$lmax)),],20)
#conc : silhouettes data are too noisy for blade length estimation : we do not use them
#
# hcol vesus hins (sheath + length)
hcnot <- Hcol$Tremie12
hinot <- Hins$Tremie12
hinot$Lins <- hinot$L
hinot <- hinot[,-grep('L$',colnames(hinot))]
comp <- merge(hinot,hcnot)
plot(comp$Lins,comp$L, xlim=c(0,80), ylim=c(0,80))
head(comp[order(-abs(comp$Lins-comp$L)),],20)
# conc : for some plants on 9/5/12, Hcol notation seems to be Hcol + leaf length : discard data
who <- comp[abs(comp$Lins-comp$L) > 5, c('Source','N','nff','organ','rank','L')]
Hcol$Tremie12 <- rbind(who,Hcol$Tremie12)
Hcol$Tremie12 <- Hcol$Tremie12[!duplicated(Hcol$Tremie12),][-(1:nrow(who)),]
#check
hcnot <- Hcol$Tremie12
comp <- merge(hinot,hcnot)
plot(comp$Lins,comp$L, xlim=c(0,80), ylim=c(0,80))
#
# hcol versus silhouettes
sil <- hclbc$Tremie12
comp <- merge(hcnot,sil)
plot(comp$L,comp$Hcol, xlim=c(0,80), ylim=c(0,80))
head(comp[order(-abs(comp$Hcol-comp$L)),],20)
#
# hins versus silhouettes
comp <- merge(hinot,sil)
plot(comp$Lins,comp$Hcol, xlim=c(0,80), ylim=c(0,80))
head(comp[order(-abs(comp$Lins-comp$Hcol)),],20)
# conc: Hins data seems better, then Hcol, then silhouette
#
# Compil Hcol data
Hcdb <- Hins
# add hcol data
Hcdb$Tremie13 <- Hcol$Tremie13
hcol <- Hcol$Tremie12
hcol$Lhcol <- hcol$L
hcol <- hcol[,-grep('L$',colnames(hcol))]
compil <- merge(Hcdb$Tremie12,hcol,all=TRUE)
compil$L <- ifelse(is.na(compil$L),compil$Lhcol,compil$L)
Hcdb$Tremie12 <- compil
# add silhouette data
for (g in names(hclbc)) {
  sil <- hclbc[[g]][,c('Source','N','rank','Hcol')]
  sil$organ <- 'col'
  compil <- merge(Hcdb[[g]], sil, all=TRUE)
  compil$L <- ifelse(is.na(compil$L),compil$Hcol,compil$L)
  Hcdb[[g]] <- compil
}
#inspect data
view_notdim(Hcdb,c(0,80))
#
# Compilation of sampled dimensions
#
# initialise with filtered scan data

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
