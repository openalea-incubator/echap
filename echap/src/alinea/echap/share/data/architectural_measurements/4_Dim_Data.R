#
#
#           Compilations donees dimensions
#
#
# extraction raw dimensions per rank from notations files
#
Lsheath <- lapply(notdb, function(x) dim_notations(x,'sheath_length'))
Linternode <- lapply(notdb, function(x) dim_notations(x,'internode_length'))
#Hcol are directe notation oof collar heights
Hcol <- lapply(notdb, function(x) dim_notations(x,'Hcol_'))
#Hins are internode + sheath lengths from notation of internode+sheath dimensions
Hins <- lapply(notdb, function(x) dim_notations(x,'Hins_'))
Lbnot <- lapply(notdb, function(x) dim_notations(x,'blade_length'))
# extraction plant level notations main stem
msdb <- lapply(notdb, function(x) plant_notations(x))
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
stem <- curvdb$Tremie13$ranktop == 0
curvdb$Tremie13$rank <- HS - curvdb$Tremie13$ranktop + 1
curvdb$Tremie13$rank[stem] <- 0
hclbc <- lapply(curvdb, HcLbCurv)
#
view_dim(hclbc,'Hcol',ylim=c(0,80))
view_dim(hclbc,'Lb')
#
# extraction of main stem dimension data from scandb
scandim <- lapply(scandb, function(x) {
  cols <- c('Source', 'N','id_Axe', 'rank','stat','lmax','wmax','A_bl')
  x[x$id_Axe=='MB',cols[cols%in%colnames(x)]]})
#
# Consolidation of Tremie12 data by crossing of sources
#
# scanned blade length versus notation blade length
comp <- merge(Lbnot$Tremie12,scandim$Tremie12)[,c('Source','N','rank','L','lmax')]
plot(comp$lmax,comp$L,xlim=c(0,30), ylim=c(0,30))
head(comp[order(-abs(comp$L-comp$lmax)),],20)
#conc : scan data okay for blade data, Notation are redundant : no merge done
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
who <- comp[abs(comp$Lins-comp$L) > 5, c('Source','N','rank','L')]
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
# Dimensions tagged plants
#
# initiatlise with blade length from dynamic notations
dimtdb <- lapply(tagged, function(x) {df <- dimTagged(x);df$Lb <- df$L;df[,-grep('L$',colnames(df))]})
# sheath, internode, col from notations on final sampling (always occcured after full expansion)
for (g in genos) {
  dimtdb[[g]]$Source <- 'tagged_plants'
  dimtdb[[g]]$Wb <- NA
  dimtdb[[g]]$Ab <- NA
  dimtdb[[g]] <- add_dimt(dimtdb[[g]], Lsheath[[g]], 'Ls')
  dimtdb[[g]] <- add_dimt(dimtdb[[g]], Linternode[[g]], 'Li')
  dimtdb[[g]] <- add_dimt(dimtdb[[g]], Hcdb[[g]], 'Hc')
  #remove data from sampled dimensions objects
  for (w in c('Lsheath','Linternode','Hcdb')) {
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
  lapply(split(dim,dim$N), function(x) points(x$rank,x$Hc/2,col=x$N,pch=1,cex=0.7,type='p'))
})
#
# Export
#
cols <- c('label', 'Source', 'N', 'nff', 'rank', 'Lb','Wb','Ab','Ls','Li','Hc')
write.csv(do.call('rbind', sapply(genos, function(g) {dim <- dimtdb[[g]]; dim$label = g; dim[,cols]}, simplify=FALSE)), 'Compil_Dim_treated_archi_tagged.csv',row.names=FALSE)
#
# Dimension data from other detructive samplings
# ----------------------------------------------
#
# Compilation of sampled dimensions
#
# initialise with filtered scan data
dimsdb <- lapply(scandim,function(dim) {
  dat <- dim[,c('Source','N','rank')]
  dat$Lb <- NA
  if ('lmax' %in% colnames(dim))
    dat$Lb <- ifelse(dim$stat < 3, dim$lmax, NA)
  dat$Wb <- NA
  if ('wmax' %in% colnames(dim))
    dat$Wb <- ifelse(dim$stat < 2, dim$wmax, NA)
  dat$Ab <- NA
  if ('A_bl' %in% colnames(dim))
    dat$Ab <- ifelse(dim$stat < 2, dim$A_bl, NA)
  dat
})
# add organ dimensions
dimsdb <- add_dims(dimsdb, Lsheath, 'Ls')
dimsdb <- add_dims(dimsdb, Linternode, 'Li')
dimsdb <- add_dims(dimsdb, Hcdb, 'Hc')
# add plant notations (nff, Nflig) and estimate HS
for (g in names(dimsdb)) {
  dimsdb[[g]] <- merge(dimsdb[[g]],msdb[[g]][,c('Source', 'N', 'nff', 'Nflig')], all.x=TRUE)
  dates <- sapply(dimsdb[[g]]$Source, function(s) format(as.Date(strsplit(s,split='_')[[1]][3], '%d%m%y'),'%d/%m/%Y'))
  TT <- TTlin[[g]]$TT[match(dates,TTlin[[g]]$Date)]
  dimsdb[[g]]$HS <- slopeM[[g]] * (TT - TToM[[g]])
}
# filter growing organs using Nflig, HS, remove empty lines
dimsdb <- lapply(dimsdb, function(dim) {
  # leaves : keep ligulated ones (observed or estimated)
  for (w in c('Lb','Wb', 'Ab'))
    dim[[w]] <- ifelse(is.na(dim$Nflig), ifelse(dim$rank <= dim$HS, dim[[w]], NA),ifelse(dim$rank <= dim$Nflig, dim[[w]], NA))
  # sheath keep if rank + 0.3 <= HS
  dim$Ls <- ifelse((dim$rank + 0.3) <= dim$HS, dim$Ls, NA)
  # internode, col : keep if rank + 1.6 <= HS or rank < 8 
  for (w in c('Li','Hc'))
    dim[[w]] <-  ifelse(((dim$rank + 1.6) <= dim$HS) | (dim$rank < 8), dim[[w]], NA)
  # filter empty lines
  empty <- sapply(seq(nrow(dim)), function(irow) all(is.na(unlist(dim[irow,c('Lb','Wb','Ab','Ls','Li','Hc')]))))
  dim[!empty,]
})
#
# check
#
view_dim(dimsdb,'Lb')
view_dim(dimsdb,'Wb', c(0,3))
view_dim(dimsdb,'Ab', c(0,45))
view_dim(dimsdb,'Ls')
view_dim(dimsdb,'Li')
view_dim(dimsdb,'Hc', c(0,80))
#
# Export
#
cols <- c('label', 'Source', 'N', 'nff', 'rank', 'Lb','Wb','Ab','Ls','Li','Hc')
write.csv(do.call('rbind', sapply(names(dimsdb), function(g) {dim <- dimsdb[[g]]; dim$label = g; dim[,cols]}, simplify=FALSE)), 'Compil_Dim_treated_archi_sampled.csv',row.names=FALSE)
#
# plant level variables
#
# Stem diameters
# 
diams <- lapply(msdb, function(dim) dim[!is.na(dim$Daxe_mm),c('Source','N','nff','Daxe_mm')])
write.csv(do.call('rbind', sapply(names(diams), function(g) {dim <- diams[[g]]; dim$d=dim$Daxe_mm / 10;dim$label = g; dim[,c('label','Source','N','nff','d')]}, simplify=FALSE)), 'Compil_Diameter.csv',row.names=FALSE)
#
# ToDo : cumul area/cumul length par plant from scanned data
